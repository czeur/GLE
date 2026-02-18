# Single-sample trial: generate data, run glasso_path, plot F1 vs lambda
# Use this to explore lambda ranges for different DGPs before multi-sample simulations

source("R/load_all.R")

# --- Configuration ---
d <- 20
n <- 5000
q_threshold <- 0.9
graph_type <- "BA"
m <- 1  # edges per step in BA model
distribution <- "mpareto"
lambda_range <- 10^seq(-2, 0, by = 0.05)
seed <- 42

# --- Output directory ---
output_dir <- file.path("output", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(output_dir, recursive = TRUE)

# --- Generate DGP ---
set.seed(seed)
graph_info <- generate_graph(graph_type, d, m = m)
data <- generate_data(distribution, n, d, graph_info$Gamma)

cat("Graph:", vcount(graph_info$graph), "nodes,", ecount(graph_info$graph), "edges\n")
cat("Data:", nrow(data), "x", ncol(data), "\n")
cat("Running glasso_path across", length(lambda_range), "lambdas...\n")

# --- Estimate across all lambdas ---
t0 <- proc.time()
path <- glasso_path(data, lambda_range, q_threshold)
elapsed <- (proc.time() - t0)[3]
cat("Elapsed:", round(elapsed, 1), "s\n")

# --- Compute F1 at each lambda ---
f1_scores <- sapply(path$results, function(res) {
  F1score(graph_info$Theta, res$Theta_hat)
})

# Number of estimated edges at each lambda
n_edges <- sapply(path$results, function(res) {
  adj <- abs(res$Theta_hat) > 0
  diag(adj) <- FALSE
  sum(adj) / 2
})

true_edges <- ecount(graph_info$graph)

# --- Print summary ---
best_idx <- which.max(f1_scores)
cat("\nBest F1:", round(f1_scores[best_idx], 3),
    "at lambda =", round(lambda_range[best_idx], 4),
    "(log10 =", round(log10(lambda_range[best_idx]), 2), ")\n")
cat("True edges:", true_edges, " Estimated edges at best:", n_edges[best_idx], "\n")

# --- Plot F1 vs log10(lambda) ---
pdf(file.path(output_dir, "f1_vs_lambda.pdf"), width = 8, height = 5)
par(mfrow = c(1, 2))

# F1 score
plot(log10(lambda_range), f1_scores, type = "l", lwd = 2,
     xlab = expression(log[10](lambda)), ylab = "F1 score",
     main = paste0("F1 vs lambda (d=", d, ", n=", n, ", ", graph_type, ")"))
abline(v = log10(lambda_range[best_idx]), col = "red", lty = 2)
points(log10(lambda_range[best_idx]), f1_scores[best_idx],
       col = "red", pch = 19, cex = 1.5)

# Number of edges
plot(log10(lambda_range), n_edges, type = "l", lwd = 2,
     xlab = expression(log[10](lambda)), ylab = "Number of edges",
     main = "Estimated edges vs lambda")
abline(h = true_edges, col = "blue", lty = 2)
legend("topright", legend = paste("True:", true_edges), col = "blue", lty = 2)
dev.off()

# --- Save results ---
trial_results <- list(
  config = list(d = d, n = n, q_threshold = q_threshold,
                graph_type = graph_type, m = m, distribution = distribution,
                seed = seed),
  lambda_range = lambda_range,
  f1_scores = f1_scores,
  n_edges = n_edges,
  true_edges = true_edges,
  best_lambda = lambda_range[best_idx],
  best_f1 = f1_scores[best_idx],
  elapsed = elapsed,
  Theta_true = graph_info$Theta,
  path = path
)
saveRDS(trial_results, file.path(output_dir, "trial_results.rds"))

cat("\nResults saved to:", output_dir, "\n")
