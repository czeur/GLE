# Compare glasso_c vs EGLearn
# Same data, same DGP. Compare F1 vs tuning parameter and timing.

source("R/load_all.R")

# --- Configuration ---
d <- 20
n <- 5000
q_threshold <- 0.9
graph_type <- "BA"
m <- 1
distribution <- "mpareto"
seed <- 42

# Tuning ranges (adjust manually to cover the F1 peak)
lambda_range <- 10^seq(-1.5, 0, by = 0.15)  # ~10 values for glasso_c
rho_range <- seq(0.05, 0.5, length.out = 10) # for EGLearn

# --- Generate data (shared) ---
set.seed(seed)
graph_info <- generate_graph(graph_type, d, m = m)
data <- generate_data(distribution, n, d, graph_info$Gamma)
true_graph <- graph_info$graph
true_edges <- ecount(true_graph)

cat("d =", d, " n =", n, " true edges =", true_edges, "\n\n")

# --- F1 from igraph objects ---
F1_graph <- function(g_true, g_est) {
  n_inter <- ecount(igraph::intersection(g_true, g_est))
  n_fp <- ecount(g_est) - n_inter
  n_fn <- ecount(g_true) - n_inter
  prec <- if (n_inter + n_fp > 0) n_inter / (n_inter + n_fp) else 0
  rec  <- if (n_inter + n_fn > 0) n_inter / (n_inter + n_fn) else 0
  if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
}

# Convert glasso_c result to igraph for consistent comparison
theta_to_graph <- function(Theta_hat, d) {
  adj <- Theta_hat != 0
  diag(adj) <- FALSE
  graph_from_adjacency_matrix(adj * 1, mode = "undirected", diag = FALSE)
}

# ============================
# Method 1: glasso_c
# ============================
cat("Running glasso_c (", length(lambda_range), " lambdas)...\n")
t_glasso <- proc.time()
path <- glasso_path(data, lambda_range, q_threshold)
t_glasso <- (proc.time() - t_glasso)[3]

f1_glasso <- sapply(path$results, function(res) {
  g_est <- theta_to_graph(res$Theta_hat, d)
  F1_graph(true_graph, g_est)
})

edges_glasso <- sapply(path$results, function(res) {
  sum(res$Theta_hat[upper.tri(res$Theta_hat)] != 0)
})

cat(sprintf("  Time: %.2f s\n", t_glasso))
best_g <- which.max(f1_glasso)
cat(sprintf("  Best F1: %.3f at lambda = %.4f (edges = %d)\n\n",
            f1_glasso[best_g], lambda_range[best_g], edges_glasso[best_g]))

# ============================
# Method 2: EGLearn (neighborhood selection)
# ============================
cat("Running EGLearn (", length(rho_range), " rhos)...\n")
t_eglearn <- proc.time()
fit_eg <- eglearn(data, p = q_threshold, rholist = rho_range, reg_method = "ns")
t_eglearn <- (proc.time() - t_eglearn)[3]

f1_eglearn <- sapply(fit_eg$graph, function(g) F1_graph(true_graph, g))

edges_eglearn <- sapply(fit_eg$graph, ecount)

cat(sprintf("  Time: %.2f s\n", t_eglearn))
best_e <- which.max(f1_eglearn)
cat(sprintf("  Best F1: %.3f at rho = %.4f (edges = %d)\n\n",
            f1_eglearn[best_e], rho_range[best_e], edges_eglearn[best_e]))

# ============================
# Summary
# ============================
cat("=== Summary ===\n")
cat(sprintf("%-12s  Best F1   Time\n", "Method"))
cat(sprintf("%-12s  %.3f     %.2f s\n", "glasso_c", max(f1_glasso), t_glasso))
cat(sprintf("%-12s  %.3f     %.2f s\n", "EGLearn", max(f1_eglearn), t_eglearn))
cat(sprintf("Speedup: %.1fx\n", t_eglearn / t_glasso))

# ============================
# Plot
# ============================
dir.create("output", showWarnings = FALSE)
pdf("output/compare_glasso_eglearn.pdf", width = 12, height = 5)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))

# Panel 1: F1 vs lambda (glasso_c)
plot(log10(lambda_range), f1_glasso, type = "b", lwd = 2, pch = 19,
     xlab = expression(log[10](lambda)), ylab = "F1 score",
     main = "glasso_c", ylim = c(0, 1), col = "blue")
abline(v = log10(lambda_range[best_g]), col = "red", lty = 2)
text(log10(lambda_range[best_g]), f1_glasso[best_g] + 0.05,
     sprintf("F1=%.3f", f1_glasso[best_g]), col = "red", cex = 0.8)
mtext(sprintf("Time: %.2f s", t_glasso), side = 3, line = 0, cex = 0.7)

# Panel 2: F1 vs rho (EGLearn)
plot(rho_range, f1_eglearn, type = "b", lwd = 2, pch = 19,
     xlab = expression(rho), ylab = "F1 score",
     main = "EGLearn (NS)", ylim = c(0, 1), col = "darkgreen")
abline(v = rho_range[best_e], col = "red", lty = 2)
text(rho_range[best_e], f1_eglearn[best_e] + 0.05,
     sprintf("F1=%.3f", f1_eglearn[best_e]), col = "red", cex = 0.8)
mtext(sprintf("Time: %.2f s", t_eglearn), side = 3, line = 0, cex = 0.7)

# Panel 3: Edge count comparison
plot(log10(lambda_range), edges_glasso, type = "b", lwd = 2, pch = 19,
     xlab = "Tuning parameter index", ylab = "Number of edges",
     main = "Edge count", ylim = range(c(edges_glasso, edges_eglearn, true_edges)),
     col = "blue", xaxt = "n")
axis(1, at = log10(lambda_range), labels = seq_along(lambda_range))
points(log10(lambda_range[1]) + (seq_along(rho_range) - 1) *
         diff(range(log10(lambda_range))) / (length(rho_range) - 1),
       edges_eglearn, type = "b", lwd = 2, pch = 17, col = "darkgreen")
abline(h = true_edges, col = "red", lty = 2)
legend("topright", legend = c("glasso_c", "EGLearn", paste("True:", true_edges)),
       col = c("blue", "darkgreen", "red"), pch = c(19, 17, NA),
       lty = c(1, 1, 2), lwd = 2, cex = 0.8)

dev.off()
cat("\nPlot saved to output/compare_glasso_eglearn.pdf\n")

# ============================
# Save results
# ============================
comparison <- list(
  config = list(d = d, n = n, q_threshold = q_threshold,
                graph_type = graph_type, m = m, seed = seed),
  glasso_c = list(lambda_range = lambda_range, f1 = f1_glasso,
                  edges = edges_glasso, time = t_glasso),
  eglearn = list(rho_range = rho_range, f1 = f1_eglearn,
                 edges = edges_eglearn, time = t_eglearn),
  true_edges = true_edges
)
saveRDS(comparison, "output/compare_glasso_eglearn.rds")
cat("Results saved to output/compare_glasso_eglearn.rds\n")
