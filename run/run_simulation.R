# Multi-sample parallel simulation
# Runs nsim replications for a fixed DGP and lambda range, computes F1 scores

source("R/load_all.R")
source("R/parallel.R")
source("R/output.R")

# --- Configuration ---
config <- list(
  d = 20,
  n = 5000,
  q_threshold = 0.9,
  graph_type = "BA",
  m = 1,
  distribution = "mpareto",
  nsim = 100,
  ncores = 10,
  seed = 100,
  iter_max = 1000
)
lambda_range <- 10^seq(-2, 0, by = 0.1)

# --- Output ---
output_dir <- create_output_dir(
  prefix = paste0("sim_", config$graph_type, "_d", config$d))

save_config(c(config, list(
  lambda_range_log10 = paste(range(log10(lambda_range)), collapse = " to "),
  nlambda = length(lambda_range)
)), output_dir)

cat("Output directory:", output_dir, "\n")

# --- Generate graph (shared across replications) ---
set.seed(config$seed)
graph_info <- generate_graph(config$graph_type, config$d, m = config$m)
Theta_true <- graph_info$Theta
Gamma <- graph_info$Gamma
true_edges <- ecount(graph_info$graph)

cat("Graph:", vcount(graph_info$graph), "nodes,", true_edges, "edges\n")
cat("Running", config$nsim, "replications on", config$ncores, "cores...\n")

# --- Define single-replication function ---
# (this will be called on each worker)
one_replication <- function(rep_id) {
  data <- generate_data(config$distribution, config$n, config$d, Gamma)
  path <- glasso_path(data, lambda_range, config$q_threshold,
                      config$iter_max)

  f1 <- sapply(path$results, function(res) F1score(Theta_true, res$Theta_hat))
  n_edges <- sapply(path$results, function(res) {
    adj <- abs(res$Theta_hat) > 0
    diag(adj) <- FALSE
    sum(adj) / 2
  })

  list(f1 = f1, n_edges = n_edges)
}

# --- Setup function for workers ---
worker_setup <- function() {
  source("R/load_all.R")
}

# --- Run parallel simulation ---
t0 <- proc.time()
results <- run_parallel(
  sim_fn = one_replication,
  nsim = config$nsim,
  ncores = config$ncores,
  seed = config$seed,
  setup_fn = worker_setup,
  export_vars = c("config", "Gamma", "Theta_true", "lambda_range"),
  export_env = environment()
)
elapsed <- (proc.time() - t0)[3]
cat("Total elapsed:", round(elapsed, 1), "s\n")

# --- Aggregate results ---
f1_matrix <- do.call(rbind, lapply(results, `[[`, "f1"))  # nsim x nlambda
edges_matrix <- do.call(rbind, lapply(results, `[[`, "n_edges"))

f1_mean <- colMeans(f1_matrix)
f1_sd <- apply(f1_matrix, 2, sd)
f1_q25 <- apply(f1_matrix, 2, quantile, 0.25)
f1_q75 <- apply(f1_matrix, 2, quantile, 0.75)

best_idx <- which.max(f1_mean)
cat("\nBest mean F1:", round(f1_mean[best_idx], 3),
    "+/-", round(f1_sd[best_idx], 3),
    "at lambda =", round(lambda_range[best_idx], 4), "\n")

# --- Plot ---
save_plot(function() {
  par(mfrow = c(1, 2))

  # F1 score with IQR band
  log_lam <- log10(lambda_range)
  plot(log_lam, f1_mean, type = "l", lwd = 2,
       xlab = expression(log[10](lambda)), ylab = "F1 score",
       main = paste0("F1 vs lambda (d=", config$d, ", n=", config$n,
                     ", nsim=", config$nsim, ")"),
       ylim = c(0, 1))
  polygon(c(log_lam, rev(log_lam)), c(f1_q25, rev(f1_q75)),
          col = rgb(0, 0, 1, 0.2), border = NA)
  lines(log_lam, f1_mean, lwd = 2)
  abline(v = log_lam[best_idx], col = "red", lty = 2)

  # Number of edges
  edges_mean <- colMeans(edges_matrix)
  plot(log_lam, edges_mean, type = "l", lwd = 2,
       xlab = expression(log[10](lambda)), ylab = "Mean edges",
       main = "Estimated edges vs lambda")
  abline(h = true_edges, col = "blue", lty = 2)
  legend("topright", legend = paste("True:", true_edges), col = "blue", lty = 2)
}, name = "simulation_results", output_dir = output_dir)

# --- Save ---
sim_results <- list(
  config = config,
  lambda_range = lambda_range,
  f1_matrix = f1_matrix,
  edges_matrix = edges_matrix,
  f1_mean = f1_mean,
  f1_sd = f1_sd,
  true_edges = true_edges,
  Theta_true = Theta_true,
  elapsed = elapsed
)
save_results(sim_results, "simulation_results", output_dir)

cat("Results saved to:", output_dir, "\n")
