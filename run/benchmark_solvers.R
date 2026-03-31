# Benchmark: compare C++ Gram-CD vs glmnet solver
# Runs a single setting and reports F1 + timing for both methods
#
# Usage: Rscript --vanilla run/benchmark_solvers.R

source("R/load_all.R")

cat("=== Solver Benchmark ===\n\n")

# Settings to benchmark
configs <- list(
  list(d = 20, m = 2, k_mult = 5),
  list(d = 100, m = 2, k_mult = 5)
)

lambda_range <- 10^seq(-1, 1, by = 0.1)
nsim <- 5  # small number for quick benchmarking
seed <- 42

for (cfg in configs) {
  d <- cfg$d
  m <- cfg$m
  k <- round(d * cfg$k_mult)
  n <- floor(k^(1 / 0.7))
  q_threshold <- 1 - k / n

  cat(sprintf("--- BA(%d), d=%d, k=%d, n=%d ---\n", m, d, k, n))

  # Generate graph
  set.seed(seed)
  graph_info <- generate_graph("BA", d, m = m)
  Theta_true <- graph_info$Theta
  Gamma <- graph_info$Gamma

  for (rep in seq_len(nsim)) {
    set.seed(seed + rep)
    data <- generate_data("mpareto", n, d, Gamma)

    # --- C++ Gram-CD ---
    t0 <- proc.time()
    path_cpp <- glasso_path(data, lambda_range, q_threshold, method = "cpp")
    time_cpp <- (proc.time() - t0)[3]
    f1_cpp <- sapply(path_cpp$results, function(r) F1score(Theta_true, r$Theta_hat))
    best_f1_cpp <- max(f1_cpp)

    # --- glmnet ---
    t0 <- proc.time()
    path_glm <- glasso_path(data, lambda_range, q_threshold, method = "glmnet")
    time_glm <- (proc.time() - t0)[3]
    f1_glm <- sapply(path_glm$results, function(r) F1score(Theta_true, r$Theta_hat))
    best_f1_glm <- max(f1_glm)

    cat(sprintf("  rep %d: cpp=%.2fs (F1=%.3f) | glmnet=%.2fs (F1=%.3f)\n",
                rep, time_cpp, best_f1_cpp, time_glm, best_f1_glm))
  }
  cat("\n")
}
