# Quick benchmark: 4 settings, 3 reps, comparing cpp vs glmnet vs EGLearn
# Usage: Rscript --vanilla run/quick_bench.R

setwd("/Users/czhou3/Documents/zhou/claudebot/projects/graphical/GLE")
source("R/load_all.R")
source("R/eglearn.R")

cat("=== Quick Benchmark: C++ Gram-CD vs glmnet vs EGLearn ===\n\n")

# 4 settings: (d, m, distribution)
configs <- list(
  list(d = 20,  m = 2, k_mult = 5,   dist = "mpareto"),
  list(d = 100, m = 2, k_mult = 5,   dist = "mpareto"),
  list(d = 100, m = 1, k_mult = 2.5, dist = "mpareto"),
  list(d = 100, m = 2, k_mult = 5,   dist = "maxstable")
)

lambda_range <- 10^seq(-1, 1, by = 0.2)  # 11 values (faster than 21)
rho_range <- seq(0.025, 0.525, by = 0.05) # 11 values
nsim <- 3
seed <- 42

F1_graph <- function(g_true, g_est) {
  n_inter <- ecount(igraph::intersection(g_true, g_est))
  n_fp <- ecount(g_est) - n_inter
  n_fn <- ecount(g_true) - n_inter
  prec <- if (n_inter + n_fp > 0) n_inter / (n_inter + n_fp) else 0
  rec  <- if (n_inter + n_fn > 0) n_inter / (n_inter + n_fn) else 0
  if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
}

results <- list()

for (ci in seq_along(configs)) {
  cfg <- configs[[ci]]
  d <- cfg$d; m <- cfg$m
  k <- round(d * cfg$k_mult)
  n <- floor(k^(1 / 0.7))
  q_threshold <- 1 - k / n

  tag <- sprintf("BA(%d) d=%d k=%d %s", m, d, k, cfg$dist)
  cat(sprintf("--- %s (n=%d) ---\n", tag, n))

  set.seed(seed + ci)
  graph_info <- generate_graph("BA", d, m = m)
  Theta_true <- graph_info$Theta
  Gamma <- graph_info$Gamma
  true_graph <- graph_info$graph

  for (rep in seq_len(nsim)) {
    set.seed(seed + ci * 100 + rep)
    data <- generate_data(cfg$dist, n, d, Gamma)

    # --- C++ Gram-CD ---
    t0 <- proc.time()
    path_cpp <- tryCatch(
      glasso_path(data, lambda_range, q_threshold, method = "cpp"),
      error = function(e) NULL)
    time_cpp <- (proc.time() - t0)[3]
    if (!is.null(path_cpp)) {
      f1_cpp <- sapply(path_cpp$results, function(r) F1score(Theta_true, r$Theta_hat))
      best_f1_cpp <- max(f1_cpp)
    } else {
      best_f1_cpp <- NA
    }

    # --- glmnet ---
    t0 <- proc.time()
    path_glm <- tryCatch(
      glasso_path(data, lambda_range, q_threshold, method = "glmnet"),
      error = function(e) NULL)
    time_glm <- (proc.time() - t0)[3]
    if (!is.null(path_glm)) {
      f1_glm <- sapply(path_glm$results, function(r) F1score(Theta_true, r$Theta_hat))
      best_f1_glm <- max(f1_glm)
    } else {
      best_f1_glm <- NA
    }

    # --- EGLearn-ns ---
    t0 <- proc.time()
    fit_eg <- tryCatch(
      eglearn2(data, p = q_threshold, rholist = rho_range, reg_method = "ns"),
      error = function(e) NULL)
    time_eg <- (proc.time() - t0)[3]
    if (!is.null(fit_eg)) {
      f1_eg <- sapply(fit_eg$graph, function(g) F1_graph(true_graph, g))
      best_f1_eg <- max(f1_eg)
    } else {
      best_f1_eg <- NA
    }

    cat(sprintf("  rep %d: cpp=%.1fs (F1=%.3f) | glmnet=%.1fs (F1=%.3f) | EGLearn=%.1fs (F1=%.3f)\n",
                rep, time_cpp, best_f1_cpp, time_glm, best_f1_glm, time_eg, best_f1_eg))

    results[[length(results) + 1]] <- data.frame(
      setting = tag, rep = rep,
      time_cpp = time_cpp, f1_cpp = best_f1_cpp,
      time_glm = time_glm, f1_glm = best_f1_glm,
      time_eg = time_eg, f1_eg = best_f1_eg
    )
  }
  cat("\n")
}

# Summary
df <- do.call(rbind, results)
cat("=== Summary (means across reps) ===\n")
agg <- aggregate(cbind(time_cpp, f1_cpp, time_glm, f1_glm, time_eg, f1_eg) ~ setting,
                 data = df, FUN = mean, na.rm = TRUE)
print(agg, digits = 3)
