# Quick benchmark: glmnet solver vs EGLearn (no C++ compilation needed)
# Usage: Rscript --vanilla run/quick_bench_glmnet_only.R

# Run from repo root: Rscript --vanilla run/quick_bench_glmnet_only.R

library(graphicalExtremes)
library(igraph)
library(glmnet)

source("R/glasso_glmnet.R")
source("R/f1score.R")
source("R/dgp.R")
source("R/eglearn.R")

# est_W in pure R (no C++ needed)
est_W_R <- function(data, q.threshold) {
  n <- nrow(data); d <- ncol(data)
  # Rank transform
  mdata <- apply(data, 2, function(x) 1 / (1 - rank(x) / (n + 1)))
  threshold_val <- 1 / (1 - q.threshold)

  cov_mat <- matrix(0, d, d)
  trace_sum <- 0
  for (k in 1:d) {
    keep <- which(mdata[, k] > threshold_val)
    if (length(keep) < 2) next
    w <- log(mdata[keep, -k]) - log(mdata[keep, k])
    cov_k <- cov(w)
    cov_mat[-k, -k] <- cov_mat[-k, -k] + cov_k / d
    trace_sum <- trace_sum + sum(cov_k) / (d^3)
  }
  # Trim
  Winv <- tryCatch(solve(cov_mat), error = function(e) NULL)
  if (!is.null(Winv)) {
    c_max <- 1 / sum(Winv)
    c_use <- min(trace_sum, c_max)
  } else {
    c_use <- trace_sum
  }
  trimmed <- cov_mat - c_use * matrix(1, d, d)
  # Correlation matrix
  sd_inv <- 1 / sqrt(diag(trimmed))
  result <- trimmed * (sd_inv %o% sd_inv)
  list(cov = result)
}

glasso_path_R <- function(data, lambda_range, q_threshold, iter_max = 200) {
  d <- ncol(data)
  W <- est_W_R(data, q_threshold)
  c_val <- 1 / (d * eigen(W$cov, symmetric = TRUE, only.values = TRUE)$values[1])
  results <- vector("list", length(lambda_range))
  for (i in seq_along(lambda_range)) {
    results[[i]] <- glasso_c_reest_glmnet(W$cov, lambda_range[i], c_val, iter_max)
  }
  list(W = W, c = c_val, lambda_range = lambda_range, results = results)
}

F1_graph <- function(g_true, g_est) {
  n_inter <- ecount(igraph::intersection(g_true, g_est))
  n_fp <- ecount(g_est) - n_inter
  n_fn <- ecount(g_true) - n_inter
  prec <- if (n_inter + n_fp > 0) n_inter / (n_inter + n_fp) else 0
  rec  <- if (n_inter + n_fn > 0) n_inter / (n_inter + n_fn) else 0
  if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
}

cat("=== Quick Benchmark: glmnet EGLasso vs EGLearn ===\n\n")

configs <- list(
  list(d = 20,  m = 2, k_mult = 5,   dist = "mpareto"),
  list(d = 100, m = 2, k_mult = 5,   dist = "mpareto"),
  list(d = 100, m = 1, k_mult = 2.5, dist = "mpareto"),
  list(d = 100, m = 2, k_mult = 5,   dist = "maxstable")
)

lambda_range <- 10^seq(-1, 1, by = 0.2)  # 11 values
rho_range <- seq(0.025, 0.525, by = 0.05) # 11 values
nsim <- 3
seed <- 42
results_list <- list()

for (ci in seq_along(configs)) {
  cfg <- configs[[ci]]
  d <- cfg$d; m <- cfg$m
  k <- round(d * cfg$k_mult)
  n <- floor(k^(1 / 0.7))
  q_threshold <- 1 - k / n

  tag <- sprintf("BA(%d) d=%d k=%d %s", m, d, k, cfg$dist)
  cat(sprintf("--- %s (n=%d, q=%.3f) ---\n", tag, n, q_threshold))

  set.seed(seed + ci)
  graph_info <- generate_graph("BA", d, m = m)
  Theta_true <- graph_info$Theta
  Gamma <- graph_info$Gamma
  true_graph <- graph_info$graph
  cat(sprintf("  True edges: %d\n", ecount(true_graph)))

  for (rep in seq_len(nsim)) {
    set.seed(seed + ci * 100 + rep)
    data <- generate_data(cfg$dist, n, d, Gamma)

    # --- glmnet EGLasso ---
    t0 <- proc.time()
    path_glm <- tryCatch(
      glasso_path_R(data, lambda_range, q_threshold),
      error = function(e) { cat(sprintf("    glmnet error: %s\n", e$message)); NULL })
    time_glm <- (proc.time() - t0)[3]
    if (!is.null(path_glm)) {
      f1_glm <- sapply(path_glm$results, function(r) F1score(Theta_true, r$Theta_hat))
      best_f1_glm <- max(f1_glm, na.rm = TRUE)
    } else {
      best_f1_glm <- NA
    }

    # --- EGLearn-ns ---
    t0 <- proc.time()
    fit_eg <- tryCatch(
      eglearn2(data, p = q_threshold, rholist = rho_range, reg_method = "ns"),
      error = function(e) { cat(sprintf("    EGLearn error: %s\n", e$message)); NULL })
    time_eg <- (proc.time() - t0)[3]
    if (!is.null(fit_eg)) {
      f1_eg <- sapply(fit_eg$graph, function(g) F1_graph(true_graph, g))
      best_f1_eg <- max(f1_eg, na.rm = TRUE)
    } else {
      best_f1_eg <- NA
    }

    cat(sprintf("  rep %d: glmnet=%.1fs (F1=%.3f) | EGLearn=%.1fs (F1=%.3f) | ratio=%.1fx\n",
                rep, time_glm, best_f1_glm, time_eg, best_f1_eg,
                time_eg / max(time_glm, 0.01)))

    results_list[[length(results_list) + 1]] <- data.frame(
      setting = tag, rep = rep,
      time_glm = time_glm, f1_glm = best_f1_glm,
      time_eg = time_eg, f1_eg = best_f1_eg
    )
  }
  cat("\n")
}

df <- do.call(rbind, results_list)
cat("=== Summary (means across reps) ===\n")
agg <- aggregate(cbind(time_glm, f1_glm, time_eg, f1_eg) ~ setting,
                 data = df, FUN = mean, na.rm = TRUE)
agg$speed_ratio <- agg$time_eg / agg$time_glm
print(agg, digits = 3)
