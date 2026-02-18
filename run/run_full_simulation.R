# Full simulation: glasso_c vs EGLearn across all settings
# Run on remote machine with 15 cores
#
# Output: one RDS file per setting in output/full_sim/

source("R/load_all.R")
source("R/parallel.R")

# === Configuration ===
ncores <- 15
nsim <- 100
seed_base <- 42

lambda_range <- 10^seq(-1, 1, by = 0.1)   # 21 values
rho_range <- seq(0.025, 0.525, by = 0.025) # 21 values

# Build settings grid
settings <- list()
for (m_val in c(1, 2)) {
  k_mults <- if (m_val == 1) c(0.5, 1, 2.5) else c(0.5, 1, 5)
  for (d in c(20, 100)) {
    for (k_mult in k_mults) {
      k <- round(d * k_mult)
      n <- floor(k^(1 / 0.7))
      for (distribution in c("mpareto", "maxstable")) {
        settings[[length(settings) + 1]] <- list(
          m = m_val, d = d, k = k, n = n,
          q_threshold = 1 - k / n,
          distribution = distribution
        )
      }
    }
  }
}

cat(sprintf("Total settings: %d\n", length(settings)))
cat(sprintf("Samples per setting: %d\n", nsim))
cat(sprintf("Cores: %d\n", ncores))
cat(sprintf("Lambda values: %d, Rho values: %d\n\n",
            length(lambda_range), length(rho_range)))

# === F1 from igraph objects (for EGLearn) ===
F1_graph <- function(g_true, g_est) {
  n_inter <- ecount(igraph::intersection(g_true, g_est))
  n_fp <- ecount(g_est) - n_inter
  n_fn <- ecount(g_true) - n_inter
  prec <- if (n_inter + n_fp > 0) n_inter / (n_inter + n_fp) else 0
  rec  <- if (n_inter + n_fn > 0) n_inter / (n_inter + n_fn) else 0
  if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
}

# === Output directory ===
dir.create("output/full_sim", recursive = TRUE, showWarnings = FALSE)

# === Run each setting ===
for (si in seq_along(settings)) {
  s <- settings[[si]]
  tag <- sprintf("BA%d_%s_d%d_k%d", s$m, s$distribution, s$d, s$k)
  outfile <- sprintf("output/full_sim/%s.rds", tag)

  if (file.exists(outfile)) {
    cat(sprintf("[%d/%d] %s — already exists, skipping\n", si, length(settings), tag))
    next
  }

  cat(sprintf("[%d/%d] %s (d=%d, k=%d, n=%d, q=%.4f)\n",
              si, length(settings), tag, s$d, s$k, s$n, s$q_threshold))

  # Generate graph (fixed seed per setting)
  set.seed(seed_base + si)
  graph_info <- generate_graph("BA", s$d, m = s$m)
  true_graph <- graph_info$graph
  true_edges <- ecount(true_graph)
  Gamma <- graph_info$Gamma
  Theta_true <- graph_info$Theta

  cat(sprintf("  True edges: %d\n", true_edges))

  # One replication function
  nlambda <- length(lambda_range)
  nrho <- length(rho_range)

  one_rep <- function(rep_id) {
    # Generate data
    data <- generate_data(s$distribution, s$n, s$d, Gamma)

    # --- glasso_c ---
    res_g <- tryCatch({
      t0 <- proc.time()
      path <- glasso_path(data, lambda_range, s$q_threshold)
      time_glasso <- (proc.time() - t0)[3]
      f1_g <- sapply(path$results, function(r) F1score(Theta_true, r$Theta_hat))
      edges_g <- sapply(path$results, function(r) {
        adj <- r$Theta_hat != 0; diag(adj) <- FALSE; sum(adj) / 2
      })
      list(f1_g = f1_g, edges_g = edges_g, time_glasso = time_glasso)
    }, error = function(e) {
      list(f1_g = rep(NA, nlambda), edges_g = rep(NA, nlambda), time_glasso = NA)
    })

    # --- EGLearn ---
    res_e <- tryCatch({
      t0 <- proc.time()
      fit_eg <- eglearn(data, p = s$q_threshold, rholist = rho_range, reg_method = "ns")
      time_eglearn <- (proc.time() - t0)[3]
      f1_e <- sapply(fit_eg$graph, function(g) F1_graph(true_graph, g))
      edges_e <- sapply(fit_eg$graph, ecount)
      list(f1_e = f1_e, edges_e = edges_e, time_eglearn = time_eglearn)
    }, error = function(e) {
      list(f1_e = rep(NA, nrho), edges_e = rep(NA, nrho), time_eglearn = NA)
    })

    c(res_g, res_e)
  }

  # Worker setup
  worker_setup <- function() {
    source("R/load_all.R")
  }

  # Run in parallel
  t_total <- proc.time()
  results <- run_parallel(
    sim_fn = one_rep,
    nsim = nsim,
    ncores = ncores,
    seed = seed_base + si,
    setup_fn = worker_setup,
    export_vars = c("s", "Gamma", "Theta_true", "true_graph",
                    "lambda_range", "rho_range", "nlambda", "nrho", "F1_graph"),
    export_env = environment()
  )
  elapsed <- (proc.time() - t_total)[3]

  # Aggregate
  f1_glasso   <- do.call(rbind, lapply(results, `[[`, "f1_g"))
  f1_eglearn  <- do.call(rbind, lapply(results, `[[`, "f1_e"))
  edges_glasso  <- do.call(rbind, lapply(results, `[[`, "edges_g"))
  edges_eglearn <- do.call(rbind, lapply(results, `[[`, "edges_e"))
  time_glasso  <- sapply(results, `[[`, "time_glasso")
  time_eglearn <- sapply(results, `[[`, "time_eglearn")

  n_fail_g <- sum(is.na(time_glasso))
  n_fail_e <- sum(is.na(time_eglearn))

  best_g <- which.max(colMeans(f1_glasso, na.rm = TRUE))
  best_e <- which.max(colMeans(f1_eglearn, na.rm = TRUE))
  cat(sprintf("  glasso_c: mean F1=%.3f at log10(lam)=%.1f | mean time=%.2fs/sample",
              mean(f1_glasso[, best_g], na.rm = TRUE), log10(lambda_range[best_g]),
              mean(time_glasso, na.rm = TRUE)))
  if (n_fail_g > 0) cat(sprintf(" | %d failed", n_fail_g))
  cat("\n")
  cat(sprintf("  EGLearn:  mean F1=%.3f at rho=%.3f     | mean time=%.2fs/sample",
              mean(f1_eglearn[, best_e], na.rm = TRUE), rho_range[best_e],
              mean(time_eglearn, na.rm = TRUE)))
  if (n_fail_e > 0) cat(sprintf(" | %d failed", n_fail_e))
  cat("\n")
  cat(sprintf("  Total elapsed: %.0fs\n\n", elapsed))

  # Save
  sim_result <- list(
    config = list(
      m = s$m, d = s$d, k = s$k, n = s$n,
      q_threshold = s$q_threshold,
      distribution = s$distribution,
      nsim = nsim, ncores = ncores,
      seed_base = seed_base, setting_index = si
    ),
    lambda_range = lambda_range,
    rho_range = rho_range,
    f1_glasso = f1_glasso,       # nsim x 21
    f1_eglearn = f1_eglearn,     # nsim x 21
    edges_glasso = edges_glasso, # nsim x 21
    edges_eglearn = edges_eglearn,
    time_glasso = time_glasso,   # length nsim
    time_eglearn = time_eglearn,
    true_edges = true_edges,
    Theta_true = Theta_true,
    elapsed = elapsed
  )
  saveRDS(sim_result, outfile)
  cat(sprintf("  Saved: %s\n\n", outfile))
}

cat("All settings complete.\n")
