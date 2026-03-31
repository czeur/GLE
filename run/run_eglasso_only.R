# EGLasso-only simulation (rerun after M fix)
# EGLearn results from previous run (output/full_sim_3103/) are still valid
#
# Usage: nohup Rscript --vanilla run/run_eglasso_only.R > output_eglasso.log 2>&1 &
#        tail -f output_eglasso.log
#
# Output: one RDS file per setting in output/full_sim_0401/

log <- function(...) {
  cat(...); flush.console()
}

source("R/load_all.R")

# === Configuration ===
ncores <- 15
nsim <- 100
seed_base <- 42

lambda_range <- 10^seq(-1.2, 0, by = 0.1)  # 13 values

# Build settings grid (same as run_full_simulation.R)
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

log(sprintf("Total settings: %d\n", length(settings)))
log(sprintf("Samples per setting: %d\n", nsim))
log(sprintf("Cores: %d\n", ncores))
log(sprintf("Lambda values: %d\n\n", length(lambda_range)))

# === Output directory ===
dir.create("output/full_sim_0401", recursive = TRUE, showWarnings = FALSE)

# === Run each setting ===
for (si in seq_along(settings)) {
  s <- settings[[si]]
  tag <- sprintf("BA%d_%s_d%d_k%d", s$m, s$distribution, s$d, s$k)
  outfile <- sprintf("output/full_sim_0401/%s.rds", tag)

  if (file.exists(outfile)) {
    log(sprintf("[%d/%d] %s — already exists, skipping\n", si, length(settings), tag))
    next
  }

  log(sprintf("[%d/%d] %s (d=%d, k=%d, n=%d, q=%.4f)\n",
              si, length(settings), tag, s$d, s$k, s$n, s$q_threshold))

  # Generate graph (fixed seed per setting — same as full simulation)
  set.seed(seed_base + si)
  graph_info <- generate_graph("BA", s$d, m = s$m)
  true_graph <- graph_info$graph
  true_edges <- ecount(true_graph)
  Gamma <- graph_info$Gamma
  Theta_true <- graph_info$Theta

  log(sprintf("  True edges: %d\n", true_edges))

  nlambda <- length(lambda_range)

  one_rep <- function(rep_id) {
    data <- generate_data(s$distribution, s$n, s$d, Gamma)

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

    res_g
  }

  worker_setup <- function() {
    source("R/load_all.R")
  }

  # Run in parallel
  source("R/parallel.R")
  t_total <- proc.time()
  results <- run_parallel(
    sim_fn = one_rep,
    nsim = nsim,
    ncores = ncores,
    seed = seed_base + si,
    setup_fn = worker_setup,
    export_vars = c("s", "Gamma", "Theta_true", "true_graph",
                    "lambda_range", "nlambda"),
    export_env = environment()
  )
  elapsed <- (proc.time() - t_total)[3]

  # Aggregate
  f1_glasso   <- do.call(rbind, lapply(results, `[[`, "f1_g"))
  edges_glasso  <- do.call(rbind, lapply(results, `[[`, "edges_g"))
  time_glasso  <- sapply(results, `[[`, "time_glasso")

  n_fail_g <- sum(is.na(time_glasso))

  best_g <- which.max(colMeans(f1_glasso, na.rm = TRUE))
  log(sprintf("  EGLasso: mean F1=%.3f at log10(lam)=%.1f | mean time=%.2fs/sample",
              mean(f1_glasso[, best_g], na.rm = TRUE), log10(lambda_range[best_g]),
              mean(time_glasso, na.rm = TRUE)))
  if (n_fail_g > 0) log(sprintf(" | %d failed", n_fail_g))
  log("\n")
  log(sprintf("  Total elapsed: %.0fs\n\n", elapsed))

  # Save (same structure as full sim, but without EGLearn fields)
  sim_result <- list(
    config = list(
      m = s$m, d = s$d, k = s$k, n = s$n,
      q_threshold = s$q_threshold,
      distribution = s$distribution,
      nsim = nsim, ncores = ncores,
      seed_base = seed_base, setting_index = si
    ),
    lambda_range = lambda_range,
    f1_glasso = f1_glasso,
    edges_glasso = edges_glasso,
    time_glasso = time_glasso,
    true_edges = true_edges,
    Theta_true = Theta_true,
    elapsed = elapsed
  )
  saveRDS(sim_result, outfile)
  log(sprintf("  Saved: %s\n\n", outfile))
}

log("All settings complete.\n")
