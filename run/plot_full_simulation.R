# Plot results from run_full_simulation.R
# Run locally after syncing output/full_sim/*.rds
#
# Produces:
#   8 F1 figures: one per (graph, d, model) combination
#     x-axis = k/d ratio, boxplots of best F1 per sample for glasso_c vs EGLearn
#   8 time figures: same layout, boxplots of estimation time per sample

# === Load all results ===
rds_dir <- "output/full_sim"
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(rds_files) == 0) stop("No RDS files found in ", rds_dir)

all_results <- lapply(rds_files, readRDS)
names(all_results) <- gsub("\\.rds$", "", basename(rds_files))

cat(sprintf("Loaded %d settings\n\n", length(all_results)))

# === Organize by (graph, d, model) ===
# Each group has 3 k/d ratios

groups <- list()
for (nm in names(all_results)) {
  res <- all_results[[nm]]
  cfg <- res$config
  key <- sprintf("BA(%d)_d%d_%s", cfg$m, cfg$d, cfg$distribution)
  if (is.null(groups[[key]])) groups[[key]] <- list()
  groups[[key]][[nm]] <- res
}

cat(sprintf("Groups: %d\n", length(groups)))
for (key in sort(names(groups))) {
  cat(sprintf("  %s: %d k-values\n", key, length(groups[[key]])))
}

# === Colors ===
col_glasso <- "steelblue"
col_eglearn <- "seagreen"

# === Plot ===
dir.create("output/full_sim/plots", recursive = TRUE, showWarnings = FALSE)

# === Helper: plot one F1 boxplot figure ===
plot_f1 <- function(grp, key) {
  k_over_d <- sapply(grp, function(r) r$config$k / r$config$d)
  grp <- grp[order(k_over_d)]
  k_over_d <- sort(k_over_d)
  n_ratios <- length(k_over_d)

  cfg1 <- grp[[1]]$config
  title_str <- sprintf("BA(%d), d=%d, %s", cfg1$m, cfg1$d, cfg1$distribution)

  best_f1_glasso <- lapply(grp, function(r) apply(r$f1_glasso, 1, max))
  best_f1_eglearn <- lapply(grp, function(r) apply(r$f1_eglearn, 1, max))

  bp_data <- list()
  bp_cols <- c()
  bp_at <- c()
  for (i in seq_len(n_ratios)) {
    bp_data[[2 * i - 1]] <- best_f1_glasso[[i]]
    bp_data[[2 * i]]     <- best_f1_eglearn[[i]]
    bp_cols <- c(bp_cols, col_glasso, col_eglearn)
    bp_at <- c(bp_at, i * 3 - 1.5, i * 3 - 0.5)
  }

  pdf(sprintf("output/full_sim/plots/f1_%s.pdf", key), width = 7, height = 5)
  par(mar = c(5.5, 4.5, 3, 1))
  boxplot(bp_data, at = bp_at, col = bp_cols, names = rep("", length(bp_data)),
          ylim = c(0, 1), ylab = "Best F1 score (per sample)",
          main = title_str, outline = FALSE, boxwex = 0.8)

  label_pos <- sapply(seq_len(n_ratios), function(i) mean(bp_at[c(2*i-1, 2*i)]))
  for (i in seq_len(n_ratios)) {
    mtext(sprintf("k/d=%.1f\nk=%d, n=%d", k_over_d[i],
                  grp[[i]]$config$k, grp[[i]]$config$n),
          side = 1, at = label_pos[i], line = 2, cex = 0.75)
  }

  legend("bottomleft", legend = c("glasso_c", "EGLearn"),
         fill = c(col_glasso, col_eglearn), cex = 0.9, bg = "white")
  dev.off()
  cat(sprintf("  Saved: output/full_sim/plots/f1_%s.pdf\n", key))
}

# === Helper: plot one time boxplot figure ===
plot_time <- function(grp, key) {
  k_over_d <- sapply(grp, function(r) r$config$k / r$config$d)
  grp <- grp[order(k_over_d)]
  k_over_d <- sort(k_over_d)
  n_ratios <- length(k_over_d)

  cfg1 <- grp[[1]]$config
  title_str <- sprintf("BA(%d), d=%d, %s", cfg1$m, cfg1$d, cfg1$distribution)

  bp_data <- list()
  bp_cols <- c()
  bp_at <- c()
  for (i in seq_len(n_ratios)) {
    bp_data[[2 * i - 1]] <- grp[[i]]$time_glasso
    bp_data[[2 * i]]     <- grp[[i]]$time_eglearn
    bp_cols <- c(bp_cols, col_glasso, col_eglearn)
    bp_at <- c(bp_at, i * 3 - 1.5, i * 3 - 0.5)
  }

  pdf(sprintf("output/full_sim/plots/time_%s.pdf", key), width = 7, height = 5)
  par(mar = c(5.5, 4.5, 3, 1))
  boxplot(bp_data, at = bp_at, col = bp_cols, names = rep("", length(bp_data)),
          ylab = "Estimation time per sample (s)",
          main = title_str, outline = FALSE, boxwex = 0.8)

  label_pos <- sapply(seq_len(n_ratios), function(i) mean(bp_at[c(2*i-1, 2*i)]))
  for (i in seq_len(n_ratios)) {
    mtext(sprintf("k/d=%.1f\nk=%d, n=%d", k_over_d[i],
                  grp[[i]]$config$k, grp[[i]]$config$n),
          side = 1, at = label_pos[i], line = 2, cex = 0.75)
  }

  legend("topleft", legend = c("glasso_c", "EGLearn"),
         fill = c(col_glasso, col_eglearn), cex = 0.9, bg = "white")
  dev.off()
  cat(sprintf("  Saved: output/full_sim/plots/time_%s.pdf\n", key))
}

# === Generate all figures ===
cat("\nGenerating F1 figures...\n")
for (key in sort(names(groups))) plot_f1(groups[[key]], key)

cat("\nGenerating time figures...\n")
for (key in sort(names(groups))) plot_time(groups[[key]], key)

# === Summary table ===
summary_rows <- lapply(all_results, function(res) {
  cfg <- res$config
  best_f1_g <- apply(res$f1_glasso, 1, max)
  best_f1_e <- apply(res$f1_eglearn, 1, max)

  data.frame(
    graph = sprintf("BA(%d)", cfg$m),
    model = cfg$distribution,
    d = cfg$d,
    k = cfg$k,
    k_over_d = cfg$k / cfg$d,
    n = cfg$n,
    true_edges = res$true_edges,
    f1_glasso_mean = mean(best_f1_g),
    f1_glasso_sd = sd(best_f1_g),
    f1_eglearn_mean = mean(best_f1_e),
    f1_eglearn_sd = sd(best_f1_e),
    time_glasso_mean = mean(res$time_glasso),
    time_eglearn_mean = mean(res$time_eglearn),
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[order(summary_df$graph, summary_df$d, summary_df$k, summary_df$model), ]
rownames(summary_df) <- NULL

cat("\n=== Summary ===\n")
print(summary_df, digits = 3)

write.csv(summary_df, "output/full_sim/plots/summary.csv", row.names = FALSE)
cat("\nSaved: output/full_sim/plots/summary.csv\n")
