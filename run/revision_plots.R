# Generate publication-quality plots for GLE paper revision
# Input: [indir]/*.rds (24 settings, 100 reps each)
# Output: output/revision_plots/*.pdf
#
# Usage: Rscript run/revision_plots.R [indir]
#        Default indir: output/full_sim

# === Load all results ===
args <- commandArgs(trailingOnly = TRUE)
rds_dir <- if (length(args) >= 1) args[1] else "output/full_sim"
rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) stop("No RDS files found in ", rds_dir)

all_results <- lapply(rds_files, readRDS)
names(all_results) <- gsub("\\.rds$", "", basename(rds_files))
cat(sprintf("Loaded %d settings\n", length(all_results)))

# === Organize by (graph, d, distribution) ===
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
  kvals <- sapply(groups[[key]], function(r) r$config$k)
  cat(sprintf("  %s: k = %s\n", key, paste(sort(kvals), collapse = ", ")))
}

# === Output directory ===
outdir <- "output/revision_plots"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# === Colors ===
col_glasso  <- "steelblue"
col_eglearn <- "seagreen"

# === Helper: compute best F1 per sample with NA handling ===
best_f1 <- function(f1_matrix) {
  vals <- apply(f1_matrix, 1, function(row) {
    finite <- row[is.finite(row)]
    if (length(finite) == 0) return(NA_real_)
    max(finite)
  })
  vals
}

# === Helper: 4-panel F1 boxplot for one distribution ===
plot_f1_4panel <- function(groups_list, dist, filename) {
  # Panel order: BA(1) d=20, BA(1) d=100, BA(2) d=20, BA(2) d=100
  panel_keys <- sprintf("BA(%d)_d%d_%s", rep(1:2, each = 2), rep(c(20, 100), 2), dist)
  panel_labels <- c("BA(1), d = 20", "BA(1), d = 100",
                     "BA(2), d = 20", "BA(2), d = 100")

  pdf(file.path(outdir, filename), width = 7, height = 6.5)
  par(mfrow = c(2, 2), mar = c(4.5, 4, 2.5, 0.8), oma = c(0, 0, 0, 0))

  for (p in seq_along(panel_keys)) {
    key <- panel_keys[p]
    grp <- groups_list[[key]]
    if (is.null(grp)) { plot.new(); next }

    # Sort by k/d
    k_over_d <- sapply(grp, function(r) r$config$k / r$config$d)
    grp <- grp[order(k_over_d)]
    k_over_d <- sort(k_over_d)
    n_ratios <- length(k_over_d)

    bf1_glasso  <- lapply(grp, function(r) best_f1(r$f1_glasso))
    bf1_eglearn <- lapply(grp, function(r) best_f1(r$f1_eglearn))

    bp_data <- list()
    bp_cols <- c()
    bp_at   <- c()
    for (i in seq_len(n_ratios)) {
      bp_data[[2 * i - 1]] <- bf1_glasso[[i]]
      bp_data[[2 * i]]     <- bf1_eglearn[[i]]
      bp_cols <- c(bp_cols, col_glasso, col_eglearn)
      bp_at   <- c(bp_at, i * 3 - 1.5, i * 3 - 0.5)
    }

    boxplot(bp_data, at = bp_at, col = bp_cols,
            names = rep("", length(bp_data)),
            ylim = c(0, 1), ylab = "Best F1 score",
            main = panel_labels[p], outline = FALSE, boxwex = 0.8,
            cex.main = 1.1, cex.lab = 1, cex.axis = 0.9)

    # x-axis labels
    label_pos <- sapply(seq_len(n_ratios), function(i) mean(bp_at[c(2*i-1, 2*i)]))
    for (i in seq_len(n_ratios)) {
      kd_label <- k_over_d[i]
      # Format nicely: show as fraction if clean, else decimal
      kd_str <- if (kd_label == round(kd_label)) sprintf("%.0f", kd_label) else sprintf("%.1f", kd_label)
      mtext(bquote(italic(k/d) == .(kd_str)),
            side = 1, at = label_pos[i], line = 2, cex = 0.7)
    }

    # Legend in first panel only
    if (p == 1) {
      legend("bottomright", legend = c("EGLasso", "EGLearn"),
             fill = c(col_glasso, col_eglearn), cex = 0.8, bg = "white",
             border = NA, box.lwd = 0)
    }
  }

  dev.off()
  cat(sprintf("Saved: %s\n", file.path(outdir, filename)))
}

# === Helper: 4-panel timing boxplot for one distribution ===
plot_time_4panel <- function(groups_list, dist, filename) {
  panel_keys <- sprintf("BA(%d)_d%d_%s", rep(1:2, each = 2), rep(c(20, 100), 2), dist)
  panel_labels <- c("BA(1), d = 20", "BA(1), d = 100",
                     "BA(2), d = 20", "BA(2), d = 100")

  pdf(file.path(outdir, filename), width = 7, height = 6.5)
  par(mfrow = c(2, 2), mar = c(4.5, 4, 2.5, 0.8), oma = c(0, 0, 0, 0))

  for (p in seq_along(panel_keys)) {
    key <- panel_keys[p]
    grp <- groups_list[[key]]
    if (is.null(grp)) { plot.new(); next }

    # Sort by k/d
    k_over_d <- sapply(grp, function(r) r$config$k / r$config$d)
    grp <- grp[order(k_over_d)]
    k_over_d <- sort(k_over_d)
    n_ratios <- length(k_over_d)

    bp_data <- list()
    bp_cols <- c()
    bp_at   <- c()
    for (i in seq_len(n_ratios)) {
      tg <- grp[[i]]$time_glasso
      te <- grp[[i]]$time_eglearn
      # Replace non-finite with NA
      tg[!is.finite(tg)] <- NA
      te[!is.finite(te)] <- NA
      bp_data[[2 * i - 1]] <- tg
      bp_data[[2 * i]]     <- te
      bp_cols <- c(bp_cols, col_glasso, col_eglearn)
      bp_at   <- c(bp_at, i * 3 - 1.5, i * 3 - 0.5)
    }

    # Compute y-range across all data (log scale)
    all_times <- unlist(bp_data)
    all_times <- all_times[is.finite(all_times) & all_times > 0]
    if (length(all_times) == 0) { plot.new(); next }
    ylim <- range(all_times)

    boxplot(bp_data, at = bp_at, col = bp_cols,
            names = rep("", length(bp_data)),
            ylab = "Time per sample (seconds)", log = "y",
            main = panel_labels[p], outline = FALSE, boxwex = 0.8,
            cex.main = 1.1, cex.lab = 1, cex.axis = 0.9)

    # x-axis labels
    label_pos <- sapply(seq_len(n_ratios), function(i) mean(bp_at[c(2*i-1, 2*i)]))
    for (i in seq_len(n_ratios)) {
      kd_label <- k_over_d[i]
      kd_str <- if (kd_label == round(kd_label)) sprintf("%.0f", kd_label) else sprintf("%.1f", kd_label)
      mtext(bquote(italic(k/d) == .(kd_str)),
            side = 1, at = label_pos[i], line = 2, cex = 0.7)
    }

    # Legend in first panel only
    if (p == 1) {
      legend("topleft", legend = c("EGLasso", "EGLearn"),
             fill = c(col_glasso, col_eglearn), cex = 0.8, bg = "white",
             border = NA, box.lwd = 0)
    }
  }

  dev.off()
  cat(sprintf("Saved: %s\n", file.path(outdir, filename)))
}

# === Generate main text figures (maxstable) ===
cat("\n--- Main text figures (maxstable) ---\n")
plot_f1_4panel(groups, "maxstable", "f1_maxstable.pdf")
plot_time_4panel(groups, "maxstable", "time_maxstable.pdf")

# === Generate online appendix figures (mpareto) ===
cat("\n--- Online appendix figures (mpareto) ---\n")
plot_f1_4panel(groups, "mpareto", "f1_mpareto.pdf")
plot_time_4panel(groups, "mpareto", "time_mpareto.pdf")

# === Summary table ===
cat("\n--- Summary ---\n")
summary_rows <- lapply(all_results, function(res) {
  cfg <- res$config
  bf1_g <- best_f1(res$f1_glasso)
  bf1_e <- best_f1(res$f1_eglearn)
  data.frame(
    graph = sprintf("BA(%d)", cfg$m),
    dist  = cfg$distribution,
    d     = cfg$d,
    k     = cfg$k,
    kd    = round(cfg$k / cfg$d, 1),
    f1_eglasso  = round(mean(bf1_g, na.rm = TRUE), 3),
    f1_eglearn  = round(mean(bf1_e, na.rm = TRUE), 3),
    time_eglasso = round(mean(res$time_glasso, na.rm = TRUE), 2),
    time_eglearn = round(mean(res$time_eglearn, na.rm = TRUE), 2),
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[order(summary_df$graph, summary_df$dist, summary_df$d, summary_df$k), ]
rownames(summary_df) <- NULL
print(summary_df, row.names = FALSE)

write.csv(summary_df, file.path(outdir, "summary.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %s\n", file.path(outdir, "summary.csv")))

cat("\nAll revision plots generated.\n")
