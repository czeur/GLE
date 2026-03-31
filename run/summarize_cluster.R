# Summarize cluster results from output/full_sim_3103
rds_dir <- "output/full_sim_3103"
files <- list.files(rds_dir, pattern = "[.]rds$", full.names = TRUE)

rows <- list()
for (f in files) {
  res <- readRDS(f)
  cfg <- res$config
  best_f1_g <- apply(res$f1_glasso, 1, max, na.rm = TRUE)
  best_f1_e <- apply(res$f1_eglearn, 1, max, na.rm = TRUE)
  rows[[length(rows) + 1]] <- data.frame(
    graph = sprintf("BA(%d)", cfg$m),
    dist = cfg$distribution,
    d = cfg$d,
    k = cfg$k,
    kd = round(cfg$k / cfg$d, 1),
    f1_eglasso = mean(best_f1_g, na.rm = TRUE),
    f1_eglearn = mean(best_f1_e, na.rm = TRUE),
    time_eglasso = mean(res$time_glasso, na.rm = TRUE),
    time_eglearn = mean(res$time_eglearn, na.rm = TRUE)
  )
}
df <- do.call(rbind, rows)
df <- df[order(df$graph, df$d, df$k, df$dist), ]
df$speed_ratio <- df$time_eglearn / df$time_eglasso
df$f1_gap <- df$f1_eglearn - df$f1_eglasso

cat("=== Cluster results (100 reps each) ===\n\n")
print(df, digits = 3, row.names = FALSE)

cat("\n\n=== d=100 settings only ===\n\n")
print(df[df$d == 100, c("graph","dist","k","kd","f1_eglasso","f1_eglearn","f1_gap",
                          "time_eglasso","time_eglearn","speed_ratio")],
      digits = 3, row.names = FALSE)
