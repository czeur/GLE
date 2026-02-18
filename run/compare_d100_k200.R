# Compare glasso_c vs EGLearn
# d=100, k=200, n=floor(200^(1/0.7))=1937, maxstable
# Both BA(1) and BA(2) graphs

source("R/load_all.R")

# --- Configuration ---
d <- 100
k <- 500
n <- floor(k^(1/0.7))  # 7172
q_threshold <- 1 - k / n
distribution <- "maxstable"
seed <- 42

lambda_range <- 10^seq(-1.5, 0.5, by = 0.1)  # wider range for d=100
rho_range <- seq(0.05, 0.5, length.out = 10)

cat(sprintf("d=%d, k=%d, n=%d, q=%.4f, dist=%s\n\n",
            d, k, n, q_threshold, distribution))

# --- F1 from igraph objects ---
F1_graph <- function(g_true, g_est) {
  n_inter <- ecount(igraph::intersection(g_true, g_est))
  n_fp <- ecount(g_est) - n_inter
  n_fn <- ecount(g_true) - n_inter
  prec <- if (n_inter + n_fp > 0) n_inter / (n_inter + n_fp) else 0
  rec  <- if (n_inter + n_fn > 0) n_inter / (n_inter + n_fn) else 0
  if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
}

theta_to_graph <- function(Theta_hat, d) {
  adj <- Theta_hat != 0
  diag(adj) <- FALSE
  graph_from_adjacency_matrix(adj * 1, mode = "undirected", diag = FALSE)
}

run_comparison <- function(m_val) {
  cat(sprintf("========== BA(%d) ==========\n", m_val))

  set.seed(seed)
  graph_info <- generate_graph("BA", d, m = m_val)
  data <- generate_data(distribution, n, d, graph_info$Gamma)
  true_graph <- graph_info$graph
  true_edges <- ecount(true_graph)
  cat(sprintf("True edges: %d\n", true_edges))

  # glasso_c
  cat(sprintf("Running glasso_c (%d lambdas)...\n", length(lambda_range)))
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

  best_g <- which.max(f1_glasso)
  cat(sprintf("  Time: %.2f s | Best F1: %.3f at lambda=%.4f (edges=%d)\n",
              t_glasso, f1_glasso[best_g], lambda_range[best_g], edges_glasso[best_g]))

  # EGLearn
  cat(sprintf("Running EGLearn (%d rhos)...\n", length(rho_range)))
  t_eglearn <- proc.time()
  fit_eg <- eglearn(data, p = q_threshold, rholist = rho_range, reg_method = "ns")
  t_eglearn <- (proc.time() - t_eglearn)[3]

  f1_eglearn <- sapply(fit_eg$graph, function(g) F1_graph(true_graph, g))
  edges_eglearn <- sapply(fit_eg$graph, ecount)

  best_e <- which.max(f1_eglearn)
  cat(sprintf("  Time: %.2f s | Best F1: %.3f at rho=%.4f (edges=%d)\n\n",
              t_eglearn, f1_eglearn[best_e], rho_range[best_e], edges_eglearn[best_e]))

  list(
    m = m_val, true_edges = true_edges,
    glasso_c = list(f1 = f1_glasso, edges = edges_glasso, time = t_glasso, best = best_g),
    eglearn = list(f1 = f1_eglearn, edges = edges_eglearn, time = t_eglearn, best = best_e)
  )
}

# Run both
res_ba1 <- run_comparison(1)
res_ba2 <- run_comparison(2)

# --- Plot ---
dir.create("output", showWarnings = FALSE)
pdf("output/compare_d100_k500_ba12.pdf", width = 12, height = 10)
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1))

for (res in list(res_ba1, res_ba2)) {
  m_val <- res$m
  bg <- res$glasso_c$best
  be <- res$eglearn$best

  # Panel 1: F1 vs lambda
  plot(log10(lambda_range), res$glasso_c$f1, type = "b", lwd = 2, pch = 19,
       xlab = expression(log[10](lambda)), ylab = "F1 score",
       main = sprintf("glasso_c — BA(%d)", m_val), ylim = c(0, 1), col = "blue")
  abline(v = log10(lambda_range[bg]), col = "red", lty = 2)
  text(log10(lambda_range[bg]), res$glasso_c$f1[bg] + 0.05,
       sprintf("F1=%.3f", res$glasso_c$f1[bg]), col = "red", cex = 0.8)
  mtext(sprintf("Time: %.1f s", res$glasso_c$time), side = 3, line = 0, cex = 0.7)

  # Panel 2: F1 vs rho
  plot(rho_range, res$eglearn$f1, type = "b", lwd = 2, pch = 19,
       xlab = expression(rho), ylab = "F1 score",
       main = sprintf("EGLearn — BA(%d)", m_val), ylim = c(0, 1), col = "darkgreen")
  abline(v = rho_range[be], col = "red", lty = 2)
  text(rho_range[be], res$eglearn$f1[be] + 0.05,
       sprintf("F1=%.3f", res$eglearn$f1[be]), col = "red", cex = 0.8)
  mtext(sprintf("Time: %.1f s", res$eglearn$time), side = 3, line = 0, cex = 0.7)

  # Panel 3: Edge count
  plot(log10(lambda_range), res$glasso_c$edges, type = "b", lwd = 2, pch = 19,
       xlab = "Tuning parameter index", ylab = "Number of edges",
       main = sprintf("Edges — BA(%d)", m_val),
       ylim = range(c(res$glasso_c$edges, res$eglearn$edges, res$true_edges)),
       col = "blue", xaxt = "n")
  axis(1, at = log10(lambda_range), labels = seq_along(lambda_range))
  points(log10(lambda_range[1]) + (seq_along(rho_range) - 1) *
           diff(range(log10(lambda_range))) / (length(rho_range) - 1),
         res$eglearn$edges, type = "b", lwd = 2, pch = 17, col = "darkgreen")
  abline(h = res$true_edges, col = "red", lty = 2)
  legend("topright", legend = c("glasso_c", "EGLearn", paste("True:", res$true_edges)),
         col = c("blue", "darkgreen", "red"), pch = c(19, 17, NA),
         lty = c(1, 1, 2), lwd = 2, cex = 0.8)
}

dev.off()
cat("Plot saved to output/compare_d100_k500_ba12.pdf\n")

# Save results
comparison <- list(
  config = list(d = d, k = k, n = n, q_threshold = q_threshold,
                distribution = distribution, seed = seed,
                lambda_range = lambda_range, rho_range = rho_range),
  BA1 = res_ba1, BA2 = res_ba2
)
saveRDS(comparison, "output/compare_d100_k500_ba12.rds")
cat("Results saved to output/compare_d100_k500_ba12.rds\n")
