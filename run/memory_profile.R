# Memory profiling: EGLasso vs EGLearn
# Measures peak memory usage for a single fit at each dimension
#
# Usage: Rscript --vanilla run/memory_profile.R

# Run from repo root: Rscript --vanilla run/memory_profile.R
source("R/load_all.R")
source("R/eglearn.R")

# Use pryr::mem_change or gc() to track memory
measure_memory <- function(expr) {
  gc(reset = TRUE)
  gc_before <- gc(verbose = FALSE)
  mem_before <- gc_before[2, 2]  # Vcells used (MB)

  result <- eval(expr)

  gc_after <- gc(verbose = FALSE)
  mem_after <- gc_after[2, 2]
  peak <- gc_after[2, 6]  # max used (MB)

  list(result = result, mem_used = mem_after - mem_before, peak = peak)
}

cat("=== Memory Profile: EGLasso vs EGLearn ===\n\n")

lambda_range <- 10^seq(-1.2, 0, by = 0.1)
rho_range <- seq(0.025, 0.525, by = 0.025)
seed <- 42

for (d in c(20, 50, 100, 200)) {
  m <- 2
  k <- round(d * 5)
  if (d >= 200) k <- round(d * 2.5)  # keep k/d reasonable for large d
  n <- floor(k^(1 / 0.7))
  q_threshold <- 1 - k / n

  cat(sprintf("--- d=%d, BA(2), k=%d, n=%d ---\n", d, k, n))

  set.seed(seed)
  gi <- generate_graph("BA", d, m = m)
  set.seed(seed + 1)
  data <- generate_data("mpareto", n, d, gi$Gamma)

  # EGLasso
  gc(reset = TRUE)
  gc(verbose = FALSE)
  t0 <- proc.time()
  mem_before_g <- sum(gc(verbose = FALSE)[, 2])
  path_g <- tryCatch(
    glasso_path(data, lambda_range, q_threshold, method = "cpp"),
    error = function(e) NULL)
  mem_after_g <- sum(gc(verbose = FALSE)[, 2])
  time_g <- (proc.time() - t0)[3]
  peak_g <- sum(gc(verbose = FALSE)[, 6])

  # EGLearn
  gc(reset = TRUE)
  gc(verbose = FALSE)
  t0 <- proc.time()
  mem_before_e <- sum(gc(verbose = FALSE)[, 2])
  fit_e <- tryCatch(
    eglearn2(data, p = q_threshold, rholist = rho_range, reg_method = "ns"),
    error = function(e) NULL)
  mem_after_e <- sum(gc(verbose = FALSE)[, 2])
  time_e <- (proc.time() - t0)[3]
  peak_e <- sum(gc(verbose = FALSE)[, 6])

  cat(sprintf("  EGLasso:  time=%6.2fs  mem_delta=%6.1fMB  peak=%6.1fMB\n",
              time_g, mem_after_g - mem_before_g, peak_g))
  cat(sprintf("  EGLearn:  time=%6.2fs  mem_delta=%6.1fMB  peak=%6.1fMB\n",
              time_e, mem_after_e - mem_before_e, peak_e))
  cat("\n")
}
