# Parallel simulation utilities

library(parallel)

#' Run a function in parallel over nsim replications
#'
#' @param sim_fn Function to call for each replication. Must accept a single
#'   integer argument (the replication index) and return arbitrary results.
#' @param nsim Number of replications
#' @param ncores Number of parallel cores (default: detectCores() - 1)
#' @param batch_size Number of jobs per batch (default: ncores)
#' @param seed Base seed for reproducibility. Each replication gets a
#'   deterministic seed derived from this.
#' @param setup_fn Optional function called once per worker to set up
#'   the environment (e.g., load libraries, source files).
#'   Receives no arguments.
#' @param export_vars Character vector of variable names to export to workers
#' @param export_env Environment from which to export variables (default: parent)
#' @return List of length nsim with results from each replication
run_parallel <- function(sim_fn, nsim, ncores = NULL, batch_size = NULL,
                         seed = 42, setup_fn = NULL,
                         export_vars = NULL, export_env = parent.frame()) {
  if (is.null(ncores)) ncores <- max(1, detectCores() - 1)
  if (is.null(batch_size)) batch_size <- ncores

  # Generate deterministic seeds for each replication
  set.seed(seed)
  rep_seeds <- sample.int(1e7, nsim)

  cl <- makeCluster(ncores)
  on.exit(stopCluster(cl), add = TRUE)

  # Export sim_fn and user variables
  clusterExport(cl, "sim_fn", envir = environment())
  if (!is.null(export_vars)) {
    clusterExport(cl, export_vars, envir = export_env)
  }

  # Run setup function on each worker
  if (!is.null(setup_fn)) {
    clusterCall(cl, setup_fn)
  }

  # Run in batches
  results <- vector("list", nsim)
  for (start in seq(1, nsim, by = batch_size)) {
    end <- min(start + batch_size - 1, nsim)
    batch_idx <- start:end
    batch_seeds <- rep_seeds[batch_idx]

    clusterExport(cl, "batch_seeds", envir = environment())

    batch_results <- parLapply(cl, seq_along(batch_idx), function(j) {
      set.seed(batch_seeds[j])
      sim_fn(batch_idx[j])
    })

    results[batch_idx] <- batch_results

    cat(sprintf("  Completed %d / %d replications\n", end, nsim))
  }

  results
}
