# Run graphical LASSO across a range of lambda values
# No model selection -- returns all results for post-hoc evaluation

#' Estimate graphical model across a lambda range
#'
#' Given raw data, computes est_W once and runs glasso_c_reest at each lambda.
#' Returns all estimated precision matrices and sparsity patterns.
#'
#' @param data Data matrix (n x d)
#' @param lambda_range Vector of lambda values (should be sorted ascending)
#' @param q_threshold Quantile threshold for est_W
#' @param iter_max Maximum GLASSO iterations (default 1000)
#' @return List with:
#'   - W: the estimated W matrix (output of est_W)
#'   - c: the constant used
#'   - lambda_range: the lambdas used
#'   - results: list of length(lambda_range), each with Theta_hat and graph
glasso_path <- function(data, lambda_range, q_threshold,
                        iter_max = 1000) {
  d <- ncol(data)

  # Step 1: Compute W matrix once
  W <- est_W(data, q_threshold)

  # Step 2: Compute constant c
  c_val <- 1 / (d * eigen(W$cov)$values[1])

  # Step 3: Run glasso_c_reest at each lambda
  results <- vector("list", length(lambda_range))
  for (i in seq_along(lambda_range)) {
    results[[i]] <- glasso_c_reest(W$cov, lambda_range[i], c_val, iter_max)
  }

  list(
    W = W,
    c = c_val,
    lambda_range = lambda_range,
    results = results
  )
}
