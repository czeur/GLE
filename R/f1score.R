# F1 score for graph recovery evaluation

#' Compute F1 score between true and estimated precision matrices
#'
#' Compares the sparsity patterns (off-diagonal zero/nonzero) of two
#' precision matrices and returns the F1 score.
#'
#' @param Theta True precision matrix (d x d)
#' @param Theta_hat Estimated precision matrix (d x d)
#' @return F1 score (scalar in [0, 1])
F1score <- function(Theta, Theta_hat) {
  d <- nrow(Theta)
  if (d != nrow(Theta_hat)) stop("Dimension mismatch")

  true <- Theta[upper.tri(Theta)] != 0
  est <- Theta_hat[upper.tri(Theta_hat)] != 0

  tp <- sum(true & est)
  fp <- sum(!true & est)
  fn <- sum(true & !est)

  precision <- if (tp + fp > 0) tp / (tp + fp) else 0
  recall <- if (tp + fn > 0) tp / (tp + fn) else 0

  if (precision + recall > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    0
  }
}
