# BIC-based lambda selection (optional utility)
# Requires: igraph, graphicalExtremes, est_W, glasso_c_reest

#' Select lambda via BIC/MBIC
#'
#' @param data Data matrix (n x d)
#' @param lambda.range Vector of lambda values (0 = auto-generate)
#' @param q.threshold Quantile threshold
#' @param iter.max Maximum GLASSO iterations
#' @param BIC_method One of "mean", "max", "MBIC", "maxMBIC", "medianMBIC"
#' @return List with lambda (selected), Theta (estimated precision matrix)
glasso_bic <- function(data, lambda.range = 0, q.threshold = 0.9,
                       iter.max = 1000, BIC_method = "maxMBIC") {

  if (length(lambda.range) == 1 && lambda.range == 0) {
    lambda.range <- 10 ^ seq(from = -2, to = 0, by = 0.1)
  }

  n <- nrow(data)
  d <- ncol(data)
  k <- n * (1 - q.threshold)

  nlambda <- length(lambda.range)
  IC <- rep(Inf, nlambda)

  W <- est_W(data, q.threshold)
  ev <- eigen(W$cov, only.values = TRUE)$values
  c_val <- 1 / (d * min(ev[ev > 1e-10]))  # smallest positive eigenvalue

  for (i in 1:nlambda) {
    lambda <- lambda.range[i]
    res <- glasso_c_reest(W$cov, lambda, c_val, iter.max)
    graph <- igraph::graph_from_adjacency_matrix(
      !(res$graph), diag = FALSE, mode = "undirected")

    if (!igraph::is.connected(graph)) {
      break
    }

    det_Theta_part <- numeric(d)
    for (j in 1:d) det_Theta_part[j] <- det(res$Theta_hat[-j, -j])

    if (sum(det_Theta_part <= 0) > 0) {
      Gamma_hat <- graphicalExtremes::complete_Gamma(
        Gamma = graphicalExtremes::Sigma2Gamma(W$cov), graph = graph)
      Theta <- graphicalExtremes::Gamma2Theta(Gamma_hat)
    } else {
      Theta <- res$Theta_hat
    }

    BIC <- rep(NA, d)
    for (j in 1:d) {
      Theta_part <- Theta[-j, -j]
      missing <- sum(res$graph[-j, -j])
      p <- ((d - 1) * d / 2 - missing / 2)

      if (BIC_method %in% c("mean", "max")) {
        BIC[j] <- -log(det(Theta_part)) +
          sum(diag(W$subcovlist[[j]] %*% Theta_part)) + p * log(k) / k
      } else if (BIC_method %in% c("MBIC", "maxMBIC", "medianMBIC")) {
        BIC[j] <- -log(det(Theta_part)) +
          sum(diag(W$subcovlist[[j]] %*% Theta_part)) +
          p * log(k) * log(log(d - 1)) / k
      }
    }

    if (BIC_method %in% c("max", "maxMBIC")) {
      IC[i] <- max(BIC)
    } else if (BIC_method %in% c("mean", "MBIC")) {
      IC[i] <- mean(BIC)
    } else if (BIC_method == "medianMBIC") {
      IC[i] <- quantile(BIC, 0.5)
    }
  }

  if (min(IC) < Inf) {
    lambda.select <- lambda.range[which.min(IC)]
  } else {
    lambda.select <- min(lambda.range)
  }

  if (which.min(IC) == 1) {
    cat("\n Warning: best lambda is at the lower bound!\n")
  } else if (which.min(IC) == nlambda) {
    cat("\n Warning: best lambda is at the upper bound!\n")
  }

  res <- glasso_c_reest(W$cov, lambda.select, c_val, thres, iter.max)
  Theta <- res$Theta_hat
  Theta[res$graph] <- 0
  list(lambda = lambda.select, Theta = Theta)
}
