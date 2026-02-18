# Data generating processes: graph generators + data generators

library(graphicalExtremes)
library(igraph)

#' Generate a graph and corresponding Gamma/Theta matrices
#'
#' @param type Graph type: "BA" (Barabasi-Albert) or "tree" (random tree)
#' @param d Dimension (number of nodes)
#' @param m Number of edges per new node (BA model, default 1 = tree)
#' @param weight_range Range for edge weights (default c(2, 5))
#' @return List with Theta, Gamma, graph (igraph object)
generate_graph <- function(type = c("BA", "tree"), d, m = 1,
                           weight_range = c(2, 5)) {
  type <- match.arg(type)

  if (type == "BA") {
    graph <- sample_pa(n = d, m = m, zero.appeal = 1, directed = FALSE)
    W_mat <- as_adj(graph, sparse = FALSE) *
      matrix(runif(d^2, weight_range[1], weight_range[2]), nrow = d)
    W_mat[lower.tri(W_mat)] <- t(W_mat)[lower.tri(W_mat)]
    Theta <- diag(rowSums(W_mat)) - W_mat
    Gamma <- Theta2Gamma(Theta)

  } else if (type == "tree") {
    # Random tree: each new node attaches to a random existing node
    adj_mat <- matrix(0, d, d)
    for (k in 2:d) {
      node_to <- sample.int(k - 1, 1)
      adj_mat[k, node_to] <- 1
      adj_mat[node_to, k] <- 1
    }
    graph <- graph_from_adjacency_matrix(adj_mat, diag = FALSE, mode = "undirected")

    # Generate edge weights and complete the variogram
    G_vec <- runif(ecount(graph), min = weight_range[1] / 4, max = weight_range[2] / 4)
    G <- matrix(0, d, d)
    edge_index <- as.matrix(ends(graph, E(graph)))
    for (i in 1:nrow(edge_index)) {
      G[edge_index[i, 1], edge_index[i, 2]] <- G_vec[i]
      G[edge_index[i, 2], edge_index[i, 1]] <- G_vec[i]
    }
    Gamma <- complete_Gamma(graph = graph, Gamma = G)

    # Compute Theta from Gamma
    Theta <- matrix(0, d, d)
    for (k in 1:d) {
      Theta[-k, -k] <- solve(Gamma2Sigma(Gamma, k = k))
    }
  }

  list(Theta = Theta, Gamma = Gamma, graph = graph)
}

#' Generate data from a multivariate extreme value model
#'
#' @param distribution "mpareto" or "maxstable"
#' @param n Sample size
#' @param d Dimension
#' @param Gamma Variogram matrix
#' @return Data matrix (n x d)
generate_data <- function(distribution = c("mpareto", "maxstable"),
                          n, d, Gamma) {
  distribution <- match.arg(distribution)

  if (distribution == "mpareto") {
    rmpareto(n = n, model = "HR", d = d, par = Gamma)
  } else if (distribution == "maxstable") {
    rmstable(n = n, model = "HR", d = d, par = Gamma)
  }
}
