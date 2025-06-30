library(graphicalExtremes)
library(igraph)
library(parallel)

source("functions/sim_graph.R")

set.seed(100)

# Set dimensions, 20 or 100
# For each d value, there is corresponding
# n: sample size
# q.threshold: choice of threshold
# q: q=1 or q=2 in the BA(d, q) model
# nsim: number of replication

d <- 20
qset <- c(1, 2)
nsim <- 100

n <- 500
q.threshold <- 0.85

lambda.range <- seq(from = -2, to = 0, by = 0.1)
lambda.range <- 10 ^ lambda.range

# Simulate a BA(d,q) model as the "true graph"

for (q in qset){
  graph <- sample_pa(n = d, m = q, zero.appeal = 1, directed = FALSE)
  W_mat <- as_adj(graph, sparse = FALSE) * matrix(runif(d^2, 2, 5), nrow = d)
  W_mat[lower.tri(W_mat)] <- t(W_mat)[lower.tri(W_mat)]
  Theta <- diag(rowSums(W_mat)) - W_mat
  Gamma <- Theta2Gamma(Theta)
  
  cat("\n Finish setting up Theta! \n")

  all_objects <- ls()
  print(all_objects)

  num_cores <- 10  # Use one less core than available
  cl <- makeCluster(num_cores)
  clusterExport(cl, all_objects)
  clusterEvalQ(cl, library(graphicalExtremes))
  clusterEvalQ(cl, library(glmnet))
  clusterEvalQ(cl, library(igraph))

  cat("\n Read to start simulation! \n")

  start_time <- Sys.time()
  # Apply the function in parallel
  results <- parLapply(cl, 1:nsim, function(i) {
    sim_graphs_multilambda(Theta, n, d, q.threshold, lambda.range, thres = 1e-3, iter.max = 1e3, BIC_method = "maxMBIC")
  })

  # Stop the cluster
  stopCluster(cl)

  cat("\n Time used for simulation: ", Sys.time() - start_time)
  cat("\n")

  fn =  paste("results/Multiple tuning_d_", d, "_q_", q, ".Rdata", sep = "")
  save(results, Theta, lambda.range, file = fn)
}
