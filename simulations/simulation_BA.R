library(graphicalExtremes)
library(igraph)
library(parallel)

source("functions/sim_graph.R")

# Path to store figures
figpath <- "figures/"

set.seed(100)

# Set dimensions, 20 or 50
# For each d value, there is corresponding
# n: sample size
# q.threshold: choice of threshold
# nsim: number of replication

# d <- 100
d <- 20
nsim <- 100

if (d == 20) {
  n <- 5000
  q.threshold <- 0.9
  lambda.range <- 10 ^ seq(from = -2, to = 0, by = 0.1)
#   m <- 100
} else if (d == 100){
  n <- 10000
  q.threshold <- 0.95
  lambda.range <- 10 ^ seq(from = -2.4, to = -1.4, by = 0.2)
} else{
  stop("\n One needs to specify the sample size and threshold! \n")
}

# Simulate a graph with the BA(d, q) model.
q <- 1
graph <- sample_pa(n = d, m = q, zero.appeal = 1, directed = FALSE)
graph.true <- graph
W_mat <- as_adj(graph, sparse = FALSE) * matrix(runif(d^2, 2, 5), nrow = d)
W_mat[lower.tri(W_mat)] <- t(W_mat)[lower.tri(W_mat)]
Theta <- diag(rowSums(W_mat)) - W_mat
Gamma <- Theta2Gamma(Theta)

cat("\n Finish setting up Theta! \n")
##################################

# Single function for simulation


all_objects <- ls()
print(all_objects)

num_cores <- 10  # Use one less core than available
cl <- makeCluster(num_cores)
clusterExport(cl, all_objects)
clusterEvalQ(cl, library(graphicalExtremes))
clusterEvalQ(cl, library(glmnet))
clusterEvalQ(cl, library(igraph))

cat("\n Ready to start simulation! \n")

start_time <- Sys.time()
# Apply the function in parallel

total_jobs <- nsim
batch_size <- 10  # Number of jobs to run at once

results <- list()

for (i in seq(1, total_jobs, by = batch_size)) {
  jobs_to_run <- min(batch_size, total_jobs - i + 1)
  random_seeds <- sample.int(1e6, jobs_to_run)
  clusterExport(cl, "random_seeds")
  results <- c(results, parLapply(cl, 1:jobs_to_run, function(j) {
      set.seed(random_seeds[j])
      sim_graph(Theta, Gamma, n, d, lambda.range = lambda.range, q.threshold, thres = 1e-3, iter.max = 1e3, BIC_method = "maxMBIC")
  }))
}

# Stop the cluster
stopCluster(cl)

cat("\n Time used for simulation: ", Sys.time() - start_time)
cat("\n")

fn =  paste("results/Results_BA_", q, "_d_", d, "_n_", n, "_sim_", nsim, ".Rdata", sep = "")
save(results, Theta, file = fn)
