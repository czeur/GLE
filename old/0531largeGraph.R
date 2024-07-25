library(graphicalExtremes)
library(igraph)
library(glasso)

set.seed(34)
d <- 200
n <- 100000
q.threshold <- 0.95

# d <- 20
# n <- 5000
# q.threshold <- 0.9

# Simulate a simple tree
# Start from nodes 2->1, every time a new node randomly attaches to an old node

adj.mat <- matrix(0,d,d)
for (k in 2:d){
  node.to <- floor(runif(1,1,k-0.5))
  adj.mat[k,node.to] <- 1
  adj.mat[node.to,k] <- 1
}
graph.true <- graph_from_adjacency_matrix(adj.mat,diag = FALSE,mode = "undirected")
plot(graph.true,main='True',edge.width=2,vertex.size=1,vertex.label=NA)

# Now generate variogram Gamma

G.vec <- runif(ecount(graph.true), min = 0.5, max = 1)
G <- matrix(0,d,d)
edge.index <- matrix(as.numeric(igraph::ends(graph.true,igraph::E(graph.true))),ncol=2)
for (i in 1:nrow(edge.index)){
  G[edge.index[i,1],edge.index[i,2]] <- G.vec[i]
  G[edge.index[i,2],edge.index[i,1]] <- G.vec[i]
}
Gamma <- complete_Gamma(graph = graph.true, Gamma = G) 

# Calculate the true precision matrix Theta
Theta <- matrix(0,d,d)
for (k in 1:d){
  Theta[-k,-k] <- solve(Gamma2Sigma(Gamma,k=k)) # function Gamma2Sigma is used to calculate each Theta^{(k)}
}

##################################

# Simulate data
start_time <- Sys.time()

data <- rmpareto(n=n, model="HR", d=d, par=Gamma) 

# Rank transformation
for (i in 1:d){
  x <- data[,i]
  rx <- rank(x)/(n+1)
  # data[,i] <- -1/(log(rx))
  data[,i] <- 1/(1-rx)
}

finish_sim <- Sys.time()

# Compute S_0, the estimate for Sigma_0 = Sigma + M*1^T%*%1, manuscript proposition 3.5
cov.mat <- matrix(0,d,d) # this is S_0
cov.k.list <- list()
nk.list <- rep(NA,d)
# Set a threshold to truncate each dimension
for (k in 1:d){
  data.truncated <- data[data[,k]>quantile(data[,k],q.threshold),]
  # data.truncated <- data[data[,k]>1,]
  w <- log(data.truncated[,-k]) - matrix(log(data.truncated[,k]),nrow(data.truncated),d-1)
  cov.k <- cov(w)
  cov.mat[-k,-k] <- cov.mat[-k,-k] + cov.k/d
  cov.k.list[[k]] <- cov.k
  nk.list[k] <- nrow(data.truncated)
}
finish_est <- Sys.time()


###########################

# GLASSO

rho <- 1.91
fit <- glasso(cov.mat + matrix(0,d,d), rho, penalize.diagonal=F)
Theta_hat <- fit$wi
Theta_graph <- solve(cov.mat + matrix(0,d,d))
Theta_graph <- Theta_graph - mean(Theta_graph)
Theta_graph[Theta_hat==0] <- 0
Theta_hat <- Theta_graph

###########################
finish_graph <- Sys.time()

# Plot results
par(mfrow=c(1,2))
l1 <- layout.fruchterman.reingold(graph.true)
plot(graph.true,layout=l1,main='True',edge.width=2,vertex.size=1,vertex.label=NA)

adj.est <- abs(Theta_hat)
diag(adj.est) <- rep(0,d)
adj.est.index <- adj.est
adj.est.index[adj.est>1e-4] <- 1
# adj.est <- adj.est/max(adj.est[adj.est>0])
graph_est <- graph_from_adjacency_matrix(adj.est.index,diag = FALSE,mode = "undirected")
plot(graph_est,
     layout=l1,main='Estimated',edge.width=1.5,vertex.size=1,vertex.label=NA)

cat("\n Time used for simulation: ", finish_sim - start_time)
cat("\n Time used for estimation: ", finish_est - finish_sim)
cat("\n Time used for graph: ", finish_graph - finish_est)
cat("\n")
