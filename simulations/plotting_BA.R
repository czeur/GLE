library(graphicalExtremes)
library(igraph)

# Path to store figures
figpath <- "figures/"

# Path storing results
respath <- "results/"


# Set dimensions, 100 or 20
# For each d value, there is corresponding sample size
# nsim: number of replication

# d <- 50
d <- 20

if (d == 20) {
  n <- 5000
  q.threshold <- 0.9
#   m <- 100
} else if (d == 100){
  n <- 10000
  q.threshold <- 0.95
} else{
  stop("\n One needs to specify the sample size and threshold! \n")
}

nsim <- 100
q <- 1

# Load results
resfn <- paste("results/Results_BA_", q, "_d_", d, "_n_", n, "_sim_", nsim, ".Rdata", sep = "")
load(resfn)

# Plot the true graph
adj.mat <- matrix(0, d, d)
adj.mat[Theta != 0] <- 1

graph.true <- graph_from_adjacency_matrix(adj.mat, diag = FALSE,mode = "undirected")

# true_file <-paste(figpath, "true_graph_d_", d, ".eps", sep = "")
# postscript(true_file, width = 6, height = 8)
# plot(graph.true, edge.width = 3, vertex.size = 2, vertex.label = NA)
# title(main = "True", cex.main = 2)
# dev.off()

res.mat <- matrix(0, d, d)
for (i in 1:nsim){
  res.mat <- res.mat + (results[[i]]$Theta != 0)
}


fig_file <-paste(figpath, "Graph_BA_q_", q, "_d_", d, "_n_", n, "_sim_", nsim, ".pdf", sep = "")
pdf(fig_file, width = 12, height = 8)

par(mfrow=c(1,2))

l1 <- layout.fruchterman.reingold(graph.true)
plot(graph.true, edge_colors <- 0, layout = l1, edge.width = 3,vertex.size = 2, vertex.label = NA)
title(main = "True", cex.main = 2)

graph_est <- graph_from_adjacency_matrix(res.mat, diag = FALSE, mode = "undirected", weighted = TRUE)

edge_weights <- E(graph_est)$weight

if(max(edge_weights) != nsim) stop("maximum edge weigths is not correct!")

# normalized_weights <- edge_weights / max(edge_weights)
# edge_colors <- gray(1 - normalized_weights)

# edges <- as_data_frame(graph_est, what = "edges")
# layout <- layout_nicely(graph_est)
edge_order <- order(E(graph_est)$weight, decreasing = FALSE)

normalized_weights <- E(graph_est)$weight / max(E(graph_est)$weight)

normalized_weights[normalized_weights < 0.5] <- 0

edge_colors <- gray(1 - normalized_weights)

edge_colors[normalized_weights == 0] <- rgb(0, 0, 0, alpha = 0.3)

# Plot the graph with edge width proportional to weight and higher weights on top
# plot(g_sorted, edge.color = gray(1 - normalized_weights), edge.width = 3*normalized_weights, vertex.size = 2, vertex.label = NA)

plot(graph_est, layout = l1, edge.color = edge_colors, edge.width = 3*normalized_weights, edge.order = edge_order, vertex.size = 2, vertex.label = NA)
title(main = "Estimated (aggregated)", cex.main = 2)
dev.off()
