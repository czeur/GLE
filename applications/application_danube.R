library(graphicalExtremes)
library(igraph)
# library(parallel)

source("functions/BIC.R")

source("applications/functions_paper.R")

load("applications/coords_danube.Rdata")
# Path to store figures
figpath <- "figures/"

# ids <- (1:31)[-(23:27)]
ids <- (1:31)
X <- danube$data_clustered[,ids]
d <- ncol(X)
# flow_edges_26 <- danube$flow_edges[1:25,]
# flow_edges_26 <- apply(flow_edges_26, 1:2, function(x){
#   if(x <= 22) return(x)
#   else return(x-5)
# })
flow_edges_26 <- danube$flow_edges
g <- graph_from_edgelist(flow_edges_26)
# g <- graph_from_edgelist(danube$flow_edges)
g <- set_graph_parameters(g)
coords_danube <- coords_danube[ids,]

n <- nrow(X)
p <- 1 - floor(n^.7)/n
Y <- data2mpareto(X, p)
rholist <- seq(.005, .5, length.out=100)

eglearn_fit <- eglearn2(Y, rholist=rholist, reg_method="ns", ic="hr")
connected_tmp <- sapply(eglearn_fit$graph, is_connected)
eglearn_fit$graph <- lapply(eglearn_fit$graph, set_graph_parameters)
eglearn_fit$graph_ic <- lapply(eglearn_fit$graph_ic, set_graph_parameters)


eglearn_fit_mbic <- eglearn_fit$graph_ic$mbic
eglearn_fit_select1 <- eglearn_fit$graph[[15]]
eglearn_fit_select2 <- eglearn_fit$graph[[40]]
eglearn_fit_select3 <- eglearn_fit$graph[[81]]

glasso_fit_mbic <- glasso_bic(X, lambda.range = 10^seq(from = -2, to = 0, by = 0.1), q.threshold = p, thres = 0.01, iter.max = 1000, BIC_method = "maxMBIC")


glasso_graph_mbic <- graph_from_adjacency_matrix((glasso_fit_mbic$Theta != 0), diag = FALSE, mode = "undirected")


# vertex_attr(glasso_graph_mbic) <- list(name = colnames(X))

# graph_est <- graph_from_adjacency_matrix((glasso_fit$Theta != 0), diag = FALSE, mode = "undirected")



fig_file <-paste(figpath, "application_danube.pdf", sep = "")
pdf(fig_file, width = 12, height = 8)

par(mfrow=c(1,2))

plot(eglearn_fit_mbic, layout = coords_danube)
title(main = "EGLearn", cex.main = 2)
plot(glasso_graph_mbic, layout = coords_danube)
title(main = "EGLasso", cex.main = 2)
dev.off()



glasso_fit_select1 <- glasso_bic(X, lambda.range = 0.4, q.threshold = p, thres = 1e-3, iter.max = 1000, BIC_method = "maxMBIC")

glasso_fit_select2 <- glasso_bic(X, lambda.range = 0.44, q.threshold = p, thres = 1e-3, iter.max = 1000, BIC_method = "maxMBIC")

glasso_fit_select3 <- glasso_bic(X, lambda.range = 0.445, q.threshold = p, thres = 1e-3, iter.max = 1000, BIC_method = "maxMBIC")


glasso_graph_select1 <- graph_from_adjacency_matrix((glasso_fit_select1$Theta != 0), diag = FALSE, mode = "undirected")

glasso_graph_select2 <- graph_from_adjacency_matrix((glasso_fit_select2$Theta != 0), diag = FALSE, mode = "undirected")

glasso_graph_select3 <- graph_from_adjacency_matrix((glasso_fit_select3$Theta != 0), diag = FALSE, mode = "undirected")

# vertex_attr(glasso_graph_select) <- list(name = colnames(X))

fig_file <-paste(figpath, "application_danube_select.pdf", sep = "")
pdf(fig_file, width = 12, height = 16)

par(mfrow = c(3, 2), mar = c(2, 4, 2, 2) )

plot(eglearn_fit_select1, layout = coords_danube)
title(main = "EGLearn (rho = 0.075)", cex.main = 2)
plot(glasso_graph_select1, layout = coords_danube)
title(main = "EGLasso (gamma = 0.4)", cex.main = 2)

plot(eglearn_fit_select2, layout = coords_danube)
title(main = "EGLearn (rho = 0.2)", cex.main = 2)
plot(glasso_graph_select2, layout = coords_danube)
title(main = "EGLasso (gamma = 0.44)", cex.main = 2)


plot(eglearn_fit_select3, layout = coords_danube)
title(main = "EGLearn (rho = 0.405)", cex.main = 2)
plot(glasso_graph_select3, layout = coords_danube)
title(main = "EGLasso (gamma = 0.445)" , cex.main = 2)

dev.off()
