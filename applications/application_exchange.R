library(graphicalExtremes)
library(igraph)
# library(parallel)

source("functions/BIC.R")

source("applications/functions_paper.R")

load("applications/coords_exchange.Rdata")
load("applications/currency.Rdata")
# Path to store figures
figpath <- "figures/"

n <- nrow(X)
p <- 1 - floor(n^.7)/n
Y <- data2mpareto(X, p)
rholist <- seq(.005, 1, length.out=200)


eglearn_fit <- eglearn2(Y, rholist=rholist, reg_method="ns", ic="hr")
connected_tmp <- sapply(eglearn_fit$graph, is_connected)
eglearn_fit$graph <- lapply(eglearn_fit$graph, function(gr) {
  vertex_attr(gr) <- list(name = colnames(X))
  return(set_graph_parameters(gr))
})
eglearn_fit$graph_ic <- lapply(eglearn_fit$graph_ic, function(gr) {
  vertex_attr(gr) <- list(name = colnames(X))
  return(set_graph_parameters(gr))
})

eglearn_fit_mbic <- eglearn_fit$graph_ic$mbic
eglearn_fit_select <- eglearn_fit$graph[[106]]


glasso_fit_mbic <- glasso_bic(X, lambda.range = 10^seq(from = -2, to = 0, by = 0.1), q.threshold = p, thres = 1e-3, iter.max = 1000, BIC_method = "maxMBIC")


glasso_graph_mbic <- graph_from_adjacency_matrix((glasso_fit_mbic$Theta != 0), diag = FALSE, mode = "undirected")

vertex_attr(glasso_graph_mbic) <- list(name = colnames(X))

fig_file <-paste(figpath, "application_exchange.pdf", sep = "")
pdf(fig_file, width = 12, height = 8)

par(mfrow=c(1,2))

plot(eglearn_fit_mbic, layout = coords_exchange)
title(main = "EGLearn", cex.main = 2)
plot(glasso_graph_mbic, layout = coords_exchange)
title(main = "EGLasso", cex.main = 2)
dev.off()



glasso_fit_select <- glasso_bic(X, lambda.range = 0.7, q.threshold = p, thres = 1e-3, iter.max = 1000, BIC_method = "maxMBIC")


glasso_graph_select <- graph_from_adjacency_matrix((glasso_fit_select$Theta != 0), diag = FALSE, mode = "undirected")

vertex_attr(glasso_graph_select) <- list(name = colnames(X))

fig_file <-paste(figpath, "application_exchange_select.pdf", sep = "")
pdf(fig_file, width = 12, height = 8)

par(mfrow=c(1,2))

plot(eglearn_fit_select, layout = coords_exchange)
title(main = "EGLearn", cex.main = 2)
plot(glasso_graph_select, layout = coords_exchange)
title(main = "EGLasso", cex.main = 2)
dev.off()
