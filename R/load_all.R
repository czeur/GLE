# Load all project functions in correct dependency order
# Usage: source("R/load_all.R")

# Core packages
library(graphicalExtremes)
library(igraph)

# C++ implementations (compiled once, cached by Rcpp)
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("src/est_W.cpp")
Rcpp::sourceCpp("src/glasso_c.cpp")

# R wrappers (depend on C++ being compiled)
source("R/est_W.R")
source("R/glasso_c.R")

# Core utilities (no internal dependencies)
source("R/f1score.R")

# Estimation functions (depend on glasso_c, est_W)
source("R/glasso_path.R")
source("R/bic.R")

# DGP module (depends on graphicalExtremes, igraph)
source("R/dgp.R")

# Optional: extended methods + plotting
# source("R/eglearn.R")
# source("R/incoherence.R")
# source("R/plotting.R")
