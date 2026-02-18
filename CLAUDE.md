# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

R codebase for the paper "Graphical lasso for extremes" (Wan & Zhou, 2025, https://arxiv.org/abs/2307.15004). Implements sparse precision matrix estimation under extreme value models using a modified graphical LASSO approach.

## Requirements

- R version 4.5.0+
- Core: `Rcpp`, `RcppArmadillo`, `graphicalExtremes`, `igraph`
- Optional: `glmnet`, `ggplot2`, `egg`, `knitr`, `parallel`

## Architecture

```
src/                C++ implementations (RcppArmadillo)
├── glasso_c.cpp    Core graphical LASSO with built-in LASSO solver (no glmnet dependency)
└── est_W.cpp       W matrix estimation via rank transformation

R/                  Consolidated R functions
├── load_all.R      Sources everything in correct order — use: source("R/load_all.R")
├── glasso_c.R      Thin R wrappers for C++ glasso_c, glasso_c_reest
├── est_W.R         Thin R wrapper for C++ est_W
├── glasso_path.R   Run glasso_c across a range of lambdas (main estimation workflow)
├── bic.R           BIC-based lambda selection (optional utility)
├── f1score.R       F1 score for graph recovery evaluation
├── dgp.R           Graph generators (BA, tree) + data generators (mpareto, maxstable)
├── eglearn.R       Extended eglearn2(), glasso_mb2(), mychol() (from applications)
├── incoherence.R   GLNSi(), Gamma2Inc() incoherence measures
├── plotting.R      ggplot theme, graph visualization, save_myplot()
├── parallel.R      Parallel simulation wrapper (run_parallel)
└── output.R        Experiment tracking: timestamped dirs, config saving

run/                Runnable scripts
├── trial_single_sample.R   Single-sample exploration of F1 vs lambda
└── run_simulation.R        Multi-sample parallel simulation

functions/          Original R implementations (kept for reference)
simulations/        Original simulation scripts (Figures 3-5)
applications/       Original real-data applications (Figures 6-7)
```

### Data pipeline (new workflow)

1. Generate graph + variogram → `generate_graph()` (R/dgp.R)
2. Generate multivariate Pareto data → `generate_data()` (R/dgp.R)
3. Estimate across lambda range → `glasso_path()` (R/glasso_path.R)
4. Evaluate → `F1score()` at each lambda (R/f1score.R)

### Running

```r
# Single-sample trial (explore lambda ranges)
source("run/trial_single_sample.R")

# Multi-sample parallel simulation
source("run/run_simulation.R")

# Or load functions interactively
source("R/load_all.R")
```

## Key Parameters

- Dimension `d`: 10–200
- Sample size `n`: 1,000–100,000
- Quantile threshold `q_threshold`: 0.85–0.95
- Lambda range: typically `10^seq(-2, 0, by = 0.05)`
- Coordinate descent max iterations: 1000

## Original Figure Reproduction

| Figure | Simulation Script | Plotting Script |
|--------|-------------------|-----------------|
| 3 | `simulations/simulation_BA_tuning.R` | `simulations/plotting_tuning.R` |
| 4 | `simulations/simulation_BA.R` | `simulations/plotting_BA.R` |
| 5 | `simulations/simulation_BA_F1.R` | `simulations/plotting_F1.R` |
| 6 | `applications/application_exchange.R` | (self-contained) |
| 7 | `applications/application_danube.R` | (self-contained) |
