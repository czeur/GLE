# Graphical Lasso for Extremes

R implementation for the paper: Wan, P. & Zhou, C. (2025). *Graphical lasso for extremes*. [arXiv:2307.15004](https://arxiv.org/abs/2307.15004)

## Requirements

- R >= 4.5.0
- **Core** (needed on remote): `Rcpp`, `RcppArmadillo`, `graphicalExtremes`, `igraph`
- **Optional** (for plotting/applications): `glmnet`, `ggplot2`, `egg`, `knitr`, `parallel`

Install core packages:

```r
install.packages(c("Rcpp", "RcppArmadillo", "graphicalExtremes", "igraph"))
```

## File Structure

```
src/                    C++ source (RcppArmadillo)
├── glasso_c.cpp        Core graphical LASSO + built-in LASSO solver
└── est_W.cpp           W matrix estimation via rank transformation

R/                      R functions (sourced, not installed)
├── load_all.R          Entry point: source("R/load_all.R")
├── glasso_c.R          R wrappers → glasso_c(), glasso_c_reest()
├── est_W.R             R wrapper  → est_W()
├── glasso_path.R       glasso_path(): run estimation across lambda range
├── f1score.R           F1score(): graph recovery evaluation
├── dgp.R              generate_graph(), generate_data()
├── bic.R               glasso_bic(): BIC-based lambda selection (optional)
├── parallel.R          run_parallel(): parallel simulation wrapper
├── output.R            create_output_dir(), save_config(), save_results()
├── eglearn.R           eglearn2(), glasso_mb2() (for applications)
├── incoherence.R       GLNSi(), Gamma2Inc() (incoherence measures)
└── plotting.R          ggplot theme, graph viz, save_myplot()

run/                    Executable scripts
├── trial_single_sample.R   Single-sample F1 vs lambda exploration
└── run_simulation.R        Multi-sample parallel simulation

output/                 Generated results (not tracked in git)

functions/              Original R code (kept for reference)
simulations/            Original simulation scripts (Figures 3-5)
applications/           Original application scripts (Figures 6-7)
```

## Dependency Graph

```
src/glasso_c.cpp ──┐
src/est_W.cpp ─────┤
                   ▼
            R/load_all.R  ← sources everything below
                   │
    ┌──────────────┼──────────────┐
    ▼              ▼              ▼
R/est_W.R    R/glasso_c.R   R/f1score.R
    │              │
    └──────┬───────┘
           ▼
    R/glasso_path.R          R/dgp.R
           │                    │
           └────────┬───────────┘
                    ▼
         run/trial_single_sample.R  (single sample)
         run/run_simulation.R       (multi-sample, also uses R/parallel.R + R/output.R)
```

## Workflow

The workflow has two phases: **estimation** (compute-heavy, run on remote) and **plotting** (lightweight, run locally). The boundary is the `.rds` file.

### Step 1: Single-sample trial (explore lambda range)

Run on remote or locally (fast for small d):

```r
source("run/trial_single_sample.R")
```

Edit the configuration at the top of the script:

```r
d <- 20; n <- 5000; q_threshold <- 0.9
graph_type <- "BA"; m <- 1
lambda_range <- 10^seq(-2, 0, by = 0.05)
```

**Output**: `output/<timestamp>/trial_results.rds` + `f1_vs_lambda.pdf`

Use this to identify a good lambda range for each (d, n) setting before running full simulations.

### Step 2: Multi-sample simulation (remote)

Edit `run/run_simulation.R` config, then run on remote:

```r
source("run/run_simulation.R")
```

**Output**: `output/<prefix>_<timestamp>/simulation_results.rds` + `config.yaml` + plots

### Step 3: Plot from saved results (local)

Copy the `output/` folder (or just the `.rds` files) to your local machine. Then:

```r
res <- readRDS("output/<folder>/simulation_results.rds")

# res contains:
#   res$config        — simulation parameters
#   res$lambda_range   — lambda values used
#   res$f1_matrix      — nsim x nlambda matrix of F1 scores
#   res$edges_matrix   — nsim x nlambda matrix of edge counts
#   res$f1_mean, f1_sd — summary statistics
#   res$true_edges     — number of true edges
#   res$Theta_true     — true precision matrix
#   res$elapsed        — total time in seconds

# Example: plot F1 vs lambda
plot(log10(res$lambda_range), res$f1_mean, type = "l",
     xlab = "log10(lambda)", ylab = "Mean F1")
```

For single-sample trials:

```r
res <- readRDS("output/<folder>/trial_results.rds")

# res contains:
#   res$config         — DGP parameters
#   res$lambda_range   — lambda values
#   res$f1_scores      — F1 at each lambda
#   res$n_edges        — estimated edges at each lambda
#   res$true_edges     — true edge count
#   res$path$results   — list of Theta_hat and graph at each lambda
```

### Remote vs Local Summary

| Step | Where | What | Output |
|------|-------|------|--------|
| 1. Explore lambda | Remote (or local if small d) | `trial_single_sample.R` | `trial_results.rds` |
| 2. Simulate | Remote | `run_simulation.R` | `simulation_results.rds` |
| 3. Plot | Local | `readRDS()` + custom plotting | Figures |

**What to transfer**: only the `output/` folder (`.rds` files + `config.yaml`). No code or C++ compilation needed on the local machine for plotting.

## Interactive Use

```r
source("R/load_all.R")

# Generate a graph
g <- generate_graph("BA", d = 20, m = 1)

# Generate data
data <- generate_data("mpareto", n = 5000, d = 20, Gamma = g$Gamma)

# Estimate across lambdas
path <- glasso_path(data, lambda_range = 10^seq(-2, 0, by = 0.1), q_threshold = 0.9)

# Evaluate
f1 <- sapply(path$results, function(r) F1score(g$Theta, r$Theta_hat))
```

## Key Parameters

| Parameter | Typical range | Notes |
|-----------|--------------|-------|
| `d` | 10–200 | Dimension. d=200 takes ~45s per lambda |
| `n` | 1,000–100,000 | Sample size |
| `q_threshold` | 0.85–0.95 | Quantile threshold for exceedances |
| `lambda_range` | `10^seq(-1, 1, ...)` | Penalty. Larger d may need higher lambdas |
| `m` | 1–2 | BA model: m=1 gives a tree |

## Performance

Measured on a single core:

| Setting | DGP time | Per-lambda | 41 lambdas |
|---------|----------|------------|------------|
| d=20, n=5000 | 0.3s | 11ms | 0.4s |
| d=200, n=1938 | 14s | 45s | ~16 min |

The C++ implementation is ~20x faster than the original R+glmnet version.

## Original Paper Figures

The original simulation and application scripts are preserved in `simulations/` and `applications/`:

| Figure | Simulation | Plotting |
|--------|-----------|----------|
| 3 | `simulations/simulation_BA_tuning.R` | `simulations/plotting_tuning.R` |
| 4 | `simulations/simulation_BA.R` | `simulations/plotting_BA.R` |
| 5 | `simulations/simulation_BA_F1.R` | `simulations/plotting_F1.R` |
| 6 | `applications/application_exchange.R` | (self-contained) |
| 7 | `applications/application_danube.R` | (self-contained) |
