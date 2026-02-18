# Code Review Checklist

Review order follows the dependency chain (bottom-up). Files grouped by how much scrutiny they need.

## Priority 1: C++ core (everything depends on these)

- [ ] `src/est_W.cpp` — Full C++ rewrite of `functions/est_W.R`. Check rank transformation, covariance loop, trim_matrix logic.
- [ ] `src/glasso_c.cpp` — Full C++ rewrite of `functions/glasso_c_function.R` + built-in LASSO solver (replaces glmnet). Check: `lasso_cd()` coordinate descent, `glasso_c_cpp()` main loop, `glasso_c_reest_cpp()` zeroing of Theta_hat entries.

## Priority 2: R wrappers for C++

- [ ] `R/est_W.R` — Thin wrapper, calls `est_W_cpp()`. New file.
- [ ] `R/glasso_c.R` — Thin wrapper, calls `glasso_c_cpp()` and `glasso_c_reest_cpp()`. Converts C++ `umat` (numeric 0/1) to logical. New file.

## Priority 3: New files (no original counterpart)

- [ ] `R/glasso_path.R` — Entirely new. Core estimation workflow: computes W once, loops glasso_c_reest over lambdas. Short (~40 lines).
- [ ] `R/parallel.R` — Rewrite of inline parallel code from `simulations/simulation_BA.R`. Generalized `run_parallel()` with configurable ncores, batch_size, seed, setup_fn, export_vars.
- [ ] `R/output.R` — Entirely new. `create_output_dir()`, `save_config()`, `save_results()`, `save_plot()`.
- [ ] `R/load_all.R` — Entirely new. Sources everything in dependency order. Check the ordering is correct.

## Priority 4: Run scripts (top-level, use everything)

- [ ] `run/trial_single_sample.R` — Entirely new. Single-sample F1 vs lambda exploration. Saves .rds + .pdf.
- [ ] `run/run_simulation.R` — Entirely new. Multi-sample parallel simulation. Saves .rds + config.yaml + plots.

## Priority 5: Moderate changes (extracted + parametrized)

- [ ] `R/dgp.R` — Extracted from `simulations/simulation_BA.R` (BA graph) and `simulations/highd_simulation.R` (tree). Wrapped into `generate_graph()` and `generate_data()`. Check parameter handling matches originals.
- [ ] `R/plotting.R` — Extracted from `applications/functions_paper.R`. Changed `size` to `linewidth` in `theme_fct()`. Subset of original (only `set_graph_parameters`, `save_myplot`, `my_palette`, `theme_fct`).

## Priority 6: Trivial changes (quick scan)

- [ ] `R/f1score.R` — Near-identical to `functions/F1score.R`. Minor formatting.
- [ ] `R/bic.R` — Near-identical to `functions/BIC.R`. Removed `library(caret)`, style cleanup.
- [ ] `R/eglearn.R` — Verbatim extraction from `applications/functions_paper.R` (`mychol`, `glasso_mb2`, `eglearn2`).
- [ ] `R/incoherence.R` — Verbatim extraction from `applications/functions_paper.R` (`GLNSi`, `Gamma2Inc`, `incoherence`, `incoherence_pess`).

## Priority 7: Documentation

- [ ] `CLAUDE.md`
- [ ] `README.md`
