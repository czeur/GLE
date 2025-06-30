# R Code for "Graphical lasso for extremes"

This repository contains the R code used in the paper:

Wan, P., & Zhou, C. (2025). Graphical lasso for extremes. https://arxiv.org/abs/2307.15004

## Description

This code was developed to conduct the simulation analyses and applications in this paper. It reproduces the  figures presented in the paper: Figures 3-7.

## File Structure

- `functions/` – Basic functions used for the algorithm
- `simulations/` – R scripts for simulations and visualization
- `applications/` – R scripts for applications
- `README.md` – This file

## Details corresponding to the paper

### Simulation graphs

- `simulations/simulation_BA_turning.R' – produces all simulations results needed for Figure 3
- `simulations/plotting_tuning.R' – uses the simulated results above to produce Figure 3

- `simulations/simulation_BA.R' – produces all simulations results needed for Figure 4
- `simulations/plotting_BA.R' – uses the simulated results above to produce Figure 4

- `simulations/simulation_BA_F1.R' – produces all simulations results needed for Figure 5
- `simulations/plotting_F1.R' – uses the simulated results above to produce Figure 5

### Application graphs

- `applications/application_exchange.R' - produces Figure 6
- `applications/application_danube.R' – produces Figure 7

## Requirements

- R version 4.5.0
