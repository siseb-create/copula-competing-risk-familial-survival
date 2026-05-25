# Copula-Based Competing Risk Model for Familial Survival Data

This repository contains the R code used to reproduce the simulation study associated with the second project of the thesis.

The project focuses on a competing risk survival model for familial data, where the first event may correspond either to breast cancer or ovarian cancer. Familial dependence is modeled using a Gaussian copula whose correlation structure is parameterized through kinship coefficients.

## Repository structure

- `R/`: R functions used for data generation, marginal competing risk functions, copula calculations, likelihood evaluation, optimization, and variance estimation.
- `reproduce_one_simulation.R`: main script used to generate one reproducible simulated data set and estimate the model.
- `reproduce_B_simulations.R`: optional script for repeated simulations.
- `output/`: folder where results from one simulation are saved locally.
- `output_B/`: folder where results from repeated simulations are saved locally.
- `figures/`: folder where generated figures are saved locally.
- `sessionInfo.txt`: information about the R session used to run the code.

## Two sample-size reproducibility script

The file `reproduce_two_sample_sizes.R` runs one synthetic simulation for \(I = 400\) families and one synthetic simulation for \(I = 800\) families.

For each sample size, the script generates a synthetic family-structured survival data set, estimates the proposed model, computes the robust variance-covariance matrix, and saves the parameter estimates, variance diagonal, and standard errors.

The script is intended to provide reviewers with a direct demonstration that the proposed estimation procedure runs for two different sample sizes.

The results are saved in the `output_two_sample_sizes/` folder.

## Main script

The main reproducible script is:

```r
source("reproduce_two_sample_sizes.R")

