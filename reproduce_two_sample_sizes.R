############################################################
## Reproducible simulations for two sample sizes
##
## This script runs one synthetic simulation for I = 400 and
## one synthetic simulation for I = 800.
##
## For each sample size, the script:
## 1. generates one synthetic family-structured data set;
## 2. prepares the univariate and bivariate likelihood data;
## 3. estimates the proposed model parameters;
## 4. computes the robust variance-covariance matrix;
## 5. saves the estimated parameters, variance diagonal,
##    and standard errors.
##
## The goal is to provide reviewers with a direct reproducible
## demonstration that the proposed estimation procedure runs.
############################################################

rm(list = ls())

############################################################
## 0. Required packages
############################################################

library(boot)
library(survival)
library(cubature)
library(VineCopula)
library(numDeriv)
library(condMVNorm)
library(pbivnorm)
library(mvtnorm)

############################################################
## 1. Source functions
############################################################

source("R/Stat_Des.R")
source("R/Proba.R")
source("R/Generate.Covariates.R")
source("R/Generate.Survival.R")
source("R/kinship.R")
source("R/indices.R")
source("R/Fonctions.R")

source("R/Neg_loglik_Biv_bis_vec.R")
source("R/Neg_loglik_uni_vec.R")

source("R/Neg_loglik_Biv_bis.R")
source("R/Neg_loglik_uni.R")
source("R/D2.R")

source("R/Optimisation_bi.R")
source("R/Optimisation_Uni.R")
source("R/Optimisation_2Optim.R")

source("R/Poids_Biv.R")
source("R/Poids_Uni.R")

source("R/F_Bsurv.R")
source("R/DerivCop.R")

############################################################
## 2. Output folders
############################################################

dir.create("output_two_sample_sizes", recursive = TRUE, showWarnings = FALSE)
dir.create("output_two_sample_sizes/data", recursive = TRUE, showWarnings = FALSE)
dir.create("output_two_sample_sizes/estimates", recursive = TRUE, showWarnings = FALSE)
dir.create("output_two_sample_sizes/variances", recursive = TRUE, showWarnings = FALSE)
dir.create("output_two_sample_sizes/summary", recursive = TRUE, showWarnings = FALSE)
dir.create("output_two_sample_sizes/errors", recursive = TRUE, showWarnings = FALSE)

############################################################
## 3. Reproducibility seed
############################################################

set.seed(12345)

############################################################
## 4. True marginal model parameters
############################################################

lambda1 <- 1 / 1000
alpha1  <- 4
beta11  <- 1
beta12  <- -1

theta10 <- c(
  log(alpha1),
  log(lambda1),
  beta11,
  beta12
)

lambda2 <- 1 / 970
alpha2  <- 3.55
beta21  <- 1
beta22  <- -1

theta20 <- c(
  log(alpha2),
  log(lambda2),
  beta21,
  beta22
)

## In this reproducible example, both causes are generated
## using the same marginal parameter vector.
theta20 <- theta10

params0 <- c(theta10, theta20)

############################################################
## 5. True dependence parameters
############################################################

h11 <- 0.4
h12 <- 0.2
h22 <- 0.4

H <- c(h11, h12, h22)

############################################################
## 6. Simulation settings
############################################################

m.a <- 4
v.a <- 2

p <- 0.5

min.c <- 5
max.c <- 10

############################################################
## 7. Utility function: safe standard errors
############################################################

compute_standard_errors <- function(vcov_mat) {
  
  if (!is.matrix(vcov_mat)) {
    stop("Var_asymp must return a variance-covariance matrix.")
  }
  
  variance_diag <- diag(vcov_mat)
  
  if (any(variance_diag < 0, na.rm = TRUE)) {
    warning("Some diagonal variance estimates are negative.")
  }
  
  se_hat <- sqrt(pmax(variance_diag, 0))
  
  return(
    list(
      variance_diag = variance_diag,
      standard_errors = se_hat
    )
  )
}

############################################################
## 8. Function to run one simulation for a given I
############################################################

run_one_simulation <- function(I_value) {
  
  cat("\n========================================\n")
  cat("Running one simulation for I =", I_value, "\n")
  cat("========================================\n")
  
  ##########################################################
  ## 8.1 Create indices
  ##########################################################
  
  indices.proband <- creer.indices.proband(I_value)
  
  indices.non.proband <- creer.indices.non.proband(I_value)
  
  indices.bivariee <- creer.indices.bivariee(I_value)
  
  ##########################################################
  ## 8.2 Generate synthetic family-structured data
  ##########################################################
  
  data <- Generate.data(
    I = I_value,
    theta10 = theta10,
    theta20 = theta20,
    kinship = kinship,
    h11 = h11,
    h22 = h22,
    h12 = h12,
    m.a = m.a,
    v.a = v.a,
    p = p,
    min.c = min.c,
    max.c = max.c
  )
  
  saveRDS(
    data,
    file = paste0(
      "output_two_sample_sizes/data/simulated_data_I",
      I_value,
      ".rds"
    )
  )
  
  ##########################################################
  ## 8.3 Prepare univariate data list
  ##########################################################
  
  Data_Uni <- prepare_uni_datalist(
    data = data,
    indices.proband = indices.proband,
    indices.non.proband = indices.non.proband,
    kin = kin,
    I = I_value
  )
  
  ##########################################################
  ## 8.4 Prepare bivariate data list
  ##########################################################
  
  Data_Biv <- create.data.bivarie(
    data = data,
    indices.bivariee = indices.bivariee,
    kinJK = kinJK,
    I = I_value
  )
  
  D <- data_bivariate_list(Data_Biv)
  
  ##########################################################
  ## 8.5 Estimate model parameters
  ##########################################################
  
  result <- compute.mle(
    params = params0,
    Data_Uni = Data_Uni,
    D = D,
    H = H
  )
  
  theta1_hat <- result[1:4]
  theta2_hat <- result[5:8]
  params_hat <- result[1:8]
  H_hat      <- result[9:11]
  
  ##########################################################
  ## 8.6 Generate survival probability blocks
  ##########################################################
  
  B_surv <- generer_B_surv(
    D = D,
    t1 = theta1_hat,
    t2 = theta2_hat
  )
  
  ##########################################################
  ## 8.7 Compute bivariate likelihood weights
  ##########################################################
  
  WB <- list(
    B1 = Poid.Biv101_E11(
      params = params_hat,
      D$B101$E11,
      h11 = H_hat[1],
      h12 = H_hat[2]
    ),
    
    B2 = Poid.Biv101_E21(
      params = params_hat,
      D$B101$E21,
      h11 = H_hat[1],
      h12 = H_hat[2],
      h22 = H_hat[3]
    ),
    
    B3 = Poid.Biv101_E12(
      params = params_hat,
      D$B101$E12,
      h11 = H_hat[1],
      h12 = H_hat[2],
      h22 = H_hat[3]
    ),
    
    B4 = Poid.Biv101_E22(
      params = params_hat,
      D$B101$E22,
      h12 = H_hat[2],
      h22 = H_hat[3]
    ),
    
    B5 = Poid.Biv110_E11(
      params = params_hat,
      D$B110$E11,
      h11 = H_hat[1],
      h12 = H_hat[2]
    ),
    
    B6 = Poid.Biv110_E21(
      params = params_hat,
      D$B110$E21,
      h11 = H_hat[1],
      h12 = H_hat[2],
      h22 = H_hat[3]
    ),
    
    B7 = Poid.Biv110_E12(
      params = params_hat,
      D$B110$E12,
      h11 = H_hat[1],
      h12 = H_hat[2],
      h22 = H_hat[3]
    ),
    
    B8 = Poid.Biv110_E22(
      params = params_hat,
      D$B110$E22,
      h12 = H_hat[2],
      h22 = H_hat[3]
    ),
    
    B9 = Poid.Biv100_E1(
      params = params_hat,
      D$B100$E1,
      h11 = H_hat[1],
      h12 = H_hat[2],
      h22 = H_hat[3]
    ),
    
    B10 = Poid.Biv100_E2(
      params = params_hat,
      D$B100$E2,
      h11 = H_hat[1],
      h12 = H_hat[2],
      h22 = H_hat[3]
    )
  )
  
  ##########################################################
  ## 8.8 Compute univariate likelihood weights
  ##########################################################
  
  WU.1 <- Uni.poids.delta0.E1(
    params = params_hat,
    data.non.proband.delta.0.E.proband1 = Data_Uni$d0_P1,
    h11 = H_hat[1],
    h12 = H_hat[2]
  )
  
  WU.2 <- Uni.poids.delta0.E2(
    params = params_hat,
    data.non.proband.delta.0.E.proband2 = Data_Uni$d0_P2,
    h12 = H_hat[2],
    h22 = H_hat[3]
  )
  
  ##########################################################
  ## 8.9 Robust variance-covariance matrix
  ##########################################################
  
  vcov_hat <- Var_asymp(
    Theta = result,
    Data_Uni = Data_Uni,
    D = D,
    B_surv = B_surv,
    I = I_value,
    WU.1 = WU.1,
    WU.2 = WU.2,
    WB = WB
  )
  
  variance_outputs <- compute_standard_errors(vcov_hat)
  
  variance_diag <- variance_outputs$variance_diag
  
  se_hat <- variance_outputs$standard_errors
  
  ##########################################################
  ## 8.10 Save outputs
  ##########################################################
  
  saveRDS(
    result,
    file = paste0(
      "output_two_sample_sizes/estimates/parameter_estimates_I",
      I_value,
      ".rds"
    )
  )
  
  saveRDS(
    vcov_hat,
    file = paste0(
      "output_two_sample_sizes/variances/variance_covariance_I",
      I_value,
      ".rds"
    )
  )
  
  saveRDS(
    variance_diag,
    file = paste0(
      "output_two_sample_sizes/variances/variance_diagonal_I",
      I_value,
      ".rds"
    )
  )
  
  saveRDS(
    se_hat,
    file = paste0(
      "output_two_sample_sizes/variances/standard_errors_I",
      I_value,
      ".rds"
    )
  )
  
  ##########################################################
  ## 8.11 Build and save summary table
  ##########################################################
  
  summary_table <- data.frame(
    parameter = paste0("theta_", seq_along(result)),
    estimate = as.numeric(result),
    variance = as.numeric(variance_diag),
    standard_error = as.numeric(se_hat)
  )
  
  saveRDS(
    summary_table,
    file = paste0(
      "output_two_sample_sizes/summary/summary_table_I",
      I_value,
      ".rds"
    )
  )
  
  write.csv(
    summary_table,
    file = paste0(
      "output_two_sample_sizes/summary/summary_table_I",
      I_value,
      ".csv"
    ),
    row.names = FALSE
  )
  
  ##########################################################
  ## 8.12 Return results
  ##########################################################
  
  cat("Simulation for I =", I_value, "completed successfully.\n")
  
  return(
    list(
      I = I_value,
      parameter_estimates = result,
      variance_covariance = vcov_hat,
      variance_diagonal = variance_diag,
      standard_errors = se_hat,
      summary_table = summary_table
    )
  )
}

############################################################
## 9. Run simulations for I = 400 and I = 800
############################################################

error_log <- list()

results_I400 <- tryCatch(
  run_one_simulation(I_value = 400),
  error = function(e) {
    error_log$I400 <<- conditionMessage(e)
    NULL
  }
)

results_I800 <- tryCatch(
  run_one_simulation(I_value = 800),
  error = function(e) {
    error_log$I800 <<- conditionMessage(e)
    NULL
  }
)

############################################################
## 10. Save combined results
############################################################

combined_results <- list(
  I400 = results_I400,
  I800 = results_I800
)

saveRDS(
  combined_results,
  file = "output_two_sample_sizes/summary/combined_results_I400_I800.rds"
)

saveRDS(
  error_log,
  file = "output_two_sample_sizes/errors/error_log.rds"
)

############################################################
## 11. Save session information
############################################################

writeLines(
  capture.output(sessionInfo()),
  con = "output_two_sample_sizes/sessionInfo_two_sample_sizes.txt"
)

############################################################
## 12. End
############################################################

cat("\nTwo sample-size simulations completed.\n")
cat("Results are saved in output_two_sample_sizes/.\n")
cat("Session information saved in output_two_sample_sizes/sessionInfo_two_sample_sizes.txt.\n")

if (length(error_log) > 0) {
  cat("\nSome simulations failed. See output_two_sample_sizes/errors/error_log.rds.\n")
}