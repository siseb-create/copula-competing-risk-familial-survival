############################################################
## Optional B-simulation script
##
## NOT RUN BY DEFAULT.
##
## This script runs repeated simulations with B = 30
## replications by default.
##
## For each sample size, I = 400 and I = 800, the script
## generates B synthetic family-structured survival data sets,
## estimates the proposed model, computes the robust
## variance-covariance matrix, and saves the parameter
## estimates, variance diagonals, and standard errors.
##
## This script can be computationally expensive and may take
## substantial time to run.
##
## If executed, all outputs are saved locally in the output_B/
## folder, which is created automatically.
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

dir.create("output_B", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/data", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/estimates", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/variances", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/summary", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/errors", recursive = TRUE, showWarnings = FALSE)

############################################################
## 3. Reproducibility seed
############################################################

set.seed(678910)

############################################################
## 4. B-simulation settings
############################################################

B <- 30

sample_sizes <- c(400, 800)

############################################################
## 5. True marginal model parameters
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

## In this reproducible example, both causes are generated
## using the same marginal parameter vector.
theta20 <- theta10

params0 <- c(theta10, theta20)

############################################################
## 6. True dependence parameters
############################################################

h11 <- 0.4
h12 <- 0.2
h22 <- 0.4

H <- c(h11, h12, h22)

############################################################
## 7. Simulation settings
############################################################

m.a <- 4
v.a <- 2

p <- 0.5

min.c <- 5
max.c <- 10

############################################################
## 8. Utility function: compute variance diagonal and SE
############################################################

compute_standard_errors <- function(vcov_mat) {
  
  if (!is.matrix(vcov_mat)) {
    stop("Var_asymp must return a variance-covariance matrix.")
  }
  
  variance_diag <- diag(vcov_mat)
  
  if (any(variance_diag < 0, na.rm = TRUE)) {
    warning("Some diagonal variance estimates are negative. Negative values are truncated to 0 before taking square roots.")
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
## 9. Function to run one simulation for a given I
############################################################

run_one_simulation <- function(I_value, b_value) {
  
  cat("\n----------------------------------------\n")
  cat("Running simulation b =", b_value, "for I =", I_value, "\n")
  cat("----------------------------------------\n")
  
  ##########################################################
  ## 9.1 Create indices
  ##########################################################
  
  indices.proband <- creer.indices.proband(I_value)
  
  indices.non.proband <- creer.indices.non.proband(I_value)
  
  indices.bivariee <- creer.indices.bivariee(I_value)
  
  ##########################################################
  ## 9.2 Generate synthetic family-structured data
  ##########################################################
  
  data <- Generate.data(
    I = I_value,
    theta1 = theta10,
    theta2 = theta20,
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
      "output_B/data/simulated_data_I",
      I_value,
      "_b",
      b_value,
      ".rds"
    )
  )
  
  ##########################################################
  ## 9.3 Prepare univariate data list
  ##########################################################
  
  Data_Uni <- prepare_uni_datalist(
    data = data,
    indices.proband = indices.proband,
    indices.non.proband = indices.non.proband,
    kin = kin,
    I = I_value
  )
  
  ##########################################################
  ## 9.4 Prepare bivariate data list
  ##########################################################
  
  Data_Biv <- create.data.bivarie(
    data = data,
    indices.bivariee = indices.bivariee,
    kinJK = kinJK,
    I = I_value
  )
  
  D <- data_bivariate_list(Data_Biv)
  
  ##########################################################
  ## 9.5 Estimate model parameters
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
  ## 9.6 Generate survival probability blocks
  ##########################################################
  
  B_surv <- generer_B_surv(
    D = D,
    t1 = theta1_hat,
    t2 = theta2_hat
  )
  
  ##########################################################
  ## 9.7 Compute bivariate likelihood weights
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
  ## 9.8 Compute univariate likelihood weights
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
  ## 9.9 Robust variance-covariance matrix
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
  ## 9.10 Summary table for this replication
  ##########################################################
  
  summary_table <- data.frame(
    parameter = paste0("theta_", seq_along(result)),
    estimate = as.numeric(result),
    variance = as.numeric(variance_diag),
    standard_error = as.numeric(se_hat)
  )
  
  ##########################################################
  ## 9.11 Save replication-level outputs
  ##########################################################
  
  saveRDS(
    result,
    file = paste0(
      "output_B/estimates/parameter_estimates_I",
      I_value,
      "_b",
      b_value,
      ".rds"
    )
  )
  
  saveRDS(
    vcov_hat,
    file = paste0(
      "output_B/variances/variance_covariance_I",
      I_value,
      "_b",
      b_value,
      ".rds"
    )
  )
  
  saveRDS(
    variance_diag,
    file = paste0(
      "output_B/variances/variance_diagonal_I",
      I_value,
      "_b",
      b_value,
      ".rds"
    )
  )
  
  saveRDS(
    se_hat,
    file = paste0(
      "output_B/variances/standard_errors_I",
      I_value,
      "_b",
      b_value,
      ".rds"
    )
  )
  
  saveRDS(
    summary_table,
    file = paste0(
      "output_B/summary/summary_table_I",
      I_value,
      "_b",
      b_value,
      ".rds"
    )
  )
  
  write.csv(
    summary_table,
    file = paste0(
      "output_B/summary/summary_table_I",
      I_value,
      "_b",
      b_value,
      ".csv"
    ),
    row.names = FALSE
  )
  
  cat("Simulation b =", b_value, "for I =", I_value, "completed successfully.\n")
  
  return(
    list(
      I = I_value,
      b = b_value,
      parameter_estimates = result,
      variance_covariance = vcov_hat,
      variance_diagonal = variance_diag,
      standard_errors = se_hat,
      summary_table = summary_table
    )
  )
}

############################################################
## 10. Run B simulations for each sample size
############################################################

all_results <- list()

error_log <- list()

for (I_value in sample_sizes) {
  
  cat("\n========================================\n")
  cat("Starting B simulations for I =", I_value, "\n")
  cat("Number of replications B =", B, "\n")
  cat("========================================\n")
  
  results_I <- vector("list", B)
  
  parameter_matrix <- NULL
  variance_diag_matrix <- NULL
  standard_error_matrix <- NULL
  
  for (b in seq_len(B)) {
    
    fit_b <- tryCatch(
      run_one_simulation(
        I_value = I_value,
        b_value = b
      ),
      error = function(e) {
        
        error_message <- paste0(
          "Failure for I = ",
          I_value,
          ", b = ",
          b,
          ": ",
          conditionMessage(e)
        )
        
        cat(error_message, "\n")
        
        error_log[[paste0("I", I_value, "_b", b)]] <<- error_message
        
        saveRDS(
          error_log,
          file = "output_B/errors/error_log.rds"
        )
        
        return(NULL)
      }
    )
    
    if (is.null(fit_b)) {
      next
    }
    
    results_I[[b]] <- fit_b
    
    parameter_matrix <- rbind(
      parameter_matrix,
      fit_b$parameter_estimates
    )
    
    variance_diag_matrix <- rbind(
      variance_diag_matrix,
      fit_b$variance_diagonal
    )
    
    standard_error_matrix <- rbind(
      standard_error_matrix,
      fit_b$standard_errors
    )
    
    ########################################################
    ## Save running results after each successful replication
    ########################################################
    
    saveRDS(
      results_I,
      file = paste0(
        "output_B/summary/results_I",
        I_value,
        ".rds"
      )
    )
    
    saveRDS(
      parameter_matrix,
      file = paste0(
        "output_B/estimates/parameter_estimates_I",
        I_value,
        ".rds"
      )
    )
    
    saveRDS(
      variance_diag_matrix,
      file = paste0(
        "output_B/variances/variance_diagonal_I",
        I_value,
        ".rds"
      )
    )
    
    saveRDS(
      standard_error_matrix,
      file = paste0(
        "output_B/variances/standard_errors_I",
        I_value,
        ".rds"
      )
    )
  }
  
  ############################################################
  ## 11. Monte Carlo summaries for this sample size
  ############################################################
  
  if (!is.null(parameter_matrix)) {
    
    mean_parameter_estimates <- colMeans(
      parameter_matrix,
      na.rm = TRUE
    )
    
    mean_variance_diagonal <- colMeans(
      variance_diag_matrix,
      na.rm = TRUE
    )
    
    mean_standard_errors <- colMeans(
      standard_error_matrix,
      na.rm = TRUE
    )
    
    summary_table <- data.frame(
      parameter = paste0("theta_", seq_along(mean_parameter_estimates)),
      mean_estimate = as.numeric(mean_parameter_estimates),
      mean_variance = as.numeric(mean_variance_diagonal),
      mean_standard_error = as.numeric(mean_standard_errors)
    )
    
    saveRDS(
      mean_parameter_estimates,
      file = paste0(
        "output_B/summary/mean_parameter_estimates_I",
        I_value,
        ".rds"
      )
    )
    
    saveRDS(
      mean_variance_diagonal,
      file = paste0(
        "output_B/summary/mean_variance_diagonal_I",
        I_value,
        ".rds"
      )
    )
    
    saveRDS(
      mean_standard_errors,
      file = paste0(
        "output_B/summary/mean_standard_errors_I",
        I_value,
        ".rds"
      )
    )
    
    saveRDS(
      summary_table,
      file = paste0(
        "output_B/summary/summary_table_I",
        I_value,
        ".rds"
      )
    )
    
    write.csv(
      summary_table,
      file = paste0(
        "output_B/summary/summary_table_I",
        I_value,
        ".csv"
      ),
      row.names = FALSE
    )
    
  } else {
    
    summary_table <- NULL
    
    error_log[[paste0("I", I_value, "_summary")]] <- paste0(
      "No successful simulations for I = ",
      I_value
    )
  }
  
  all_results[[paste0("I", I_value)]] <- list(
    results = results_I,
    parameter_estimates = parameter_matrix,
    variance_diagonal = variance_diag_matrix,
    standard_errors = standard_error_matrix,
    summary_table = summary_table
  )
  
  saveRDS(
    all_results,
    file = "output_B/summary/all_results_B_simulations.rds"
  )
  
  saveRDS(
    error_log,
    file = "output_B/errors/error_log.rds"
  )
}

############################################################
## 12. Save session information
############################################################

writeLines(
  capture.output(sessionInfo()),
  con = "output_B/sessionInfo_B_simulations.txt"
)

############################################################
## 13. End
############################################################

cat("\nB simulations completed.\n")
cat("Default number of replications B =", B, "\n")
cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
cat("Results are saved in output_B/.\n")
cat("Session information saved in output_B/sessionInfo_B_simulations.txt.\n")

if (length(error_log) > 0) {
  cat("Some simulations failed. See output_B/errors/error_log.rds.\n")
}
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

dir.create("output_B", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/data", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/estimates", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/variances", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/summary", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/errors", recursive = TRUE, showWarnings = FALSE)

############################################################
## 3. Reproducibility seed
############################################################

set.seed(12345)

############################################################
## 4. B-simulation settings
############################################################

B <- 30

sample_sizes <- c(400, 800)

############################################################
## 5. True marginal model parameters
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


## In this reproducible example, both causes are generated
## using the same marginal parameter vector.
theta20 <- theta10

params0 <- c(theta10, theta20)

############################################################
## 6. True dependence parameters
############################################################

h11 <- 0.4
h12 <- 0.2
h22 <- 0.4

H <- c(h11, h12, h22)

############################################################
## 7. Simulation settings
############################################################

m.a <- 4
v.a <- 2

p <- 0.5

min.c <- 5
max.c <- 10

############################################################
## 8. Utility function: compute variance diagonal and SE
############################################################

compute_standard_errors <- function(vcov_mat) {
  
  if (!is.matrix(vcov_mat)) {
    stop("Var_asymp must return a variance-covariance matrix.")
  }
  
  variance_diag <- diag(vcov_mat)
  
  if (any(variance_diag < 0, na.rm = TRUE)) {
    warning("Some diagonal variance estimates are negative. Negative values are truncated to 0 before taking square roots.")
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
## 9. Function to run one simulation for a given I
############################################################

run_one_simulation <- function(I_value, b_value) {
  
  cat("\n----------------------------------------\n")
  cat("Running simulation b =", b_value, "for I =", I_value, "\n")
  cat("----------------------------------------\n")
  
  ##########################################################
  ## 9.1 Create indices
  ##########################################################
  
  indices.proband <- creer.indices.proband(I_value)
  
  indices.non.proband <- creer.indices.non.proband(I_value)
  
  indices.bivariee <- creer.indices.bivariee(I_value)
  
  ##########################################################
  ## 9.2 Generate synthetic family-structured data
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
      "output_B/data/simulated_data_I",
      I_value,
      "_b",
      b_value,
      ".rds"
    )
  )
  
  ##########################################################
  ## 9.3 Prepare univariate data list
  ##########################################################
  
  Data_Uni <- prepare_uni_datalist(
    data = data,
    indices.proband = indices.proband,
    indices.non.proband = indices.non.proband,
    kin = kin,
    I = I_value
  )
  
  ##########################################################
  ## 9.4 Prepare bivariate data list
  ##########################################################
  
  Data_Biv <- create.data.bivarie(
    data = data,
    indices.bivariee = indices.bivariee,
    kinJK = kinJK,
    I = I_value
  )
  
  D <- data_bivariate_list(Data_Biv)
  
  ##########################################################
  ## 9.5 Estimate model parameters
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
  ## 9.6 Generate survival probability blocks
  ##########################################################
  
  B_surv <- generer_B_surv(
    D = D,
    t1 = theta1_hat,
    t2 = theta2_hat
  )
  
  ##########################################################
  ## 9.7 Compute bivariate likelihood weights
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
  ## 9.8 Compute univariate likelihood weights
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
  ## 9.9 Robust variance-covariance matrix
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
  ## 9.10 Summary table for this replication
  ##########################################################
  
  summary_table <- data.frame(
    parameter = paste0("theta_", seq_along(result)),
    estimate = as.numeric(result),
    variance = as.numeric(variance_diag),
    standard_error = as.numeric(se_hat)
  )
  
  ##########################################################
  ## 9.11 Save replication-level outputs
  ##########################################################
  
  saveRDS(
    result,
    file = paste0(
      "output_B/estimates/parameter_estimates_I",
      I_value,
      "_b",
      b_value,
      ".rds"
    )
  )
  
  saveRDS(
    vcov_hat,
    file = paste0(
      "output_B/variances/variance_covariance_I",
      I_value,
      "_b",
      b_value,
      ".rds"
    )
  )
  
  saveRDS(
    variance_diag,
    file = paste0(
      "output_B/variances/variance_diagonal_I",
      I_value,
      "_b",
      b_value,
      ".rds"
    )
  )
  
  saveRDS(
    se_hat,
    file = paste0(
      "output_B/variances/standard_errors_I",
      I_value,
      "_b",
      b_value,
      ".rds"
    )
  )
  
  saveRDS(
    summary_table,
    file = paste0(
      "output_B/summary/summary_table_I",
      I_value,
      "_b",
      b_value,
      ".rds"
    )
  )
  
  write.csv(
    summary_table,
    file = paste0(
      "output_B/summary/summary_table_I",
      I_value,
      "_b",
      b_value,
      ".csv"
    ),
    row.names = FALSE
  )
  
  cat("Simulation b =", b_value, "for I =", I_value, "completed successfully.\n")
  
  return(
    list(
      I = I_value,
      b = b_value,
      parameter_estimates = result,
      variance_covariance = vcov_hat,
      variance_diagonal = variance_diag,
      standard_errors = se_hat,
      summary_table = summary_table
    )
  )
}

############################################################
## 10. Run B simulations for each sample size
############################################################

all_results <- list()

error_log <- list()

for (I_value in sample_sizes) {
  
  cat("\n========================================\n")
  cat("Starting B simulations for I =", I_value, "\n")
  cat("Number of replications B =", B, "\n")
  cat("========================================\n")
  
  results_I <- vector("list", B)
  
  parameter_matrix <- NULL
  variance_diag_matrix <- NULL
  standard_error_matrix <- NULL
  
  for (b in seq_len(B)) {
    
    fit_b <- tryCatch(
      run_one_simulation(
        I_value = I_value,
        b_value = b
      ),
      error = function(e) {
        
        error_message <- paste0(
          "Failure for I = ",
          I_value,
          ", b = ",
          b,
          ": ",
          conditionMessage(e)
        )
        
        cat(error_message, "\n")
        
        error_log[[paste0("I", I_value, "_b", b)]] <<- error_message
        
        saveRDS(
          error_log,
          file = "output_B/errors/error_log.rds"
        )
        
        return(NULL)
      }
    )
    
    if (is.null(fit_b)) {
      next
    }
    
    results_I[[b]] <- fit_b
    
    parameter_matrix <- rbind(
      parameter_matrix,
      fit_b$parameter_estimates
    )
    
    variance_diag_matrix <- rbind(
      variance_diag_matrix,
      fit_b$variance_diagonal
    )
    
    standard_error_matrix <- rbind(
      standard_error_matrix,
      fit_b$standard_errors
    )
    
    ########################################################
    ## Save running results after each successful replication
    ########################################################
    
    saveRDS(
      results_I,
      file = paste0(
        "output_B/summary/results_I",
        I_value,
        ".rds"
      )
    )
    
    saveRDS(
      parameter_matrix,
      file = paste0(
        "output_B/estimates/parameter_estimates_I",
        I_value,
        ".rds"
      )
    )
    
    saveRDS(
      variance_diag_matrix,
      file = paste0(
        "output_B/variances/variance_diagonal_I",
        I_value,
        ".rds"
      )
    )
    
    saveRDS(
      standard_error_matrix,
      file = paste0(
        "output_B/variances/standard_errors_I",
        I_value,
        ".rds"
      )
    )
  }
  
  ############################################################
  ## 11. Monte Carlo summaries for this sample size
  ############################################################
  
  if (!is.null(parameter_matrix)) {
    
    mean_parameter_estimates <- colMeans(
      parameter_matrix,
      na.rm = TRUE
    )
    
    mean_variance_diagonal <- colMeans(
      variance_diag_matrix,
      na.rm = TRUE
    )
    
    mean_standard_errors <- colMeans(
      standard_error_matrix,
      na.rm = TRUE
    )
    
    summary_table <- data.frame(
      parameter = paste0("theta_", seq_along(mean_parameter_estimates)),
      mean_estimate = as.numeric(mean_parameter_estimates),
      mean_variance = as.numeric(mean_variance_diagonal),
      mean_standard_error = as.numeric(mean_standard_errors)
    )
    
    saveRDS(
      mean_parameter_estimates,
      file = paste0(
        "output_B/summary/mean_parameter_estimates_I",
        I_value,
        ".rds"
      )
    )
    
    saveRDS(
      mean_variance_diagonal,
      file = paste0(
        "output_B/summary/mean_variance_diagonal_I",
        I_value,
        ".rds"
      )
    )
    
    saveRDS(
      mean_standard_errors,
      file = paste0(
        "output_B/summary/mean_standard_errors_I",
        I_value,
        ".rds"
      )
    )
    
    saveRDS(
      summary_table,
      file = paste0(
        "output_B/summary/summary_table_I",
        I_value,
        ".rds"
      )
    )
    
    write.csv(
      summary_table,
      file = paste0(
        "output_B/summary/summary_table_I",
        I_value,
        ".csv"
      ),
      row.names = FALSE
    )
    
  } else {
    
    summary_table <- NULL
    
    error_log[[paste0("I", I_value, "_summary")]] <- paste0(
      "No successful simulations for I = ",
      I_value
    )
  }
  
  all_results[[paste0("I", I_value)]] <- list(
    results = results_I,
    parameter_estimates = parameter_matrix,
    variance_diagonal = variance_diag_matrix,
    standard_errors = standard_error_matrix,
    summary_table = summary_table
  )
  
  saveRDS(
    all_results,
    file = "output_B/summary/all_results_B_simulations.rds"
  )
  
  saveRDS(
    error_log,
    file = "output_B/errors/error_log.rds"
  )
}

############################################################
## 12. Save session information
############################################################

writeLines(
  capture.output(sessionInfo()),
  con = "output_B/sessionInfo_B_simulations.txt"
)

############################################################
## 13. End
############################################################

cat("\nB simulations completed.\n")
cat("Default number of replications B =", B, "\n")
cat("Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
cat("Results are saved in output_B/.\n")
cat("Session information saved in output_B/sessionInfo_B_simulations.txt.\n")

if (length(error_log) > 0) {
  cat("Some simulations failed. See output_B/errors/error_log.rds.\n")
}
