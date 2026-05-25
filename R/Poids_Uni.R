############################################################
## Weight functions for the univariate likelihood
##
## This file contains the functions used to compute the
## univariate likelihood weights for the different observation
## patterns.
##
## These weights enter the univariate likelihood contributions
## used to estimate the marginal parameters of the proposed
## survival model, conditionally on the current intra-family
## dependence parameters.
############################################################

# =========================================================
# OUTILS
# =========================================================
clip01 <- function(u, eps = 1e-10) {
  pmin(pmax(u, eps), 1 - eps)
}

clip_rho <- function(rho, eps = 1e-8) {
  pmin(pmax(rho, -1 + eps), 1 - eps)
}

# log(a + b) stable
logspace_add <- function(logx, logy) {
  m <- pmax(logx, logy)
  m + log(exp(logx - m) + exp(logy - m))
}

# =========================================================
# COPULE GAUSSIENNE BIVARIEE
# =========================================================
gauss_cop_logpdf <- function(u1, u2, rho, eps_u = 1e-10, eps_rho = 1e-8) {
  u1  <- clip01(u1, eps_u)
  u2  <- clip01(u2, eps_u)
  rho <- clip_rho(rho, eps_rho)
  
  z1 <- qnorm(u1)
  z2 <- qnorm(u2)
  
  den <- 1 - rho^2
  
  -0.5 * log(den) +
    (2 * rho * z1 * z2 - rho^2 * (z1^2 + z2^2)) / (2 * den)
}

gauss_cop_pdf <- function(u1, u2, rho, eps_u = 1e-10, eps_rho = 1e-8) {
  exp(gauss_cop_logpdf(u1, u2, rho, eps_u, eps_rho))
}

gauss_cop_loghfunc2 <- function(u1, u2, rho, eps_u = 1e-10, eps_rho = 1e-8) {
  u1  <- clip01(u1, eps_u)
  u2  <- clip01(u2, eps_u)
  rho <- clip_rho(rho, eps_rho)
  
  z1 <- qnorm(u1)
  z2 <- qnorm(u2)
  
  arg <- (z1 - rho * z2) / sqrt(1 - rho^2)
  
  pnorm(arg, log.p = TRUE)
}

# =========================================================
# POIDS DELTA0 / E1
# =========================================================
Uni.poids.delta0.E1 <- function(params, data.non.proband.delta.0.E.proband1, h11, h12) {
  theta1 <- params[1:4]
  theta2 <- params[5:8]
  
  # log P(epsilon = 1 | X) et log P(epsilon = 2 | X)
  logterme1 <- l.F1_inf_vec(
    theta1, theta2,
    data.non.proband.delta.0.E.proband1$X1.non.proband,
    data.non.proband.delta.0.E.proband1$X2.non.proband
  )
  
  logterme2 <- l.F2_inf_vec(
    theta1, theta2,
    data.non.proband.delta.0.E.proband1$X1.non.proband,
    data.non.proband.delta.0.E.proband1$X2.non.proband
  )
  
  # composante copule pour la branche 1
  u11 <- survie_vec1(
    data.non.proband.delta.0.E.proband1$Y.non.proband,
    theta1, theta2,
    data.non.proband.delta.0.E.proband1$X1.non.proband,
    data.non.proband.delta.0.E.proband1$X2.non.proband
  )
  
  u12 <- survie_vec1(
    data.non.proband.delta.0.E.proband1$Y.proband,
    theta1, theta2,
    data.non.proband.delta.0.E.proband1$X1.proband,
    data.non.proband.delta.0.E.proband1$X2.proband
  )
  
  logterme3 <- gauss_cop_loghfunc2(
    u11, u12,
    h11 * data.non.proband.delta.0.E.proband1$kin
  )
  
  # composante copule pour la branche 2
  u21 <- survie_vec2(
    data.non.proband.delta.0.E.proband1$Y.non.proband,
    theta1, theta2,
    data.non.proband.delta.0.E.proband1$X1.non.proband,
    data.non.proband.delta.0.E.proband1$X2.non.proband
  )
  
  u22 <- survie_vec1(
    data.non.proband.delta.0.E.proband1$Y.proband,
    theta1, theta2,
    data.non.proband.delta.0.E.proband1$X1.proband,
    data.non.proband.delta.0.E.proband1$X2.proband
  )
  
  logterme4 <- gauss_cop_loghfunc2(
    u21, u22,
    h12 * data.non.proband.delta.0.E.proband1$kin
  )
  
  # poids stables numériquement
  logA <- logterme1 + logterme3 #+ logterme5
  logB <- logterme2 + logterme4 #+ logterme6
  logD <- logspace_add(logA, logB)
  
  W1 <- exp(logA - logD)
  W2 <- exp(logB - logD)
  
  list(W1, W2)
}


# =========================================================
# POIDS DELTA0 / E2
# =========================================================
Uni.poids.delta0.E2 <- function(params, data.non.proband.delta.0.E.proband2, h12, h22) {
  theta1 <- params[1:4]
  theta2 <- params[5:8]
  
  # log P(epsilon = 1 | X) et log P(epsilon = 2 | X)
  logterme1 <- l.F1_inf_vec(
    theta1, theta2,
    data.non.proband.delta.0.E.proband2$X1.non.proband,
    data.non.proband.delta.0.E.proband2$X2.non.proband
  )
  
  logterme2 <- l.F2_inf_vec(
    theta1, theta2,
    data.non.proband.delta.0.E.proband2$X1.non.proband,
    data.non.proband.delta.0.E.proband2$X2.non.proband
  )
  
  # composante copule pour la branche 1
  u11 <- survie_vec1(
    data.non.proband.delta.0.E.proband2$Y.non.proband,
    theta1, theta2,
    data.non.proband.delta.0.E.proband2$X1.non.proband,
    data.non.proband.delta.0.E.proband2$X2.non.proband
  )
  
  u12 <- survie_vec2(
    data.non.proband.delta.0.E.proband2$Y.proband,
    theta1, theta2,
    data.non.proband.delta.0.E.proband2$X1.proband,
    data.non.proband.delta.0.E.proband2$X2.proband
  )
  
  logterme3 <- gauss_cop_loghfunc2(
    u11, u12,
    h12 * data.non.proband.delta.0.E.proband2$kin
  )
  
  # composante copule pour la branche 2
  u21 <- survie_vec2(
    data.non.proband.delta.0.E.proband2$Y.non.proband,
    theta1, theta2,
    data.non.proband.delta.0.E.proband2$X1.non.proband,
    data.non.proband.delta.0.E.proband2$X2.non.proband
  )
  
  u22 <- survie_vec2(
    data.non.proband.delta.0.E.proband2$Y.proband,
    theta1, theta2,
    data.non.proband.delta.0.E.proband2$X1.proband,
    data.non.proband.delta.0.E.proband2$X2.proband
  )
  
  logterme4 <- gauss_cop_loghfunc2(
    u21, u22,
    h22 * data.non.proband.delta.0.E.proband2$kin
  )
  
  # poids stables numériquement
  logA <- logterme1 + logterme3 #+ logterme5
  logB <- logterme2 + logterme4 #+ logterme6
  logD <- logspace_add(logA, logB)
  
  W1 <- exp(logA - logD)
  W2 <- exp(logB - logD)
  
  list(W1, W2)
}
