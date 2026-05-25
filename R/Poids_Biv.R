############################################################
## Weight functions for the bivariate likelihood
##
## This file contains the functions used to compute the
## bivariate likelihood weights for the different observation
## patterns. These weights enter the bivariate likelihood
## contributions used to estimate the intra-family dependence
## parameters H = (h11, h12, h22).
############################################################


# ==============================================================================
# 1. OUTILS ET STABILITÉ NUMÉRIQUE
# ==============================================================================

clip01 <- function(x, eps = 1e-10) {
  pmin(pmax(x, eps), 1 - eps)
}

clip_rho <- function(r, eps = 1e-8) {
  pmin(pmax(unlist(r), -1 + eps), 1 - eps)
}

safe_log <- function(x, eps = 1e-300) {
  log(pmax(x, eps))
}

safe_sqrt <- function(x, eps = 1e-12) {
  sqrt(pmax(x, eps))
}

logspace_add2 <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

logspace_add4 <- function(a, b, c, d) {
  m <- pmax(pmax(a, b), pmax(c, d))
  m + log(exp(a - m) + exp(b - m) + exp(c - m) + exp(d - m))
}

empty_weights <- function(k) {
  replicate(k, numeric(0), simplify = FALSE)
}


# ==============================================================================
# 3. CALCUL DES POIDS (E-STEP)
# ==============================================================================

# --- POIDS BIV101 ---
Poid.Biv101_E11 <- function(params, Data, h11, h12) {
  if (nrow(Data) == 0) return(empty_weights(2))
  theta1 <- params[1:4]; theta2 <- params[5:8]
  
  F1<-l.F1_inf( theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)
  F2<-l.F2_inf( theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)
  
  logt1 <- F1 + log.C011(u1=clip01(survie_vec1(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), u2=clip01(survie_vec1(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), u3=clip01(survie_vec1(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), r12=h11*Data$kin1J, r13=h11*Data$kinJK, r23=h11*Data$kin1K)
  logt2 <- F2 + log.C011(u1=clip01(survie_vec2(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), u2=clip01(survie_vec1(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), u3=clip01(survie_vec1(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), r12=h12*Data$kin1J, r13=h12*Data$kinJK, r23=h11*Data$kin1K)
  logden <- logspace_add2(logt1, logt2)
  list(exp(logt1 - logden), exp(logt2 - logden))
}

Poid.Biv101_E21 <- function(params, Data, h11, h12, h22) {
  if (nrow(Data) == 0) return(empty_weights(2))
  theta1 <- params[1:4]; theta2 <- params[5:8]
  
  F1<-l.F1_inf( theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)
  F2<-l.F2_inf( theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)
  
  logt1 <- F1 + log.C011(u1=clip01(survie_vec1(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), u2=clip01(survie_vec2(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), u3=clip01(survie_vec1(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), r12=h12*Data$kin1J, r13=h11*Data$kinJK, r23=h12*Data$kin1K)
  logt2 <- F2 + log.C011(u1=clip01(survie_vec2(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), u2=clip01(survie_vec2(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), u3=clip01(survie_vec1(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), r12=h22*Data$kin1J, r13=h12*Data$kinJK, r23=h12*Data$kin1K)
  logden <- logspace_add2(logt1, logt2)
  list(exp(logt1 - logden), exp(logt2 - logden))
}

Poid.Biv101_E12 <- function(params, Data, h11, h12, h22) {
  if (nrow(Data) == 0) return(empty_weights(2))
  theta1 <- params[1:4]; theta2 <- params[5:8]
  
  F1<-l.F1_inf( theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)
  F2<-l.F2_inf( theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)
  
  logt1 <- F1 + log.C011(u1=clip01(survie_vec1(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), u2=clip01(survie_vec1(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), u3=clip01(survie_vec2(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), r12=h11*Data$kin1J, r13=h12*Data$kinJK, r23=h12*Data$kin1K)
  logt2 <- F2 + log.C011(u1=clip01(survie_vec2(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), u2=clip01(survie_vec1(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), u3=clip01(survie_vec2(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), r12=h12*Data$kin1J, r13=h22*Data$kinJK, r23=h12*Data$kin1K)
  logden <- logspace_add2(logt1, logt2)
  list(exp(logt1 - logden), exp(logt2 - logden))
}

Poid.Biv101_E22 <- function(params, Data, h12, h22) {
  if (nrow(Data) == 0) return(empty_weights(2))
  theta1 <- params[1:4]; theta2 <- params[5:8]
  
  F1<-l.F1_inf( theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)
  F2<-l.F2_inf( theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)
  
  logt1 <- F1 + log.C011(u1=clip01(survie_vec1(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), u2=clip01(survie_vec2(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), u3=clip01(survie_vec2(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), r12=h12*Data$kin1J, r13=h12*Data$kinJK, r23=h22*Data$kin1K)
  logt2 <- F2 + log.C011(u1=clip01(survie_vec2(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), u2=clip01(survie_vec2(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), u3=clip01(survie_vec2(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), r12=h22*Data$kin1J, r13=h22*Data$kinJK, r23=h22*Data$kin1K)
  logden <- logspace_add2(logt1, logt2)
  list(exp(logt1 - logden), exp(logt2 - logden))
}

# --- POIDS BIV110 ---
Poid.Biv110_E11 <- function(params, Data, h11, h12) {
  if (nrow(Data) == 0) return(empty_weights(2))
  theta1 <- params[1:4]; theta2 <- params[5:8]
  
  F1<-l.F1_inf(theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)
  F2<-l.F2_inf(theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)
  
  logt1 <- F1 + log.C011(clip01(survie_vec1(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), clip01(survie_vec1(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), clip01(survie_vec1(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), h11*Data$kin1K, h11*Data$kinJK, h11*Data$kin1J)
  logt2 <- F2 + log.C011(clip01(survie_vec2(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), clip01(survie_vec1(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), clip01(survie_vec1(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), h12*Data$kin1K, h12*Data$kinJK, h11*Data$kin1J)
  logden <- logspace_add2(logt1, logt2)
  list(exp(logt1 - logden), exp(logt2 - logden))
}

Poid.Biv110_E21 <- function(params, Data, h11, h12, h22) {
  if (nrow(Data) == 0) return(empty_weights(2))
  theta1 <- params[1:4]; theta2 <- params[5:8]
  
  F1<-l.F1_inf(theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)
  F2<-l.F2_inf(theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)
  
  logt1 <- F1 + log.C011(clip01(survie_vec1(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), clip01(survie_vec2(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), clip01(survie_vec1(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), h12*Data$kin1K, h11*Data$kinJK, h12*Data$kin1J)
  logt2 <- F2 + log.C011(clip01(survie_vec2(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), clip01(survie_vec2(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), clip01(survie_vec1(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), h22*Data$kin1K, h12*Data$kinJK, h12*Data$kin1J)
  logden <- logspace_add2(logt1, logt2)
  list(exp(logt1 - logden), exp(logt2 - logden))
}

Poid.Biv110_E12 <- function(params, Data, h11, h12, h22) {
  if (nrow(Data) == 0) return(empty_weights(2))
  theta1 <- params[1:4]; theta2 <- params[5:8]
  
  F1<-l.F1_inf(theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)
  F2<-l.F2_inf(theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)
  
  logt1 <- F1 + log.C011(clip01(survie_vec1(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), clip01(survie_vec1(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), clip01(survie_vec2(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), h11*Data$kin1K, h12*Data$kinJK, h12*Data$kin1J)
  logt2 <- F2 + log.C011(clip01(survie_vec2(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), clip01(survie_vec1(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), clip01(survie_vec2(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), h12*Data$kin1K, h22*Data$kinJK, h12*Data$kin1J)
  logden <- logspace_add2(logt1, logt2)
  list(exp(logt1 - logden), exp(logt2 - logden))
}

Poid.Biv110_E22 <- function(params, Data, h12, h22) {
  if (nrow(Data) == 0) return(empty_weights(2))
  theta1 <- params[1:4]; theta2 <- params[5:8]
  
  F1<-l.F1_inf(theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)
  F2<-l.F2_inf(theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)
  
  logt1 <- F1 + log.C011(clip01(survie_vec1(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), clip01(survie_vec2(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), clip01(survie_vec2(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), h12*Data$kin1K, h12*Data$kinJK, h22*Data$kin1J)
  logt2 <- F2 + log.C011(clip01(survie_vec2(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), clip01(survie_vec2(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband)), clip01(survie_vec2(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), h22*Data$kin1K, h22*Data$kinJK, h22*Data$kin1J)
  logden <- logspace_add2(logt1, logt2)
  list(exp(logt1 - logden), exp(logt2 - logden))
}

# --- POIDS BIV100 ---
Poid.Biv100_E1 <- function(params, Data, h11, h12, h22) {
  if (nrow(Data) == 0) return(empty_weights(4))
  theta1 <- params[1:4]; theta2 <- params[5:8]
  u_prob <- clip01(survie_vec1(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband))
  
  FJ1<-l.F1_inf( theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)
  FJ2<-l.F2_inf( theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)
  
  FK1<-l.F1_inf(theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)
  FK2<-l.F2_inf(theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)
  
  logt1 <- FJ1 + FK1 + log.C001(clip01(survie_vec1(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), clip01(survie_vec1(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), u_prob, h11*Data$kinJK, h11*Data$kin1J, h11*Data$kin1K)
  logt2 <- FJ2 + FK1 + log.C001(clip01(survie_vec2(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), clip01(survie_vec1(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), u_prob, h12*Data$kinJK, h11*Data$kin1J, h12*Data$kin1K)
  logt3 <- FK1 + FJ1 + log.C001(clip01(survie_vec1(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), clip01(survie_vec2(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), u_prob, h11*Data$kinJK, h12*Data$kin1J, h12*Data$kin1K)
  logt4 <- FK2 + FJ2 + log.C001(clip01(survie_vec2(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), clip01(survie_vec2(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), u_prob, h12*Data$kinJK, h12*Data$kin1J, h22*Data$kin1K)
  
  logden <- logspace_add4(logt1, logt2, logt3, logt4)
  list(exp(logt1 - logden), exp(logt2 - logden), exp(logt3 - logden), exp(logt4 - logden))
}

Poid.Biv100_E2 <- function(params, Data, h11, h12, h22) {
  if (nrow(Data) == 0) return(empty_weights(4))
  theta1 <- params[1:4]; theta2 <- params[5:8]
  u_prob <- clip01(survie_vec2(Data$Y.proband, theta1, theta2, Data$X1.proband, Data$X2.proband))
  
  FJ1<-l.F1_inf( theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)
  FJ2<-l.F2_inf( theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)
  
  FK1<-l.F1_inf(theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)
  FK2<-l.F2_inf(theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)
  
  logt1 <- FJ1 + FK1 + log.C001(clip01(survie_vec1(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), clip01(survie_vec1(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), u_prob, h11*Data$kinJK, h11*Data$kin1J, h12*Data$kin1K)
  logt2 <- FJ2 + FK1 + log.C001(clip01(survie_vec2(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), clip01(survie_vec1(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), u_prob, h12*Data$kinJK, h12*Data$kin1J, h22*Data$kin1K)
  logt3 <- FK1 + FJ1 + log.C001(clip01(survie_vec1(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), clip01(survie_vec2(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), u_prob, h12*Data$kinJK, h12*Data$kin1J, h22*Data$kin1K)
  logt4 <- FK2 + FJ2 + log.C001(clip01(survie_vec2(Data$Y.non.probandIJ, theta1, theta2, Data$X1.non.probandIJ, Data$X2.non.probandIJ)), clip01(survie_vec2(Data$Y.non.probandIK, theta1, theta2, Data$X1.non.probandIK, Data$X2.non.probandIK)), u_prob, h22*Data$kinJK, h22*Data$kin1J, h22*Data$kin1K)
  
  logden <- logspace_add4(logt1, logt2, logt3, logt4)
  list(exp(logt1 - logden), exp(logt2 - logden), exp(logt3 - logden), exp(logt4 - logden))
}
