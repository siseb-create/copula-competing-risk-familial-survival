# =========================================================
# 1. OUTILS NUMERIQUES & COPULES
# =========================================================
fu <- function(u, eps = 1e-10) pmin(pmax(u, eps), 1 - eps)

safe_log <- function(x, eps = 1e-300) log(pmax(x, eps))

gauss_cop_logpdf <- function(u1, u2, rho, eps_u = 1e-10, eps_rho = 1e-8) {
  u1  <- clip01(u1, eps_u)
  u2  <- clip01(u2, eps_u)
  rho <- clip_rho(rho, eps_rho)
  z1 <- qnorm(u1)
  z2 <- qnorm(u2)
  den <- 1 - rho^2
  -0.5 * log(den) + (2 * rho * z1 * z2 - rho^2 * (z1^2 + z2^2)) / (2 * den)
}



# --- PROBANDS ---
NegLogLik.proband.E1_vec <- function(params, data.proband.E1) {
  theta1 <- params[1:4]; theta2 <- params[5:8]
  terme1 <- l.densite1(data.proband.E1$Y, theta1, theta2, data.proband.E1$X1, data.proband.E1$X2)
  terme2 <- safe_log(1 - survie.T(data.proband.E1$a, theta1, theta2, data.proband.E1$X1, data.proband.E1$X2))
  v1<-cbind(-(terme1 - terme2),data.proband.E1$ID)
  return(v1)
}

NegLogLik.proband.E2_vec <- function(params, data.proband.E2) {
  theta1 <- params[1:4]; theta2 <- params[5:8]
  terme1 <- l.densite2(data.proband.E2$Y, theta1, theta2, data.proband.E2$X1, data.proband.E2$X2)
  terme2 <- safe_log(1 - survie.T(data.proband.E2$a, theta1, theta2, data.proband.E2$X1, data.proband.E2$X2))
  v1<-cbind(-(terme1 - terme2),data.proband.E2$ID)
  return(v1)
}


# --- NON-PROBANDS DELTA = 1 ---
# Note : Les fonctions NegLogLik.non.proband.delta1.E... suivent la structure :


NegLogLik.non.proband.delta1.E1.E1_vec <- function(params, data, h11) {
  theta1 <- params[1:4]; theta2 <- params[5:8]
  terme1 <- l.densite1(data$Y.non.proband, theta1, theta2, data$X1.non.proband, data$X2.non.proband)
  u1 <- fu(survie_vec1(data$Y.non.proband, theta1, theta2, data$X1.non.proband, data$X2.non.proband))
  u2 <- fu(survie_vec1(data$Y.proband, theta1, theta2, data$X1.proband, data$X2.proband))
  rho <- h11 * data$kin
  v1<-cbind(-(terme1 + gauss_cop_logpdf(u1, u2, rho)),data$ID)
  return(v1)
}

NegLogLik.non.proband.delta1.E2.E1_vec <- function(params, data, h12) {
  theta1 <- params[1:4]; theta2 <- params[5:8]
  terme1 <- l.densite1(data$Y.non.proband, theta1, theta2, data$X1.non.proband, data$X2.non.proband)
  u1 <- fu(survie_vec1(data$Y.non.proband, theta1, theta2, data$X1.non.proband, data$X2.non.proband))
  u2 <- fu(survie_vec2(data$Y.proband, theta1, theta2, data$X1.proband, data$X2.proband))
  rho <- h12 * data$kin
  v1<-cbind(-(terme1 + gauss_cop_logpdf(u1, u2, rho)),data$ID)
  return(v1)
}

NegLogLik.non.proband.delta1.E1.E2_vec <- function(params, data, h12) {
  theta1 <- params[1:4]; theta2 <- params[5:8]
  terme1 <- l.densite2(data$Y.non.proband, theta1, theta2, data$X1.non.proband, data$X2.non.proband)
  u1 <- fu(survie_vec2(data$Y.non.proband, theta1, theta2, data$X1.non.proband, data$X2.non.proband))
  u2 <- fu(survie_vec1(data$Y.proband, theta1, theta2, data$X1.proband, data$X2.proband))
  rho <- h12 * data$kin
  v1<-cbind(-(terme1 + gauss_cop_logpdf(u1, u2, rho)),data$ID)
  return(v1)
}

NegLogLik.non.proband.delta1.E2.E2_vec <- function(params, data, h22) {
  theta1 <- params[1:4]; theta2 <- params[5:8]
  terme1 <- l.densite2(data$Y.non.proband, theta1, theta2, data$X1.non.proband, data$X2.non.proband)
  u1 <- fu(survie_vec2(data$Y.non.proband, theta1, theta2, data$X1.non.proband, data$X2.non.proband))
  u2 <- fu(survie_vec2(data$Y.proband, theta1, theta2, data$X1.proband, data$X2.proband))
  rho <- h22 * data$kin
  v1<-cbind(-(terme1 + gauss_cop_logpdf(u1, u2, rho)),data$ID)
  return(v1)
}

# --- NON-PROBANDS DELTA = 0 ---
NegLogLik.non.proband.delta0.E1_vec <- function(params, data, h11, h12, WU.1) {
  theta1 <- params[1:4]; theta2 <- params[5:8]
  #data=Data_Uni$d0_P1
  # Cas où NP appartient à E1
  u1_1 <- fu(survie_vec1(data$Y.non.proband, theta1, theta2, data$X1.non.proband, data$X2.non.proband))
  u2_1 <- fu(survie_vec1(data$Y.proband, theta1, theta2, data$X1.proband, data$X2.proband))
  t12 <- gauss_cop_logpdf(u1_1, u2_1, h11 * data$kin)
  # Cas où NP appartient à E2
  u1_2 <- fu(survie_vec2(data$Y.non.proband, theta1, theta2, data$X1.non.proband, data$X2.non.proband))
  t22 <- gauss_cop_logpdf(u1_2, u2_1, h12 * data$kin)
  -sum(WU.1[[1]] * t12, na.rm = TRUE) - sum(WU.1[[2]] * t22, na.rm = TRUE)
  v1<-cbind(-((WU.1[[1]] * t12) + (WU.1[[2]] * t22)),data$ID)
  return(v1)
}

NegLogLik.non.proband.delta0.E2_vec <- function(params, data, h12, h22, WU.2) {
  theta1 <- params[1:4]; theta2 <- params[5:8]
  
  
  # Cas où NP appartient à E1
  u1_1 <- fu(survie_vec1(data$Y.non.proband, theta1, theta2, data$X1.non.proband, data$X2.non.proband))
  u2_2 <- fu(survie_vec2(data$Y.proband, theta1, theta2, data$X1.proband, data$X2.proband))
  t12 <- gauss_cop_logpdf(u1_1, u2_2, h12 * data$kin)
  # Cas où NP appartient à E2
  u1_2 <- fu(survie_vec2(data$Y.non.proband, theta1, theta2, data$X1.non.proband, data$X2.non.proband))
  t22 <- gauss_cop_logpdf(u1_2, u2_2, h22 * data$kin)

  v1<-cbind(-((WU.2[[1]] * t12) + (WU.2[[2]] * t22)),data$ID)
  
  return(v1)
}

# =========================================================
# 4. LOG-VRAISEMBLANCE COMPLETE & WRAPPERS
# =========================================================

Negloglik.uni_vec <- function(params, Data_Uni, H, WU.1, WU.2) {
  
  h11 <- H[1]; h12 <- H[2]; h22 <- H[3]
  
  #val <- 0
  
  if (nrow(Data_Uni$p_E1) > 0) val1 <-   NegLogLik.proband.E1_vec(params, Data_Uni$p_E1)
  if (nrow(Data_Uni$p_E2) > 0) val2 <-  NegLogLik.proband.E2_vec(params, Data_Uni$p_E2)
  
  if (nrow(Data_Uni$d1_P1_NP1) > 0) val3 <-  NegLogLik.non.proband.delta1.E1.E1_vec(params, Data_Uni$d1_P1_NP1, h11)
  if (nrow(Data_Uni$d1_P2_NP1) > 0) val4 <-  NegLogLik.non.proband.delta1.E2.E1_vec(params, Data_Uni$d1_P2_NP1, h12)
  if (nrow(Data_Uni$d1_P1_NP2) > 0) val5 <-  NegLogLik.non.proband.delta1.E1.E2_vec(params, Data_Uni$d1_P1_NP2, h12)
  if (nrow(Data_Uni$d1_P2_NP2) > 0) val6 <- NegLogLik.non.proband.delta1.E2.E2_vec(params, Data_Uni$d1_P2_NP2, h22)
  
  if (nrow(Data_Uni$d0_P1) > 0) val7 <-NegLogLik.non.proband.delta0.E1_vec(params, Data_Uni$d0_P1, h11, h12, WU.1)
  if (nrow(Data_Uni$d0_P2) > 0) val8 <-NegLogLik.non.proband.delta0.E2_vec(params, Data_Uni$d0_P2, h12, h22, WU.2)
  
  val<-rbind(val1,val2,val3,val4,val5,val6,val7, val8)
  
  val_mean <- aggregate(val[, 1] ~ val[, 2], FUN = sum, na.rm=TRUE)
  #print(val)
  return(val_mean[, 2])
  
}


Negloglik.uni_full_vec <- function(params, Data_Uni, H,WU.1, WU.2) {
  

  # 3. Appel de la fonction de vraisemblance univariée avec les poids calculés
  val <- Negloglik.uni_vec(params=params, Data_Uni=Data_Uni, H=H, WU.1 = WU.1,WU.2 = WU.2)
  #print(val)
  return(val)
}



# Version robuste pour optim()
Negloglik.uni_safe <- function(params, Data_Uni, H,WU.1, WU.2) {
  val <- Negloglik.uni_full(params = params, Data_Uni=Data_Uni, H=H,WU.1=WU.1, WU.2=WU.2)
  if (!is.finite(val)) return(1e12)
  return(val)
}

# Gradient numérique via numDeriv
Negloglik.uni_grad <- function(params, Data_Uni, H,WU.1, WU.2) {
  numDeriv::grad(func = Negloglik.uni_safe, x = params, method = "Richardson",
                 method.args = list(eps = 1e-4, d = 0.0005, r = 6), Data_Uni=Data_Uni, H=H,WU.1=WU.1, WU.2=WU.2)
}

