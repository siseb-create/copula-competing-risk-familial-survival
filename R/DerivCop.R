
# =========================================================
# OUTILS DE STABILITÉ NUMÉRIQUE
# =========================================================
clamp01 <- function(x, eps = 1e-12) {
  pmin(pmax(unlist(x), eps), 1 - eps)
}

clip_rho <- function(r, eps = 1e-8) {
  pmin(pmax(unlist(r), -1 + eps), 1 - eps)
}

# =========================================================
# C011 = P(U1 <= u1 | U2 = u2, U3 = u3)
# ENTIÈREMENT VECTORISÉE (LOI NORMALE CONDITIONNELLE)
# =========================================================
log.C011 <- function(u1, u2, u3, r12, r13, r23, eps = 1e-12) {
  # Clipping et transformation Z
  u1  <- clamp01(u1, eps); u2 <- clamp01(u2, eps); u3 <- clamp01(u3, eps)
  r12 <- clip_rho(r12); r13 <- clip_rho(r13); r23 <- clip_rho(r23)
  
  Z1 <- qnorm(u1); Z2 <- qnorm(u2); Z3 <- qnorm(u3)
  
  # Déterminant partiel pour la conditionnelle
  D <- pmax(1 - r23^2, eps)
  
  # Moments de la distribution conditionnelle : Z1 | (Z2, Z3)
  mu  <- (Z2 * (r12 - r13 * r23) + Z3 * (r13 - r12 * r23)) / D
  var <- pmax(1 - (r12^2 + r13^2 - 2 * r12 * r13 * r23) / D, eps)
  sdv <- sqrt(var)
  
  pnorm(Z1, mean = mu, sd = sdv, log.p = TRUE)
}

C011 <- function(u1, u2, u3, r12, r13, r23, eps = 1e-12) {
  exp(log.C011(u1, u2, u3, r12, r13, r23, eps = eps))
}

# =========================================================
# C001 = P(U1 <= u1, U2 <= u2 | U3 = u3)
# OPTION A : pbivnorm (VECTORISATION SANS BOUCLE)
# =========================================================
C001.bis <- function(u1, u2, u3, r12, r13, r23, eps = 1e-12) {
  # Transformation et clipping
  u1  <- clamp01(u1, eps); u2 <- clamp01(u2, eps); u3 <- clamp01(u3, eps)
  r12 <- clip_rho(r12); r13 <- clip_rho(r13); r23 <- clip_rho(r23)
  
  Z1 <- qnorm(u1); Z2 <- qnorm(u2); Z3 <- qnorm(u3)
  
  # Moments conditionnels de la bivariée (Z1, Z2 | Z3)
  M1 <- Z3 * r13
  M2 <- Z3 * r23
  
  V11 <- pmax(1 - r13^2, eps)
  V22 <- pmax(1 - r23^2, eps)
  V12 <- r12 - r13 * r23
  
  # Normalisation pour pbivnorm (transformation en normale standard corrélée)
  sd1 <- sqrt(V11)
  sd2 <- sqrt(V22)
  
  rho_cond <- pmin(pmax(V12 / (sd1 * sd2), -1 + eps), 1 - eps)
  upper1   <- (Z1 - M1) / sd1
  upper2   <- (Z2 - M2) / sd2
  
  # Appel vectorisé à pbivnorm (C++ sous le capot)
  pbivnorm(upper1, upper2, rho = rho_cond)
}

# Version log sécurisée pour les poids
log.C001<- function(u1, u2, u3, r12, r13, r23) {
  val <- C001.bis(u1, u2, u3, r12, r13, r23)
  log(pmax(val, 1e-300))
}

# =========================================================
# log.C111 : DENSITÉ DE LA COPULE GAUSSIENNE 3D
# VERSION OPTIMISÉE (ALGÈBRE DIRECTE ET VECTORISÉE)
# =========================================================
log.C111 <- function(r12, r13, r23, u1, u2, u3, eps = 1e-12) {
  u1 <- clamp01(u1, eps); u2 <- clamp01(u2, eps); u3 <- clamp01(u3, eps)
  r12 <- clip_rho(r12); r13 <- clip_rho(r13); r23 <- clip_rho(r23)
  
  z1 <- qnorm(u1); z2 <- qnorm(u2); z3 <- qnorm(u3)
  
  # Déterminant 3x3 direct
  detS <- 1 - r12^2 - r13^2 - r23^2 + 2 * r12 * r13 * r23
  detS <- pmax(detS, eps)
  
  # Éléments de la matrice de précision (Inverse de Sigma)
  invDet <- 1 / detS
  a11_minus_1 <- ((1 - r23^2) * invDet) - 1
  a22_minus_1 <- ((1 - r13^2) * invDet) - 1
  a33_minus_1 <- ((1 - r12^2) * invDet) - 1
  
  a12 <- (r13 * r23 - r12) * invDet
  a13 <- (r12 * r23 - r13) * invDet
  a23 <- (r12 * r13 - r23) * invDet
  
  # Forme quadratique optimisée (z' (Sigma^-1 - I) z)
  quad <- a11_minus_1 * z1^2 + a22_minus_1 * z2^2 + a33_minus_1 * z3^2 + 
    2 * (a12 * z1 * z2 + a13 * z1 * z3 + a23 * z2 * z3)
  
  -0.5 * (log(detS) + quad)
}

C111 <- function(r12, r13, r23, u1, u2, u3, eps = 1e-12) {
  exp(log.C111(r12, r13, r23, u1, u2, u3, eps = eps))
}
