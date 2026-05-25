
# =========================================================
# 0. OUTIL
# =========================================================
clip01 <- function(v, eps = 1e-12) {
  pmin(pmax(v, eps), 1 - eps)
}

if (!requireNamespace("statmod", quietly = TRUE)) {
  stop("Le package 'statmod' doit être installé.")
}

# =========================================================
# 1. QUADRATURES PRÉ-CALCULÉES
# =========================================================
# Gauss-Legendre sur [0, t]
GL <- statmod::gauss.quad(64, kind = "legendre")
GL_nodes <- 0.5 * (GL$nodes + 1)
GL_wts   <- 0.5 * GL$weights

# Gauss-Laguerre sur [0, +Inf)
GL_inf <- statmod::gauss.quad(64, kind = "laguerre")
GL_inf_nodes <- GL_inf$nodes
GL_inf_wts   <- GL_inf$weights

# =========================================================
# 2. EXTRACTION DES PARAMÈTRES WEIBULL
# =========================================================
.par_weibull <- function(theta, X1, X2) {
  alpha  <- exp(theta[1])
  lambda <- exp(theta[2])
  beta1  <- theta[3]
  beta2  <- theta[4]
  z <- lambda * exp(beta1 * X1 + beta2 * X2)
  list(alpha = alpha, lambda = lambda, beta1 = beta1, beta2 = beta2, z = z)
}

# =========================================================
# 3. FONCTIONS DE BASE
# =========================================================
hasard <- function(t, theta, X1, X2) {
  p <- .par_weibull(theta, X1, X2)
  p$alpha * p$z * pmax(t, 1e-12)^(p$alpha - 1)
}

hasard.cum <- function(t, theta, X1, X2) {
  p <- .par_weibull(theta, X1, X2)
  p$z * pmax(t, 0)^p$alpha
}

survie.T <- function(t, theta1, theta2, X1, X2) {
  exp(-(hasard.cum(t, theta1, X1, X2) + hasard.cum(t, theta2, X1, X2)))
}

CIF.a.integrer <- function(t, gam1, theta1, theta2, X1, X2) {
  hasard(t, gam1, X1, X2) * survie.T(t, theta1, theta2, X1, X2)
}

# =========================================================
# 4. CIF : F_l(t|X)
# =========================================================
CIF <- function(t, gam1, theta1, theta2, X1, X2) {
  n <- max(length(t), length(X1), length(X2))
  t  <- rep_len(t,  n)
  X1 <- rep_len(X1, n)
  X2 <- rep_len(X2, n)
  
  out <- numeric(n)
  id <- is.finite(t) & (t > 0)
  
  if (!any(id)) return(out)
  
  tt <- t[id]
  x1 <- X1[id]
  x2 <- X2[id]
  
  pg <- .par_weibull(gam1,   x1, x2)
  p1 <- .par_weibull(theta1, x1, x2)
  p2 <- .par_weibull(theta2, x1, x2)
  
  # Transformation u = t * s, s in [0,1]
  U <- outer(tt, GL_nodes)
  
  zg <- matrix(pg$z, nrow = length(tt), ncol = length(GL_nodes))
  z1 <- matrix(p1$z, nrow = length(tt), ncol = length(GL_nodes))
  z2 <- matrix(p2$z, nrow = length(tt), ncol = length(GL_nodes))
  
  logFU <- log(pg$alpha) +
    log(zg) +
    (pg$alpha - 1) * log(pmax(U, 1e-12)) -
    z1 * U^p1$alpha -
    z2 * U^p2$alpha
  
  FU <- exp(logFU)
  out[id] <- tt * drop(FU %*% GL_wts)
  
  pmax(out, 0)
}

# =========================================================
# 5. F_l(inf|X)
# =========================================================
# Proba1(theta1, theta2, ...) = P(epsilon = 1 | X)
Proba1 <- function(theta1, theta2, X1, X2) {
  n  <- max(length(X1), length(X2))
  X1 <- rep_len(X1, n)
  X2 <- rep_len(X2, n)
  
  p1 <- .par_weibull(theta1, X1, X2)
  p2 <- .par_weibull(theta2, X1, X2)
  
  m <- length(GL_inf_nodes)
  U <- matrix(GL_inf_nodes, nrow = n, ncol = m, byrow = TRUE)
  
  z1m <- matrix(p1$z, nrow = n, ncol = m)
  z2m <- matrix(p2$z, nrow = n, ncol = m)
  
  # intégrande brute :
  # alpha1*z1*u^(alpha1-1)*exp(-z1*u^a1 - z2*u^a2)
  # en Laguerre, on multiplie par exp(u)
  logFU <- log(p1$alpha) +
    log(z1m) +
    (p1$alpha - 1) * log(pmax(U, 1e-12)) -
    z1m * U^p1$alpha -
    z2m * U^p2$alpha +
    U
  
  FU <- exp(logFU)
  out <- drop(FU %*% GL_inf_wts)
  
  pmin(pmax(out, 0), 1)
}

# =========================================================
# 6. P(epsilon=1 | T<=a, X)
# =========================================================
PProba <- function(a, theta1, theta2, X1, X2) {
  p1 <- CIF(a, theta1, theta1, theta2, X1, X2)
  p2 <- CIF(a, theta2, theta1, theta2, X1, X2)
  clip01(p1 / pmax(p1 + p2, 1e-12))
}

# =========================================================
# 7. SURVIES CONDITIONNELLES PAR CAUSE
# =========================================================
# proband : conditionnement sur T<=a et epsilon=l
survie.cond.proband <- function(t, a, gam1, theta1, theta2, X1, X2, res = 0) {
  A <- CIF(t, gam1, theta1, theta2, X1, X2)
  B <- CIF(a, gam1, theta1, theta2, X1, X2)
  1 - A / pmax(B, 1e-12) - res
}

# conditionnement sur epsilon=l
survie.cond <- function(t, gam1, theta1, theta2, X1, X2, res = 0) {
  A <- CIF(t, gam1, theta1, theta2, X1, X2)
  
  if (isTRUE(all.equal(gam1, theta1, tolerance = 0))) {
    B <- Proba1(theta1, theta2, X1, X2)
  } else if (isTRUE(all.equal(gam1, theta2, tolerance = 0))) {
    B <- Proba1(theta2, theta1, X1, X2)
  } else {
    stop("gam1 doit être égal à theta1 ou theta2.")
  }
  
  1 - A / pmax(B, 1e-12) - res
}

# =========================================================
# 8. SURVIES CONDITIONNELLES PAR CAUSE VECTEURS
# P(T>t | epsilon=l, X)
# =========================================================
survie_vec1 <- function(t, theta1, theta2, X1, X2) {
  num <- pmax(
    F1_t_sup(t, theta1, theta2, X1, X2),
    0
  )
  den <- pmax(
    F1_inf(theta1, theta2, X1, X2),
    1e-12
  )
  num / den
}

survie_vec2 <- function(t, theta1, theta2, X1, X2) {
  num <- pmax(
    F2_t_sup(t, theta1, theta2, X1, X2),
    0
  )
  den <- pmax(
    F2_inf(theta1, theta2, X1, X2),
    1e-12
  )
  num / den
}

# =========================================================
# 9. MATRICE DE DÉPENDANCE
# =========================================================
Calculer.V <- function(epsilon, kinship, h11, h12, h22) {
  epsilon <- as.numeric(epsilon)
  
  e1 <- 2 - epsilon
  e2 <- epsilon - 1
  
  V11 <- h11 * kinship * tcrossprod(e1)
  V22 <- h22 * kinship * tcrossprod(e2)
  V12 <- h12 * kinship * (e1 %o% e2)
  V21 <- h12 * kinship * (e2 %o% e1)
  
  corr <- diag(1 - h11 * e1 - h22 * e2)
  
  V <- corr + V11 + V22 + V12 + V21
  V <- (V + t(V)) / 2
  V
}

# =========================================================
# 10. FONCTIONS ADDITIONNELLES
# =========================================================
densite1 <- function(t, theta1, theta2, X1, X2) {
  hasard(t, theta = theta1, X1, X2) * survie.T(t, theta1, theta2, X1, X2)
}

densite2 <- function(t, theta1, theta2, X1, X2) {
  hasard(t, theta = theta2, X1, X2) * survie.T(t, theta1, theta2, X1, X2)
}

l.densite1 <- function(t, theta1, theta2, X1, X2) {
  log(pmax(hasard(t = t, theta = theta1, X1 = X1, X2 = X2), 1e-300)) +
    log(pmax(survie.T(t, theta1, theta2, X1, X2), 1e-300))
}

l.densite2 <- function(t, theta1, theta2, X1, X2) {
  log(pmax(hasard(t, theta = theta2, X1 = X1, X2 = X2), 1e-300)) +
    log(pmax(survie.T(t, theta1, theta2, X1, X2), 1e-300))
}

# =========================================================
# 11. F_l(inf|X)
# =========================================================
F1_inf <- function(theta1, theta2, X1, X2) {
  Proba1(theta1, theta2, X1, X2)
}

l.F1_inf <- function(theta1, theta2, X1, X2) {
  log(pmax(F1_inf(theta1, theta2, X1, X2), 1e-300))
}

l.F1_inf_vec <- function(theta1, theta2, X1, X2) {
  log(pmax(F1_inf(theta1, theta2, X1, X2), 1e-300))
}

F2_inf <- function(theta1, theta2, X1, X2) {
  Proba1(theta2, theta1, X1, X2)
}

l.F2_inf <- function(theta1, theta2, X1, X2) {
  log(pmax(F2_inf(theta1, theta2, X1, X2), 1e-300))
}

l.F2_inf_vec <- function(theta1, theta2, X1, X2) {
  log(pmax(F2_inf(theta1, theta2, X1, X2), 1e-300))
}

# =========================================================
# 12. F_l(t,inf|X) NON NORMALISÉES
# F_l(inf|X) - F_l(t|X) = P(T>t, epsilon=l | X)
# =========================================================
F1_t_sup <- function(t, theta1, theta2, X1, X2) {
  pmax(
    F1_inf(theta1, theta2, X1, X2) -
      CIF(t, theta1, theta1, theta2, X1, X2),
    0
  )
}

F2_t_sup <- function(t, theta1, theta2, X1, X2) {
  pmax(
    F2_inf(theta1, theta2, X1, X2) -
      CIF(t, theta2, theta1, theta2, X1, X2),
    0
  )
}