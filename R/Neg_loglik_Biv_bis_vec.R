############################################################
## Vectorized bivariate likelihood contributions for B-hat
##
## This file contains vectorized versions of the family-specific
## bivariate likelihood contributions. They are not used for
## parameter estimation directly, but only for computing the
## empirical covariance component B-hat in the robust asymptotic
## variance calculation.
############################################################

# =========================================================
# BLOC 111 : Trivarié complet (8 termes)
# =========================================================
Neg.loglik.Biv111_vec <- function(D, h11, h12, h22, b) {
  
  t1 <- if(nrow(D$E111)>0) (log.C111(r12=h11*D$E111$kin1J, r13=h11*D$E111$kin1K, r23=h11*D$E111$kinJK, u1=b[1], u2=b[2], u3=b[3])) else 0
  t2 <- if(nrow(D$E211)>0) (log.C111(r12=h12*D$E211$kin1J, r13=h12*D$E211$kin1K, r23=h11*D$E211$kinJK, u1=b[4], u2=b[5], u3=b[6])) else 0
  t3 <- if(nrow(D$E121)>0) (log.C111(r12=h12*D$E121$kin1J, r13=h11*D$E121$kin1K, r23=h12*D$E121$kinJK, u1=b[7], u2=b[8], u3=b[9])) else 0
  t4 <- if(nrow(D$E221)>0) (log.C111(r12=h22*D$E221$kin1J, r13=h12*D$E221$kin1K, r23=h12*D$E221$kinJK, u1=b[10],u2=b[11],u3=b[12])) else 0
  t5 <- if(nrow(D$E112)>0) (log.C111(r12=h11*D$E112$kin1J, r13=h12*D$E112$kin1K, r23=h12*D$E112$kinJK, u1=b[13],u2=b[14],u3=b[15])) else 0
  t6 <- if(nrow(D$E212)>0) (log.C111(r12=h12*D$E212$kin1J, r13=h22*D$E212$kin1K, r23=h12*D$E212$kinJK, u1=b[16],u2=b[17],u3=b[18])) else 0
  t7 <- if(nrow(D$E122)>0) (log.C111(r12=h12*D$E122$kin1J, r13=h12*D$E122$kin1K, r23=h22*D$E122$kinJK, u1=b[19],u2=b[20],u3=b[21])) else 0
  t8 <- if(nrow(D$E222)>0) (log.C111(r12=h22*D$E222$kin1J, r13=h22*D$E222$kin1K, r23=h22*D$E222$kinJK, u1=b[22],u2=b[23],u3=b[24])) else 0
  res<-rbind(
    cbind(t1,D$E111$ID),
    cbind(t2,D$E211$ID),
    cbind(t3,D$E121$ID),
    cbind(t4,D$E221$ID),
    cbind(t5,D$E112$ID),
    cbind(t6,D$E212$ID),
    cbind(t7,D$E122$ID),
    cbind(t8,D$E222$ID)
  )
  return(res)
}

# =========================================================
# BLOC 101 : Bivarié type 1 (4 sous-blocs)
# =========================================================
Neg.loglik.Biv101_vec <- function(D, WB, h11, h12, h22, b) {
  b1 <- if(nrow(D$E11)>0) (WB$B1[[1]]*log.C011(r12=h11*D$E11$kin1J, r13=h11*D$E11$kinJK, r23=h11*D$E11$kin1K, u1=b[1], u2=b[2], u3=b[3]) + 
                                WB$B1[[2]]*log.C011(r12=h12*D$E11$kin1J, r13=h12*D$E11$kinJK, r23=h11*D$E11$kin1K, u1=b[4], u2=b[5], u3=b[6])) else 0
  b2 <- if(nrow(D$E21)>0) (WB$B2[[1]]*log.C011(r12=h12*D$E21$kin1J, r13=h11*D$E21$kinJK, r23=h12*D$E21$kin1K, u1=b[7], u2=b[8], u3=b[9]) + 
                                WB$B2[[2]]*log.C011(r12=h22*D$E21$kin1J, r13=h12*D$E21$kinJK, r23=h12*D$E21$kin1K, u1=b[10],u2=b[11],u3=b[12])) else 0
  b3 <- if(nrow(D$E12)>0) (WB$B3[[1]]*log.C011(r12=h11*D$E12$kin1J, r13=h12*D$E12$kinJK, r23=h12*D$E12$kin1K, u1=b[13],u2=b[14],u3=b[15]) + 
                                WB$B3[[2]]*log.C011(r12=h12*D$E12$kin1J, r13=h22*D$E12$kinJK, r23=h12*D$E12$kin1K, u1=b[16],u2=b[17],u3=b[18])) else 0
  b4 <- if(nrow(D$E22)>0) (WB$B4[[1]]*log.C011(r12=h12*D$E22$kin1J, r13=h12*D$E22$kinJK, r23=h22*D$E22$kin1K, u1=b[19],u2=b[20],u3=b[21]) + 
                                WB$B4[[2]]*log.C011(r12=h22*D$E22$kin1J, r13=h22*D$E22$kinJK, r23=h22*D$E22$kin1K, u1=b[22],u2=b[23],u3=b[24])) else 0
  res<-rbind(
    cbind(b1,D$E11$ID),
    cbind(b2,D$E21$ID),
    cbind(b3,D$E12$ID),
    cbind(b4,D$E22$ID)
  )
  return(res)
}

# =========================================================
# BLOC 110 : Bivarié type 2 (4 sous-blocs)
# =========================================================
# =========================================================
Neg.loglik.Biv110_vec <- function(D, WB, h11, h12, h22, b) {
  b1 <- if(nrow(D$E11)>0) (WB$B5[[1]]*log.C011(r12=h11*D$E11$kin1K, r13=h11*D$E11$kinJK, r23=h11*D$E11$kin1J, u1=b[1], u2=b[2], u3=b[3]) + 
                                WB$B5[[2]]*log.C011(r12=h12*D$E11$kin1K, r13=h12*D$E11$kinJK, r23=h11*D$E11$kin1J, u1=b[4], u2=b[5], u3=b[6])) else 0
  b2 <- if(nrow(D$E21)>0) (WB$B6[[1]]*log.C011(r12=h12*D$E21$kin1K, r13=h11*D$E21$kinJK, r23=h12*D$E21$kin1J, u1=b[7], u2=b[8], u3=b[9]) + 
                                WB$B6[[2]]*log.C011(r12=h22*D$E21$kin1K, r13=h12*D$E21$kinJK, r23=h12*D$E21$kin1J, u1=b[10],u2=b[11],u3=b[12])) else 0
  b3 <- if(nrow(D$E12)>0) (WB$B7[[1]]*log.C011(r12=h11*D$E12$kin1K, r13=h12*D$E12$kinJK, r23=h12*D$E12$kin1J, u1=b[13],u2=b[14],u3=b[15]) + 
                                WB$B7[[2]]*log.C011(r12=h12*D$E12$kin1K, r13=h22*D$E12$kinJK, r23=h12*D$E12$kin1J, u1=b[16],u2=b[17],u3=b[18])) else 0
  b4 <- if(nrow(D$E22)>0) (WB$B8[[1]]*log.C011(r12=h12*D$E22$kin1K, r13=h12*D$E22$kinJK, r23=h22*D$E22$kin1J, u1=b[19],u2=b[20],u3=b[21]) + 
                                WB$B8[[2]]*log.C011(r12=h22*D$E22$kin1K, r13=h22*D$E22$kinJK, r23=h22*D$E22$kin1J, u1=b[22],u2=b[23],u3=b[24])) else 0
  res<-rbind(
    cbind(b1,D$E11$ID),
    cbind(b2,D$E21$ID),
    cbind(b3,D$E12$ID),
    cbind(b4,D$E22$ID)
  )
  return(res)
}

# =========================================================
# BLOC 100 : Univarié (2 sous-blocs avec mélanges de 4)
# =========================================================
Neg.loglik.Biv100_vec <- function(D, WB, h11, h12, h22, b) {
  
  b1 <- if(nrow(D$E1)>0) {
    s1 <- log.C001(r12=h11*D$E1$kin1J, r13=h11*D$E1$kinJK, r23=h11*D$E1$kin1K, u1=b[1], u2=b[3], u3=b[2])
    s2 <- log.C001(r12=h12*D$E1$kin1J, r13=h12*D$E1$kinJK, r23=h11*D$E1$kin1K, u1=b[4], u2=b[6], u3=b[5])
    s3 <- log.C001(r12=h12*D$E1$kin1J, r13=h11*D$E1$kinJK, r23=h12*D$E1$kin1K, u1=b[7], u2=b[9], u3=b[8])
    s4 <- log.C001(r12=h22*D$E1$kin1J, r13=h12*D$E1$kinJK, r23=h12*D$E1$kin1K, u1=b[10],u2=b[12],u3=b[11])
    (WB$B9[[1]]*s1 + WB$B9[[2]]*s2 + WB$B9[[3]]*s3 + WB$B9[[4]]*s4)
  } else 0
  
  b2 <- if(nrow(D$E2)>0) {
    s1 <- log.C001(r12=h11*D$E2$kin1J, r13=h12*D$E2$kinJK, r23=h12*D$E2$kin1K, u1=b[13],u2=b[15],u3=b[14])
    s2 <- log.C001(r12=h12*D$E2$kin1J, r13=h22*D$E2$kinJK, r23=h12*D$E2$kin1K, u1=b[16],u2=b[18],u3=b[17])
    s3 <- log.C001(r12=h12*D$E2$kin1J, r13=h12*D$E2$kinJK, r23=h22*D$E2$kin1K, u1=b[19],u2=b[21],u3=b[20])
    s4 <- log.C001(r12=h22*D$E2$kin1J, r13=h22*D$E2$kinJK, r23=h22*D$E2$kin1K, u1=b[22],u2=b[24],u3=b[23])
    (WB$B10[[1]]*s1 + WB$B10[[2]]*s2 + WB$B10[[3]]*s3 + WB$B10[[4]]*s4)
  } else 0
  res<-rbind(
    cbind(b1,D$E1$ID),
    cbind(b2,D$E2$ID))
  return(res)
}


# =========================================================
# VRAISEMBLANCE NÉGATIVE GLOBALE
# =========================================================
Negloglik.Biv_vec <- function(H, D, B_surv, WB) {
  
  #print(H)
  
  h11 <- H[1]; h12 <- H[2]; h22 <- H[3]
  
  # Somme des vraisemblances par bloc
  L111 <- Neg.loglik.Biv111_vec(D=D$B111, h11=h11, h12=h12, h22=h22, B_surv[1:24])
  L101 <- Neg.loglik.Biv101_vec(D=D$B101, WB=WB, h11=h11, h12=h12, h22=h22, B_surv[25:48])
  L110 <- Neg.loglik.Biv110_vec(D=D$B110, WB=WB, h11=h11, h12=h12, h22=h22, B_surv[49:72])
  L100 <- Neg.loglik.Biv100_vec(D=D$B100, WB=WB, h11=h11, h12=h12, h22=h22, B_surv[73:96])
  
  # Inversion du signe pour minimisation (NLL)
  val <- rbind(L111, L101, L110, L100)
  
  val_mean <- aggregate(val[, 1] ~ val[, 2], FUN = sum, na.rm=TRUE)
  #print(val)
  return(val_mean[, 2])
}



