
# =========================================================
# 2. PRÉPARATION DES BASES DE DONNÉES (DATA MANAGEMENT)
# =========================================================
# Cette étape prépare les 8 segments nécessaires à la vraisemblance
prepare_uni_datalist <- function(data, indices.proband, indices.non.proband, kin, I) {
  
  # Probands
  data.proband <- create.data.proband(data, indices.proband)
  names(data.proband)[names(data.proband) == "famID"] <- "ID"
  
  # Non-Probands
  data.non.proband <- create.data.non.proband(data, indices.non.proband, kin, I)
  names(data.non.proband)[names(data.non.proband) == "famID"] <- "ID"
  
  list(
    p_E1 = data.proband[data.proband$epsilon == 1, ],
    p_E2 = data.proband[data.proband$epsilon == 2, ],
    
    d1_P1_NP1 = data.non.proband[data.non.proband$delta.non.proband == 1 & 
                                   data.non.proband$Epsilon.proband == 1 & 
                                   data.non.proband$Epsilon.non.proband == 1, -6],
    
    d1_P2_NP1 = data.non.proband[data.non.proband$delta.non.proband == 1 & 
                                   data.non.proband$Epsilon.proband == 2 & 
                                   data.non.proband$Epsilon.non.proband == 1, -6],
    
    d1_P1_NP2 = data.non.proband[data.non.proband$delta.non.proband == 1 & 
                                   data.non.proband$Epsilon.proband == 1 & 
                                   data.non.proband$Epsilon.non.proband == 2, -6],
    
    d1_P2_NP2 = data.non.proband[data.non.proband$delta.non.proband == 1 & 
                                   data.non.proband$Epsilon.proband == 2 & 
                                   data.non.proband$Epsilon.non.proband == 2, -6],
    
    d0_P1 = data.non.proband[data.non.proband$delta.non.proband == 0 & 
                               data.non.proband$Epsilon.proband == 1, -6],
    
    d0_P2 = data.non.proband[data.non.proband$delta.non.proband == 0 & 
                               data.non.proband$Epsilon.proband == 2, -6]
  )
}


# Fonction de regroupement de la base bivariée

data_bivariate_list<-function(Data_Biv){
  # --- Création de la structure de données regroupée ---
  DataList <- list(
    B111 = list(
      E111 = Data_Biv[Data_Biv$delta.non.probandIJ==1 & Data_Biv$delta.non.probandIK==1 & Data_Biv$Epsilon.proband==1 & Data_Biv$Epsilon.non.probandIJ==1 & Data_Biv$Epsilon.non.probandIK==1, ],
      E211 = Data_Biv[Data_Biv$delta.non.probandIJ==1 & Data_Biv$delta.non.probandIK==1 & Data_Biv$Epsilon.proband==2 & Data_Biv$Epsilon.non.probandIJ==1 & Data_Biv$Epsilon.non.probandIK==1, ],
      E121 = Data_Biv[Data_Biv$delta.non.probandIJ==1 & Data_Biv$delta.non.probandIK==1 & Data_Biv$Epsilon.proband==1 & Data_Biv$Epsilon.non.probandIJ==2 & Data_Biv$Epsilon.non.probandIK==1, ],
      E221 = Data_Biv[Data_Biv$delta.non.probandIJ==1 & Data_Biv$delta.non.probandIK==1 & Data_Biv$Epsilon.proband==2 & Data_Biv$Epsilon.non.probandIJ==2 & Data_Biv$Epsilon.non.probandIK==1, ],
      E112 = Data_Biv[Data_Biv$delta.non.probandIJ==1 & Data_Biv$delta.non.probandIK==1 & Data_Biv$Epsilon.proband==1 & Data_Biv$Epsilon.non.probandIJ==1 & Data_Biv$Epsilon.non.probandIK==2, ],
      E212 = Data_Biv[Data_Biv$delta.non.probandIJ==1 & Data_Biv$delta.non.probandIK==1 & Data_Biv$Epsilon.proband==2 & Data_Biv$Epsilon.non.probandIJ==1 & Data_Biv$Epsilon.non.probandIK==2, ],
      E122 = Data_Biv[Data_Biv$delta.non.probandIJ==1 & Data_Biv$delta.non.probandIK==1 & Data_Biv$Epsilon.proband==1 & Data_Biv$Epsilon.non.probandIJ==2 & Data_Biv$Epsilon.non.probandIK==2, ],
      E222 = Data_Biv[Data_Biv$delta.non.probandIJ==1 & Data_Biv$delta.non.probandIK==1 & Data_Biv$Epsilon.proband==2 & Data_Biv$Epsilon.non.probandIJ==2 & Data_Biv$Epsilon.non.probandIK==2, ]
    ),
    B101 = list(
      E11 = Data_Biv[Data_Biv$delta.non.probandIJ==0 & Data_Biv$delta.non.probandIK==1 & Data_Biv$Epsilon.proband==1 & Data_Biv$Epsilon.non.probandIK==1, ],
      E21 = Data_Biv[Data_Biv$delta.non.probandIJ==0 & Data_Biv$delta.non.probandIK==1 & Data_Biv$Epsilon.proband==2 & Data_Biv$Epsilon.non.probandIK==1, ],
      E12 = Data_Biv[Data_Biv$delta.non.probandIJ==0 & Data_Biv$delta.non.probandIK==1 & Data_Biv$Epsilon.proband==1 & Data_Biv$Epsilon.non.probandIK==2, ],
      E22 = Data_Biv[Data_Biv$delta.non.probandIJ==0 & Data_Biv$delta.non.probandIK==1 & Data_Biv$Epsilon.proband==2 & Data_Biv$Epsilon.non.probandIK==2, ]
    ),
    B110 = list(
      E11 = Data_Biv[Data_Biv$delta.non.probandIK==0 & Data_Biv$delta.non.probandIJ==1 & Data_Biv$Epsilon.proband==1 & Data_Biv$Epsilon.non.probandIJ==1, ],
      E21 = Data_Biv[Data_Biv$delta.non.probandIK==0 & Data_Biv$delta.non.probandIJ==1 & Data_Biv$Epsilon.proband==2 & Data_Biv$Epsilon.non.probandIJ==1, ],
      E12 = Data_Biv[Data_Biv$delta.non.probandIK==0 & Data_Biv$delta.non.probandIJ==1 & Data_Biv$Epsilon.proband==1 & Data_Biv$Epsilon.non.probandIJ==2, ],
      E22 = Data_Biv[Data_Biv$delta.non.probandIK==0 & Data_Biv$delta.non.probandIJ==1 & Data_Biv$Epsilon.proband==2 & Data_Biv$Epsilon.non.probandIJ==2, ]
    ),
    B100 = list(
      E1 = Data_Biv[Data_Biv$delta.non.probandIJ==0 & Data_Biv$delta.non.probandIK==0 & Data_Biv$Epsilon.proband==1, ],
      E2 = Data_Biv[Data_Biv$delta.non.probandIJ==0 & Data_Biv$delta.non.probandIK==0 & Data_Biv$Epsilon.proband==2, ]
    )
  )
}


# =========================================================
# GÉNÉRATION DU VECTEUR B_SURV (96 ÉLÉMENTS)
# =========================================================

generer_B_surv <- function(D, t1, t2) {
  
  # --- BLOC 111 (b1 à b24) ---
  b1  <- survie_vec1(D$B111$E111$Y.proband, t1, t2, D$B111$E111$X1.proband, D$B111$E111$X2.proband)
  b2  <- survie_vec1(D$B111$E111$Y.non.probandIJ, t1, t2, D$B111$E111$X1.non.probandIJ, D$B111$E111$X2.non.probandIJ)
  b3  <- survie_vec1(D$B111$E111$Y.non.probandIK, t1, t2, D$B111$E111$X1.non.probandIK, D$B111$E111$X2.non.probandIK)
  
  b4  <- survie_vec2(D$B111$E211$Y.proband, t1, t2, D$B111$E211$X1.proband, D$B111$E211$X2.proband)
  b5  <- survie_vec1(D$B111$E211$Y.non.probandIJ, t1, t2, D$B111$E211$X1.non.probandIJ, D$B111$E211$X2.non.probandIJ)
  b6  <- survie_vec1(D$B111$E211$Y.non.probandIK, t1, t2, D$B111$E211$X1.non.probandIK, D$B111$E211$X2.non.probandIK)
  
  b7  <- survie_vec1(D$B111$E121$Y.proband, t1, t2, D$B111$E121$X1.proband, D$B111$E121$X2.proband)
  b8  <- survie_vec2(D$B111$E121$Y.non.probandIJ, t1, t2, D$B111$E121$X1.non.probandIJ, D$B111$E121$X2.non.probandIJ)
  b9  <- survie_vec1(D$B111$E121$Y.non.probandIK, t1, t2, D$B111$E121$X1.non.probandIK, D$B111$E121$X2.non.probandIK)
  
  b10 <- survie_vec2(D$B111$E221$Y.proband, t1, t2, D$B111$E221$X1.proband, D$B111$E221$X2.proband)
  b11 <- survie_vec2(D$B111$E221$Y.non.probandIJ, t1, t2, D$B111$E221$X1.non.probandIJ, D$B111$E221$X2.non.probandIJ)
  b12 <- survie_vec1(D$B111$E221$Y.non.probandIK, t1, t2, D$B111$E221$X1.non.probandIK, D$B111$E221$X2.non.probandIK)
  
  b13 <- survie_vec1(D$B111$E112$Y.proband, t1, t2, D$B111$E112$X1.proband, D$B111$E112$X2.proband)
  b14 <- survie_vec1(D$B111$E112$Y.non.probandIJ, t1, t2, D$B111$E112$X1.non.probandIJ, D$B111$E112$X2.non.probandIJ)
  b15 <- survie_vec2(D$B111$E112$Y.non.probandIK, t1, t2, D$B111$E112$X1.non.probandIK, D$B111$E112$X2.non.probandIK)
  
  b16 <- survie_vec2(D$B111$E212$Y.proband, t1, t2, D$B111$E212$X1.proband, D$B111$E212$X2.proband)
  b17 <- survie_vec1(D$B111$E212$Y.non.probandIJ, t1, t2, D$B111$E212$X1.non.probandIJ, D$B111$E212$X2.non.probandIJ)
  b18 <- survie_vec2(D$B111$E212$Y.non.probandIK, t1, t2, D$B111$E212$X1.non.probandIK, D$B111$E212$X2.non.probandIK)
  
  b19 <- survie_vec1(D$B111$E122$Y.proband, t1, t2, D$B111$E122$X1.proband, D$B111$E122$X2.proband)
  b20 <- survie_vec2(D$B111$E122$Y.non.probandIJ, t1, t2, D$B111$E122$X1.non.probandIJ, D$B111$E122$X2.non.probandIJ)
  b21 <- survie_vec2(D$B111$E122$Y.non.probandIK, t1, t2, D$B111$E122$X1.non.probandIK, D$B111$E122$X2.non.probandIK)
  
  b22 <- survie_vec2(D$B111$E222$Y.proband, t1, t2, D$B111$E222$X1.proband, D$B111$E222$X2.proband)
  b23 <- survie_vec2(D$B111$E222$Y.non.probandIJ, t1, t2, D$B111$E222$X1.non.probandIJ, D$B111$E222$X2.non.probandIJ)
  b24 <- survie_vec2(D$B111$E222$Y.non.probandIK, t1, t2, D$B111$E222$X1.non.probandIK, D$B111$E222$X2.non.probandIK)
  
  # --- BLOC 101 (b25 à b48) ---
  b25 <- survie_vec1(D$B101$E11$Y.non.probandIJ, t1, t2, D$B101$E11$X1.non.probandIJ, D$B101$E11$X2.non.probandIJ)
  b26 <- survie_vec1(D$B101$E11$Y.proband, t1, t2, D$B101$E11$X1.proband, D$B101$E11$X2.proband)
  b27 <- survie_vec1(D$B101$E11$Y.non.probandIK, t1, t2, D$B101$E11$X1.non.probandIK, D$B101$E11$X2.non.probandIK)
  
  b28 <- survie_vec2(D$B101$E11$Y.non.probandIJ, t1, t2, D$B101$E11$X1.non.probandIJ, D$B101$E11$X2.non.probandIJ)
  b29 <- survie_vec1(D$B101$E11$Y.proband, t1, t2, D$B101$E11$X1.proband, D$B101$E11$X2.proband)
  b30 <- survie_vec1(D$B101$E11$Y.non.probandIK, t1, t2, D$B101$E11$X1.non.probandIK, D$B101$E11$X2.non.probandIK)
  
  b31 <- survie_vec1(D$B101$E21$Y.non.probandIJ, t1, t2, D$B101$E21$X1.non.probandIJ, D$B101$E21$X2.non.probandIJ)
  b32 <- survie_vec2(D$B101$E21$Y.proband, t1, t2, D$B101$E21$X1.proband, D$B101$E21$X2.proband)
  b33 <- survie_vec1(D$B101$E21$Y.non.probandIK, t1, t2, D$B101$E21$X1.non.probandIK, D$B101$E21$X2.non.probandIK)
  
  b34 <- survie_vec2(D$B101$E21$Y.non.probandIJ, t1, t2, D$B101$E21$X1.non.probandIJ, D$B101$E21$X2.non.probandIJ)
  b35 <- survie_vec2(D$B101$E21$Y.proband, t1, t2, D$B101$E21$X1.proband, D$B101$E21$X2.proband)
  b36 <- survie_vec1(D$B101$E21$Y.non.probandIK, t1, t2, D$B101$E21$X1.non.probandIK, D$B101$E21$X2.non.probandIK)
  
  b37 <- survie_vec1(D$B101$E12$Y.non.probandIJ, t1, t2, D$B101$E12$X1.non.probandIJ, D$B101$E12$X2.non.probandIJ)
  b38 <- survie_vec1(D$B101$E12$Y.proband, t1, t2, D$B101$E12$X1.proband, D$B101$E12$X2.proband)
  b39 <- survie_vec2(D$B101$E12$Y.non.probandIK, t1, t2, D$B101$E12$X1.non.probandIK, D$B101$E12$X2.non.probandIK)
  
  b40 <- survie_vec2(D$B101$E12$Y.non.probandIJ, t1, t2, D$B101$E12$X1.non.probandIJ, D$B101$E12$X2.non.probandIJ)
  b41 <- survie_vec1(D$B101$E12$Y.proband, t1, t2, D$B101$E12$X1.proband, D$B101$E12$X2.proband)
  b42 <- survie_vec2(D$B101$E12$Y.non.probandIK, t1, t2, D$B101$E12$X1.non.probandIK, D$B101$E12$X2.non.probandIK)
  
  b43 <- survie_vec1(D$B101$E22$Y.non.probandIJ, t1, t2, D$B101$E22$X1.non.probandIJ, D$B101$E22$X2.non.probandIJ)
  b44 <- survie_vec2(D$B101$E22$Y.proband, t1, t2, D$B101$E22$X1.proband, D$B101$E22$X2.proband)
  b45 <- survie_vec2(D$B101$E22$Y.non.probandIK, t1, t2, D$B101$E22$X1.non.probandIK, D$B101$E22$X2.non.probandIK)
  
  b46 <- survie_vec2(D$B101$E22$Y.non.probandIJ, t1, t2, D$B101$E22$X1.non.probandIJ, D$B101$E22$X2.non.probandIJ)
  b47 <- survie_vec2(D$B101$E22$Y.proband, t1, t2, D$B101$E22$X1.proband, D$B101$E22$X2.proband)
  b48 <- survie_vec2(D$B101$E22$Y.non.probandIK, t1, t2, D$B101$E22$X1.non.probandIK, D$B101$E22$X2.non.probandIK)
  
  # --- BLOC 110 (b49 à b72) ---
  b49 <- survie_vec1(D$B110$E11$Y.non.probandIK, t1, t2, D$B110$E11$X1.non.probandIK, D$B110$E11$X2.non.probandIK)
  b50 <- survie_vec1(D$B110$E11$Y.proband, t1, t2, D$B110$E11$X1.proband, D$B110$E11$X2.proband)
  b51 <- survie_vec1(D$B110$E11$Y.non.probandIJ, t1, t2, D$B110$E11$X1.non.probandIJ, D$B110$E11$X2.non.probandIJ)
  
  b52 <- survie_vec2(D$B110$E11$Y.non.probandIK, t1, t2, D$B110$E11$X1.non.probandIK, D$B110$E11$X2.non.probandIK)
  b53 <- survie_vec1(D$B110$E11$Y.proband, t1, t2, D$B110$E11$X1.proband, D$B110$E11$X2.proband)
  b54 <- survie_vec1(D$B110$E11$Y.non.probandIJ, t1, t2, D$B110$E11$X1.non.probandIJ, D$B110$E11$X2.non.probandIJ)
  
  b55 <- survie_vec1(D$B110$E21$Y.non.probandIK, t1, t2, D$B110$E21$X1.non.probandIK, D$B110$E21$X2.non.probandIK)
  b56 <- survie_vec2(D$B110$E21$Y.proband, t1, t2, D$B110$E21$X1.proband, D$B110$E21$X2.proband)
  b57 <- survie_vec1(D$B110$E21$Y.non.probandIJ, t1, t2, D$B110$E21$X1.non.probandIJ, D$B110$E21$X2.non.probandIJ)
  
  b58 <- survie_vec2(D$B110$E21$Y.non.probandIK, t1, t2, D$B110$E21$X1.non.probandIK, D$B110$E21$X2.non.probandIK)
  b59 <- survie_vec2(D$B110$E21$Y.proband, t1, t2, D$B110$E21$X1.proband, D$B110$E21$X2.proband)
  b60 <- survie_vec1(D$B110$E21$Y.non.probandIJ, t1, t2, D$B110$E21$X1.non.probandIJ, D$B110$E21$X2.non.probandIJ)
  
  b61 <- survie_vec1(D$B110$E12$Y.non.probandIK, t1, t2, D$B110$E12$X1.non.probandIK, D$B110$E12$X2.non.probandIK)
  b62 <- survie_vec1(D$B110$E12$Y.proband, t1, t2, D$B110$E12$X1.proband, D$B110$E12$X2.proband)
  b63 <- survie_vec2(D$B110$E12$Y.non.probandIJ, t1, t2, D$B110$E12$X1.non.probandIJ, D$B110$E12$X2.non.probandIJ)
  
  b64 <- survie_vec2(D$B110$E12$Y.non.probandIK, t1, t2, D$B110$E12$X1.non.probandIK, D$B110$E12$X2.non.probandIK)
  b65 <- survie_vec1(D$B110$E12$Y.proband, t1, t2, D$B110$E12$X1.proband, D$B110$E12$X2.proband)
  b66 <- survie_vec2(D$B110$E12$Y.non.probandIJ, t1, t2, D$B110$E12$X1.non.probandIJ, D$B110$E12$X2.non.probandIJ)
  
  b67 <- survie_vec1(D$B110$E22$Y.non.probandIK, t1, t2, D$B110$E22$X1.non.probandIK, D$B110$E22$X2.non.probandIK)
  b68 <- survie_vec2(D$B110$E22$Y.proband, t1, t2, D$B110$E22$X1.proband, D$B110$E22$X2.proband)
  b69 <- survie_vec2(D$B110$E22$Y.non.probandIJ, t1, t2, D$B110$E22$X1.non.probandIJ, D$B110$E22$X2.non.probandIJ)
  
  b70 <- survie_vec2(D$B110$E22$Y.non.probandIK, t1, t2, D$B110$E22$X1.non.probandIK, D$B110$E22$X2.non.probandIK)
  b71 <- survie_vec2(D$B110$E22$Y.proband, t1, t2, D$B110$E22$X1.proband, D$B110$E22$X2.proband)
  b72 <- survie_vec2(D$B110$E22$Y.non.probandIJ, t1, t2, D$B110$E22$X1.non.probandIJ, D$B110$E22$X2.non.probandIJ)
  
  # --- BLOC 100 (b73 à b96) ---
  b73 <- survie_vec1(D$B100$E1$Y.non.probandIJ, t1, t2, D$B100$E1$X1.non.probandIJ, D$B100$E1$X2.non.probandIJ)
  b74 <- survie_vec1(D$B100$E1$Y.non.probandIK, t1, t2, D$B100$E1$X1.non.probandIK, D$B100$E1$X2.non.probandIK)
  b75 <- survie_vec1(D$B100$E1$Y.proband, t1, t2, D$B100$E1$X1.proband, D$B100$E1$X2.proband)
  
  b76 <- survie_vec2(D$B100$E1$Y.non.probandIJ, t1, t2, D$B100$E1$X1.non.probandIJ, D$B100$E1$X2.non.probandIJ)
  b77 <- survie_vec1(D$B100$E1$Y.non.probandIK, t1, t2, D$B100$E1$X1.non.probandIK, D$B100$E1$X2.non.probandIK)
  b78 <- survie_vec1(D$B100$E1$Y.proband, t1, t2, D$B100$E1$X1.proband, D$B100$E1$X2.proband)
  
  b79 <- survie_vec1(D$B100$E1$Y.non.probandIJ, t1, t2, D$B100$E1$X1.non.probandIJ, D$B100$E1$X2.non.probandIJ)
  b80 <- survie_vec2(D$B100$E1$Y.non.probandIK, t1, t2, D$B100$E1$X1.non.probandIK, D$B100$E1$X2.non.probandIK)
  b81 <- survie_vec1(D$B100$E1$Y.proband, t1, t2, D$B100$E1$X1.proband, D$B100$E1$X2.proband)
  
  b82 <- survie_vec2(D$B100$E1$Y.non.probandIJ, t1, t2, D$B100$E1$X1.non.probandIJ, D$B100$E1$X2.non.probandIJ)
  b83 <- survie_vec2(D$B100$E1$Y.non.probandIK, t1, t2, D$B100$E1$X1.non.probandIK, D$B100$E1$X2.non.probandIK)
  b84 <- survie_vec1(D$B100$E1$Y.proband, t1, t2, D$B100$E1$X1.proband, D$B100$E1$X2.proband)
  
  b85 <- survie_vec1(D$B100$E2$Y.non.probandIJ, t1, t2, D$B100$E2$X1.non.probandIJ, D$B100$E2$X2.non.probandIJ)
  b86 <- survie_vec1(D$B100$E2$Y.non.probandIK, t1, t2, D$B100$E2$X1.non.probandIK, D$B100$E2$X2.non.probandIK)
  b87 <- survie_vec2(D$B100$E2$Y.proband, t1, t2, D$B100$E2$X1.proband, D$B100$E2$X2.proband)
  
  b88 <- survie_vec2(D$B100$E2$Y.non.probandIJ, t1, t2, D$B100$E2$X1.non.probandIJ, D$B100$E2$X2.non.probandIJ)
  b89 <- survie_vec1(D$B100$E2$Y.non.probandIK, t1, t2, D$B100$E2$X1.non.probandIK, D$B100$E2$X2.non.probandIK)
  b90 <- survie_vec2(D$B100$E2$Y.proband, t1, t2, D$B100$E2$X1.proband, D$B100$E2$X2.proband)
  
  b91 <- survie_vec1(D$B100$E2$Y.non.probandIJ, t1, t2, D$B100$E2$X1.non.probandIJ, D$B100$E2$X2.non.probandIJ)
  b92 <- survie_vec2(D$B100$E2$Y.non.probandIK, t1, t2, D$B100$E2$X1.non.probandIK, D$B100$E2$X2.non.probandIK)
  b93 <- survie_vec2(D$B100$E2$Y.proband, t1, t2, D$B100$E2$X1.proband, D$B100$E2$X2.proband)
  
  b94 <- survie_vec2(D$B100$E2$Y.non.probandIJ, t1, t2, D$B100$E2$X1.non.probandIJ, D$B100$E2$X2.non.probandIJ)
  b95 <- survie_vec2(D$B100$E2$Y.non.probandIK, t1, t2, D$B100$E2$X1.non.probandIK, D$B100$E2$X2.non.probandIK)
  b96 <- survie_vec2(D$B100$E2$Y.proband, t1, t2, D$B100$E2$X1.proband, D$B100$E2$X2.proband)
  
  # On assemble tout dans une liste (ou un vecteur avec unlist)
  res <- list(
    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12,
    b13, b14, b15, b16, b17, b18, b19, b20, b21, b22, b23, b24,
    b25, b26, b27, b28, b29, b30, b31, b32, b33, b34, b35, b36,
    b37, b38, b39, b40, b41, b42, b43, b44, b45, b46, b47, b48,
    b49, b50, b51, b52, b53, b54, b55, b56, b57, b58, b59, b60,
    b61, b62, b63, b64, b65, b66, b67, b68, b69, b70, b71, b72,
    b73, b74, b75, b76, b77, b78, b79, b80, b81, b82, b83, b84,
    b85, b86, b87, b88, b89, b90, b91, b92, b93, b94, b95, b96
  )
  
  return(res)
}