# =========================================================
# Fonction qui extrait les statistiques d'une base Data_Biv
# =========================================================
get_stats_biv <- function(Data_Biv) {
  
  groupes <- list(
    "Data_Biv111_E111" = with(Data_Biv, delta.non.probandIJ == 1 & delta.non.probandIK == 1 &
                                Epsilon.proband == 1 & Epsilon.non.probandIJ == 1 & Epsilon.non.probandIK == 1),
    "Data_Biv111_E211" = with(Data_Biv, delta.non.probandIJ == 1 & delta.non.probandIK == 1 &
                                Epsilon.proband == 2 & Epsilon.non.probandIJ == 1 & Epsilon.non.probandIK == 1),
    
    "Data_Biv111_E121" = with(Data_Biv, delta.non.probandIJ == 1 & delta.non.probandIK == 1 &
                                Epsilon.proband == 1 & Epsilon.non.probandIJ == 2 & Epsilon.non.probandIK == 1),
    "Data_Biv111_E221" = with(Data_Biv, delta.non.probandIJ == 1 & delta.non.probandIK == 1 &
                                Epsilon.proband == 2 & Epsilon.non.probandIJ == 2 & Epsilon.non.probandIK == 1),
    
    "Data_Biv111_E112" = with(Data_Biv, delta.non.probandIJ == 1 & delta.non.probandIK == 1 &
                                Epsilon.proband == 1 & Epsilon.non.probandIJ == 1 & Epsilon.non.probandIK == 2),
    "Data_Biv111_E212" = with(Data_Biv, delta.non.probandIJ == 1 & delta.non.probandIK == 1 &
                                Epsilon.proband == 2 & Epsilon.non.probandIJ == 1 & Epsilon.non.probandIK == 2),
    
    "Data_Biv111_E122" = with(Data_Biv, delta.non.probandIJ == 1 & delta.non.probandIK == 1 &
                                Epsilon.proband == 1 & Epsilon.non.probandIJ == 2 & Epsilon.non.probandIK == 2),
    "Data_Biv111_E222" = with(Data_Biv, delta.non.probandIJ == 1 & delta.non.probandIK == 1 &
                                Epsilon.proband == 2 & Epsilon.non.probandIJ == 2 & Epsilon.non.probandIK == 2),
    
    "Data_Biv101_E11"  = with(Data_Biv, delta.non.probandIJ == 0 & delta.non.probandIK == 1 &
                                Epsilon.proband == 1 & Epsilon.non.probandIK == 1),
    "Data_Biv101_E21"  = with(Data_Biv, delta.non.probandIJ == 0 & delta.non.probandIK == 1 &
                                Epsilon.proband == 2 & Epsilon.non.probandIK == 1),
    "Data_Biv101_E12"  = with(Data_Biv, delta.non.probandIJ == 0 & delta.non.probandIK == 1 &
                                Epsilon.proband == 1 & Epsilon.non.probandIK == 2),
    "Data_Biv101_E22"  = with(Data_Biv, delta.non.probandIJ == 0 & delta.non.probandIK == 1 &
                                Epsilon.proband == 2 & Epsilon.non.probandIK == 2),
    
    "Data_Biv110_E11"  = with(Data_Biv, delta.non.probandIK == 0 & delta.non.probandIJ == 1 &
                                Epsilon.proband == 1 & Epsilon.non.probandIJ == 1),
    "Data_Biv110_E21"  = with(Data_Biv, delta.non.probandIK == 0 & delta.non.probandIJ == 1 &
                                Epsilon.proband == 2 & Epsilon.non.probandIJ == 1),
    "Data_Biv110_E12"  = with(Data_Biv, delta.non.probandIK == 0 & delta.non.probandIJ == 1 &
                                Epsilon.proband == 1 & Epsilon.non.probandIJ == 2),
    "Data_Biv110_E22"  = with(Data_Biv, delta.non.probandIK == 0 & delta.non.probandIJ == 1 &
                                Epsilon.proband == 2 & Epsilon.non.probandIJ == 2),
    
    "Data_Biv100_E1"   = with(Data_Biv, delta.non.probandIJ == 0 & delta.non.probandIK == 0 &
                                Epsilon.proband == 1),
    "Data_Biv100_E2"   = with(Data_Biv, delta.non.probandIJ == 0 & delta.non.probandIK == 0 &
                                Epsilon.proband == 2)
  )
  
  effectifs <- sapply(groupes, sum, na.rm = TRUE)
  
  nb_familles <- sapply(groupes, function(ind) {
    tmp <- Data_Biv[ind, , drop = FALSE]
    if (nrow(tmp) == 0) return(0)
    length(unique(tmp$ID))
  })
  
  list(
    effectifs = effectifs,
    nb_familles = nb_familles
  )
}


get_stats_brut <- function(data) {
  
  effectifs <- c(
    Censure = sum(data$delta == 0, na.rm = TRUE),
    Cause1  = sum(data$delta == 1 & data$epsilon == 1, na.rm = TRUE),
    Cause2  = sum(data$delta == 1 & data$epsilon == 2, na.rm = TRUE)
  )
  
  nb_familles <- c(
    Censure = length(unique(data$famID[data$delta == 0])),
    Cause1  = length(unique(data$famID[data$delta == 1 & data$epsilon == 1])),
    Cause2  = length(unique(data$famID[data$delta == 1 & data$epsilon == 2]))
  )
  
  list(
    effectifs = effectifs,
    nb_familles = nb_familles
  )
}
