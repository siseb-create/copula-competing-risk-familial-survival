############################################################
## Iterative maximum likelihood estimation procedure
##
## This file contains the main function `compute.mle`, which
## implements the iterative optimization procedure for the
## proposed copula-based familial survival model.
##
## The algorithm updates the marginal model parameters,
## and the intra-family dependence parameters h.
############################################################

# Optimisation itt?ractive

compute.mle<-function(params,Data_Uni,D,H)
{

  Nb_itt=0
  
  repeat {
    
  Theta0 <-c(params,H)
  
  params_op<-compute.mle_uni(params=params,Data_Uni=Data_Uni,H=H)
  params<-params_op
  
  H.op <-  compute.mle_Bi(params=params,D=D,H=H)
  Theta <- c(params,H.op)
 
  #print(Theta)
 
    Nb_itt=Nb_itt+1
    
  #DeltaTheta <- sqrt(sum((Theta-Theta0)^2)) / sqrt(length(Theta))
    
  if (max(abs(Theta-Theta0))<1e-3){break}
    
    H<-H.op
    
    }
  return(Theta)
} 
