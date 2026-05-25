
############################################################
## Marginal parameter update
##
## This function updates the marginal model parameters
## conditionally on the current value of the intra-family
## dependence parameters H = (h11, h12, h22).
##
## The function first computes the univariate likelihood weights
## associated with censored non-proband observations, and then
## minimizes the negative univariate log-likelihood using
## numerical optimization.
############################################################

compute.mle_uni<-function(params,Data_Uni,H)
{
  #print(params)
  
  h11<-H[1]#exp(l.H[1])
  h12<-H[2]#exp(l.H[2])
  h22<-H[3]#exp(l.H[3])
  
  WU.1<-Uni.poids.delta0.E1(params=params,data.non.proband.delta.0.E.proband1=Data_Uni$d0_P1,h11=h11,h12=h12)
  WU.2<-Uni.poids.delta0.E2(params=params,data.non.proband.delta.0.E.proband2=Data_Uni$d0_P2,h12=h12,h22=h22)

  params.op <-  optim(par = params,fn = Negloglik.uni_safe,gr=Negloglik.uni_grad,
                      Data_Uni=Data_Uni,H=H,WU.1=WU.1,WU.2=WU.2,
                      control = list(maxit=5000,reltol=1e-6)#,method ="BFGS" #avec reltol=1e-3
                      )
return(params.op$par)
} 


#1.3830420 -6.6989233  0.8521414 -0.9138038  1.4100580 -6.9143266  0.8119666 -0.9335378
