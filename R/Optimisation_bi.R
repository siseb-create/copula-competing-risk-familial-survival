
compute.mle_Bi<-function(params,D,H)
{
  
  h11 <- H[1]; h12 <- H[2]; h22 <- H[3]
  #params=params0
  theta1=params[1:4]
  theta2=params[5:8]

  B_surv<-generer_B_surv(D, t1=theta1, t2=theta2)
  # Calcul dynamique des poids WB (dépendent de H)
  WB <- list(
    B1 = Poid.Biv101_E11(params, D$B101$E11, h11, h12),
    B2 = Poid.Biv101_E21(params, D$B101$E21, h11, h12, h22),
    B3 = Poid.Biv101_E12(params, D$B101$E12, h11, h12, h22),
    B4 = Poid.Biv101_E22(params, D$B101$E22, h12, h22),
    B5 = Poid.Biv110_E11(params, D$B110$E11, h11, h12),
    B6 = Poid.Biv110_E21(params, D$B110$E21, h11, h12, h22),
    B7 = Poid.Biv110_E12(params, D$B110$E12, h11, h12, h22),
    B8 = Poid.Biv110_E22(params, D$B110$E22, h12, h22),
    B9 = Poid.Biv100_E1(params, D$B100$E1, h11, h12, h22),
    B10= Poid.Biv100_E2(params, D$B100$E2, h11, h12, h22)
  )
  
  eps <- 1e-8
  
  # H = (h11, h12, h22)
  ui <- rbind(
    c( 1,  0,  0),   #  h11 >=  eps;
    c( 0,  1,  0),   #  h12 >=  eps
    c( 0,  0,  1),   #  h22 >=  eps
    c(-1,  0,  0),   # -h11 >= -(1-eps)  => h11 <= 1-eps
    c( 0, -1,  0),   # -h12 >= -(1-eps)  => h12 <= 1-eps
    c( 0,  0, -1),   # -h22 >= -(1-eps)  => h22 <= 1-eps
    c( 1, -1,  0),   #  h11 - h12 >= 0   => h12 <= h11
   c( 0, -1,  1)    # -h12 + h22 >= 0   => h12 <= h22
  )
  
  ci <- c(eps, eps, eps, -(1 - eps), -(1 - eps), -(1 - eps), 0, 0
          )

  h.op <-  constrOptim( theta = H,  f= Negloglik.Biv,grad=Negloglik.Biv_grad, 
                        mu = 1e-04,
                        ui= ui, ci = ci, D = D,B_surv= B_surv,WB=WB, 
                        control = list(reltol = 1e-6, maxit = 5000)#,method ="BFGS"#"Nelder-Mead"
                        )
  return(h.op$par)
} 


#h.op <-  optim(par = H,  fn= Negloglik.Biv_full,#grad= Negloglik.Biv_grad, 
#                      lower= c(0,0,0)+eps, upper = c(1,1,1)-eps, D = D,B_surv= B_surv,  params = params,WB=WB, 
#                      control = list(factr = 1e4,  # Remplace reltol
#                                     pgtol = 0, maxit = 5000),method = "L-BFGS-B"
#)
