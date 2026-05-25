############################################################
## Conditional functions and data generation functions
##
## This file contains the conditional functions and data
## generation procedures required to simulate synthetic
## family-structured survival data under the proposed model.
############################################################


Generate.proband.epsilon <- function(a,theta1,theta2,X1,X2)
{ifelse(runif(1)<PProba(a,theta1,theta2,X1,X2),1,2)}


Generate.non.proband.epsilon <- function(theta1,theta2,X1,X2)
{ifelse(runif(1)<Proba1(theta1,theta2,X1,X2),1,2)}



Generate.survival.proband <- function(a,gam1,theta1,theta2,X1,X2)
{
T.proband <- uniroot(survie.cond.proband,lower=0,upper=a,a=a,
                     gam1=gam1,theta1=theta1,theta2=theta2,X1=X1,X2=X2,res=runif(1))$root
Z.proband <- qnorm(survie.cond(t=T.proband,gam1=gam1,theta1=theta1,theta2=theta2,X1=X1,X2=X2))
return(data.frame(T.proband=T.proband,Z.proband=Z.proband))
}


Generate.survival.non.proband <- function(Z.proband,X1.non.proband,X2.non.proband,epsilon.non.proband,theta1,theta2,V)
{
Z.non.proband <- rcmvnorm(n=1,mean=rep(0,6),sigma=V,dependent.ind=(2:6),given.ind=1,X.given=Z.proband)[1,]
U.non.proband <- pnorm(Z.non.proband)
T.non.proband <- rep(-999,5)
for (i in 1:5)
{
gam1 <- (2-epsilon.non.proband[i])*theta1+(epsilon.non.proband[i]-1)*theta2
#gam2 <- (epsilon.non.proband[i]-1)*theta1+(2-epsilon.non.proband[i])*theta2
X1 <- X1.non.proband[i]
X2 <- X2.non.proband[i]
res <- U.non.proband[i]
T.non.proband[i] <- uniroot(survie.cond,lower=0,upper=20,gam1=gam1,theta1=theta1,theta2=theta2,X1=X1,X2=X2,res=res)$root
}
T.non.proband
}


Generate.survival.family <- function(theta1,theta2,kinship,h11,h22,h12,m.a,v.a,p,min.c,max.c)
{
a <- Generate.a(m.a,v.a)
C <- c(a,Generate.C(min.c,max.c))
X1.proband <- Generate.X1(p)
X2.proband <- Generate.X2() 
epsilon.proband <- Generate.proband.epsilon(a,theta1,theta2,X1.proband,X2.proband)
gam1.epsilon <- (2-epsilon.proband)*theta1+(epsilon.proband-1)*theta2
#gam2.epsilon <- (epsilon.proband-1)*theta1+(2-epsilon.proband)*theta2
survival.proband <- Generate.survival.proband(a=a,gam1=gam1.epsilon,theta1=theta1,theta2=theta2,X1=X1.proband,X2=X2.proband)
T.proband <- survival.proband$T.proband
Z.proband <- survival.proband$Z.proband

X1.non.proband <- rep(-999,5)
X2.non.proband <- rep(-999,5)
epsilon.non.proband <- rep(-999,5)

for(i in 1:5)
{
X1.non.proband[i] <- Generate.X1(p)
X2.non.proband[i] <- Generate.X2()
epsilon.non.proband[i] <- Generate.non.proband.epsilon(theta1,theta2,X1.non.proband[i],X2.non.proband[i])
}

epsilon <- c(epsilon.proband,epsilon.non.proband)
V <- Calculer.V(epsilon,kinship,h11,h12,h22)
T.non.proband <- Generate.survival.non.proband(Z.proband,X1.non.proband,X2.non.proband,epsilon.non.proband,theta1,theta2,V)

X1 <- c(X1.proband,X1.non.proband)
X2 <- c(X2.proband,X2.non.proband)
T <- c(T.proband,T.non.proband)

Y <- pmin(T,C)
delta <- as.numeric(T<C)

data.frame(a=a,Y=Y,delta=delta,X1=X1,X2=X2,epsilon=epsilon)

}

Generate.data <- function(I,theta1,theta2,kinship,h11,h22,h12,m.a,v.a,p,min.c,max.c)
{
  data <- NULL
  for(i in 1:I)
  {
    data.family <- Generate.survival.family(theta1,theta2,kinship,h11,h22,h12,m.a,v.a,p,min.c,max.c)
    data <- rbind(data,data.family)
  }
  data$ID<-rep(1:I, each = 6)
  data
}



