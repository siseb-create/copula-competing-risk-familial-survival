source("Proba.R")

Generate.X1 <- function(p) {sample(x=c(-1,1),size=1,replace=TRUE,prob=c(1-p,p))}

#Generate.X2 <- function(p2=0.6) {sample(x=c(-1,1),size=1,replace=TRUE,prob=c(1-p2,p2))}
Generate.X2 <- function() {runif(1,min=-1,max=1)}

Generate.a <- function(m,v)
{
ab <- calculer.a.b.from.m.v(m,v)
shape <- ab$a
scale <- ab$b
rweibull(1,shape=shape,scale=scale)
} 

Generate.C <- function(min,max) {runif(5,min=min.c,max=max.c)} 
