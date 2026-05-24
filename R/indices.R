creer.indices.proband <- function(I)
{
seq(from=1,to=6*I-5,by=6)
}

creer.indices.non.proband <- function(I)
{
indices.non.proband <- NULL
for(i in 1:I)
{
for(j in 2:6)
{
indices.non.proband <- rbind(indices.non.proband,c(6*i-5,6*i+j-6))
}
}
indices.non.proband
}

creer.indices.bivariee<-function(I)
{
indices.bivariee <- NULL
for (i in 1:I) 
{
elements <- (6*i - 4):(6*i)
combs <- t(combn(elements, 2))
indices.bivariee <- rbind(indices.bivariee, cbind(6*i - 5, combs))
}
indices.bivariee
}


create.data.proband <- function(data,indices.proband) {data[indices.proband,-3]}

create.data.non.proband <- function(data,indices.non.proband,kin,I)
{
data.non.proband <- cbind(data[indices.non.proband[,1],-c(1,3,7)],data[indices.non.proband[,2],-1])
data.non.proband <- cbind(data.non.proband,rep(kin,I))
colnames(data.non.proband) <- c("Y.proband","X1.proband","X2.proband","Epsilon.proband",
                                "Y.non.proband","delta.non.proband","X1.non.proband","X2.non.proband","Epsilon.non.proband","ID","kin")
data.non.proband
}



create.data.bivarie<-function(data, indices.bivariee,kinJK,I)
{
data.bivarie<- cbind(data[indices.bivariee[,1],c(-1,-3,-7)],
                     data[indices.bivariee[,2], c(-1,-7)],
                     data[indices.bivariee[,3], -1],
                     do.call(rbind, replicate(I, kinJK, simplify = FALSE))
                     )
colnames(data.bivarie)<-c("Y.proband","X1.proband","X2.proband","Epsilon.proband",
                          "Y.non.probandIJ","delta.non.probandIJ","X1.non.probandIJ","X2.non.probandIJ","Epsilon.non.probandIJ",
                          "Y.non.probandIK","delta.non.probandIK","X1.non.probandIK","X2.non.probandIK","Epsilon.non.probandIK","ID",
                          "kin1J","kin1K","kinJK")
data.bivarie
}

#Data_Biv<-create.data.bivarie(data, indices.bivariee,kinJK,I)
