rm(list=ls())
library(bio3d)

setwd("~/privateerpython/project_alliance/glycampdbfiles/omannose/IGgHCd/")

man7  <- read.pdb("cluster1.pdb")
IGgHCd <- read.pdb("AF-P01880-F1-model_v1.pdb")

n <- length(man7$atom[,1])
O1 <- as.numeric(man7$atom[1,9:11]) #O1
N <- as.numeric(IGgHCd$atom[max(which(IGgHCd$atom[,7] ==225))-1,9:11])
O1N<- N - O1
O1Nx <- rep(O1N[1],n)
O1Ny <- rep(O1N[2],n)
O1Nz <- rep(O1N[3],n)
O1N <- cbind(O1Nx,O1Ny,O1Nz)

man7$atom[,9:11] <- man7$atom[,9:11] + O1N
x <- numeric(0)
x[1:3]<- as.vector(t(man7$atom[1,9:11]))
for (i in 2:n){
x <- c(x,as.vector(t(man7$atom[i,9:11])))
} 
man7$xyz <- x      
write.pdb(man7,"tglycan.pdb")    
