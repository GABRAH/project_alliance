rm(list=ls())
  library(bio3d)

#
#This is a script which changed the anomers for all of the glycans *** DO NOT RUN THIS AGAIN *** I didn't change the directories for a reason. Don't do it.   
#


O1switch <- function(glycan){
  man5  <- read.pdb(glycan)
  O1 <- as.numeric(man5$atom[4,9:11]) #O1
  C1 <- as.numeric(man5$atom[3,9:11]) #C1
  H1 <- as.numeric(man5$atom[2,9:11]) #H1
  
  C1O1 <- O1 - C1
  x <- sqrt(sum(C1O1*C1O1))
  nC1O1 <- C1O1 / x
  O1 <- C1 + 1.4*nC1O1
  man5$atom[2,9:11] <- t(O1)
  man5$xyz[4:6] <- O1
  
  C1H1 <- H1 - C1
  x<-sqrt(sum(C1H1*C1H1))
  nC1H1 <- C1H1 / x
  H1 <- C1 + 1.09*nC1H1
  man5$atom[4,9:11] <- t(H1)
  man5$xyz[10:12] <- H1
  
  man5$atom <- man5$atom[-1,] 
  man5$xyz <- man5$xyz[-c(1:3)]
  man5$atom[which(man5$atom[,5] =="ROH" | man5$atom[,5] =="4YA" ),5] <- "4YB"
  man5$atom[which(man5$atom[,7] == 1),7] <- 2
  write.pdb(man5,file = glycan)  
}

setwd("~/privateerpython/project_alliance/glycampdbfiles/omannosePrivateer/man9/")
O1switch("cluster1_1.pdb")
O1switch("cluster2_4.pdb")
O1switch("cluster3_4.pdb")
O1switch("cluster4_1.pdb")
O1switch("cluster4_2.pdb")
O1switch("cluster5_1.pdb")
O1switch("cluster5_1.pdb")
O1switch("cluster6_1.pdb")
O1switch("cluster6_2.pdb")
O1switch("cluster6_3.pdb")
O1switch("cluster7_3.pdb")
setwd("../man8_3/")
O1switch("cluster5_2.pdb")
setwd("../man8_2/")
O1switch("cluster1_1.pdb")
setwd("../man5")
O1switch("cluster1_1.pdb")

setwd("~/privateerpython/project_alliance/glycampdbfiles/glycanfrags/chito/")
O1switch("chito2.pdb")


setwd("/mnt/7825866d-3ecf-419f-bfb2-731c265f392f/Data/piliwork/Glycans/tria/1/")

#setwd("/mnt/7825866d-3ecf-419f-bfb2-731c265f392f/Data/highmanose/omannose/man5/")
#O1switch("cluster1.pdb")
#O1switch("cluster2.pdb")
#O1switch("cluster3.pdb")
#O1switch("cluster4.pdb")
#setwd("/mnt/7825866d-3ecf-419f-bfb2-731c265f392f/Data/highmanose/omannose/man6_1/")
#O1switch("cluster1.pdb")
#O1switch("cluster2.pdb")
#O1switch("cluster3.pdb")
O1switch("cluster4.pdb")
O1switch("cluster5.pdb")
setwd("/mnt/7825866d-3ecf-419f-bfb2-731c265f392f/Data/highmanose/omannose/man6_2/")
O1switch("cluster1.pdb")
O1switch("cluster2.pdb")
O1switch("cluster3.pdb")
O1switch("cluster4.pdb")
setwd("/mnt/7825866d-3ecf-419f-bfb2-731c265f392f/Data/highmanose/omannose/man6_3/")
O1switch("cluster1.pdb")
O1switch("cluster2.pdb")
O1switch("cluster3.pdb")
O1switch("cluster4.pdb")
setwd("/mnt/7825866d-3ecf-419f-bfb2-731c265f392f/Data/highmanose/omannose/man7_1/")
O1switch("cluster1.pdb")
O1switch("cluster2.pdb")
O1switch("cluster3.pdb")
O1switch("cluster4.pdb")
O1switch("cluster5.pdb")
setwd("/mnt/7825866d-3ecf-419f-bfb2-731c265f392f/Data/highmanose/omannose/man7_2/")
O1switch("cluster1.pdb")
O1switch("cluster2.pdb")
O1switch("cluster3.pdb")
O1switch("cluster4.pdb")
O1switch("cluster5.pdb")
setwd("/mnt/7825866d-3ecf-419f-bfb2-731c265f392f/Data/highmanose/omannose/man7_3/")
O1switch("cluster1.pdb")
O1switch("cluster2.pdb")
O1switch("cluster3.pdb")
O1switch("cluster4.pdb")
O1switch("cluster5.pdb")
setwd("/mnt/7825866d-3ecf-419f-bfb2-731c265f392f/Data/highmanose/omannose/man7_4/")
O1switch("cluster1.pdb")
O1switch("cluster2.pdb")
O1switch("cluster3.pdb")
O1switch("cluster4.pdb")
setwd("/mnt/7825866d-3ecf-419f-bfb2-731c265f392f/Data/highmanose/omannose/man8_1/")
O1switch("cluster1.pdb")
O1switch("cluster2.pdb")
O1switch("cluster3.pdb")
O1switch("cluster5.pdb")
setwd("/mnt/7825866d-3ecf-419f-bfb2-731c265f392f/Data/highmanose/omannose/man8_2/")
O1switch("cluster1.pdb")
O1switch("cluster2.pdb")
O1switch("cluster3.pdb")
O1switch("cluster4.pdb")
O1switch("cluster5.pdb")
setwd("/mnt/7825866d-3ecf-419f-bfb2-731c265f392f/Data/highmanose/omannose/man8_3/")
O1switch("cluster1.pdb")
O1switch("cluster2.pdb")
O1switch("cluster3.pdb")
O1switch("cluster4.pdb")
O1switch("cluster5.pdb")
setwd("/mnt/7825866d-3ecf-419f-bfb2-731c265f392f/Data/highmanose/omannose/man9/")
O1switch("cluster1.pdb")
O1switch("cluster2.pdb")
O1switch("cluster3.pdb")
O1switch("cluster4.pdb")
O1switch("cluster5.pdb")
O1switch("cluster6.pdb")






