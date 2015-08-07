



#load packages for use
library(permute)
library(lattice)
library(vegan)
library(scatterplot3d)
library(MASS)
library(nlme)
library(mgcv)
library(cluster)
library(rgl)
library(ape)
library(picante)
library(FactoMineR)
library(gdata)
library(WriteXLS)
library(plyr)
library(caret)
library(BiodiversityR)
library(gtools)
library(AppliedPredictiveModeling)
library(limma)
library(vegetarian)
library(survival)
transparentTheme(trans = 0.4)


### Start with the islands set

setwd("~/Bioinformatics Work/Meth & RNA/WGCNA_Methylation_profiles/WGCNA_meth_set")

#Concatenate
#islands
meth_list <- list()
for (i in 1:32) 
{
  input_file = paste("Tumours_binary.",i,"_processed.txt_Island.txt",sep="")
  meth_list[[i]] <- read.table(input_file, sep ="\t", header = TRUE)
}
Tumours_all_binary <- do.call(rbind,meth_list)

## Write out this one too
write.table(Tumours_all_binary, "Meth_TSS200_island.txt", sep = "\t")

#Concatenate
#shore
meth_list <- list()
for (i in 1:32) 
{
  input_file = paste("Tumours_binary.",i,"_processed.txt_Shore.txt",sep="")
  meth_list[[i]] <- read.table(input_file, sep ="\t", header = TRUE)
}
Tumours_all_binary <- do.call(rbind,meth_list)

## Write out this one too
write.table(Tumours_all_binary, "Meth_TSS200_Shore.txt", sep = "\t")

