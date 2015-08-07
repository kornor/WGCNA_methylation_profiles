#####

######### This is a script to perform set up for a WGCNA on 
############ methylation results for genes with high EXP-METH correlation

#### set working directory and load libraries
setwd("~/Bioinformatics Work/Meth & RNA/WGCNA_Methylation_profiles")

library(plyr)
library(gplots)
library(ggplot2)


## load in the islands and shores files

islands <- read.table("Meth_TSS200_island.txt", sep = "\t", header = TRUE)
shores <- read.table("Meth_TSS200_Shore.txt", sep = "\t", header = TRUE)


list <- which(duplicated(islands[,1]))
islands <- islands[-list,]

list <- which(duplicated(shores[,1]))
shores <- shores[-list,]

write.table(islands, "Meth_islands.txt", sep = "\t")
write.table(shores, "Meth_shores.txt", sep = "\t")


islands <- read.table("Meth_islands.txt", sep = "\t", 
                      header = TRUE, row.names = 1)
shores <- read.table("Meth_shores.txt", sep = "\t", 
                     header = TRUE, row.names = 1)


## match the rows and write out again
cols <- intersect(colnames(islands), colnames(shores))
rows <- intersect(rownames(islands), rownames(shores))

islands <- islands[,cols]
shores <- shores[,cols]


######## these are now absolutely ready for use

islands <- read.table("Meth_islands.txt", sep = "\t", 
                      header = TRUE, row.names = 1)
shores <- read.table("Meth_shores.txt", sep = "\t", 
                     header = TRUE, row.names = 1)



#### create a new dataframe for the ave of the values in both

  meth <- matrix(nrow = ,ncol = )

# Calculate ave between
c1.r = 1:(ncol(Meth))
c1.p = 1:(ncol(Meth))
for(i in 1:(ncol(Meth))) {
  c1 = mean(Meth[,i],Exp[,i])
  c1.r[i] = c1$estimate
  c1.p[i] = c1$p.value
}




#### Now load in the clinical file and trim to basal only

clin <- read.table("Clinical_final_TCGA.txt", sep = "\t", header = TRUE)
pam50 <- as.factor(clin$PAM50Call_RNAseq)


bas_shores <- shores[pam50 == 'Basal',]
write.table(bas_shores, "Meth_shores_basal.txt", sep = "\t")


bas_islands <- islands[pam50 == 'Basal',]
write.table(bas_islands, "Meth_islands_basal.txt", sep = "\t")
