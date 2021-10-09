library(MASS)
library(gplots)
library(parallel)
library(BiocGenerics)
library(multtest)
library(compiler)
library(EMMREML)
library(ape)
library(LDheatmap)
library(scatterplot3d)
library(plotly)

#source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
#source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

library(GAPIT3)

myG <- read.delim("127_hapmap.txt", header = F)
myY <- read.delim("new_working_pheno.txt", header = T)

####
setwd("farmCPU_pc1/")

myAnalysis <- GAPIT(
  G = myG,
  Y = myY[,c(1,2,3,4,5)],
  PCA.total = 1,
  kinship.algorithm=c("EMMA"),
  Inter.Plot = TRUE,
  model = "FarmCPU",
  Multiple_analysis = TRUE,
  Major.allele.zero = T
  #SNP.MAF = 0.05
)


setwd("../farmCPU_pc2/")

myAnalysis <- GAPIT(
  G = myG,
  Y = myY[,c(1,2,3,4,5)],
  PCA.total = 2,
  kinship.algorithm=c("EMMA"),
  Inter.Plot = TRUE,
  model = "FarmCPU",
  Multiple_analysis = TRUE,
  Major.allele.zero = T
  #SNP.MAF = 0.05
)

