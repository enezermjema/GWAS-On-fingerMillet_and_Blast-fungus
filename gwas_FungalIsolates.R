getwd()
#setwd("rGWAS/pc1/okhale and ie/")

# importing dataset

myG_fmb <- read.delim("../../fmb_hmp_final_28Jan.txt", header = F)
myY_fmb <- read.delim("../../all_isolate_incremented.txt")
myY_ea <- myY_fmb[,c(1,2,3,6,7)]

# pc1

myAnalysis <- GAPIT(
  G = myG_fmb,
  Y = myY_ea[,c(1,2,3,4,5)],
  PCA.total = 1,
  kinship.algorithm=c("EMMA"),
  Inter.Plot = TRUE,
  model = "FarmCPU",
  Multiple_analysis = TRUE,
  Major.allele.zero = T
  #SNP.MAF = 0.05
)

# pc2
setwd("../pc2/")

myAnalysis <- GAPIT(
  G = myG_fmb,
  Y = myY_ea[,c(1,2,3,4,5)],
  PCA.total = 2,
  kinship.algorithm=c("EMMA"),
  Inter.Plot = TRUE,
  model = "FarmCPU",
  Multiple_analysis = TRUE,
  Major.allele.zero = T
  #SNP.MAF = 0.05
)


