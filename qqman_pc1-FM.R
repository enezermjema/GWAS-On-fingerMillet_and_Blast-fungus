library(ggplot2)
library(qqman)
library(ggpubr)
library(gridExtra)

####### plotting results for the first one PC
getwd()
setwd("GWAS/qqman_plots/pc1/")

# importing results

pc1_e43 <- read.csv("e43.csv")
pc1_k1 <- read.csv("k1.csv")
pc1_t1 <- read.csv("t1.csv")
pc1_u8 <- read.csv("u8.csv")

# converting chromosome column into numeric

pc1_e43$Chromosome <- as.numeric(pc1_e43$Chromosome)
pc1_k1$Chromosome <- as.numeric(pc1_k1$Chromosome)
pc1_t1$Chromosome <- as.numeric(pc1_t1$Chromosome)
pc1_u8$Chromosome <- as.numeric(pc1_u8$Chromosome)


# manhattan plots

manhattan(pc1_e43, main = "E43 isolate", chr = "Chromosome", bp = "Position", snp = "SNP", p = "P.value", cex.lab = 1.8, cex.main = 2, cex=1.8, cex.axis=1.5, col=c("blue4","orange3"), suggestiveline=F, genomewideline=-log10(5.008515245e-6), ylim=c(0, 12))
manhattan(pc1_k1, main = "K1 isolate", chr = "Chromosome", bp = "Position", snp = "SNP", p = "P.value", cex.lab = 1.8, cex.main = 2, cex=1.8, cex.axis=1.5, col=c("blue4","orange3"), suggestiveline=F, genomewideline=-log10(5.008515245e-6), ylim=c(0, 12))
manhattan(pc1_t1, main = "T1 isolate", chr = "Chromosome", bp = "Position", snp = "SNP", p = "P.value", cex.lab = 1.8, cex.main = 2, cex=1.8, cex.axis=1.5, col=c("blue4","orange3"), suggestiveline=F, genomewideline=-log10(5.008515245e-6), ylim=c(0, 15))
manhattan(pc1_u8, main = "U8 isolate", chr = "Chromosome", bp = "Position", snp = "SNP", p = "P.value", cex.lab = 1.8, cex.main = 2, cex=1.8, cex.axis=1.5, col=c("blue4","orange3"), suggestiveline=F, genomewideline=-log10(5.008515245e-6), ylim=c(0, 12))

# qqplots

qq(pvector = pc1_e43$P.value, cex.lab = 1.1, cex.main = 2, cex=1.5, cex.axis=1, col="blue")
qq(pvector = pc1_k1$P.value, cex.lab = 1.1, cex.main = 2, cex=1.5, cex.axis=1, col="blue")
qq(pvector = pc1_t1$P.value, cex.lab = 1.1, cex.main = 2, cex=1.5, cex.axis=1, col="blue")
qq(pvector = pc1_u8$P.value, cex.lab = 1.1, cex.main = 2, cex=1.5, cex.axis=1, col="blue")




