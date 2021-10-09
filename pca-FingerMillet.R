library(adegenet)
library(factoextra)
library("RColorBrewer")
library("magrittr")
library(ade4)
library(ggplot2)
library(ggpubr)
library(dplyr)

getwd()

geno_hmp <- read.delim(file = 'final_geno_pheno_notdiplo.hmp.txt', check.names = FALSE, header = TRUE, stringsAsFactors = FALSE) #load marker data
pop_file <- read.delim("only_ea-127.txt", header = F) # loading population file

## filtering ea

geno_127 <- select(geno_hmp, 
                   c("100002","203542","208725","215852","215863","215984","216035","216036","216039","216045",
                     "216049","228202","229731","230105","235699",
                     "236446","237443","237475","238299","238305","238460","242114","242115",
                     "242116","242617","242618","242619_1","242619_2","245085",
                     "AAUFM_02","AAUFM_04","AAUFM_06","AAUFM_07","AAUFM_10","AAUFM_15","AAUFM_17",
                     "AAUFM_19","AAUFM_28","AAUFM_32","AAUFM_42","AAUFM_44","AAUFM_50","BKFM_0005",
                     "BKFM_0006","BKFM_0007","BKFM_0009","BKFM_0010","BKFM_0011","BKFM_0014",
                     "BKFM_0018","BKFM_0020","BKFM_0023","BKFM_0028","BKFM_0029","BKFM_0034",
                     "BKFM_0038","BKFM_0042","BKFM_0043","BKFM_0047","EK_01","EK_02","EK_03",
                     "EK_04","EK_08","EX_ALUPE","IE_2335","IE_2367","IE_2464","KATFM_1","KNE_10057",
                     "KNE_626","KNE_796","MD_20","NANJALA","P_221","P_227","P_318","PI_318897",
                     "PI_321114","PI_321125","PW_001022","SA_77_SADC","TADESSE","TZ_3876",
                     "TZA_116","TZA_141","TZA_1565","TZA_1570","TZA_1585","TZA_1617","TZA_1622",
                     "TZA_1629","TZA_1636","TZA_1637","TZA_1638","TZA_164","TZA_1642","TZA_1647",
                     "TZA_1649","TZA_1650","TZA_1653","TZA_1660","TZA_1667","TZA_1688","TZA_1689",
                     "TZA_1704","TZA_206","TZA_2361","TZA_2720","TZA_2924","TZA_2975","TZA_2999",
                     "TZA_3114","TZA_3721","TZA_3876","TZA_3989","TZA_3995","TZA_4134","TZA_626",
                     "TZA_699","TZA_703","TZA_709","TZA_86","TZA_87","TZA_98","TZA709","WHITESS_5"
                   ))


hapmap <- cbind(geno_hmp[,1:11], geno_127) # genotype for gwas
write.table(hapmap, "fm_haploiFormat_1208.txt", sep = "\t", quote = F, row.names = F)
##### filtering over

geno_numeric <- hap_to_G(hapmap, 11) # make a numeric marker matrix 
row.names(geno_numeric) <- geno_numeric$`rs#` 

transposed_matrix <- t(geno_numeric[, 12:ncol(geno_numeric)])
#row.names(transposed_matrix) <- NULL

genind_obj <- df2genind(transposed_matrix, ploidy=2, type="PA", sep = "\t") # create a genind object from our matrix

pop(genind_obj) <- paste(pop_file$V2) # appending country info
indNames(genind_obj) <- paste(pop_file$V1) # appending the individual names 

#indNames(genind_obj)
#write.table(as.data.frame(transposed_matrix), "pca_matrix.txt", quote = F, sep = "\t")

x <- tab(genind_obj, NA.method="mean") # relpacing NAs for PCA analysis

pca2 <- dudi.pca(x, scale = F, scannf = FALSE, nf = 4) # pca analysis

<- get_eigenvalue(pca2)

fviz_contrib(pca2, choice = "ind", axes = c(1,2), top = 20) # contribution plot

# scree plot
scree_1 <- fviz_eig(pca2, main = "Screeplot - Eigenvalues", linecolor =  "black", ncp = 10,
                  barcolor = heat.colors(10), barfill = heat.colors(10), addlabels = T, ylim = c(0,30),
                  ggtheme = theme_classic(), font.title = c(14,"bold","black"),
                  font.x = c(14,"bold","black"),
                  font.y =  c(14,"bold","black")) +
                  theme(text = element_text(size=rel(4.5)))# plot scree plot on percentage
scree_1
# factoextra-based plots-grouping by country
## PC 1+2

#options(ggrepel.max.overlaps = 10)


groups <- as.factor(genind_obj$pop)

col <- brewer.pal(n=5, name = "Dark2") # n = number of groups

ind_pc12_1 <- fviz_pca_ind(pca2,
                         axes = c(1, 2),
                         label = "none",
                         col.ind = groups, # color by groups
                         palette = col,
                         #addEllipses = T, # Concentration ellipses
                         #ellipse.type = "confidence",
                         legend.title = "Groups",
                         title="PC 1 vs 2",
                         repel = TRUE,
                         font.title = c(12,"bold","black"),
                         font.x = c(12,"bold","black"),
                         font.y =  c(12,"bold","black")
) + labs(x = "PC1 (25.8%)", y = "PC2 (13.9%)")
ind_pc12_1


ind_pc13_1 <- fviz_pca_ind(pca2,
                         axes = c(1, 3),
                         label = "none",
                         col.ind = groups, # color by groups
                         palette = col,
                         #addEllipses = TRUE, # Concentration ellipses
                         #ellipse.type = "confidence",
                         legend.title = "Groups",
                         title="PC 1 vs 3",
                         repel = TRUE,
                         font.title = c(12,"bold","black"),
                         font.x = c(12,"bold","black"),
                         font.y =  c(12,"bold","black")
) + labs(x = "PC1 (25.8%)", y = "PC3 (6.8%)")
ind_pc13_1
## ggarrange
ggarrange(ind_pc12, ind_pc13,
          labels = "AUTO"
)



# contributions PCA
fviz_pca_ind(pca1, col.ind = "cos2",
             axes = c(1, 3),
             gradient.cols = c("blue","yellow","red"),
             repel = T
)

fviz_pca_ind(pca1, col.ind = "cos2",
             axes = c(1, 2),
             gradient.cols = c("blue","yellow","red"),
             repel = T
)

fviz_pca_ind(pca1, col.ind = "cos2",
             gradient.cols = c("blue","yellow","red"),
             repel = T
)



