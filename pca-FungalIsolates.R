setwd("fmb/")
getwd()

fmb <- read.delim("fmb_hmp_final_14April_trial.txt") 
fmb_ea <- fmb[,c(1:222)] # extracting only ea

###################
ea_only <- select(fmb,
                  c("E1","E12","E13","E14","E15","E16","E18","E19","E20","E21","E22","E24","E26","E27","E28","E29","E3",
                    "E30","E31","E32","E33","E35","E36","E37","E38","E39","E40","E40","E41","E42","E44","E45","E46","E47",
                    "E48","E49","E5","E51","E52","E53","E54","E55","E56","E57","E58","E60","E61","E62","E63","K10","K11","K12",
                    "K13","K14","K15","K16","K17","K18","K20","K21","K23","K24","K25","K26","K27","K29","K30","K31","K32","K33",
                    "K35","K36","K37","K38","K39","K4","K40","K41","K42","K43","K44","K45","K7","K8","K9","T10","T13","T14","T15",
                    "T16","T17","T18","T19","T2","T20","T21","T22","T23","T24","T25","T26","T27","T29","T3","T30","T31","T32","T33",
                    "T34","T35","T36","T37","T38","T40","T41","T42","T43","T44","T45","T46","T47","T48","T49","T5","T50","T51","T52",
                    "T53","T54","T56","T57","T58","T7","T8","T9","U1","U11","U12","U13","U14","U15","U16","U17","U18","U19","U2","U20",
                    "U21","U22","U23","U24","U25","U26","U27","U28","U29","U3","U30","U31","U32","U33","U35","U36","U37","U38","U39",
                    "U4","U40","U41","U42","U43","U44","U45","U48","U49","U5","U50","U51","U52","U53","U55","U56","U57","U58","U6","U7",
                    "U9"
                  )
  
)

fmb_ea1 <- cbind(fmb[,c(1:11)], ea_only)

#####################
#write.table(fmb_ea, "fmb_ea-only_hmp.txt", quote = F, row.names = F, sep = "\t")
target_iso <- c("E1","E12","E13","E14","E15","E16","E18","E19","E20","E21","E22","E24","E26","E27","E28","E29","E3",
                "E30","E31","E32","E33","E35","E36","E37","E38","E39","E40","E40","E41","E42","E44","E45","E46","E47",
                "E48","E49","E5","E51","E52","E53","E54","E55","E56","E57","E58","E60","E61","E62","E63","K10","K11","K12",
                "K13","K14","K15","K16","K17","K18","K20","K21","K23","K24","K25","K26","K27","K29","K30","K31","K32","K33",
                "K35","K36","K37","K38","K39","K4","K40","K41","K42","K43","K44","K45","K7","K8","K9","T10","T13","T14","T15",
                "T16","T17","T18","T19","T2","T20","T21","T22","T23","T24","T25","T26","T27","T29","T3","T30","T31","T32","T33",
                "T34","T35","T36","T37","T38","T40","T41","T42","T43","T44","T45","T46","T47","T48","T49","T5","T50","T51","T52",
                "T53","T54","T56","T57","T58","T7","T8","T9","U1","U11","U12","U13","U14","U15","U16","U17","U18","U19","U2","U20",
                "U21","U22","U23","U24","U25","U26","U27","U28","U29","U3","U30","U31","U32","U33","U35","U36","U37","U38","U39",
                "U4","U40","U41","U42","U43","U44","U45","U48","U49","U5","U50","U51","U52","U53","U55","U56","U57","U58","U6","U7",
                "U9"
)

pop_fmb <- read.delim("fmb_origin.txt", header = F) # loading population file
pop_fmb1 <- filter(pop_fmb, V1 %in% target_iso)


geno_numeric_fmb <- hap_to_G(fmb_ea1, 11) # make a numeric marker matrix 
row.names(geno_numeric_fmb) <- geno_numeric_fmb$`rs#` 

transposed_matrix <- t(geno_numeric_fmb[, 12:ncol(geno_numeric_fmb)])
#row.names(transposed_matrix) <- NULL

genind_obj <- df2genind(transposed_matrix, ploidy=2, type="PA", sep = "\t") # create a genind object from our matrix

pop(genind_obj) <- paste(pop_fmb1$V2) # appending country info
indNames(genind_obj) <- paste(pop_fmb1$V1) # appending the individual names 

#indNames(genind_obj)
#write.table(as.data.frame(transposed_matrix), "pca_matrix.txt", quote = F, sep = "\t")

x <- tab(genind_obj, NA.method="mean") # relpacing NAs for PCA analysis

pca2 <- dudi.pca(x, scale = F, scannf = FALSE, nf = 4) # pca analysis

fviz_contrib(pca2, choice = "ind", axes = c(1,2), top = 30) # contribution plot

get_eigenvalue(pca2)

# scree plot
scree_2 <- fviz_eig(pca2, main = "Screeplot - Eigenvalues", linecolor =  "black", ncp = 10,
                  barcolor = heat.colors(10), barfill = heat.colors(10), addlabels = T, ylim = c(0,45),
                  ggtheme = theme_classic(), font.title = c(14,"bold","black"),
                  font.x = c(14,"bold","black"),
                  font.y =  c(14,"bold","black")) +
  theme(text = element_text(size=rel(4.5)))# plot scree plot on percentage
scree_2
# factoextra-based plots-grouping by country
## PC 1+2
ggarrange(scree_2, scree_1,
          labels = "AUTO"
)




#options(ggrepel.max.overlaps = 10)


groups <- as.factor(genind_obj$pop)

col <- brewer.pal(n=5, name = "Dark2") # n = number of groups

ind_pc12 <- fviz_pca_ind(pca2,
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
) + labs(x = "PC1 (44.4%)", y = "PC2 (6.4%)")
ind_pc12


ind_pc13 <- fviz_pca_ind(pca2,
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
) + labs(x = "PC1 (44.4%)", y = "PC3 (5.2%)")
ind_pc13
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




