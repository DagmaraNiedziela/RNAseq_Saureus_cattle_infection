# NF-kB ####

# Easiest - take a list of genes involved in KEGG NF-kB signaling pathway 
# Merge it with logs for each strain, full joins! 
library(tidyverse)

library(readxl)
nfkB_genes <- read_excel("NFkB_genes.xlsx", sheet = 1)
nfkB_genes

library(AnnotationDbi) 
library(OrganismDbi)
library(org.Hs.eg.db) 
library("org.Bt.eg.db")

nrow(nfkB_genes) #100

# Bovine - with new logs ####
nfkB_genes_bovine <- nfkB_genes

nfkB_genes_bovine$rowname <- mapIds(org.Bt.eg.db,
                             keys=nfkB_genes_bovine$Alias,
                             column="ENSEMBL",
                             keytype="ALIAS",
                             multiVals="first")

# From code file 5_pathview 
temp <- rownames_to_column(logs_MOK023) 
temp2 <- rownames_to_column(logs_MOK124) 

nfkB_genes_bovine_logs <- left_join(nfkB_genes_bovine, temp, by = "rowname")
nfkB_genes_bovine_logs <- left_join(nfkB_genes_bovine_logs, temp2, by = "rowname", suffix = c("_MOK023", "_MOK124")) 
nrow(nfkB_genes_bovine_logs) #100

View(nfkB_genes_bovine_logs)
write.csv(nfkB_genes_bovine_logs, "nfkB_genes_bovine_logs.csv") 

# ** ECM receptor interaction - bovine ##### 

library(tidyverse)

library(readxl)
ECM_genes <- read_excel("NFkB_genes.xlsx", sheet = 3)
ECM_genes 

ECM_genes_bovine <- ECM_genes

ECM_genes_bovine$rowname <- mapIds(org.Bt.eg.db,
                                     keys=ECM_genes_bovine$Alias,
                                     column="ENSEMBL",
                                     keytype="ALIAS",
                                     multiVals="first")

# From code file 5_pathview 
temp <- rownames_to_column(logs_MOK023) 
temp2 <- rownames_to_column(logs_MOK124) 

ECM_genes_bovine_logs <- left_join(ECM_genes_bovine, temp, by = "rowname")
ECM_genes_bovine_logs <- left_join(ECM_genes_bovine_logs, temp2, by = "rowname", suffix = c("_MOK023", "_MOK124")) 
nrow(ECM_genes_bovine_logs) # 86 

View(ECM_genes_bovine_logs)
write.csv(ECM_genes_bovine_logs, "ECM_genes_bovine_logs.csv") 


# ** Tight junction - bovine ##### 

library(tidyverse)

library(readxl)
tight_genes <- read_excel("NFkB_genes.xlsx", sheet = 4)
tight_genes 

tight_genes_bovine <- tight_genes

tight_genes_bovine$rowname <- mapIds(org.Bt.eg.db,
                                   keys=tight_genes_bovine$Alias,
                                   column="ENSEMBL",
                                   keytype="ALIAS",
                                   multiVals="first")

# From code file 5_pathview 
temp <- rownames_to_column(logs_MOK023) 
temp2 <- rownames_to_column(logs_MOK124) 

tight_genes_bovine_logs <- left_join(tight_genes_bovine, temp, by = "rowname")
tight_genes_bovine_logs <- left_join(tight_genes_bovine_logs, temp2, by = "rowname", suffix = c("_MOK023", "_MOK124")) 
nrow(tight_genes_bovine_logs) #170

View(tight_genes_bovine_logs)
write.csv(tight_genes_bovine_logs, "tight_genes_bovine_logs.csv") 


# ** TNF signalling - bovine ##### 

library(tidyverse)

library(readxl)
TNF_genes <- read_excel("NFkB_genes.xlsx", sheet = 5)
TNF_genes 

TNF_genes_bovine <- TNF_genes

TNF_genes_bovine$rowname <- mapIds(org.Bt.eg.db,
                                     keys=TNF_genes_bovine$Alias,
                                     column="ENSEMBL",
                                     keytype="ALIAS",
                                     multiVals="first")

# From code file 5_pathview 
temp <- rownames_to_column(logs_MOK023) 
temp2 <- rownames_to_column(logs_MOK124) 

TNF_genes_bovine_logs <- left_join(TNF_genes_bovine, temp, by = "rowname")
TNF_genes_bovine_logs <- left_join(TNF_genes_bovine_logs, temp2, by = "rowname", suffix = c("_MOK023", "_MOK124")) 
nrow(TNF_genes_bovine_logs) #110

View(TNF_genes_bovine_logs)
write.csv(TNF_genes_bovine_logs, "TNF_genes_bovine_logs.csv") 


# ** Hippo signalling - bovine ##### 

library(tidyverse)

library(readxl)
Hippo_genes <- read_excel("NFkB_genes.xlsx", sheet = 6)
Hippo_genes 

Hippo_genes_bovine <- Hippo_genes

Hippo_genes_bovine$rowname <- mapIds(org.Bt.eg.db,
                                     keys=Hippo_genes_bovine$Alias,
                                     column="ENSEMBL",
                                     keytype="ALIAS",
                                     multiVals="first")

# From code file 5_logs 
temp <- rownames_to_column(logs_MOK023) 
temp2 <- rownames_to_column(logs_MOK124) 

Hippo_genes_bovine_logs <- left_join(Hippo_genes_bovine, temp, by = "rowname")
Hippo_genes_bovine_logs <- left_join(Hippo_genes_bovine_logs, temp2, by = "rowname", suffix = c("_MOK023", "_MOK124")) 
nrow(Hippo_genes_bovine_logs) #154

View(Hippo_genes_bovine_logs)
write.csv(Hippo_genes_bovine_logs, "Hippo_genes_bovine_logs.csv") 

# Heatmaps of genes in the pathways #####

View(nfkB_genes_bovine_logs)
View(focal_genes_bovine_logs) 
View(ECM_genes_bovine_logs) 
View(tight_genes_bovine_logs) 
View(TNF_genes_bovine_logs) 
View(Hippo_genes_bovine_logs)  

# First remove unnecessary columns 
# Then remove genes that have no expression to redue rows of the heatmap 
# Then make heatmap! 

# ** NF-kb heatmap #### 
# Remove unnecessary columns 
nfkB_genes_bovine_logs_path <- nfkB_genes_bovine_logs %>% select(-`ENTREZ ID`, -rowname, -`Gene name`) 
head(nfkB_genes_bovine_logs_path) 

# Get row names - Alias 
nfkB_genes_bovine_logs_path <- as.data.frame(nfkB_genes_bovine_logs_path)
rownames(nfkB_genes_bovine_logs_path) <- nfkB_genes_bovine_logs_path$Alias 
nfkB_genes_bovine_logs_path <- nfkB_genes_bovine_logs_path[,-1]

library(pheatmap) 

# Set column names 
colnames(nfkB_genes_bovine_logs_path) <- c("MOK023 24 hpi", "MOK023 48 hpi", "MOK023 72 hpi", "MOK023 168 hpi", "MOK124 24 hpi", "MOK124 48 hpi", "MOK124 72 hpi", "MOK124 168 hpi") 
View(nfkB_genes_bovine_logs_path)

# Remove rows that are all NA
nrow(nfkB_genes_bovine_logs_path) # 100, after removing rows with all NAs 49
nfkB_genes_bovine_logs_path <- nfkB_genes_bovine_logs_path[rowSums(is.na(nfkB_genes_bovine_logs_path)) != ncol(nfkB_genes_bovine_logs_path), ]

# Heatmap 
nfkb_genes_heatmap <- pheatmap(nfkB_genes_bovine_logs_path, na_col = "grey", 
                               cluster_cols = FALSE, cluster_rows = FALSE, 
                               cellwidth = 15, fontsize = 8) 
? pheatmap
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
# Use default color palette 

tiff("NFKB_heatmap.tiff", height = 20, width = 24, units = 'cm', 
     compression = "lzw", res = 600) 
pathway_genes_heatmap
dev.off() 

# Try to save with Cairo just in case 
install.packages("Cairo")
library(Cairo)
Cairo(file="NFKB_heatmap.png", 
      type="png",
      units="cm", 
      width=24, 
      height=20, 
      pointsize=60, 
      dpi=600)
pathway_genes_heatmap
dev.off() 
# font is very narrow even when I increase point size 
# but resolution (dpi) is good with 600) 

# Arrange NF-kB genes and TNF genes 
library(gridExtra)
library("grid")
library("lattice")
venn.plot<-grid.arrange(gTree(children=nfkb_genes_heatmap),
                        gTree(children=TNF_genes_heatmap),
                        ncol = 2,
                        widths = c(1,1),
                        heights = c(1,1))

venn.plot<-grid.arrange(nfkb_genes_heatmap,
                        TNF_genes_heatmap,
                        ncol = 2,
                        nfrow = 1,
                        widths = c(1,1),
                        heights = c(1,1)) 

# Pheatmap column 4 is a gtable that can be used to arrange the images
g <- grid.arrange(nfkb_genes_heatmap[[4]], TNF_genes_heatmap[[4]],ncol=2, top=c("A", "B"))
g 
# Can't figure out how to add A and B as labels 
?grid.arrange 

install.packages("ggpubr")
library(ggpubr)
g <- ggarrange(nfkb_genes_heatmap[[4]], TNF_genes_heatmap[[4]],ncol=2, labels=c("A", "B"))
#THIS! And Cairo 

library(Cairo)
Cairo(file="NFkB_TNF_AB.png", 
      type="png",
      units="cm", 
      width=24, 
      height=20, 
      pointsize=60, 
      dpi=600)
g
dev.off() 

# ** ECM heatmap #### 

View(ECM_genes_bovine_logs) 

# Remove unnecessary columns 
ECM_genes_bovine_logs_path <- ECM_genes_bovine_logs %>% select(-`ENTREZ ID`, -rowname, -`Gene name`) 
head(ECM_genes_bovine_logs_path) 

# Get row names - Alias 
ECM_genes_bovine_logs_path <- as.data.frame(ECM_genes_bovine_logs_path)
rownames(ECM_genes_bovine_logs_path) <- ECM_genes_bovine_logs_path$Alias 
ECM_genes_bovine_logs_path <- ECM_genes_bovine_logs_path[,-1]

library(pheatmap) 

# Set column names 
colnames(ECM_genes_bovine_logs_path) <- c("MOK023 24 hpi", "MOK023 48 hpi", "MOK023 72 hpi", "MOK023 168 hpi", "MOK124 24 hpi", "MOK124 48 hpi", "MOK124 72 hpi", "MOK124 168 hpi") 
View(ECM_genes_bovine_logs_path)

# Remove rows that are all NA
nrow(ECM_genes_bovine_logs_path) # 86, after removing rows with all NAs 30
ECM_genes_bovine_logs_path <- ECM_genes_bovine_logs_path[rowSums(is.na(ECM_genes_bovine_logs_path)) != ncol(ECM_genes_bovine_logs_path), ]

# Heatmap 
ECM_genes_heatmap <- pheatmap(ECM_genes_bovine_logs_path, na_col = "grey", 
                                  cluster_cols = FALSE, cluster_rows = FALSE, 
                                  cellwidth = 15, fontsize = 8) 
? pheatmap
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
# Use default color palette 

tiff("ECM_heatmap.tiff", height = 20, width = 24, units = 'cm', 
     compression = "lzw", res = 600) 
ECM_genes_heatmap
dev.off() 
# Might have to define cell height for all heatmaps 

# Try to save with Cairo just in case 
install.packages("Cairo")
library(Cairo)
Cairo(file="ECM_heatmap.png", 
      type="png",
      units="cm", 
      width=24, 
      height=20, 
      pointsize=60, 
      dpi=600)
pathway_genes_heatmap
dev.off() 

# ** Tight j heatmap #### 

View(tight_genes_bovine_logs)  

# Remove unnecessary columns 
tight_genes_bovine_logs_path <- tight_genes_bovine_logs %>% select(-`ENTREZ ID`, -rowname, -`Gene name`) 
head(tight_genes_bovine_logs_path) 

# Get row names - Alias 
tight_genes_bovine_logs_path <- as.data.frame(tight_genes_bovine_logs_path)
rownames(tight_genes_bovine_logs_path) <- tight_genes_bovine_logs_path$Alias 
tight_genes_bovine_logs_path <- tight_genes_bovine_logs_path[,-1]

library(pheatmap) 

# Set column names 
colnames(tight_genes_bovine_logs_path) <- c("MOK023 24 hpi", "MOK023 48 hpi", "MOK023 72 hpi", "MOK023 168 hpi", "MOK124 24 hpi", "MOK124 48 hpi", "MOK124 72 hpi", "MOK124 168 hpi") 
View(tight_genes_bovine_logs_path)

# Remove rows that are all NA
nrow(tight_genes_bovine_logs_path) # 170, after removing rows with all NAs 64
tight_genes_bovine_logs_path <- tight_genes_bovine_logs_path[rowSums(is.na(tight_genes_bovine_logs_path)) != ncol(tight_genes_bovine_logs_path), ]

# Heatmap 
tight_genes_heatmap <- pheatmap(tight_genes_bovine_logs_path, na_col = "grey", 
                                  cluster_cols = FALSE, cluster_rows = FALSE, 
                                  cellwidth = 15, fontsize = 8) 
? pheatmap
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
# Use default color palette 

tiff("tight_heatmap.tiff", height = 24, width = 24, units = 'cm', 
     compression = "lzw", res = 600) 
tight_genes_heatmap
dev.off() 

# Try to save with Cairo just in case 
install.packages("Cairo")
library(Cairo)
Cairo(file="tight_heatmap.png", 
      type="png",
      units="cm", 
      width=24, 
      height=20, 
      pointsize=60, 
      dpi=600)
pathway_genes_heatmap
dev.off() 

# ** TNF heatmap #### 

View(TNF_genes_bovine_logs) 

# Remove unnecessary columns 
TNF_genes_bovine_logs_path <- TNF_genes_bovine_logs %>% select(-`ENTREZ ID`, -rowname, -`Gene name`) 
head(TNF_genes_bovine_logs_path) 

# Get row names - Alias 
TNF_genes_bovine_logs_path <- as.data.frame(TNF_genes_bovine_logs_path)
rownames(TNF_genes_bovine_logs_path) <- TNF_genes_bovine_logs_path$Alias 
TNF_genes_bovine_logs_path <- TNF_genes_bovine_logs_path[,-1]

library(pheatmap) 

# Set column names 
colnames(TNF_genes_bovine_logs_path) <- c("MOK023 24 hpi", "MOK023 48 hpi", "MOK023 72 hpi", "MOK023 168 hpi", "MOK124 24 hpi", "MOK124 48 hpi", "MOK124 72 hpi", "MOK124 168 hpi") 
View(TNF_genes_bovine_logs_path)

# Remove rows that are all NA
nrow(TNF_genes_bovine_logs_path) # 110, after removing rows with all NAs 57
TNF_genes_bovine_logs_path <- TNF_genes_bovine_logs_path[rowSums(is.na(TNF_genes_bovine_logs_path)) != ncol(TNF_genes_bovine_logs_path), ]

# Heatmap 
TNF_genes_heatmap <- pheatmap(TNF_genes_bovine_logs_path, na_col = "grey", 
                              cluster_cols = FALSE, cluster_rows = FALSE, 
                              cellwidth = 15, fontsize = 8) 
? pheatmap
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
# Use default color palette 

tiff("TNF_heatmap.tiff", height = 20, width = 24, units = 'cm', 
     compression = "lzw", res = 600) 
pathway_genes_heatmap
dev.off() 

# Try to save with Cairo just in case 
install.packages("Cairo")
library(Cairo)
Cairo(file="TNF_heatmap.png", 
      type="png",
      units="cm", 
      width=24, 
      height=20, 
      pointsize=60, 
      dpi=600)
pathway_genes_heatmap
dev.off() 

# ** Hippo heatmap #### 

View(Hippo_genes_bovine_logs)   

# Remove unnecessary columns 
Hippo_genes_bovine_logs_path <- Hippo_genes_bovine_logs %>% select(-`ENTREZ ID`, -rowname, -`Gene name`) 
head(Hippo_genes_bovine_logs_path) 

# Get row names - Alias 
Hippo_genes_bovine_logs_path <- as.data.frame(Hippo_genes_bovine_logs_path)
rownames(Hippo_genes_bovine_logs_path) <- Hippo_genes_bovine_logs_path$Alias 
Hippo_genes_bovine_logs_path <- Hippo_genes_bovine_logs_path[,-1]

library(pheatmap) 

# Set column names 
colnames(Hippo_genes_bovine_logs_path) <- c("MOK023 24 hpi", "MOK023 48 hpi", "MOK023 72 hpi", "MOK023 168 hpi", "MOK124 24 hpi", "MOK124 48 hpi", "MOK124 72 hpi", "MOK124 168 hpi") 
View(Hippo_genes_bovine_logs_path)

# Remove rows that are all NA
nrow(Hippo_genes_bovine_logs_path) # 154, after removing rows with all NAs 47
Hippo_genes_bovine_logs_path <- Hippo_genes_bovine_logs_path[rowSums(is.na(Hippo_genes_bovine_logs_path)) != ncol(Hippo_genes_bovine_logs_path), ]

# Heatmap 
hippo_genes_heatmap <- pheatmap(Hippo_genes_bovine_logs_path, na_col = "grey", 
                                  cluster_cols = FALSE, cluster_rows = FALSE, 
                                  cellwidth = 15, fontsize = 8) 
? pheatmap
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
# Use default color palette 

tiff("Hippo_heatmap.tiff", height = 20, width = 24, units = 'cm', 
     compression = "lzw", res = 600) 
hippo_genes_heatmap
dev.off() 

# Try to save with Cairo just in case 
install.packages("Cairo")
library(Cairo)
Cairo(file="Hippo_heatmap.png", 
      type="png",
      units="cm", 
      width=24, 
      height=20, 
      pointsize=60, 
      dpi=600)
pathway_genes_heatmap
dev.off() 

# Combined for Hippo, ECM and tight junction heatmaps #### 

library(ggpubr)
h <- ggarrange(hippo_genes_heatmap[[4]], ECM_genes_heatmap[[4]],tight_genes_heatmap[[4]],ncol=3, labels=c("A", "B", "C"))
#THIS! And Cairo 
h

library(Cairo)
Cairo(file="Hippo_ECM_TJ_ABC.png", 
      type="png",
      units="cm", 
      width=24, 
      height=20, 
      pointsize=60, 
      dpi=600)
h
dev.off() 




