# To get DeSeq2 version 
metadata(dds)[["version"]] 
# 1.18.1 

#Load matrix - 10 groups ####
coldata <- read.csv("Samples_DeSeq2_Cowsremoved.csv")
countdata <- read.csv("all_HPI_GeneCounts_Cowsremoved.csv", row.names = "Genename") #row.names clause essential to get equal numbers of rows vs columns in coldata and countdata
nrow(coldata) # check if the two are equal - otherwise the dds matrix won't load 
ncol(countdata)

library("DESeq2")
coldata$HPI <- factor(coldata$HPI)
levels(coldata$HPI) 

# Load my dds matrix with both strains - this is used for data exploration 
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata, 
                              design = ~ HPI + Group + HPI:Group) 

rowData(dds) <- anno 							  
head(dds)

# Remove non protein coding genes 
biotype <- as.matrix(read.table(file="mart_export_biotype.2col.unique.txt", header=TRUE, row.names=1)) # file attached, import it  
te <- merge(countdata, biotype, by=0, all=FALSE) #merge my counts file with the protein coding assignment file 
head(te) 
rownames(te) <- te[,1] # this is making rownames  
te<-te[,-1] #removing column of ENSEML that’s in twice 
newdata <- te[ which(te$Gene_type=='protein_coding'), ] #remove genes that are not protein coding (miRNAs, pseudogenes, long ncRNAs etc)
head(newdata)
nrow(te) #24596 
nrow(newdata) #19981 
ncol(newdata)
newdata<-newdata[,-52] #### -52 is basically the last column (the one that says ID type – because ncol in countdata needs to be the same as nrow in coldata
countdata <- newdata #make newdata my countdata 

nrow(dds) # this to keep checking how many genes I have 
#[1] 19981 # protein coding genes only 
dds <- dds[ rowSums(counts(dds)) > 1, ] #removes genes with no counts
nrow(dds) 
#[1] 17334 

# Remove 616 48 hpi - wrong sample was sequenced 
colData(dds) # sample is a row name 
idx <- 49 # column number for 616 48 hpi sample 
dds <- dds[,-idx]
colnames(dds) 

###### NORMALISE AND CLUSTER PLOTS ####

# Normalise data with a vst function 
vsd <- vst(dds, blind=FALSE)

# visualse if data has been normalised compared to a regular log transformation 
library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds) 

#PCA of samples 
plotPCA(vsd, intgroup=c("Group", "HPI")) # The simplest PCA 
head(vsd)
View(vsd)
vsd

PCAdata <- plotPCA(vsd, intgroup=c("Group", "HPI", "Cow_ID"), returnData = TRUE) 

# My publication plot - strain (Group) is shape, HPI is color 
percentVar <- round(100 * attr(PCAdata, "percentVar")) 
library("ggplot2")
sp <- ggplot(PCAdata, aes(x = PC1, y = PC2, color = HPI, shape = Group)) +
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed() 
sp 

# ** My PCA ####
myPCA <- sp + scale_color_brewer(palette="Set1", direction = -1) + 
  scale_fill_brewer(palette = "Set1", direction = -1) + 
  theme(text = element_text(size = 14, family = "Calibri")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
myPCA

ggsave("myPCA.jpeg") 

# PCA and NF-kB normalised genes combined plot # 
# Normalised NF-kB plots are in code file 9_Plot_epithelial_cell_markers.r 

library(ggpubr)

ggarrange(myPCA, NFkB_plot, nrow = 2, ncol = 1, labels = c("A", "B"))
ggsave("PCA_NfkB_plot.jpeg", dpi = 600, width = 8, height = 8.5, units = "in")

#Load matrix for MOK023 #### 

coldata23 <- subset(coldata, coldata$Group=="MOK023") 
nrow(coldata23)
head(coldata23)

countdata23 <- read.csv("all_HPI_GeneCounts_MOK023.csv", row.names = "Genename") 
head(countdata23, 3)
ncol(countdata23) 

library("DESeq2") 

coldata23$HPI <- factor(coldata23$HPI) 
colData23$Cow_ID <- as.factor(colData23$Cow_ID)
levels(coldata23$HPI) # Check order of time points 

#remove non protein coding genes 
biotype <- as.matrix(read.table(file="mart_export_biotype.2col.unique.txt", header=TRUE, row.names=1)) # file attached, import it  
te <- merge(countdata23, biotype, by=0, all=FALSE) #merge my counts file with the protein coding assignment file
head(te)
rownames(te) <- te[,1] # this is making rownames  
te<-te[,-1] #removing column of ENSEML that’s in twice 
newdata <- te[ which(te$Gene_type=='protein_coding'), ] #remove genes that are not protein coding (miRNAs, pseudogenes, long ncRNAs etc)
head(newdata)
nrow(te) #24596
nrow(newdata) #19981 
ncol(newdata)
newdata<-newdata[,-26] #### remove the last column (the one that says ID type – because ncol in countdata needs to be the same as nrow in coldata}
countdata23 <- newdata #make newdata my countdata 

library("DESeq2")
dds23 <- DESeqDataSetFromMatrix(countData = countdata23,
                              colData = coldata23,
                              design = ~ HPI)

head(dds23)
colnames(dds23) 

nrow(dds23) 
#[1] 19981
dds23 <- dds23[ rowSums(counts(dds23)) > 1, ] #removes genes with no counts
nrow(dds23) 
#[1] 16928 

dds23 <- estimateSizeFactors(dds23) 

#####RUN DE GENE ANALYSIS - MOK023 #####

dds23 <- DESeq(dds23) #runs analysis 
# NOTE: This was run with an older version of DESeq2 and lfc shrinkage was performed as a default in this command
# Currently you can shrink log fold changes after performing this action 

res23 <- results(dds23)

mcols(res23, use.names = TRUE) #meaning of columns 
summary(res23) #gives info about outliers and low counts 

sum(!is.na(res23$pvalue)) #genes with a reported P value 
## 16906

sum(res23$padj < 0.1, na.rm=TRUE) #genes that pass adjusted P value threshold 
#8720 
sum(res23$padj < 0.05, na.rm=TRUE) #genes that pass adjusted P value threshold 
#7825 
res23_Sig <- subset(res23, padj < 0.1) #not done
nrow(res23_Sig) 

# Results with fold change above 2 - lfcThreshold - these are used for my analysis #####
res23_24h_LF <- results(dds23, lfcThreshold=1, contrast = c("HPI", "24", "0"))
res23_48h_LF <- results(dds23, lfcThreshold=1, contrast = c("HPI", "48", "0"))
res23_72h_LF <- results(dds23, lfcThreshold=1, contrast = c("HPI", "72", "0"))
res23_168h_LF <- results(dds23, lfcThreshold=1, contrast = c("HPI", "168", "0")) 

sum(res23_24h_LF$padj < 0.05, na.rm=TRUE) # 1278
sum(res23_48h_LF$padj < 0.05, na.rm=TRUE) # 2248 
sum(res23_72h_LF$padj < 0.05, na.rm=TRUE) # 1986
sum(res23_168h_LF$padj < 0.05, na.rm=TRUE) # 1750 

##### ANNOTATE GENES #####
library("AnnotationDbi")
library("OrganismDbi") 
library("org.Bt.eg.db")
columns(org.Bt.eg.db)
keytypes(org.Bt.eg.db)

head(keys(org.Bt.eg.db, keytype="ENSEMBLPROT")) #to check what those names look like 

# LF genes 
res23_24h_LF$Alias <- mapIds(org.Bt.eg.db,
                                keys=row.names(res23_24h_LF),
                                column="ALIAS",
                                keytype="ENSEMBL",
                                multiVals="first")
res23_24h_LF$Gene_Name <- mapIds(org.Bt.eg.db,
                              keys=row.names(res23_24h_LF),
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="first") #column previously called Name

res23_48h_LF$Alias <- mapIds(org.Bt.eg.db,
                          keys=row.names(res23_48h_LF),
                          column="ALIAS",
                          keytype="ENSEMBL",
                          multiVals="first")
res23_48h_LF$Gene_Name <- mapIds(org.Bt.eg.db,
                              keys=row.names(res23_48h_LF),
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="first") 

res23_72h_LF$Alias <- mapIds(org.Bt.eg.db,
                          keys=row.names(res23_72h_LF),
                          column="ALIAS",
                          keytype="ENSEMBL",
                          multiVals="first")
res23_72h_LF$Gene_Name <- mapIds(org.Bt.eg.db,
                              keys=row.names(res23_72h_LF),
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="first") 

res23_168h_LF$Alias <- mapIds(org.Bt.eg.db,
                          keys=row.names(res23_168h_LF),
                          column="ALIAS",
                          keytype="ENSEMBL",
                          multiVals="first")
res23_168h_LF$Gene_Name <- mapIds(org.Bt.eg.db,
                              keys=row.names(res23_168h_LF),
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="first")

res23_24h_LF_sig <- subset(res23_24h_LF, padj < 0.05)
res23_48h_LF_sig <- subset(res23_48h_LF, padj < 0.05)
res23_72h_LF_sig <- subset(res23_72h_LF, padj < 0.05)
res23_168h_LF_sig <- subset(res23_168h_LF, padj < 0.05) 
head(res23_24h_LF_sig)

res23_24h_LF_sigDF <- as.data.frame(res23_24h_LF_sig) 
res23_48h_LF_sigDF <- as.data.frame(res23_48h_LF_sig) 
res23_72h_LF_sigDF <- as.data.frame(res23_72h_LF_sig) 
res23_168h_LF_sigDF <- as.data.frame(res23_168h_LF_sig) 
head(res23_24h_LF_sigDF)

######EXPORT TABLES #####

# LFC threshold SIG DE genes 
write.csv(res23_24h_LF_sig, file = "LF_sig_genes_MOK023_24HPI.csv")
write.csv(res23_48h_LF_sig, file = "LF_sig_genes_MOK023_48HPI.csv")
write.csv(res23_72h_LF_sig, file = "LF_sig_genes_MOK023_72HPI.csv")
write.csv(res23_168h_LF_sig, file = "LF_sig_genes_MOK023_168HPI.csv") 


#Load matrix - MOK124 ####

coldata124 <- subset(coldata, coldata$Group=="MOK124")
#instead of giving it new files, make a subset of the first DDS (no interaction)
nrow(coldata124)
head(coldata124)
countdata124 <- read.csv("all_HPI_GeneCounts_MOK124.csv", row.names = "Genename") #fixed nrow doesnt equal ncol - row.names
ncol(countdata124)

coldata124$HPI <- factor(coldata124$HPI)
levels(coldata124$HPI) 
coldata124$Cow_ID <- factor(coldata124$Cow_ID)

#remove non protein coding genes - Paul's code 
biotype <- as.matrix(read.table(file="mart_export_biotype.2col.unique.txt", header=TRUE, row.names=1)) # file attached, import it  
te2 <- merge(countdata124, biotype, by=0, all=FALSE) #merge my counts file with the protein coding assignment file
head(te2)
rownames(te2) <- te2[,1] # this is making rownames  
te2 <- te2[,-1] #removing column of ENSEML that’s in twice 
newdata2 <- te2[ which(te2$Gene_type=='protein_coding'), ] #remove genes that are not protein coding (miRNAs, pseudogenes, long ncRNAs etc)
head(newdata2)
nrow(te2) #24596
nrow(newdata2) #19981 
ncol(newdata2) #27 
ncol(countdata124)
newdata2 <- newdata2[,-27] #### -52 is basically the last column (the one that says ID type – because ncol in countdata needs to be the same as nrow in coldata
countdata124 <- newdata2 #make newdata my countdata 

library("DESeq2")
dds124 <- DESeqDataSetFromMatrix(countData = countdata124,
                                 colData = coldata124,
                                 design = ~ HPI)

head(dds124)
colnames(dds124) 

#remove 616 48 hpi 
idx124 <- 24 
dds124 <- dds124[,-idx124]

nrow(dds124) 
#[1] 19981
dds124 <- dds124[ rowSums(counts(dds124)) > 1, ] #removes genes with no counts
nrow(dds124) 
#[1] 16884 

##### RUN DE GENE ANALYSIS - MOK0124 #####
dds124 <- DESeq(dds124) #runs analysis 

res124 <- results(dds124) 

mcols(res124, use.names = TRUE) #meaning of columns 
summary(res124) #gives info about outliers and low counts 

sum(!is.na(res124$pvalue)) #genes with a reported P value 
## 16836

sum(res124$padj < 0.1, na.rm=TRUE) #genes that pass adjusted P value threshold 
#7078
sum(res124$padj < 0.05, na.rm=TRUE) #genes that pass adjusted P value threshold 
#6032

# LfcThreshold - genes with fold change higher than 2 - this was used for data analysis #####

res124_24h_LF <- results(dds124, lfcThreshold=1, contrast = c("HPI", "24", "0"))
res124_48h_LF <- results(dds124, lfcThreshold=1, contrast = c("HPI", "48", "0"))
res124_72h_LF <- results(dds124, lfcThreshold=1, contrast = c("HPI", "72", "0"))
res124_168h_LF <- results(dds124, lfcThreshold=1, contrast = c("HPI", "168", "0")) 

sum(res124_24h_LF$padj < 0.05, na.rm=TRUE) # 2293
sum(res124_48h_LF$padj < 0.05, na.rm=TRUE) # 1979 
sum(res124_72h_LF$padj < 0.05, na.rm=TRUE) # 1428
sum(res124_168h_LF$padj < 0.05, na.rm=TRUE) # 1544 

##### ANNOTATE GENES #####
library("AnnotationDbi")
library("OrganismDbi") 
library("org.Bt.eg.db")
columns(org.Bt.eg.db)
keytypes(org.Bt.eg.db)
row.names(res124_24h)

head(keys(org.Bt.eg.db, keytype="ENSEMBLPROT")) #to check what those names look like

res124_24h_LF$Alias <- mapIds(org.Bt.eg.db,
                           keys=row.names(res124_24h_LF),
                           column="ALIAS",
                           keytype="ENSEMBL",
                           multiVals="first")
res124_24h_LF$Gene_Name <- mapIds(org.Bt.eg.db,
                               keys=row.names(res124_24h_LF),
                               column="GENENAME",
                               keytype="ENSEMBL",
                               multiVals="first") #column previously called Name


res124_48h_LF$Alias <- mapIds(org.Bt.eg.db,
                           keys=row.names(res124_48h_LF),
                           column="ALIAS",
                           keytype="ENSEMBL",
                           multiVals="first")
res124_48h_LF$Gene_Name <- mapIds(org.Bt.eg.db,
                               keys=row.names(res124_48h_LF),
                               column="GENENAME",
                               keytype="ENSEMBL",
                               multiVals="first") 

res124_72h_LF$Alias <- mapIds(org.Bt.eg.db,
                           keys=row.names(res124_72h_LF),
                           column="ALIAS",
                           keytype="ENSEMBL",
                           multiVals="first")
res124_72h_LF$Gene_Name <- mapIds(org.Bt.eg.db,
                               keys=row.names(res124_72h_LF),
                               column="GENENAME",
                               keytype="ENSEMBL",
                               multiVals="first") 

res124_168h_LF$Alias <- mapIds(org.Bt.eg.db,
                            keys=row.names(res124_168h_LF),
                            column="ALIAS",
                            keytype="ENSEMBL",
                            multiVals="first")
res124_168h_LF$Gene_Name <- mapIds(org.Bt.eg.db,
                                keys=row.names(res124_168h_LF),
                                column="GENENAME",
                                keytype="ENSEMBL",
                                multiVals="first")

res124_24h_LF_sig <- subset(res124_24h_LF, padj < 0.05)
res124_48h_LF_sig <- subset(res124_48h_LF, padj < 0.05)
res124_72h_LF_sig <- subset(res124_72h_LF, padj < 0.05)
res124_168h_LF_sig <- subset(res124_168h_LF, padj < 0.05) 

res124_24h_LF_sigDF <- as.data.frame(res124_24h_LF_sig) 
res124_48h_LF_sigDF <- as.data.frame(res124_48h_LF_sig) 
res124_72h_LF_sigDF <- as.data.frame(res124_72h_LF_sig) 
res124_168h_LF_sigDF <- as.data.frame(res124_168h_LF_sig)

######EXPORT TABLES #####

#Export the ones with LFC threshold and significant 
write.csv(res124_24h_LF_sig, file = "LF_sig_genes_MOK124_24HPI.csv")
write.csv(res124_48h_LF_sig, file = "LF_sig_genes_MOK124_48HPI.csv")
write.csv(res124_72h_LF_sig, file = "LF_sig_genes_MOK124_72HPI.csv")
write.csv(res124_168h_LF_sig, file = "LF_sig_genes_MOK124_168HPI.csv")

# List of all log fold changes and p values for every strain and time point - for supplementary file in paper ##### 
library(dplyr) 
View(res23_24h_LF_sig_DE)
res23_24h_LF_sig_DE <- res23_24h_LF_sigDF %>% select(log2FoldChange, padj) %>% tibble::rownames_to_column()
res23_48h_LF_sig_DE <- res23_48h_LF_sigDF %>% select(log2FoldChange, padj) %>% tibble::rownames_to_column()
res23_72h_LF_sig_DE <- res23_72h_LF_sigDF %>% select(log2FoldChange, padj) %>% tibble::rownames_to_column()
res23_168h_LF_sig_DE <- res23_168h_LF_sigDF %>% select(log2FoldChange, padj) %>% tibble::rownames_to_column()

res124_24h_LF_sig_DE <- res124_24h_LF_sigDF %>% select(log2FoldChange, padj) %>% tibble::rownames_to_column()
res124_48h_LF_sig_DE <- res124_48h_LF_sigDF %>% select(log2FoldChange, padj) %>% tibble::rownames_to_column()
res124_72h_LF_sig_DE <- res124_72h_LF_sigDF %>% select(log2FoldChange, padj) %>% tibble::rownames_to_column()
res124_168h_LF_sig_DE <- res124_168h_LF_sigDF %>% select(log2FoldChange, padj) %>% tibble::rownames_to_column()

LF_sig_MOK023 <- full_join(res23_24h_LF_sig_DE, res23_48h_LF_sig_DE, by = "rowname", suffix = c("_24hpi", "_48hpi")) 
LF_sig_MOK023 <- full_join(LF_sig_MOK023, res23_72h_LF_sig_DE, by = "rowname", suffix = c("_48hpi", "_72hpi")) 
LF_sig_MOK023 <- full_join(LF_sig_MOK023, res23_168h_LF_sig_DE, by = "rowname", suffix = c("_72hpi", "_168hpi"))

View(LF_sig_MOK023) 

LF_sig_MOK124 <- full_join(res124_24h_LF_sig_DE, res124_48h_LF_sig_DE, by = "rowname", suffix = c("_24hpi", "_48hpi")) 
LF_sig_MOK124 <- full_join(LF_sig_MOK124, res124_72h_LF_sig_DE, by = "rowname", suffix = c("_48hpi", "_72hpi")) 
LF_sig_MOK124 <- full_join(LF_sig_MOK124, res124_168h_LF_sig_DE, by = "rowname", suffix = c("_72hpi", "_168hpi"))

LF_sig_both <- full_join(LF_sig_MOK023, LF_sig_MOK124, by = "rowname", suffix = c("_MOK023", "_MOK124")) 
write.csv(LF_sig_both, "DE_genes_for_publication.csv") 

# At the end, annotate the full file 
library("AnnotationDbi")
library("OrganismDbi") 
library("org.Bt.eg.db")

# Normal genes 
LF_sig_both$Alias <- mapIds(org.Bt.eg.db,
                          keys=LF_sig_both$rowname,
                          column="ALIAS",
                          keytype="ENSEMBL",
                          multiVals="first")
LF_sig_both$Gene_Name <- mapIds(org.Bt.eg.db,
                              keys=LF_sig_both$rowname,
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="first") #column previously called Name

View(LF_sig_both)
write.csv(LF_sig_both, "DE_genes_for_publication_annotated.csv") 

# DE gene numbers are taken in "2_DE_gene_amounts_and_DE_gene_plot.R" ##### 
# Top 5 genes and milk protein genes are in "3_top5_genes_and_milk_genes.r" 
