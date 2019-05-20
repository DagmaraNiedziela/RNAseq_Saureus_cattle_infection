

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

# Remove non protein coding genes - this chunk of code was made by Paul Cormican
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

dds <- estimateSizeFactors(dds) #what is this for? dds has 6 rows? 

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

# sample distances heat map 
sampleDists <- dist(t(assay(vsd)))
library("pheatmap")
library("RColorBrewer") 
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$Cow_ID, vsd$HPI, vsd$Group, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# generate heat map of top 20 genes with highest normalised counts 
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df2 <- as.data.frame(colData(dds)[,c("Group","Cow_ID")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df2) 

#PCA of samples 
plotPCA(vsd, intgroup=c("Group", "HPI")) # The simplest PCA 

# Default, ntop = 500, top genes used for PCA, selected by highest variance
DESeq2:::plotPCA.DESeqTransform # gets source code 
?prcomp # This is the function used to do the PCA on data, from "stats" package

# MY PLOT! Strain is shape, HPI is color 
percentVar <- round(100 * attr(PCAdata, "percentVar")) 
library("ggplot2")
sp <- ggplot(PCAdata, aes(x = PC1, y = PC2, color = HPI, shape = Group)) +
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed() 
sp 
myPCA <- sp + scale_color_brewer(palette="Set1", direction = -1) + 
  scale_fill_brewer(palette = "Set1", direction = -1) + 
  theme(text = element_text(size = 14, family = "Calibri")) 
myPCA

ggsave("myPCA.jpeg") 

#Label PCA with Cow IDs - an option to check animal distribution on PCA 
PCAdata <- plotPCA(vsd, intgroup=c("Group", "HPI", "Cow_ID"), returnData = TRUE) 
head(PCAdata)
#for cow in PCAdata you need to set Cow_ID to factor
PCAdata$Cow_ID <- factor(PCAdata$Cow_ID) 
percentVar <- round(100 * attr(PCAdata, "percentVar")) 

library("ggplot2") 
sp <- ggplot(PCAdata, aes(x = PC1, y = PC2, color = Group, shape = HPI, label=RIN)) +
  geom_point(size =3) + geom_text(hjust=-0.3,vjust=0) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() 
sp 
sp + scale_color_brewer(palette="Set1") 

ggsave("PCA.jpg", scale = 1.5, dpi = 300) 

#plot PCA with each time point as separate grid - this was not used because the progression of my samples with time was interesting
# But is a useful code and also nice to see separation of groups at each time point 
# this one fixes the coordinates - makes these roughly square 
PCAplot2 <- ggplot(data = PCAdata) + 
  geom_point(mapping = aes(x = PC1, y = PC2, color = Group)) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  facet_wrap(~ HPI, nrow = 2) + 
  coord_fixed() + 
  theme(text = element_text(size = 14, family = "Calibri")) 

#here I unfix the scales in each square 
PCAplot3 <- ggplot(data = PCAdata) + 
  geom_point(mapping = aes(x = PC1, y = PC2, color = Group)) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  facet_wrap(~ HPI, nrow = 2, scales = "free") +
  theme(text = element_text(size = 14, family = "Calibri"))

# lovely printable version 
# for PCA3 - plot with free scales, I needed to make the output files of different scales 
library(Cairo)
ggsave("PCA_plot_faceted_perHPI3.pdf", # Name of the file to be saved
       plot      = PCAplot3, # The variable where you stored the plot
       device    = cairo_pdf, # We use Cairo here
       limitsize = FALSE,
       dpi       = 300, # Dots per inch, 300 is agood resolution for printing
       height    = 6,
       width     = 8,
       units = "in") # Units in inches

# a png file 
CairoPNG(filename = "PCA_plot_faceted_perHPI3.png", width = 800, height = 600)
PCAplot3
dev.off()

##### GENE CLUSTERING - WORK ON VST DATA ####

# Plot top 20 most variable genes 
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20) 
topVarGenes 

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Group","HPI")]) # anno is coldata identifiers, sample, strain, time 
pheatmap(mat, annotation_col = anno) 

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

#remove non protein coding genes - Paul's code 
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
newdata<-newdata[,-26] #### remove the last column (the one that says ID type – because ncol in countdata needs to be the same as nrow in coldata
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
  
###### NORMALISE AND CLUSTER PLOTS - these are not used for data exploration but useful to generate ####
vsd23 <- vst(dds23, blind=FALSE)
colnames(vsd23)

#visualse if data has been normalised compared to a regular log transformation 
library("dplyr")
library("ggplot2")

df <- bind_rows(
  as_data_frame(log2(counts(dds23, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd23)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

#sample distances and heat map 
sampleDists <- dist(t(assay(vsd23)))
library("pheatmap")
library("RColorBrewer") 
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd23$Cow_ID, vsd23$HPI, vsd23$Group, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#generate heat map of top 20 genes with highest normalised counts 

library("pheatmap")
select <- order(rowMeans(counts(dds23,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df2 <- as.data.frame(colData(dds23)[,c("HPI","Cow_ID")])
pheatmap(assay(vsd23)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df2) 
#for grouping samples use also cluster_cols=TRUE 

# PCA 
#Label PCA with Cow IDs and time as color - best plot version
PCAdata <- plotPCA(vsd23, intgroup=c("Group", "HPI", "Cow_ID"), returnData = TRUE) 
head(PCAdata)
#for cow in PCAdata you need to set Cow ID to factor 
PCAdata$Cow_ID <- factor(PCAdata$Cow_ID)

percentVar <- round(100 * attr(PCAdata, "percentVar")) 
library("ggplot2")
sp <- ggplot(PCAdata, aes(x = PC1, y = PC2, color = HPI, shape = Group, label=Cow_ID)) +
  geom_point(size =3) + geom_text(hjust=-0.3,vjust=0) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() 
sp 
sp + scale_color_brewer(palette="Set1") #without sp its regular plot, this changes the color scheme 

##### GENE CLUSTERING - WORK ON VST DATA ####
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd23)), decreasing = TRUE), 20) 
topVarGenes #I get row numbers, cant get 
anno #coldata identifiers, sample, strain, time
mat[,0] #this gives me column 0 - gene names 

mat  <- assay(vsd23)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd23)[, c("Group","HPI")])
pheatmap(mat, annotation_col = anno) 

#####RUN DE GENE ANALYSIS - MOK023 #####

dds23 <- DESeq(dds23) #runs analysis 

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

# Export results per time point (Wald test contrasts) control should be second in code (B vs A, 24 vs 0) 
res23_24h <- results(dds23, contrast = c("HPI", "24", "0")) 
res23_48h <- results(dds23, contrast = c("HPI", "48", "0"))
res23_72h <- results(dds23, contrast = c("HPI", "72", "0"))
res23_168h <- results(dds23, contrast = c("HPI", "168", "0")) 

sum(res23_24h$padj < 0.05, na.rm=TRUE) # 7433
sum(res23_48h$padj < 0.05, na.rm=TRUE) #8591
sum(res23_72h$padj < 0.05, na.rm=TRUE) #8135
sum(res23_168h$padj < 0.05, na.rm=TRUE) #7825 

resultsNames(dds23)
# [1] "Intercept"    "HPI_24_vs_0"  "HPI_48_vs_0"  "HPI_72_vs_0" 
#[5] "HPI_168_vs_0"  

# NEW - results with fold change above 2 - lfcThreshold - these are used for my analysis #####
res23_24h_LF <- results(dds23, lfcThreshold=1, contrast = c("HPI", "24", "0"))
res23_48h_LF <- results(dds23, lfcThreshold=1, contrast = c("HPI", "48", "0"))
res23_72h_LF <- results(dds23, lfcThreshold=1, contrast = c("HPI", "72", "0"))
res23_168h_LF <- results(dds23, lfcThreshold=1, contrast = c("HPI", "168", "0")) 

sum(res23_24h_LF$padj < 0.05, na.rm=TRUE) # 1278
sum(res23_48h_LF$padj < 0.05, na.rm=TRUE) # 2248 
sum(res23_72h_LF$padj < 0.05, na.rm=TRUE) # 1986
sum(res23_168h_LF$padj < 0.05, na.rm=TRUE) # 1750 

#PLOT GENES #####
topGene <- rownames(res23)[which.min(res$padj)]
plotCounts(dds23, gene = topGene, intgroup=c("HPI")) 

#with gg plot to make this prettier, Group as color
geneCounts <- plotCounts(dds23, gene = topGene, intgroup = c("HPI", "Group", "Cow_ID"),
                         returnData = TRUE)
# and lines connecting cows 
ggplot(geneCounts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + ggtitle(topGene)

#Top 200 genes - pdf 
#Top 200 genes by adjusted p value at 24hpi 
top200Genes_24h <- res23_24h_LF[order(res23_24h_LF$padj),][1:200,]  # default sorting is ascending - smallest p value first 

pdf("top200Geneplots_ggplot_MOK023_24H.pdf")
for (i in row.names(top200Genes_24h)) {
  geneCounts <- plotCounts(dds23, gene = i, intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE) 
  geneCounts$Cow_ID <- as.factor(geneCounts$Cow_ID) 
  plot <- ggplot(geneCounts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) +
    scale_y_log10() + geom_point(size = 3) + geom_line() 
  print(plot + ggtitle(i))
}
dev.off() 

##### ANNOTATE GENES #####
library("AnnotationDbi")
library("OrganismDbi") 
library("org.Bt.eg.db")
columns(org.Bt.eg.db)
keytypes(org.Bt.eg.db)

head(keys(org.Bt.eg.db, keytype="ENSEMBLPROT")) #to check what those names look like 

# Normal genes 
res23_24h$Alias <- mapIds(org.Bt.eg.db,
                             keys=row.names(res23_24h),
                             column="ALIAS",
                             keytype="ENSEMBL",
                             multiVals="first")
res23_24h$Gene_Name <- mapIds(org.Bt.eg.db,
                                 keys=row.names(res23_24h),
                                 column="GENENAME",
                                 keytype="ENSEMBL",
                                 multiVals="first") #column previously called Name

res23_48h$Alias <- mapIds(org.Bt.eg.db,
                             keys=row.names(res23_48h),
                             column="ALIAS",
                             keytype="ENSEMBL",
                             multiVals="first")
res23_48h$Gene_Name <- mapIds(org.Bt.eg.db,
                                 keys=row.names(res23_48h),
                                 column="GENENAME",
                                 keytype="ENSEMBL",
                                 multiVals="first") 

res23_72h$Alias <- mapIds(org.Bt.eg.db,
                             keys=row.names(res23_72h),
                             column="ALIAS",
                             keytype="ENSEMBL",
                             multiVals="first")
res23_72h$Gene_Name <- mapIds(org.Bt.eg.db,
                                 keys=row.names(res23_72h),
                                 column="GENENAME",
                                 keytype="ENSEMBL",
                                 multiVals="first") 

res23_168h$Alias <- mapIds(org.Bt.eg.db,
                              keys=row.names(res23_168h),
                              column="ALIAS",
                              keytype="ENSEMBL",
                              multiVals="first")
res23_168h$Gene_Name <- mapIds(org.Bt.eg.db,
                                  keys=row.names(res23_168h),
                                  column="GENENAME",
                                  keytype="ENSEMBL",
                                  multiVals="first")


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
#all genes that are significant 
res23_24h_Sig_DFS <- as.data.frame(subset(res23_24h, padj < 0.05))
write.csv(res23_24h_Sig_DFS, file = "DE_genes_MOK023_24HPI.csv") 

res23_48h_Sig_DFS <- as.data.frame(subset(res23_48h, padj < 0.05))
write.csv(res23_48h_Sig_DFS, file = "DE_genes_MOK023_48HPI.csv")  

res23_72h_Sig_DFS <- as.data.frame(subset(res23_72h, padj < 0.05))
write.csv(res23_72h_Sig_DFS, file = "DE_genes_MOK023_72HPI.csv")  

res23_168h_Sig_DFS <- as.data.frame(subset(res23_168h, padj < 0.05))
write.csv(res23_168h_Sig_DFS, file = "DE_genes_MOK023_168HPI.csv") 

#export tables with all genes, not just P < 0.05 
res23_24h_DFS <- as.data.frame(res23_24h)
write.csv(res23_24h_DFS, file = "All_DE_genes_MOK023_24HPI.csv") 

res23_48h_DFS <- as.data.frame(res23_48h)
write.csv(res23_48h_DFS, file = "All_DE_genes_MOK023_48HPI.csv")  

res23_72h_DFS <- as.data.frame(res23_72h)
write.csv(res23_72h_DFS, file = "All_DE_genes_MOK023_72HPI.csv")  

res23_168h_DFS <- as.data.frame(res23_168h)
write.csv(res23_168h_DFS, file = "All_DE_genes_MOK023_168HPI.csv") 

# New LFC threshold ones  
write.csv(res23_24h_LF, file = "LF_DE_genes_MOK023_24HPI.csv")
write.csv(res23_48h_LF, file = "LF_DE_genes_MOK023_48HPI.csv")
write.csv(res23_72h_LF, file = "LF_DE_genes_MOK023_72HPI.csv")
write.csv(res23_168h_LF, file = "LF_DE_genes_MOK023_168HPI.csv")

# New LFC threshold SIG ones 
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

###### NORMALISE AND CLUSTER PLOTS ####
vsd124 <- vst(dds124, blind=FALSE)
colnames(vsd124)

#visualse if data has been normalised compared to a regular log transformation 
library("dplyr")
library("ggplot2")

dds124 <- estimateSizeFactors(dds124) #what is this for? dds has 6 rows? 

df <- bind_rows(
  as_data_frame(log2(counts(dds124, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd124)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

#sample distances and heat map 
sampleDists <- dist(t(assay(vsd124)))
library("pheatmap")
library("RColorBrewer") 
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd124$Cow_ID, vsd124$HPI, vsd124$Group, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#generate heat map of top 20 highest normalised count genes 

library("pheatmap")
select <- order(rowMeans(counts(dds124,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df2 <- as.data.frame(colData(dds124)[,c("HPI","Cow_ID")])
pheatmap(assay(vsd124)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df2) 
#for grouping samples use also cluster_cols=TRUE 

#PCA of samples 
#Label PCA with Cow IDs and time as color - best plot version
PCAdata <- plotPCA(vsd124, intgroup=c("Group", "HPI", "Cow_ID"), returnData = TRUE) 
head(PCAdata)
#for cow in PCAdata you need to set Cow_Id to factor - done in dds124
PCAdata$Cow_ID <- factor(PCAdata$Cow_ID)

percentVar <- round(100 * attr(PCAdata, "percentVar")) 

library("ggplot2")
sp <- ggplot(PCAdata, aes(x = PC1, y = -PC2, color = HPI, shape = Group, label=Cow_ID)) +
  geom_point(size =3) + geom_text(hjust=-0.3,vjust=0) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() 
sp 
sp + scale_color_brewer(palette="Set1") #without sp its regular plot, this changes the color scheme 

##### GENE CLUSTERING - WORK ON VST DATA ####
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd124)), decreasing = TRUE), 20) 
topVarGenes #I get row numbers, cant get 
anno #coldata identifiers, sample, strain, time
mat[,0] #this gives me column 0 - gene names 

mat  <- assay(vsd124)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd124)[, c("Group","HPI")])
pheatmap(mat, annotation_col = anno) 

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

res124_24h <- results(dds124, contrast = c("HPI", "24", "0"))
res124_48h <- results(dds124, contrast = c("HPI", "48", "0"))
res124_72h <- results(dds124, contrast = c("HPI", "72", "0"))
res124_168h <- results(dds124, contrast = c("HPI", "168", "0")) 

sum(res124_24h$padj < 0.05, na.rm=TRUE) # 9499
sum(res124_48h$padj < 0.05, na.rm=TRUE) #7151
sum(res124_72h$padj < 0.05, na.rm=TRUE) #5711
sum(res124_168h$padj < 0.05, na.rm=TRUE) #6032

resultsNames(dds124)
# [1] "Intercept"    "HPI_24_vs_0"  "HPI_48_vs_0"  "HPI_72_vs_0" 
#[5] "HPI_168_vs_0"  

# LfcThreshold - genes with fold change higher than 2 - this was used for data analysis #####

res124_24h_LF <- results(dds124, lfcThreshold=1, contrast = c("HPI", "24", "0"))
res124_48h_LF <- results(dds124, lfcThreshold=1, contrast = c("HPI", "48", "0"))
res124_72h_LF <- results(dds124, lfcThreshold=1, contrast = c("HPI", "72", "0"))
res124_168h_LF <- results(dds124, lfcThreshold=1, contrast = c("HPI", "168", "0")) 

sum(res124_24h_LF$padj < 0.05, na.rm=TRUE) # 2293
sum(res124_48h_LF$padj < 0.05, na.rm=TRUE) # 1979 
sum(res124_72h_LF$padj < 0.05, na.rm=TRUE) # 1428
sum(res124_168h_LF$padj < 0.05, na.rm=TRUE) # 1544 

# PLOT GENES #####
topGene <- rownames(res124)[which.min(res$padj)]
plotCounts(dds124, gene = topGene, intgroup=c("HPI")) 

#with gg plot to make this prettier, Group as color 
geneCounts <- plotCounts(dds124, gene = topGene, intgroup = c("HPI", "Group", "Cow_ID"),
                         returnData = TRUE)
# and lines connecting cows
ggplot(geneCounts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + ggtitle(topGene)

#Top 200 genes by upregulation at 24hpi 
top200Genes_2_24h <- res124_24h_LF[order(res124_24h_LF$padj),][1:200,]  

pdf("top200Geneplots_ggplot_MOK124_24H.pdf")
for (i in row.names(top200Genes_2_24h)) {
  geneCounts <- plotCounts(dds124, gene = i, intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE) 
  geneCounts$Cow_ID <- as.factor(geneCounts$Cow_ID) 
  plot <- ggplot(geneCounts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) +
    scale_y_log10() + geom_point(size = 3) + geom_line() 
  print(plot + ggtitle(i))
}
dev.off() 

# Top 200 genesn In BOTH strains, sorted by padj in MOK124 24hpi 
pdf("top200Geneplots_ggplot_MOK124_24H_BOTH_groups.pdf")
for (i in row.names(top200Genes_2_24h)) {
  geneCounts <- plotCounts(dds, gene = i, intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE) 
  geneCounts$Cow_ID <- as.factor(geneCounts$Cow_ID) 
  plot <- ggplot(geneCounts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID, shape = Group)) +
    scale_y_log10() + geom_point(size =2) + geom_line() 
  print(plot + ggtitle(i))
}
dev.off() 

##### ANNOTATE GENES #####
library("AnnotationDbi")
library("OrganismDbi") 
library("org.Bt.eg.db")
columns(org.Bt.eg.db)
keytypes(org.Bt.eg.db)
row.names(res124_24h)

head(keys(org.Bt.eg.db, keytype="ENSEMBLPROT")) #to check what those names look like

res124_24h$Alias <- mapIds(org.Bt.eg.db,
                          keys=row.names(res124_24h),
                          column="ALIAS",
                          keytype="ENSEMBL",
                          multiVals="first")
res124_24h$Gene_Name <- mapIds(org.Bt.eg.db,
                              keys=row.names(res124_24h),
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="first") #column previously called Name

res124_48h$Alias <- mapIds(org.Bt.eg.db,
                          keys=row.names(res124_48h),
                          column="ALIAS",
                          keytype="ENSEMBL",
                          multiVals="first")
res124_48h$Gene_Name <- mapIds(org.Bt.eg.db,
                              keys=row.names(res124_48h),
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="first") 

res124_72h$Alias <- mapIds(org.Bt.eg.db,
                          keys=row.names(res124_72h),
                          column="ALIAS",
                          keytype="ENSEMBL",
                          multiVals="first")
res124_72h$Gene_Name <- mapIds(org.Bt.eg.db,
                              keys=row.names(res124_72h),
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="first") 

res124_168h$Alias <- mapIds(org.Bt.eg.db,
                           keys=row.names(res124_168h),
                           column="ALIAS",
                           keytype="ENSEMBL",
                           multiVals="first")
res124_168h$Gene_Name <- mapIds(org.Bt.eg.db,
                               keys=row.names(res124_168h),
                               column="GENENAME",
                               keytype="ENSEMBL",
                               multiVals="first") 

##### ANNOTATE GENES for LFC shrunk results #####
library("AnnotationDbi")
library("OrganismDbi") 
library("org.Bt.eg.db")

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
#all genes that are significant 
res124_24h_Sig_DFS <- as.data.frame(subset(res124_24h, padj < 0.05))
write.csv(res124_24h_Sig_DFS, file = "DE_genes_MOK124_24HPI.csv") 

res124_48h_Sig_DFS <- as.data.frame(subset(res124_48h, padj < 0.05))
write.csv(res124_48h_Sig_DFS, file = "DE_genes_MOK124_48HPI.csv")  

res124_72h_Sig_DFS <- as.data.frame(subset(res124_72h, padj < 0.05))
write.csv(res124_72h_Sig_DFS, file = "DE_genes_MOK124_72HPI.csv")  

res124_168h_Sig_DFS <- as.data.frame(subset(res124_168h, padj < 0.05))
write.csv(res124_168h_Sig_DFS, file = "DE_genes_MOK124_168HPI.csv") 

#export tables with all genes, not just P < 0.05 
res124_24h_DFS <- as.data.frame(res124_24h)
write.csv(res124_24h_DFS, file = "All_DE_genes_MOK124_24HPI.csv") 

res124_48h_DFS <- as.data.frame(res124_48h)
write.csv(res124_48h_DFS, file = "All_DE_genes_MOK124_48HPI.csv")  

res124_72h_DFS <- as.data.frame(res124_72h)
write.csv(res124_72h_DFS, file = "All_DE_genes_MOK124_72HPI.csv")  

res124_168h_DFS <- as.data.frame(res124_168h)
write.csv(res124_168h_DFS, file = "All_DE_genes_MOK124_168HPI.csv") 

#Export the ones with low log fold changes 
write.csv(res124_24h_LF, file = "LFC_DE_genes_MOK124_24HPI.csv")
write.csv(res124_48h_LF, file = "LFC_DE_genes_MOK124_48HPI.csv")
write.csv(res124_72h_LF, file = "LFC_DE_genes_MOK124_72HPI.csv")
write.csv(res124_168h_LF, file = "LFC_DE_genes_MOK124_168HPI.csv")

#Export the ones with LFC threshold and significant - only 24 hpi exported 
write.csv(res124_24h_LF_sig, file = "LF_sig_genes_MOK124_24HPI.csv")
write.csv(res124_48h_LF_sig, file = "LF_sig_genes_MOK124_48HPI.csv")
write.csv(res124_72h_LF_sig, file = "LF_sig_genes_MOK124_72HPI.csv")
write.csv(res124_168h_LF_sig, file = "LF_sig_genes_MOK124_168HPI.csv")

# Export a full LF sig DE gene file for both strains! This for supplementary files in paper/thesis ##### 
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
LF_check <- head(LF_sig_MOK023) 
LF_check 

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
