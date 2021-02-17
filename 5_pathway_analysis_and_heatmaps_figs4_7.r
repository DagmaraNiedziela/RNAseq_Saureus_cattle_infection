# ORDERING DATA AND FILES I WORK ON ###### 
# From DE_gene_amounts_2 - lfcThreshold results 
res23_24h_LF_sig <- subset(res23_24h_LF, padj < 0.05)
res23_48h_LF_sig <- subset(res23_48h_LF, padj < 0.05)
res23_72h_LF_sig <- subset(res23_72h_LF, padj < 0.05)
res23_168h_LF_sig <- subset(res23_168h_LF, padj < 0.05)

res23_24h_LF_sig <- res23_24h_LF_sig[order (res23_24h_LF_sig$padj),] #order results by adjusted p value 
res23_48h_LF_sig <- res23_48h_LF_sig[order (res23_48h_LF_sig$padj),] 
res23_72h_LF_sig <- res23_72h_LF_sig[order (res23_72h_LF_sig$padj),] 
res23_168h_LF_sig <- res23_168h_LF_sig[order (res23_168h_LF_sig$padj),] 

# Only as a reminder - it is also in Gprofiler MOK124 

res124_24h_LF_sig <- res124_24h_LF_sig[order (res124_24h_LF_sig$padj),] #order results
res124_48h_LF_sig <- res124_48h_LF_sig[order (res124_48h_LF_sig$padj),] 
res124_72h_LF_sig <- res124_72h_LF_sig[order (res124_72h_LF_sig$padj),] 
res124_168h_LF_sig <- res124_168h_LF_sig[order (res124_168h_LF_sig$padj),] 

head(res23_24h_LF_sig) 

# Convert to human - MOK023 #### 

# Table of bovine to human gene orthologs exported from biomart 
species_conversion <- read.table(file="mart_export_bovine_to_human.txt", header=TRUE, sep = "\t")
head(species_conversion)
nrow(species_conversion) # 24068 
species_1to1 <- subset(species_conversion, Human.homology.type == "ortholog_one2one")
nrow(species_1to1) #19098 
head(species_1to1) 

# all MOK023 gene vectors to human ###### 

# y23_24 was created in Venn diagrams code file as a data frame, as shown in example below: 
res23_24h_LF_sigDF <- as.data.frame(res23_24h_LF_sig) 
res23_24h_LF_sigDF <- tibble::rownames_to_column(res23_24h_LF_sigDF) 
y23_24 <- res23_24h_LF_sigDF %>% dplyr::select(rowname) 

y23_24_human <- inner_join(species_1to1, y23_24) 

# Joining, by = "rowname"
#Warning message:
#  Column `rowname` joining factor and character vector, coercing into character vector 
nrow(y23_24) #1278 
nrow(y23_24_human) #1235
head(y23_24)
head(y23_24_human)  
write.csv(y23_24_human, "MOK023_24H_human.csv") 
names_23_24h_human <- y23_24_human$Human.gene.stable.ID
head(names_23_24h_human)
View(names_23_24h_human)

y23_48_human <- inner_join(species_1to1, y23_48) 
write.csv(y23_48_human, "MOK023_48H_human.csv") 
names_23_48h_human <- y23_48_human$Human.gene.stable.ID 

y23_72_human <- inner_join(species_1to1, y23_72) 
write.csv(y23_72_human, "MOK023_72H_human.csv") 
names_23_72h_human <- y23_72_human$Human.gene.stable.ID 

y23_168_human <- inner_join(species_1to1, y23_168) 
write.csv(y23_168_human, "MOK023_168H_human.csv") 
names_23_168h_human <- y23_168_human$Human.gene.stable.ID  


# run gProfiler on HUMAN gene orthologs ######

gprofilerresult23_24h_human <- gprofiler(names_23_24h_human, organism = "hsapiens", ordered_query = T, hier_filtering = "strong")
gprofilerresult23_48h_human <- gprofiler(names_23_48h_human, organism = "hsapiens", hier_filtering = "strong")
gprofilerresult23_72h_human <- gprofiler(names_23_72h_human, organism = "hsapiens", hier_filtering = "strong")
gprofilerresult23_168h_human <- gprofiler(names_23_168h_human, organism = "hsapiens", hier_filtering = "strong")

write.csv(gprofilerresult23_24h_human , file="Gprofiler_MOK023_24HPI_human.csv")
write.csv(gprofilerresult23_48h_human , file="Gprofiler_MOK023_48HPI_human.csv")
write.csv(gprofilerresult23_72h_human , file="Gprofiler_MOK023_72HPI_human.csv")
write.csv(gprofilerresult23_168h_human , file="Gprofiler_MOK023_168HPI_human.csv")  

# Make a file combining all KEGG lists - human pathways #### 

# MOK023 pathways - select only kegg pathways, name and p value 
# add_column is from tibble package 
library(tidyverse)
gprofiler23_24_kegg_human <- gprofilerresult23_24h_human  %>% 
  filter(domain == "keg") %>% select(term.name, overlap.size, p.value) %>% 
  arrange(p.value) 
gprofiler23_24_kegg_human <- gprofiler23_24_kegg_human %>% add_column(time =rep("MOK023_24hpi", nrow(gprofiler23_24_kegg_human)))
gprofiler23_48_kegg_human <- gprofilerresult23_48h_human  %>% 
  filter(domain == "keg") %>% select(term.name, overlap.size, p.value) %>% 
  arrange(p.value) 
gprofiler23_48_kegg_human <- gprofiler23_48_kegg_human %>% add_column(time = rep("MOK023_48hpi", nrow(gprofiler23_48_kegg_human)))
gprofiler23_72_kegg_human <- gprofilerresult23_72h_human  %>% 
  filter(domain == "keg") %>% select(term.name, overlap.size, p.value) %>% 
  arrange(p.value) 
gprofiler23_72_kegg_human <- gprofiler23_72_kegg_human %>% add_column(time =rep("MOK023_72hpi", nrow(gprofiler23_72_kegg_human)))
gprofiler23_168_kegg_human <- gprofilerresult23_168h_human  %>% 
  filter(domain == "keg") %>% select(term.name, overlap.size, p.value) %>% 
  arrange(p.value) 
gprofiler23_168_kegg_human <- gprofiler23_168_kegg_human %>% add_column(time =rep("MOK023_168hpi", nrow(gprofiler23_168_kegg_human)))

gprofiler_MOK023_kegg_human <- rbind.data.frame(gprofiler23_24_kegg_human, gprofiler23_48_kegg_human, 
                                           gprofiler23_72_kegg_human, gprofiler23_168_kegg_human) 
View(gprofiler_MOK023_kegg_human) 
write.csv(gprofiler_MOK023_kegg_human, "gprofiler_MOK023_kegg_human.csv")

########## Gprofiler MOK124 ##### 

##run gprofiler
library(gProfileR)

# Order Lfs threshold significant results by adjusted P value 
res124_24h_LF_sig <- res124_24h_LF_sig[order (res124_24h_LF_sig$padj),] #order results
res124_48h_LF_sig <- res124_48h_LF_sig[order (res124_48h_LF_sig$padj),] 
res124_72h_LF_sig <- res124_72h_LF_sig[order (res124_72h_LF_sig$padj),] 
res124_168h_LF_sig <- res124_168h_LF_sig[order (res124_168h_LF_sig$padj),] 

head(res124_168h_LF_sig) 

# all MOK124 gene lists to human ##### 

# Creation of y124_24 is shown as an example in the MOK023 section, these come from Venn diagram script file 
y124_24_human <- inner_join(species_1to1, y124_24) 
write.csv(y124_24_human, "MOK124_24H_human.csv") 
names_124_24h_human <- y124_24_human$Human.gene.stable.ID 

y124_48_human <- inner_join(species_1to1, y124_48) 
write.csv(y124_48_human, "MOK124_48H_human.csv") 
names_124_48h_human <- y124_48_human$Human.gene.stable.ID 

y124_72_human <- inner_join(species_1to1, y124_72) 
write.csv(y124_72_human, "MOK124_72H_human.csv") 
names_124_72h_human <- y124_72_human$Human.gene.stable.ID 

y124_168_human <- inner_join(species_1to1, y124_168) 
write.csv(y124_168_human, "MOK124_168H_human.csv") 
names_124_168h_human <- y124_168_human$Human.gene.stable.ID 

# run gProfiler on HUMAN gene orthologs - MOK124 ######

names_124_24h_human <- as.vector(names_124_24h_human) 
head(names_124_24h_human) 
names_124_48h_human <- as.vector(names_124_48h_human) 
names_124_72h_human <- as.vector(names_124_72h_human) 
names_124_168h_human <- as.vector(names_124_168h_human) 

gprofilerresult124_24h_human <- gprofiler(names_124_24h_human, organism = "hsapiens", ordered_query = T, hier_filtering = "strong")
gprofilerresult124_48h_human <- gprofiler(names_124_48h_human, organism = "hsapiens", hier_filtering = "strong")
gprofilerresult124_72h_human <- gprofiler(names_124_72h_human, organism = "hsapiens", hier_filtering = "strong")
gprofilerresult124_168h_human <- gprofiler(names_124_168h_human, organism = "hsapiens", hier_filtering = "strong")

write.csv(gprofilerresult124_24h_human , file="Gprofiler_MOK124_24HPI_human.csv")
write.csv(gprofilerresult124_48h_human , file="Gprofiler_MOK124_48HPI_human.csv")
write.csv(gprofilerresult124_72h_human , file="Gprofiler_MOK124_72HPI_human.csv")
write.csv(gprofilerresult124_168h_human , file="Gprofiler_MOK124_168HPI_human.csv")

# Make a file combining all KEGG lists - human pathways #### 

# Select kegg pathways, name, gene count and p value - for a table to export to supplementary materials 
# add_column is from tibble package 
library(tidyverse)
gprofiler124_24_kegg_human <- gprofilerresult124_24h_human  %>% 
  filter(domain == "keg") %>% select(term.name, overlap.size, p.value) %>% 
  arrange(p.value) 
gprofiler124_24_kegg_human <- gprofiler124_24_kegg_human %>% add_column(time =rep("MOK124_24hpi", nrow(gprofiler124_24_kegg_human)))
gprofiler124_48_kegg_human <- gprofilerresult124_48h_human  %>% 
  filter(domain == "keg") %>% select(term.name, overlap.size, p.value) %>% 
  arrange(p.value) 
gprofiler124_48_kegg_human <- gprofiler124_48_kegg_human %>% add_column(time = rep("MOK124_48hpi", nrow(gprofiler124_48_kegg_human)))
gprofiler124_72_kegg_human <- gprofilerresult124_72h_human  %>% 
  filter(domain == "keg") %>% select(term.name, overlap.size, p.value) %>% 
  arrange(p.value) 
gprofiler124_72_kegg_human <- gprofiler124_72_kegg_human %>% add_column(time =rep("MOK124_72hpi", nrow(gprofiler124_72_kegg_human)))
gprofiler124_168_kegg_human <- gprofilerresult124_168h_human  %>% 
  filter(domain == "keg") %>% select(term.name, overlap.size, p.value) %>% 
  arrange(p.value) 
gprofiler124_168_kegg_human <- gprofiler124_168_kegg_human %>% add_column(time =rep("MOK124_168hpi", nrow(gprofiler124_168_kegg_human)))

gprofiler_MOK0124_kegg_human <- rbind.data.frame(gprofiler124_24_kegg_human, gprofiler124_48_kegg_human, 
                                                gprofiler124_72_kegg_human, gprofiler124_168_kegg_human) 
View(gprofiler_MOK0124_kegg_human) 
write.csv(gprofiler_MOK0124_kegg_human, "gprofiler_MOK124_kegg_human.csv")

# Make a file with all gProfiler data ##### 

library(tidyverse) 

colnames(gprofilerresult124_24h_human)
View(gprofilerresult124_24h_human)

# MOK023 pathways - select only kegg pathways, name and p value - for heatmaps 
gprofiler23_24h_kegg <- gprofilerresult23_24h_human  %>% 
  filter(domain == "keg") %>% select(term.name, p.value) %>% arrange(p.value)
head(gprofiler23_24h_kegg)  
gprofiler23_48h_kegg <- gprofilerresult23_48h_human  %>% 
  filter(domain == "keg") %>% select(term.name, p.value) %>% arrange(p.value)
gprofiler23_72h_kegg <- gprofilerresult23_72h_human  %>% 
  filter(domain == "keg") %>% select(term.name, p.value) %>% arrange(p.value)
gprofiler23_168h_kegg <- gprofilerresult23_168h_human  %>% 
  filter(domain == "keg") %>% select(term.name, p.value) %>% arrange(p.value)

# MOK124 pathway p value selection 
gprofiler124_24h_kegg <- gprofilerresult124_24h_human  %>% 
  filter(domain == "keg") %>% select(term.name, p.value) %>% arrange(p.value)
head(gprofiler124_24h_kegg) 
gprofiler124_48h_kegg <- gprofilerresult124_48h_human  %>% 
  filter(domain == "keg") %>% select(term.name, p.value) %>% arrange(p.value)
gprofiler124_72h_kegg <- gprofilerresult124_72h_human  %>% 
  filter(domain == "keg") %>% select(term.name, p.value) %>% arrange(p.value)
gprofiler124_168h_kegg <- gprofilerresult124_168h_human  %>% 
  filter(domain == "keg") %>% select(term.name, p.value) %>% arrange(p.value)

# Join all files into one 
KEGG_MOK023 <- full_join(gprofiler23_24h_kegg, gprofiler23_48h_kegg, by="term.name", suffix = c("_24h", "_48h"))
KEGG_MOK023 <- full_join(KEGG_MOK023, gprofiler23_72h_kegg, by = "term.name", suffix = c("_48h", "_72h"))
KEGG_MOK023 <- full_join(KEGG_MOK023, gprofiler23_168h_kegg, by = "term.name", suffix = c("_72h", "_168h"))
View(KEGG_MOK023) 
write.csv(KEGG_MOK023, "Kegg_pathways_gProfiler_human_MOK023_all_HPI.csv")

KEGG_MOK124 <- full_join(gprofiler124_24h_kegg, gprofiler124_48h_kegg, by="term.name", suffix = c("_24h", "_48h"))
KEGG_MOK124 <- full_join(KEGG_MOK124, gprofiler124_72h_kegg, by = "term.name", suffix = c("_48h", "_72h"))
KEGG_MOK124 <- full_join(KEGG_MOK124, gprofiler124_168h_kegg, by = "term.name", suffix = c("_72h", "_168h"))
write.csv(KEGG_MOK124, "Kegg_pathways_gProfiler_human_MOK124_all_HPI.csv")

KEGG_both <- full_join(KEGG_MOK023, KEGG_MOK124, by = "term.name", suffix = c("_023", "_124"))
write.csv(KEGG_both, "Kegg_pathways_gProfiler_both_strains.csv") 

# TOP 10 KEGG terms! 

gprofiler23_24h_KEGG_top10 <- gprofiler23_24h_kegg %>% top_n(-10, p.value) 
gprofiler23_48h_KEGG_top10 <- gprofiler23_48h_kegg %>% top_n(-10, p.value) 
gprofiler23_72h_KEGG_top10 <- gprofiler23_72h_kegg %>% top_n(-10, p.value) 
gprofiler23_168h_KEGG_top10 <- gprofiler23_168h_kegg %>% top_n(-10, p.value) 

gprofiler124_24h_KEGG_top10 <- gprofiler124_24h_kegg %>% top_n(-10, p.value) 
gprofiler124_48h_KEGG_top10 <- gprofiler124_48h_kegg %>% top_n(-10, p.value) 
gprofiler124_72h_KEGG_top10 <- gprofiler124_72h_kegg %>% top_n(-10, p.value) 
gprofiler124_168h_KEGG_top10 <- gprofiler124_168h_kegg %>% top_n(-10, p.value)


# Join all files into one 
KEGG_MOK023_t10 <- full_join(gprofiler23_24h_KEGG_top10, gprofiler23_48h_KEGG_top10, by="term.name", suffix = c("_24h", "_48h"))
KEGG_MOK023_t10 <- full_join(KEGG_MOK023_t10, gprofiler23_72h_KEGG_top10, by = "term.name", suffix = c("_48h", "_72h"))
KEGG_MOK023_t10 <- full_join(KEGG_MOK023_t10, gprofiler23_168h_KEGG_top10, by = "term.name", suffix = c("_72h", "_168h"))
View(KEGG_MOK023_t10) 
write.csv(KEGG_MOK023_t10, "KEGG_pathways_gProfiler_human_MOK023_all_HPI_TOP10.csv")

KEGG_MOK124_t10 <- full_join(gprofiler124_24h_KEGG_top10, gprofiler124_48h_KEGG_top10, by="term.name", suffix = c("_24h", "_48h"))
KEGG_MOK124_t10 <- full_join(KEGG_MOK124_t10, gprofiler124_72h_KEGG_top10, by = "term.name", suffix = c("_48h", "_72h"))
KEGG_MOK124_t10 <- full_join(KEGG_MOK124_t10, gprofiler124_168h_KEGG_top10, by = "term.name", suffix = c("_72h", "_168h"))
write.csv(KEGG_MOK124_t10, "KEGG_pathways_gProfiler_human_MOK124_all_HPI_TOP10.csv")

KEGG_both_t10 <- full_join(KEGG_MOK023_t10, KEGG_MOK124_t10, by = "term.name", suffix = c("_023", "_124"))
write.csv(KEGG_both_t10, "KEGG_pathways_gProfiler_both_strains_TOP10.csv") 
View(KEGG_both_t10) 


# Heatmap for top 10 KEGG pathways 
library("pheatmap")
detach(package:pheatmap) # if necessary 
library("RColorBrewer") 

KEGG_both_C <- KEGG_both_t10
KEGG_both_C[is.na(KEGG_both_C)] <- as.double("NA")
rownames(KEGG_both_C) <- KEGG_both_C$term.name
KEGG_both_C <- KEGG_both_C[,-1] # Remove term.name as a column 
View(KEGG_both_C) 

colnames(KEGG_both_C) <- c("MOK023 24 hpi", "MOK023 48 hpi", "MOK023 72 hpi", "MOK023 168 hpi", "MOK124 24 hpi", "MOK124 48 hpi", "MOK124 72 hpi", "MOK124 168 hpi")

# best version of heatmap 
KEGG_heatmap <- pheatmap(KEGG_both_C, col = colors, na_col = "grey", 
                       cluster_cols = FALSE, cluster_rows = FALSE, 
                       cellwidth = 15, fontsize = 12)

tiff("my_KEGG_heatmap_top10_CODED.tiff", height = 20, width = 24, units = 'cm', 
     compression = "lzw", res = 600) 
KEGG_heatmap 
dev.off()


# GENE ONTOLOGY files 
# MOK023 pathways - select only kegg pathways, name and p value 
gprofiler23_24h_GO <- gprofilerresult23_24h_human  %>% 
  filter(domain == "BP") %>% select(term.name, p.value) %>% arrange(p.value)
head(gprofiler23_24h_GO) 
gprofiler23_48h_GO <- gprofilerresult23_48h_human  %>% 
  filter(domain == "BP") %>% select(term.name, p.value) %>% arrange(p.value)
gprofiler23_72h_GO <- gprofilerresult23_72h_human  %>% 
  filter(domain == "BP") %>% select(term.name, p.value) %>% arrange(p.value)
gprofiler23_168h_GO <- gprofilerresult23_168h_human  %>% 
  filter(domain == "BP") %>% select(term.name, p.value) %>% arrange(p.value)

# MOK124 pathway selection 
gprofiler124_24h_GO <- gprofilerresult124_24h_human  %>% 
  filter(domain == "BP") %>% select(term.name, p.value) %>% arrange(p.value)
head(gprofiler124_24h_GO) 
gprofiler124_48h_GO <- gprofilerresult124_48h_human  %>% 
  filter(domain == "BP") %>% select(term.name, p.value) %>% arrange(p.value)
gprofiler124_72h_GO <- gprofilerresult124_72h_human  %>% 
  filter(domain == "BP") %>% select(term.name, p.value) %>% arrange(p.value)
gprofiler124_168h_GO <- gprofilerresult124_168h_human  %>% 
  filter(domain == "BP") %>% select(term.name, p.value) %>% arrange(p.value)

# Join all files into one - these are not sorted by p value 
GO_MOK023 <- full_join(gprofiler23_24h_GO, gprofiler23_48h_GO, by="term.name", suffix = c("_24h", "_48h"))
GO_MOK023 <- full_join(GO_MOK023, gprofiler23_72h_GO, by = "term.name", suffix = c("_48h", "_72h"))
GO_MOK023 <- full_join(GO_MOK023, gprofiler23_168h_GO, by = "term.name", suffix = c("_72h", "_168h"))
View(GO_MOK023) 
write.csv(GO_MOK023, "GO_pathways_gProfiler_human_MOK023_all_HPI.csv")

GO_MOK124 <- full_join(gprofiler124_24h_GO, gprofiler124_48h_GO, by="term.name", suffix = c("_24h", "_48h"))
GO_MOK124 <- full_join(GO_MOK124, gprofiler124_72h_GO, by = "term.name", suffix = c("_48h", "_72h"))
GO_MOK124 <- full_join(GO_MOK124, gprofiler124_168h_GO, by = "term.name", suffix = c("_72h", "_168h"))
write.csv(GO_MOK124, "GO_pathways_gProfiler_human_MOK124_all_HPI.csv")

GO_both <- full_join(GO_MOK023, GO_MOK124, by = "term.name", suffix = c("_023", "_124"))
write.csv(GO_both, "GO_pathways_gProfiler_both_strains.csv") 

# Now get top 10 GO terms #### 

gprofiler23_24h_GO_top10 <- gprofiler23_24h_GO %>% top_n(-10, p.value) 
gprofiler23_48h_GO_top10 <- gprofiler23_48h_GO %>% top_n(-10, p.value) 
gprofiler23_72h_GO_top10 <- gprofiler23_72h_GO %>% top_n(-10, p.value) 
gprofiler23_168h_GO_top10 <- gprofiler23_168h_GO %>% top_n(-10, p.value) 

gprofiler124_24h_GO_top10 <- gprofiler124_24h_GO %>% top_n(-10, p.value) 
gprofiler124_48h_GO_top10 <- gprofiler124_48h_GO %>% top_n(-10, p.value) 
gprofiler124_72h_GO_top10 <- gprofiler124_72h_GO %>% top_n(-10, p.value) 
gprofiler124_168h_GO_top10 <- gprofiler124_168h_GO %>% top_n(-10, p.value)


# Join all files into one 
GO_MOK023_t10 <- full_join(gprofiler23_24h_GO_top10, gprofiler23_48h_GO_top10, by="term.name", suffix = c("_24h", "_48h"))
GO_MOK023_t10 <- full_join(GO_MOK023_t10, gprofiler23_72h_GO_top10, by = "term.name", suffix = c("_48h", "_72h"))
GO_MOK023_t10 <- full_join(GO_MOK023_t10, gprofiler23_168h_GO_top10, by = "term.name", suffix = c("_72h", "_168h"))
View(GO_MOK023_t10) 
write.csv(GO_MOK023_t10, "GO_pathways_gProfiler_human_MOK023_all_HPI_TOP10.csv")

GO_MOK124_t10 <- full_join(gprofiler124_24h_GO_top10, gprofiler124_48h_GO_top10, by="term.name", suffix = c("_24h", "_48h"))
GO_MOK124_t10 <- full_join(GO_MOK124_t10, gprofiler124_72h_GO_top10, by = "term.name", suffix = c("_48h", "_72h"))
GO_MOK124_t10 <- full_join(GO_MOK124_t10, gprofiler124_168h_GO_top10, by = "term.name", suffix = c("_72h", "_168h"))
write.csv(GO_MOK124_t10, "GO_pathways_gProfiler_human_MOK124_all_HPI_TOP10.csv")

GO_both_t10 <- full_join(GO_MOK023_t10, GO_MOK124_t10, by = "term.name", suffix = c("_023", "_124"))
write.csv(GO_both_t10, "GO_pathways_gProfiler_both_strains_TOP10.csv") 
View(GO_both_t10)

# Heatmap of the GO pathways with p values #### 
library("pheatmap")
detach(package:pheatmap) # if necessary 
library("RColorBrewer") 

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(GO_both_M,
         col = colors, cluster_cols = FALSE, cluster_rows = FALSE) 

##### Heatmap with greyed out NAs - top 10 GO terms! 
GO_both_C <- GO_both_t10
GO_both_C[is.na(GO_both_C)] <- as.double("NA")
rownames(GO_both_C) <- GO_both_C$term.name
GO_both_C <- GO_both_C[,-1] # Remove term.name as a column 
View(GO_both_C) 
colnames(GO_both_C) <- c("MOK023 24 hpi", "MOK023 48 hpi", "MOK023 72 hpi", "MOK023 168 hpi", "MOK124 24 hpi", "MOK124 48 hpi", "MOK124 72 hpi", "MOK124 168 hpi")

GO_heatmap <- pheatmap(GO_both_C, col = colors, na_col = "grey", 
                       cluster_cols = FALSE, cluster_rows = FALSE, 
                       cellwidth = 15, fontsize = 12)

tiff("my_GO_heatmap_top10.tiff", height = 30, width = 24, units = 'cm', 
     compression = "lzw", res = 600) 
GO_heatmap 
dev.off()
