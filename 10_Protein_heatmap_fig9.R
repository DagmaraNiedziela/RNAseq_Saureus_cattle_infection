library(readxl)
library(tidyverse) 
library(dplyr)

# Protein data import - MFI - mean fluorescence intensity ####

Gilles_protein_MFI <- read_excel("Protein_analysis_Gilles_copy.xlsx", 2) 
Gilles_protein_MFI
colnames(Gilles_protein_MFI) 
ncol(Gilles_protein_MFI) # 18, so 15 proteins 
names <- c("IL1A","IL1B","IL1RA","IL6","IL17A","IL2","IL4","IFNy","IL8","CXCL10","CCL2","CCL3","IL10","CCL4","TNFA") 

Gilles_protein_MFI_for_heatmap <- Gilles_protein_MFI
Gilles_protein_MFI_for_heatmap$condition <- paste(Gilles_protein_MFI_for_heatmap$Group, Gilles_protein_MFI_for_heatmap$HPI, sep = "_")  
View(Gilles_protein_MFI_for_heatmap)

# Create a table with all summaries
Gilles_protein_summary_for_heatmap <- Gilles_protein_MFI_for_heatmap %>% 
  group_by(HPI, Group) %>% 
  summarize_if(is.numeric, mean, na.rm = TRUE)

Gilles_protein_summary_for_heatmap 

Gilles_protein_summary_for_heatmap$condition <- paste(Gilles_protein_summary_for_heatmap$Group, Gilles_protein_summary_for_heatmap$HPI, sep = "_")  

Gilles_protein_summary_for_heatmap2 <- Gilles_protein_summary_for_heatmap %>% 
  select(-Cow_ID, -HPI, -Group) 
Gilles_protein_summary_for_heatmap2
# HPI still there for some reason 
Gilles_protein_summary_for_heatmap2 <- Gilles_protein_summary_for_heatmap2[,-1] # remove HPI 
colnames(Gilles_protein_summary_for_heatmap2)
# "IL1A"      "IL1B"      "IL1RA"     "IL6"       "IL17A"     "IL2"       "IL4"       "IFNy"      "IL8"      
# "CXCL10"    "CCL2"      "CCL3"      "IL10"      "CCL4"      "TNFA"      "condition" 

Gilles_protein_summary_for_heatmap2$condition 
# [1] "MOK023_0"   "MOK124_0"   "MOK023_7"   "MOK124_7"   "MOK023_24"  "MOK124_24"  "MOK023_48"  "MOK124_48"  "MOK023_72" 
# [10] "MOK124_72"  "MOK023_168" "MOK124_168" "MOK023_336" "MOK124_336" "MOK023_504" "MOK124_504"

Gilles_protein_summary_for_heatmap3 <- as.data.frame(Gilles_protein_summary_for_heatmap2)
Gilles_protein_summary_for_heatmap3 <- t(Gilles_protein_summary_for_heatmap3)
colnames(Gilles_protein_summary_for_heatmap3) <- Gilles_protein_summary_for_heatmap2$condition
rownames(Gilles_protein_summary_for_heatmap3)
# [1] "IL1A"      "IL1B"      "IL1RA"     "IL6"       "IL17A"     "IL2"       "IL4"       "IFNy"      "IL8"      
# [10] "CXCL10"    "CCL2"      "CCL3"      "IL10"      "CCL4"      "TNFA"      "condition"
nrow(Gilles_protein_summary_for_heatmap3) #16 

head(Gilles_protein_summary_for_heatmap3) 
ncol(Gilles_protein_summary_for_heatmap3) #16

Gilles_protein_summary_for_heatmap3 <- Gilles_protein_summary_for_heatmap3[-16,]
Gilles_protein_summary_for_heatmap4 <- Gilles_protein_summary_for_heatmap3
class(Gilles_protein_summary_for_heatmap4) # "matrix" "array" 
Gilles_protein_summary_for_heatmap4 <- as.data.frame(Gilles_protein_summary_for_heatmap4)
Gilles_protein_summary_for_heatmap4 <- Gilles_protein_summary_for_heatmap4 %>% 
  transmute_all(function(x) (as.numeric(as.character(x))))

head(Gilles_protein_summary_for_heatmap4)
nrow(Gilles_protein_summary_for_heatmap4) #15
ncol(Gilles_protein_summary_for_heatmap4) #16
rownames(Gilles_protein_summary_for_heatmap4) <- rownames(Gilles_protein_summary_for_heatmap3)[1:15]
# colnames(Gilles_protein_summary_for_heatmap4) <- c() 

# Reorder to have MOK023 then MOK124
Gilles_protein_summary_for_heatmap4_ord <- Gilles_protein_summary_for_heatmap4[,c(1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16)]
head(Gilles_protein_summary_for_heatmap4_ord)

# Data for heatmap with time points corresponding to the RNAseq ones 
# Need to order them by group and then HPI 
Gilles_protein_summary_for_heatmap4_lessHPI <- Gilles_protein_summary_for_heatmap4_ord[,c(1,3,4,5,6,9,11,12,13,14)]
head(Gilles_protein_summary_for_heatmap4_lessHPI)
colnames(Gilles_protein_summary_for_heatmap4_lessHPI) <- c("MOK023 0 hpi","MOK023 24 hpi", "MOK023 48 hpi", "MOK023 72 hpi", 
                                                           "MOK023 168 hpi", "MOK124 0 hpi", "MOK124 24 hpi", "MOK124 48 hpi", "MOK124 72 hpi", "MOK124 168 hpi") 

# Make heatmap #### 
library(pheatmap)

library(RColorBrewer)
display.brewer.all()
colors2 <- colorRampPalette(brewer.pal(9, "Reds"))(50) 

breaks_hm <-  c(0,10,50,100,500,1000,10000,15000,20000)
less_colors <- colorRampPalette(brewer.pal(9, "Reds"))(length(breaks_hm)-1)
display.brewer.pal(n = 8, name = 'Reds') 

# Heatmap divided into 3 cytokine types #### 

View(Gilles_protein_summary_for_heatmap4_lessHPI) 
Gilles_protein_summary_for_heatmap4_lessHPI
Gilles_protein_summary_for_heatmap4_lessHPI2 <- Gilles_protein_summary_for_heatmap4_lessHPI[c(1,2,3,4,15,8,6,7,13,5,11,12,14,10,9),]
Gilles_protein_summary_for_heatmap4_lessHPI2

rownames(Gilles_protein_summary_for_heatmap4_lessHPI2) 
# [1] "IL1A"   "IL1B"   "IL1RA"  "IL6"    "TNFA"   "IFNy"   "IL2"    "IL4"    "IL10"   "IL17A"  "CCL2"   "CCL3"   
# "CCL4"   "CXCL10"
# [15] "IL8" 
# Vector of row names with cytokine proteins 
row_labels <- c("IL-1\u03b1", "IL-1\u03b2", "IL-1Ra", "IL-6", "TNF-\u03b1", "IFN-\u03b3", "IL-2", "IL-4", "IL-10", "IL-17A", 
                "CCL2", "CCL3", "CCL4", "CXCL10", "IL-8") 
'\u03b3'
'\u03b1'
'\u03b2'
row_labels

cytokine_protein_heatmap_less_HPI2 <- pheatmap(Gilles_protein_summary_for_heatmap4_lessHPI2, na_col = "grey", color = less_colors,
                                               cluster_cols = FALSE, cluster_rows = FALSE, gaps_row = c(5,10),
                                               cellwidth = 15, fontsize = 8, breaks = breaks_hm, legend_breaks = c(0,1000,5000,10000,15000,20000), 
                                               angle_col = 315, labels_row = row_labels) 
											   
cytokine_protein_heatmap_less_HPI2

tiff("cytokine_protein_heatmap_3types.tiff", height = 10, width = 20, units = 'cm', 
     compression = "lzw", res = 1200) 
cytokine_protein_heatmap_less_HPI2
dev.off() 

# Try to save with Cairo just in case 
install.packages("Cairo")
library(Cairo)
Cairo(file="cytokine_protein_heatmap_3types.png", 
      type="png",
      units="cm", 
      width=20, 
      height=10, 
      pointsize=60, 
      dpi=1200)
cytokine_protein_heatmap_less_HPI2
dev.off() 

