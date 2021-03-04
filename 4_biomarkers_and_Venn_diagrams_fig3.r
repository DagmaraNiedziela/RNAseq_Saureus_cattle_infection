#####Venn diagrams ##### 
library(VennDiagram) 

# Get gene lists #####

# Reminder of what the data is: 
res23_24h_LF_sig <- subset(res23_24h_LF, padj < 0.05)
res23_48h_LF_sig <- subset(res23_48h_LF, padj < 0.05)
res23_72h_LF_sig <- subset(res23_72h_LF, padj < 0.05)
res23_168h_LF_sig <- subset(res23_168h_LF, padj < 0.05)
res124_24h_LF_sig <- subset(res124_24h_LF, padj < 0.05)
res124_48h_LF_sig <- subset(res124_48h_LF, padj < 0.05)
res124_72h_LF_sig <- subset(res124_72h_LF, padj < 0.05)
res124_168h_LF_sig <- subset(res124_168h_LF, padj < 0.05) 
# check one to see if it's ordered, otherwise order by Padj (but this should have been done in pathway analysis file) 
head(res23_24h_LF_sig) 

# extract row names as a column 
library(tidyverse)

res23_24h_LF_sigDF <- as.data.frame(res23_24h_LF_sig) 
res23_24h_LF_sigDF <- tibble::rownames_to_column(res23_24h_LF_sigDF)
head(res124_24h_LF_sigDF) 
nrow(res23_24h_LF_sigDF)
y23_24 <- res23_24h_LF_sigDF %>% dplyr::select(rowname) 
y23_24 <- as.vector(y23_24)
head(y23_24)

res124_24h_LF_sigDF <- as.data.frame(res124_24h_LF_sig) 
res124_24h_LF_sigDF <- tibble::rownames_to_column(res124_24h_LF_sigDF)
y124_24 <- res124_24h_LF_sigDF %>% dplyr::select(rowname) 
y124_24 <- as.vector(y124_24)

y24hpi <- intersect(y23_24, y124_24) 
nrow(y24hpi) #1143 

res23_48h_LF_sigDF <- as.data.frame(res23_48h_LF_sig) 
res23_48h_LF_sigDF <- tibble::rownames_to_column(res23_48h_LF_sigDF)
y23_48 <- res23_48h_LF_sigDF %>% dplyr::select(rowname) 
y23_48 <- as.vector(y23_48) 
head(y23_48)

res23_72h_LF_sigDF <- as.data.frame(res23_72h_LF_sig) 
res23_72h_LF_sigDF <- tibble::rownames_to_column(res23_72h_LF_sigDF)
y23_72 <- res23_72h_LF_sigDF %>% dplyr::select(rowname) 
y23_72 <- as.vector(y23_72) 

res23_168h_LF_sigDF <- as.data.frame(res23_168h_LF_sig) 
res23_168h_LF_sigDF <- tibble::rownames_to_column(res23_168h_LF_sigDF)
y23_168 <- res23_168h_LF_sigDF %>% dplyr::select(rowname) 
y23_168 <- as.vector(y23_168) 

res124_48h_LF_sigDF <- as.data.frame(res124_48h_LF_sig) 
res124_48h_LF_sigDF <- tibble::rownames_to_column(res124_48h_LF_sigDF)
y124_48 <- res124_48h_LF_sigDF %>% dplyr::select(rowname) 
y124_48 <- as.vector(y124_48) 

res124_72h_LF_sigDF <- as.data.frame(res124_72h_LF_sig) 
res124_72h_LF_sigDF <- tibble::rownames_to_column(res124_72h_LF_sigDF)
y124_72 <- res124_72h_LF_sigDF %>% dplyr::select(rowname) 
y124_72 <- as.vector(y124_72) 

res124_168h_LF_sigDF <- as.data.frame(res124_168h_LF_sig) 
res124_168h_LF_sigDF <- tibble::rownames_to_column(res124_168h_LF_sigDF)
y124_168 <- res124_168h_LF_sigDF %>% dplyr::select(rowname) 
y124_168 <- as.vector(y124_168) 

###Diagrams for MOK023 - create intersect vectors ##### 
#Make Venns for: MOK023 separately, MOK124 separately, and for both at each time point. 

nrow(y23_24) #1278
nrow(y23_48) # 2248
y23_24_48 <- intersect(y23_24, y23_48) # this intersect requires the tidyverse!!!! 
nrow(y23_24_48) #1034
nrow(y23_72) #1986
nrow(y23_168) #1750 
y23_48_72 <- intersect(y23_48, y23_72)
nrow(y23_48_72) #1542
y23_24_72 <- intersect(y23_24, y23_72)
y23_72_168 <- intersect(y23_72, y23_168)
V72_168 #one is blank, in all I will have to substract 1
y23_24_168 <- intersect(y23_24, y23_168)
y23_48_168 <- intersect(y23_48, y23_168)
nrow(y23_48_168) # 1206 

library(dplyr) 
# The following objects are masked from ‘package:base’:
# intersect, setdiff, setequal, union 

y23_1234 <- Reduce(intersect, list(y23_24, y23_48, y23_72, y23_168)) 
nrow(y23_1234) # 402 

y23_123 <- Reduce(intersect, list(y23_24, y23_48, y23_72))
nrow(y23_123) # 661
y23_234 <- Reduce(intersect, list(y23_48, y23_72, y23_168))
nrow(y23_234) # 1164 
y23_124 <- Reduce(intersect, list(y23_24, y23_48, y23_168))
nrow(y23_124) # 409
y23_134 <- Reduce(intersect, list(y23_24, y23_72, y23_168))
nrow(y23_134) # 405

###Draw Venn diagrams - MOK023 ######

# 1 Venn diagram with 4 circles 
grid.newpage()
venn.plot23 <- draw.quad.venn(area1 = nrow(y23_24), area2 = nrow(y23_48), area3 = nrow(y23_72), area4 = nrow(y23_168), n12 = nrow(y23_24_48), n13 = nrow(y23_24_72), n14 = nrow(y23_24_168), n23 = nrow(y23_48_72), n24 = nrow(y23_48_168), 
                              n34 = nrow(y23_72_168), n123 = nrow(y23_123), n124 = nrow(y23_124), n134 = nrow(y23_134), n234 = nrow(y23_234), n1234 = nrow(y23_1234), category = c("24 HPI", "48 HPI", "72 HPI", "168 HPI"), lty = "solid", 
                              fill            = c("#c6e1ef",
                                                  "#a1d7f4",
                                                  "#2aa7ea",
                                                  "#0c15cc"),
                              euler.d=FALSE,
                              scaled=FALSE,
                              sub=substitute( paste(bolditalic('M. bovis'))),
                              sub.fontfamily = "Arial",
                              sub.cex=1.2,
                              sub.pos=c(0.5,1),
                              lwd=1,
                              alpha           = rep(0.50,4),
                              label.col       = "#003333",
                              cex             = 1.5,
                              fontfamily      = "Arial",            
                              cat.pos=c(-5,10,10,0),
                              cat.col         = "black",
                              cat.cex         = 1.5,
                              cat.fontfamily  = "Arial",
                              cat.fontface    = 2,
                              rotation.degree = 360,
                              margin          = 0,
                              height          = 10,
                              width           = 4,
                              units           = 'cm',
                              compression     = 'lzw',
                              resolution      = 1200
)
grid.draw(venn.plot23) 
grid.newpage()

# 2 Add Venn plot title 
require(gridExtra)
venn.plot23 <- grid.arrange(gTree(children=venn.plot23), top=textGrob("MOK023 DE genes", gp=gpar(fontsize=20)))

# Writing to file 
tiff(filename = "Quad_Venn_diagram_MOK023.tiff", compression = "lzw");
grid.draw(venn.plot23);
dev.off(); 

# 3 pretty file 
Cairo(file="Quad_Venn_diagram_MOK023.png", 
      type="png",
      units="in", 
      width=10, 
      height=10, 
      pointsize=12, 
      dpi=72)
grid.draw(venn.plot23)
dev.off()

### Diagrams for MOK124 - create intersect vectors ##### 
nrow(y124_24) #2293
nrow(y124_48) # 1979
y124_24_48 <- intersect(y124_24, y124_48)
nrow(y124_24_48) # 722
nrow(y124_72) # 1428
nrow(y124_168) # 1544 
y124_48_72 <- intersect(y124_48, y124_72)
nrow(y124_48_72) # 1331
y124_24_72 <- intersect(y124_24, y124_72)
y124_72_168 <- intersect(y124_72, y124_168)
y124_24_168 <- intersect(y124_24, y124_168)
y124_48_168 <- intersect(y124_48, y124_168)
nrow(y124_48_168) # 1380

y124_1234 <- Reduce(intersect, list(y124_24, y124_48, y124_72, y124_168))
nrow(y124_1234) # 362

y124_123 <- Reduce(intersect, list(y124_24, y124_48, y124_72))
nrow(y124_123) # 382
y124_234 <- Reduce(intersect, list(y124_48, y124_72, y124_168))
nrow(y124_234) # 1181
y124_124 <- Reduce(intersect, list(y124_24, y124_48, y124_168))
nrow(y124_124) # 493
y124_134 <- Reduce(intersect, list(y124_24, y124_72, y124_168))
nrow(y124_134) # 372

###Draw Venn diagrams - MOK124 ######

# 1 Venn diagram with 4 circles 
grid.newpage() 
venn.plot124 <- draw.quad.venn(area1 = nrow(y124_24), area2 = nrow(y124_48), area3 = nrow(y124_72), area4 = nrow(y124_168), n12 = nrow(y124_24_48), n13 = nrow(y124_24_72), n14 = nrow(y124_24_168), n23 = nrow(y124_48_72), n24 = nrow(y124_48_168), 
                               n34 = nrow(y124_72_168), n123 = nrow(y124_123), n124 = nrow(y124_124), n134 = nrow(y124_134), n234 = nrow(y124_234), n1234 = nrow(y124_1234), category = c("24 HPI", "48 HPI", "72 HPI", "168 HPI"), lty = "solid", 
                               fill            = c("#eabbc4",
                                                   "#ed8e9f",
                                                   "#ef3456",
                                                   "#96031d"),
                               euler.d=FALSE,
                               scaled=FALSE,
                               lwd=1,
                               alpha           = rep(0.50,4),
                               sub=substitute( paste(bolditalic('M. tuberculosis'))),
                               sub.fontfamily = "Arial",
                               sub.cex=1.2,
                               sub.pos=c(0.5,1),
                               label.col       = "#003333",
                               cex             = 1.5,
                               fontfamily      = "Arial",
                               cat.pos=c(-5,10,10,0),
                               cat.col         = "black",
                               cat.cex         = 1.5,
                               cat.fontfamily  = "Arial",
                               cat.fontface    = 2,
                               rotation.degree = 360,
                               margin          = 0,
                               height          = 10,
                               width           = 4,
                               units           = 'cm',
                               compression     = 'lzw',
                               resolution      = 1200)		
grid.draw(venn.plot124) 

# 2 Add Venn plot title 
require(gridExtra)
venn.plot124 <- grid.arrange(gTree(children=venn.plot124), top=textGrob("MOK124 DE genes", gp=gpar(fontsize=20)))

# Writing to file - not done 
tiff(filename = "Quad_Venn_diagram_MOK124.tiff", compression = "lzw");
grid.draw(venn.plot124);
dev.off(); 

# 3 pretty plot to file 
Cairo(file="Quad_Venn_diagram_MOK124.png", 
      type="png",
      units="in", 
      width=10, 
      height=10, 
      pointsize=12, 
      dpi=72)
grid.draw(venn.plot124)
dev.off()

# Merge MOK023 and MOK124 plots into 1 panel ##### 
library("gridExtra")
library("grid")
library("lattice") 

venn.plot<-grid.arrange(venn.plot23, venn.plot124, 
                        ncol = 2, 
                        widths = c(1,1),
                        heights = 1) 

grid.newpage() 
tiff(filename = "Quad_Venn_diagram_both_strains.tiff", compression = "lzw");
grid.draw(venn.plot);
dev.off(); 

# Save to file 
Cairo(file="Quad_Venn_diagram_both_strains.png", 
      type="png",
      units="in", 
      width=10, 
      height=5, 
      pointsize=12, 
      dpi=72)
grid.draw(venn.plot)
dev.off()

### Diagram for overlaps in both strains ##### 

y_overlaps <- intersect(y23_1234, y124_1234) 
nrow(y_overlaps) # 130

#overlaps 1234 
grid.newpage()
venn_overlaps <- draw.pairwise.venn(nrow(y23_1234), nrow(y124_1234), nrow(y_overlaps), category = c("MOK023", "MOK124"), lty = rep("solid", 2), 
                   fill = c("#2aa7ea", "#ef3456"), alpha = rep(0.5, 2), cat.pos = c(0,0), 
                   euler.d=FALSE,
                   scaled=FALSE,
                   lwd=1,
                   sub=substitute( paste(bolditalic('24 HPI'))),
                   sub.fontfamily = "Arial",
                   sub.cex=1.2,
                   sub.pos=c(0.5,1),
                   label.col       = "#003333",
                   cex             = 1.75,
                   fontfamily      = "Arial",
                   cat.col         = "black",
                   cat.cex         = 1.75,
                   cat.fontfamily  = "Arial",
                   cat.fontface    = 2,
                   rotation.degree = 360,
                   margin          = 0,
                   height          = 10,
                   width           = 4,
                   units           = 'cm',
                   compression     = 'lzw',
                   resolution      = 1200, 
                   cat.just=list(c(1,0) , c(0,0))) 

venn_overlaps <- grid.arrange(gTree(children=venn_overlaps), top=textGrob("All HPI overlaps", gp=gpar(fontsize=25)))

# Writing to file 
tiff(filename = "Double_Venn_overlaps.tiff", compression = "lzw");
grid.draw(venn_overlaps);
dev.off(); 

# Pretty file 
Cairo(file="Double_Venn_overlaps.png", 
      type="png",
      units="in", 
      width=5, 
      height=5, 
      pointsize=12, 
      dpi=72)
grid.draw(venn_overlaps)
dev.off()