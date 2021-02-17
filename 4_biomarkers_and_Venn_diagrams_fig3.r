#####Venn diagrams ##### 
library(VennDiagram) 
# For a cool graphically nice plot, check Kerri Malone's code in AlvMac folder on github

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

###Diagrams for MOK023 - create vectors ##### 
#Make Venns for: MOK023 separately, MOK124 separately, and for both at each time point? 
#Definitely export overlap lists? 

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

#for more lists than 2
Reduce(intersect, list(a,b,c))

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

#Venn diagram with 3 circles 
# ** adapted Kerri's pretty plot - works ##### 
grid.newpage()
venn.plot23 <- draw.quad.venn(area1 = nrow(y23_24), area2 = nrow(y23_48), area3 = nrow(y23_72), area4 = nrow(y23_168), n12 = nrow(y23_24_48), n13 = nrow(y23_24_72), n14 = nrow(y23_24_168), n23 = nrow(y23_48_72), n24 = nrow(y23_48_168), 
                              n34 = nrow(y23_72_168), n123 = nrow(y23_123), n124 = nrow(y23_124), n134 = nrow(y23_134), n234 = nrow(y23_234), n1234 = nrow(y23_1234), category = c("24 HPI", "48 HPI", "72 HPI", "168 HPI"), lty = "blank", 
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

###Diagrams for MOK124 - create vectors - redo this every time, same vector names##### 
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
V72_168 #one is blank, in all I will have to substract 1
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

#Venn diagram with 3 circles 
# ** Kerri's pretty plot #### 
grid.newpage() 
venn.plot124 <- draw.quad.venn(area1 = nrow(y124_24), area2 = nrow(y124_48), area3 = nrow(y124_72), area4 = nrow(y124_168), n12 = nrow(y124_24_48), n13 = nrow(y124_24_72), n14 = nrow(y124_24_168), n23 = nrow(y124_48_72), n24 = nrow(y124_48_168), 
                               n34 = nrow(y124_72_168), n123 = nrow(y124_123), n124 = nrow(y124_124), n134 = nrow(y124_134), n234 = nrow(y124_234), n1234 = nrow(y124_1234), category = c("24 HPI", "48 HPI", "72 HPI", "168 HPI"), lty = "blank", 
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

# Merge MOK023 and MOK124 plots, pretty ##### 
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

###Diagrams for both strains at each time point - create vectors ##### 

y24H <- intersect(y23_24, y124_24)
y48H <- intersect(y23_48, y124_48)
y72H <- intersect(y23_72, y124_72)
y168H <- intersect(y23_168, y124_168) 
nrow(y24H) #1143
nrow(y48H) #1106
nrow(y72H) #940
nrow(y168H) #1043

y_overlaps <- intersect(y23_1234, y124_1234) 
nrow(y_overlaps) # 130

###Draw Venn diagrams - Times ######

#Venn diagram with 2 circles 
#24HPI
grid.newpage()
venn24hpi <- draw.pairwise.venn(nrow(y23_24), nrow(y124_24), nrow(y24H), category = c("MOK023", "MOK124"), lty = rep("blank", 
                    2), fill = c("#2aa7ea", "#ef3456"), alpha = rep(0.5, 2),  euler.d=FALSE,
                    scaled=FALSE,
                    lwd=1,
                    sub=substitute( paste(bolditalic('24 HPI'))),
                    sub.fontfamily = "Arial",
                    sub.cex=1.2,
                    sub.pos=c(0.5,1),
                    label.col       = "#003333",
                    cex             = 1.75,
                    fontfamily      = "Arial",
                    cat.pos=c(0,0),
                    cat.col         = "black",
                    cat.cex         = 1.75,
                    cat.fontfamily  = "Arial",
                    cat.fontface    = 2,
                    rotation.degree = 180,
                    margin          = 0,
                    height          = 10,
                    width           = 4,
                    units           = 'cm',
                    compression     = 'lzw',
                    resolution      = 1200, 
                    cat.just=list(c(1,0) , c(0,0))) 

venn24hpi <- grid.arrange(gTree(children=venn24hpi), top=textGrob("24 HPI", gp=gpar(fontsize=25)))
venn24hpi

# Writing to file 
tiff(filename = "Double_Venn_diagram_24hpi.tiff", compression = "lzw");
grid.draw(venn24hpi);
dev.off();

#48HPI
grid.newpage()
venn48hpi <- draw.pairwise.venn(nrow(y23_48), nrow(y124_48), nrow(y48H), category = c("MOK023", "MOK124"), lty = rep("blank", 
                    2), fill = c("#2aa7ea", "#ef3456") , alpha = rep(0.5, 2), cat.pos = c(0,0), 
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
                    cat.just=list(c(1,0) , c(0,0))) #no need to flip, MOK023 has more and is on left
#cat.dist - changed from rep (2 of the same), to vector, to specify that MOK124 should be further 
venn48hpi <- grid.arrange(gTree(children=venn48hpi), top=textGrob("48 HPI", gp=gpar(fontsize=25)))

# Writing to file 
tiff(filename = "Double_Venn_diagram_48hpi.tiff", compression = "lzw");
grid.draw(venn48hpi);
dev.off();

#72HPI
grid.newpage()
venn72hpi <- draw.pairwise.venn(nrow(y23_72), nrow(y124_72), nrow(y72H), category = c("MOK023", "MOK124"), lty = rep("blank", 2) , 
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
#cex - font size of numbers in circles 
#cat.cex - font size of titles 
#fontface = "bold" (for cex)
#cat.fontface = "bold" (for cat.cex) 
venn72hpi <- grid.arrange(gTree(children=venn72hpi), top=textGrob("72 HPI", gp=gpar(fontsize=25)))

# Writing to file 
tiff(filename = "Double_Venn_diagram_72hpi.tiff", compression = "lzw");
grid.draw(venn72hpi);
dev.off();

#168HPI
grid.newpage()
venn168hpi <- draw.pairwise.venn(nrow(y23_168), nrow(y124_168), nrow(y168H), category = c("MOK023", "MOK124"), lty = rep("blank", 2), 
                   fill = c("#2aa7ea", "#ef3456"), alpha = rep(0.5, 2), cat.pos = c(0.5,0.5), 
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
venn168hpi <- grid.arrange(gTree(children=venn168hpi), top=textGrob("168 HPI", gp=gpar(fontsize=25)))

# Writing to file 
tiff(filename = "Double_Venn_diagram_168hpi.tiff", compression = "lzw");
grid.draw(venn168hpi);
dev.off();

# ** Merge all 4 plots into 1 figure #### 
venn.plot2<-grid.arrange(venn24hpi, venn48hpi, venn72hpi, venn168hpi, 
                        ncol = 2, 
                        nrow = 2, 
                        widths = c(1,1),
                        heights = c(1,1)) 


grid.newpage() 
tiff(filename = "Double_Venn_diagrams_allHPI.tiff");
grid.draw(venn.plot2);
dev.off(); 

library(Cairo)
Cairo(file="Double_Venn_diagrams_allHPI.png", 
      type="png",
      units="in", 
      width=10, 
      height=10, 
      pointsize=12, 
      dpi=72)
grid.draw(venn.plot2)
dev.off()

#overlaps 1234 
grid.newpage()
venn_overlaps <- draw.pairwise.venn(nrow(y23_1234), nrow(y124_1234), nrow(y_overlaps), category = c("MOK023", "MOK124"), lty = rep("blank", 2), 
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

#Annotate overlapped genes - done for v_overlaps #####
library("AnnotationDbi")
library("OrganismDbi") 
library("org.Bt.eg.db")
columns(org.Bt.eg.db)
keytypes(org.Bt.eg.db)

typeof(y23_1234$rowname) #"character"

y23_1234$Genename <- mapIds(org.Bt.eg.db,
                                 keys=y23_1234$rowname,
                                 column="GENENAME",
                                 keytype="ENSEMBL",
                                 multiVals="first") 

y23_1234$Alias <- mapIds(org.Bt.eg.db,
                           keys=y23_1234$rowname,
                           column="ALIAS",
                           keytype="ENSEMBL",
                           multiVals="first") 

head(y23_1234)

write.csv(y23_1234, "MOK023_overlaps.csv") 
write.table(y23_1234, "MOK023_overlaps.txt", row.names = y23_1234$rowname, sep = "\t")

y124_1234$Genename <- mapIds(org.Bt.eg.db,
                            keys=y124_1234$rowname,
                            column="GENENAME",
                            keytype="ENSEMBL",
                            multiVals="first") 

y124_1234$Alias <- mapIds(org.Bt.eg.db,
                         keys=y124_1234$rowname,
                         column="ALIAS",
                         keytype="ENSEMBL",
                         multiVals="first") 

head(y124_1234)

write.csv(y124_1234, "MOK124_overlaps.csv") 

# BOth overlaps 
y_overlaps$Genename <- mapIds(org.Bt.eg.db,
                             keys=y_overlaps$rowname,
                             column="GENENAME",
                             keytype="ENSEMBL",
                             multiVals="first") 

y_overlaps$Alias <- mapIds(org.Bt.eg.db,
                          keys=y_overlaps$rowname,
                          column="ALIAS",
                          keytype="ENSEMBL",
                          multiVals="first") 

head(y_overlaps)

write.csv(y_overlaps, "Both_overlaps.csv") 

# From strain-specific pathway analysis file: Total overlaps strain specific --> these are not considered biomarkers! 
y23_total <- setdiff(y23_1234, y_overlaps)
y124_total <- setdiff(y124_1234, y_overlaps) 

# MOK023 specific 1234 

y23_total$Genename <- mapIds(org.Bt.eg.db,
                              keys=y23_total$rowname,
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="first") 

y23_total$Alias <- mapIds(org.Bt.eg.db,
                           keys=y23_total$rowname,
                           column="ALIAS",
                           keytype="ENSEMBL",
                           multiVals="first") 

head(y23_total)

write.csv(y23_total, "MOK023_allHPI_overlaps_specific.csv") 


# MOK124 specific 1234 

y124_total$Genename <- mapIds(org.Bt.eg.db,
                             keys=y124_total$rowname,
                             column="GENENAME",
                             keytype="ENSEMBL",
                             multiVals="first") 

y124_total$Alias <- mapIds(org.Bt.eg.db,
                          keys=y124_total$rowname,
                          column="ALIAS",
                          keytype="ENSEMBL",
                          multiVals="first") 

head(y124_total)

write.csv(y124_total, "MOK124_allHPI_overlaps_specific.csv") 

# Make quad Venns of strain-specific overlaps ##### 

# Data - vector reminders form file: 
# "Pathway_analysis_strain_specific_LF_NEW.R" 

# Strain specific for MOK023 
y23_24S <- setdiff(y23_24, y24H) 
nrow(y23_24S) # 135
y23_48S <- setdiff(y23_48, y48H)
y23_72S <- setdiff(y23_72, y72H)
y23_168S <- setdiff(y23_168, y168H)

# strain specific for MOK124 
y124_24S <- setdiff(y124_24, y24H) 
y124_48S <- setdiff(y124_48, y48H)
y124_72S <- setdiff(y124_72, y72H)
y124_168S <- setdiff(y124_168, y168H)

# Total overlaps strain specific 
y23_total <- setdiff(y23_1234, y_overlaps)
y124_total <- setdiff(y124_1234, y_overlaps) 

# MOK023 
head(y23_24S) 

# Intersect all combinations! 
nrow(y23_24S) # 135
nrow(y23_48S) # 1142
y23_24_48S <- intersect(y23_24S, y23_48S) # this intersect requires the tidyverse!!!! 
nrow(y23_24_48S) # 84
nrow(y23_72S) # 1046 
nrow(y23_168S) # 707
y23_48_72S <- intersect(y23_48S, y23_72S)
nrow(y23_48_72S) # 567 
y23_24_72S <- intersect(y23_24S, y23_72S)
y23_72_168S <- intersect(y23_72S, y23_168S)
y23_24_168S <- intersect(y23_24S, y23_168S)
y23_48_168S <- intersect(y23_48S, y23_168S)
nrow(y23_48_168S) # 304 

#for more lists than 2
Reduce(intersect, list(a,b,c))

library(dplyr) 
# The following objects are masked from ‘package:base’:
# intersect, setdiff, setequal, union 

y23_1234S <- Reduce(intersect, list(y23_24S, y23_48S, y23_72S, y23_168S)) 
nrow(y23_1234S) #  33 
# In the Venn diagram of 1234 overlaps 23 vs 124 there is 272 genes on the MOK023 side - WHY? 

y23_123S <- Reduce(intersect, list(y23_24S, y23_48S, y23_72S))
nrow(y23_123S) # 55 
y23_234S <- Reduce(intersect, list(y23_48S, y23_72S, y23_168S))
nrow(y23_234S) # 277 
y23_124S <- Reduce(intersect, list(y23_24S, y23_48S, y23_168S))
nrow(y23_124S) # 35 
y23_134S <- Reduce(intersect, list(y23_24S, y23_72S, y23_168S))
nrow(y23_134S) # 37

# ** 1 Make quad Venn #####
#pink - this is the prettiest 
grid.newpage()
venn.plot23S <- draw.quad.venn(area1 = nrow(y23_24S), area2 = nrow(y23_48S), area3 = nrow(y23_72S), area4 = nrow(y23_168S), n12 = nrow(y23_24_48S), n13 = nrow(y23_24_72S), n14 = nrow(y23_24_168S), n23 = nrow(y23_48_72S), n24 = nrow(y23_48_168S), 
               n34 = nrow(y23_72_168S), n123 = nrow(y23_123S), n124 = nrow(y23_124S), n134 = nrow(y23_134S), n234 = nrow(y23_234S), n1234 = nrow(y23_1234S), category = c("24 HPI", "48 HPI", "72 HPI", "168 HPI"), lty = "blank", 
               fill = c("blueviolet", "cyan3", "hotpink1", "lightpink1"), 
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

# 2 Add Venn plot title 
require(gridExtra)
venn.plot23S <- grid.arrange(gTree(children=venn.plot23S), top=textGrob("MOK023 specific genes", gp=gpar(fontsize=20)))

# 3 Writing to file 
tiff(filename = "Quad_Venn_MOK023S.tiff", compression = "lzw");
grid.draw(venn.plot2);
dev.off(); 

Cairo(file="Quad_Venn_MOK023S.png", 
      type="png",
      units="in", 
      width=10, 
      height=10, 
      pointsize=12, 
      dpi=72)
grid.draw(venn.plot23S)
dev.off()

# MOK124 

# Intersect all combinations! 
nrow(y124_24S) # 1150
nrow(y124_48S) # 873
y124_24_48S <- intersect(y124_24S, y124_48S) # this intersect requires the tidyverse!!!! 
nrow(y124_24_48S) # 82
nrow(y124_72S) # 488
nrow(y124_168S) # 501 
y124_48_72S <- intersect(y124_48S, y124_72S)
nrow(y124_48_72S) # 359 
y124_24_72S <- intersect(y124_24S, y124_72S)
y124_72_168S <- intersect(y124_72S, y124_168S)
y124_24_168S <- intersect(y124_24S, y124_168S)
y124_48_168S <- intersect(y124_48S, y124_168S)
nrow(y124_48_168S) # 270 

#for more lists than 2
Reduce(intersect, list(a,b,c))

library(dplyr) 
# The following objects are masked from ‘package:base’:
# intersect, setdiff, setequal, union 

y124_1234S <- Reduce(intersect, list(y124_24S, y124_48S, y124_72S, y124_168S)) 
nrow(y124_1234S) #  18
# In the Venn diagram of 1234 overlaps 23 vs 124 there is 272 genes on the MOK023 side - WHY? 

y124_123S <- Reduce(intersect, list(y124_24S, y124_48S, y124_72S))
nrow(y124_123S) # 25
y124_234S <- Reduce(intersect, list(y124_48S, y124_72S, y124_168S))
nrow(y124_234S) # 198  
y124_124S <- Reduce(intersect, list(y124_24S, y124_48S, y124_168S))
nrow(y124_124S) # 23
y124_134S <- Reduce(intersect, list(y124_24S, y124_72S, y124_168S))
nrow(y124_134S) # 26

# 1 Make quad Venn 
#pink - this is the prettiest 
grid.newpage()
venn.plot124S <- draw.quad.venn(area1 = nrow(y124_24S), area2 = nrow(y124_48S), area3 = nrow(y124_72S), area4 = nrow(y124_168S), n12 = nrow(y124_24_48S), n13 = nrow(y124_24_72S), n14 = nrow(y124_24_168S), n23 = nrow(y124_48_72S), n24 = nrow(y124_48_168S), 
                            n34 = nrow(y124_72_168S), n123 = nrow(y124_123S), n124 = nrow(y124_124S), n134 = nrow(y124_134S), n234 = nrow(y124_234S), n1234 = nrow(y124_1234S), category = c("24 HPI", "48 HPI", "72 HPI", "168 HPI"), lty = "blank", 
                            fill = c("blueviolet", "cyan3", "hotpink1", "lightpink1"), 
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

# 2 Add Venn plot title 
require(gridExtra)
venn.plot124S <- grid.arrange(gTree(children=venn.plot124S), top=textGrob("MOK124 specific genes", gp=gpar(fontsize=20)))

# 3 Writing to file 
tiff(filename = "Quad_Venn_MOK124S.tiff", compression = "lzw");
grid.draw(venn.plot2);
dev.off(); 

Cairo(file="Quad_Venn_MOK124S.png", 
      type="png",
      units="in", 
      width=10, 
      height=10, 
      pointsize=12, 
      dpi=72)
grid.draw(venn.plot124S)
dev.off() 

# Both strain specific plots on a grid 
venn.plotS<-grid.arrange(venn.plot23S, venn.plot124S, 
                        ncol = 2, 
                        widths = c(1,1),
                        heights = 1)  

grid.newpage() 
tiff(filename = "Quad_Venn_diagram_both_strains.tiff", compression = "lzw");
grid.draw(venn.plot);
dev.off(); 

Cairo(file="Quad_Venn_diagram_specific_both_strains.png", 
      type="png",
      units="in", 
      width=10, 
      height=5, 
      pointsize=12, 
      dpi=72)
grid.draw(venn.plotS)
dev.off()


# Annotate specific overlaps from the Venns above ####

library("AnnotationDbi")
library("OrganismDbi") 
library("org.Bt.eg.db")
columns(org.Bt.eg.db)
keytypes(org.Bt.eg.db)

typeof(y23_1234$rowname) #"character"

y23_1234S$Genename <- mapIds(org.Bt.eg.db,
                            keys=y23_1234S$rowname,
                            column="GENENAME",
                            keytype="ENSEMBL",
                            multiVals="first") 

y23_1234S$Alias <- mapIds(org.Bt.eg.db,
                         keys=y23_1234S$rowname,
                         column="ALIAS",
                         keytype="ENSEMBL",
                         multiVals="first") 

head(y23_1234S)

write.csv(y23_1234S, "MOK023_overlaps_biomarkers.csv") 
write.table(y23_1234, "MOK023_overlaps.txt", row.names = y23_1234$rowname, sep = "\t") 

y124_1234S$Genename <- mapIds(org.Bt.eg.db,
                             keys=y124_1234S$rowname,
                             column="GENENAME",
                             keytype="ENSEMBL",
                             multiVals="first") 

y124_1234S$Alias <- mapIds(org.Bt.eg.db,
                          keys=y124_1234S$rowname,
                          column="ALIAS",
                          keytype="ENSEMBL",
                          multiVals="first") 

head(y124_1234S)

write.csv(y124_1234S, "MOK124_overlaps_biomarkers.csv") 

# COMMON BIOMARKERS ####
# Biomarkers would be a list of genes that are common for both strains at all HPI 
biomarkers <- Reduce(intersect, list(y23_24, y23_48, y23_72, y23_168, y124_24, y124_48, y124_72, y124_168)) 
nrow(biomarkers) 
# 130 

# Annotate this list 
library("AnnotationDbi")
library("OrganismDbi") 
library("org.Bt.eg.db")

typeof(biomarkers$rowname) #"character"

biomarkers$Genename <- mapIds(org.Bt.eg.db,
                            keys=biomarkers$rowname,
                            column="GENENAME",
                            keytype="ENSEMBL",
                            multiVals="first") 

biomarkers$Alias <- mapIds(org.Bt.eg.db,
                         keys=biomarkers$rowname,
                         column="ALIAS",
                         keytype="ENSEMBL",
                         multiVals="first") 

head(biomarkers)

write.csv(biomarkers, "Biomarkers.csv") 
nrow(intersect(biomarkers$rowname,y_overlaps$rowname)) # NULL 
# They are the same! 

# Table of biomarkers with all log2foldchanges 
head(logs_MOK023)
head(logs_MOK124)
colnames(logs_MOK023) 
# rowname has to go to column again 
logs_MOK023_rname <- tibble::rownames_to_column(logs_MOK023)
logs_MOK124_rname <- tibble::rownames_to_column(logs_MOK124)
colnames(logs_MOK023_rname) 
colnames(biomarkers)
biomarkers_logs <- inner_join(biomarkers, logs_MOK023_rname, by = "rowname", suffix = c("_biomark", "_MOK023"))
# No need for suffixes above cause only rowname is the same column name 
biomarkers_logs <- inner_join(biomarkers_logs, logs_MOK124_rname, by = "rowname", suffix = c("_MOK023", "_MOK124"))
nrow(biomarkers_logs) 
colnames(biomarkers_logs) 
write.csv(biomarkers_logs, "Biomarkers_logs.csv") 

# How does intersect work ? 
aa <- c(1,2,3,4,5) # Using aa cause a is already taken 
bb <- c(4,5,6,7)
cc <- c(5,4,3,2,1)
dd <- c(7,6,5,4) 
aa
intersect(aa,bb)
intersect(cc,dd) 
# Intersected in the order in which they appear 
# Because all vectors are ordered by pvalue, overlapped vectors are intersected by significance (sort of) 

# Add log2 fold changes to strain specific biomarker data #### 
# Make a logs file of strain specific biomarkers by merging with LF_sig_both - list of all DE gene log2 fold changes and adjusted p values 
head(y23_1234S) 
nrow(y23_1234S) 

S23_biomarkers_logs <- inner_join(y23_1234S, LF_sig_both, by = "rowname") 
View(S23_biomarkers_logs) 
write.csv(S23_biomarkers_logs, "MOK023_biomarkers_logs.csv")

head(y124_1234S) 
nrow(y124_1234S) 

S124_biomarkers_logs <- inner_join(y124_1234S, LF_sig_both, by = "rowname") 
View(S124_biomarkers_logs) 
write.csv(S124_biomarkers_logs, "MOK124_biomarkers_logs.csv")

# Plotting normalised gene counts for biomarkers ####

# Common 130  
# In MOK023 
library(DESeq2) 
library(ggplot2)

nrow(y_overlaps) 
head(y_overlaps) 


geneCounts <- plotCounts(dds, gene = "ENSBTAG00000018156", intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE) 
geneCounts$Cow_ID <- as.factor(geneCounts$Cow_ID) 
plot <- ggplot(geneCounts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID, shape = Group)) +
  scale_y_log10() + geom_point(size = 2) + geom_line() 
plot + ggtitle(paste("ENSBTAG00000018156", y_overlaps$Alias[which(grepl("ENSBTAG00000018156", y_overlaps$rowname))], sep = " "))
plot 

which(grepl("ENSBTAG00000018156", y_overlaps$rowname)) 
?grep

pdf("Common_biomarkers_normalised counts_MOK023.pdf")
for (i in y_overlaps$rowname) {
  geneCounts <- plotCounts(dds23, gene = i, intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE) 
  geneCounts$Cow_ID <- as.factor(geneCounts$Cow_ID) 
  plot <- ggplot(geneCounts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) +
    scale_y_log10() + geom_point(size = 3) + geom_line() 
    print(plot + ggtitle(paste(i, y_overlaps$Alias[which(grepl(i, y_overlaps$rowname))], sep = " ")))
}
dev.off() 

# In MOK124 
pdf("Common_biomarkers_normalised counts_MOK124.pdf")
for (i in y_overlaps$rowname) {
  geneCounts <- plotCounts(dds124, gene = i, intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE) 
  geneCounts$Cow_ID <- as.factor(geneCounts$Cow_ID) 
  plot <- ggplot(geneCounts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) +
    scale_y_log10() + geom_point(size = 3) + geom_line() 
  print(plot + ggtitle(paste(i, y_overlaps$Alias[which(grepl(i, y_overlaps$rowname))], sep = " ")))
}
dev.off() 


# In BOTH!!!! 
pdf("Common_biomarkers_normalised counts_BOTH.pdf")
for (i in y_overlaps$rowname) {
  geneCounts <- plotCounts(dds, gene = i, intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE) 
  geneCounts$Cow_ID <- as.factor(geneCounts$Cow_ID) 
  plot <- ggplot(geneCounts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID, shape = Group)) +
    scale_y_log10() + geom_point(size =2) + geom_line() 
  print(plot + ggtitle(paste(i, y_overlaps$Alias[which(grepl(i, y_overlaps$rowname))], sep = " ")))
}
dev.off() 

# MOK023 specific 
head(y23_1234S) 

pdf("MOK023_biomarkers_normalised counts.pdf")
for (i in y23_1234S$rowname) {
  geneCounts <- plotCounts(dds23, gene = i, intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE) 
  geneCounts$Cow_ID <- as.factor(geneCounts$Cow_ID) 
  plot <- ggplot(geneCounts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) +
    scale_y_log10() + geom_point(size = 3) + geom_line() 
  print(plot + ggtitle(paste(i, y23_1234S$Alias[which(grepl(i, y23_1234S$rowname))], sep = " ")))
}
dev.off() 

# MOK124 specific 
pdf("MOK124_biomarkers_normalised counts.pdf")
for (i in y124_1234S$rowname) {
  geneCounts <- plotCounts(dds124, gene = i, intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE) 
  geneCounts$Cow_ID <- as.factor(geneCounts$Cow_ID) 
  plot <- ggplot(geneCounts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) +
    scale_y_log10() + geom_point(size = 3) + geom_line() 
  print(plot + ggtitle(paste(i, y124_1234S$Alias[which(grepl(i, y124_1234S$rowname))], sep = " ")))
}
dev.off() 
