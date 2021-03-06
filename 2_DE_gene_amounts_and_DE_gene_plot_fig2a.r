 
# Aim of this file - to create a data table of numbers of up and downregulated DE genes and then plot it 

#MOK023 #####

# Reminder of extracted values 
res23_24h_LF <- results(dds23, lfcThreshold=1, contrast = c("HPI", "0", "24"))
res23_48h_LF <- results(dds23, lfcThreshold=1, contrast = c("HPI", "0", "48"))
res23_72h_LF <- results(dds23, lfcThreshold=1, contrast = c("HPI", "0", "72"))
res23_168h_LF <- results(dds23, lfcThreshold=1, contrast = c("HPI", "0", "168")) 

sum(res23_24h_LF$padj < 0.05, na.rm=TRUE) # 1278
sum(res23_48h_LF$padj < 0.05, na.rm=TRUE) # 2248 
sum(res23_72h_LF$padj < 0.05, na.rm=TRUE) # 1986
sum(res23_168h_LF$padj < 0.05, na.rm=TRUE) # 1750 

res23_24h_LF_sig <- subset(res23_24h_LF, padj < 0.05)
res23_48h_LF_sig <- subset(res23_48h_LF, padj < 0.05)
res23_72h_LF_sig <- subset(res23_72h_LF, padj < 0.05)
res23_168h_LF_sig <- subset(res23_168h_LF, padj < 0.05)

# Get DE numbers! 
#All DE genes 

a <- nrow(subset(res23_24h_LF_sig, log2FoldChange > 0))
b <- nrow(subset(res23_24h_LF_sig, log2FoldChange < 0))
c <- nrow(subset(res23_48h_LF_sig, log2FoldChange > 0))
d <- nrow(subset(res23_48h_LF_sig, log2FoldChange < 0))
e <- nrow(subset(res23_72h_LF_sig, log2FoldChange > 0))
f <- nrow(subset(res23_72h_LF_sig, log2FoldChange < 0))
g <- nrow(subset(res23_168h_LF_sig, log2FoldChange > 0))
h <- nrow(subset(res23_168h_LF_sig, log2FoldChange < 0))
All_numbers23 <- c(a,b,c,d,e,f,g,h)
All_numbers23 
All_labels <- c("24_up", "24_down", "48_up", "48_down", "72_up", "72_down", "168_up", "168_down") 
LFnumbers_23 <- cbind(All_labels, All_numbers23, HPI, Group2, Direction)
HPI <- c(24, 24, 48, 48, 72, 72, 168, 168)
Group2 <- c("MOK023")
Group2 <- rep(Group2, 8)
Direction <- c("Up", "Down")
Direction <- rep(Direction, 4)
LFnumbers_23 
Group2

#MOK124 ##### 

res124_24h_LF <- results(dds124, lfcThreshold=1, contrast = c("HPI", "0", "24"))
res124_48h_LF <- results(dds124, lfcThreshold=1, contrast = c("HPI", "0", "48"))
res124_72h_LF <- results(dds124, lfcThreshold=1, contrast = c("HPI", "0", "72"))
res124_168h_LF <- results(dds124, lfcThreshold=1, contrast = c("HPI", "0", "168")) 

sum(res124_24h_LF$padj < 0.05, na.rm=TRUE) # 2293
sum(res124_48h_LF$padj < 0.05, na.rm=TRUE) # 1979 
sum(res124_72h_LF$padj < 0.05, na.rm=TRUE) # 1428
sum(res124_168h_LF$padj < 0.05, na.rm=TRUE) # 1544 

res124_24h_LF_sig <- subset(res124_24h_LF, padj < 0.05)
res124_48h_LF_sig <- subset(res124_48h_LF, padj < 0.05)
res124_72h_LF_sig <- subset(res124_72h_LF, padj < 0.05)
res124_168h_LF_sig <- subset(res124_168h_LF, padj < 0.05) 

# Get DE numbers! 
#All DE genes 

a2 <- nrow(subset(res124_24h_LF_sig, log2FoldChange > 0))
b2 <- nrow(subset(res124_24h_LF_sig, log2FoldChange < 0))
c2 <- nrow(subset(res124_48h_LF_sig, log2FoldChange > 0))
d2 <- nrow(subset(res124_48h_LF_sig, log2FoldChange < 0))
e2 <- nrow(subset(res124_72h_LF_sig, log2FoldChange > 0))
f2 <- nrow(subset(res124_72h_LF_sig, log2FoldChange < 0))
g2 <- nrow(subset(res124_168h_LF_sig, log2FoldChange > 0))
h2 <- nrow(subset(res124_168h_LF_sig, log2FoldChange < 0))
All_numbers <- c(a2,b2,c2,d2,e2,f2,g2,h2)
All_numbers 
All_labels <- c("24_up", "24_down", "48_up", "48_down", "72_up", "72_down", "168_up", "168_down") 
LFnumbers_124 <- cbind(All_labels, All_numbers, HPI, Group, Direction)
HPI <- c(24, 24, 48, 48, 72, 72, 168, 168)
Group <- c("MOK124")
Group <- rep(Group, 8)
Direction <- c("Up", "Down")
Direction <- rep(Direction, 4) 
a2+b2

# Plot DE gene amounts #####

LFnumbers <- rbind(LFnumbers_23, LFnumbers_124)
colnames(LFnumbers) <- c("Label", "Count", "HPI", "Strain", "Direction")
write.csv(LFnumbers, "LFnumbers.csv") 

LFnumbers2 <- read.csv("LFnumbers.csv")

library(ggplot2) 

library(ggplot2) 
library(plyr)
library(RColorBrewer) 

# Up and downregulated gene subsets 
LFnumbers2_Up <- subset(LFnumbers2, Direction == "Up") 
LFnumbers2_Down <- subset(LFnumbers2, Direction == "Down") 

# ** My pet graph - DE gene numbers for paper ####  
ggplot(data = LFnumbers2_Up) + 
  geom_col(data = LFnumbers2_Up, 
           mapping = aes(x = HPI, y = Count, fill = Direction), position = "dodge") + 
  geom_text(data = LFnumbers2_Up, 
            mapping = aes(x = HPI, y = Count, group = Strain, label = Count), vjust = 1.2) + 
  facet_wrap(LFnumbers2_Up$Strain) + 
  geom_col(data = LFnumbers2_Down, 
           mapping = aes(x = HPI, y = -Count, fill = Direction, group = Strain), position = "dodge") + 
  geom_text(data = LFnumbers2_Down, 
          mapping = aes(x = HPI, y = -Count, group = Strain, label = Count), vjust = -0.5) + 
  facet_wrap(LFnumbers2_Up$Strain) + 
  scale_fill_brewer(palette = "Set1", 
                      name="Direction",
                      breaks=c("Up", "Down"),
                      labels=c("Upregulated", "Downregulated")) + 
  xlab("Hours post infection (hpi)") + 
  ylab("Number of genes") + 
  theme(text = element_text(size = 14, family = "Calibri")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("DE_gene_for_paper.jpeg")

