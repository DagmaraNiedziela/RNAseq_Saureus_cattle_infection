# Top 5 genes up and down genes from each time point - MOK023 ####
head(res23_24h_LF_sig) 
head(res23_24h_LF_sigDF) 
res23_24h_LF_sig_top10_up 
?top_n

library(tidyverse) 
library(RColorBrewer) 
library(ggrepel) 

res23_24h_LF_sig_top5_up <- res23_24h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK023", 6), HPI = rep("24", 6), Direction = rep("Up", 6))
res23_24h_LF_sig_top5_down <- res23_24h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK023", 5), HPI = rep("24", 5), Direction = rep("Down", 5))

res23_48h_LF_sig_top5_up <- res23_48h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK023", 5), HPI = rep("48", 5), Direction = rep("Up", 5))
res23_48h_LF_sig_top5_down <- res23_48h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK023", 5), HPI = rep("48", 5), Direction = rep("Down", 5)) 

res23_72h_LF_sig_top5_up <- res23_72h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK023", 5), HPI = rep("72", 5), Direction = rep("Up", 5))
res23_72h_LF_sig_top5_down <- res23_72h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK023", 5), HPI = rep("72", 5), Direction = rep("Down", 5)) 

res23_168h_LF_sig_top5_up <- res23_168h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK023", 5), HPI = rep("168", 5), Direction = rep("Up", 5))
res23_168h_LF_sig_top5_down <- res23_168h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK023", 5), HPI = rep("168", 5), Direction = rep("Down", 5))

res23_top5 <- rbind(res23_24h_LF_sig_top5_up, res23_24h_LF_sig_top5_down, res23_48h_LF_sig_top5_up, res23_48h_LF_sig_top5_down, 
                     res23_72h_LF_sig_top5_up, res23_72h_LF_sig_top5_down, res23_168h_LF_sig_top5_up, res23_168h_LF_sig_top5_down)
View(res23_top5) 
lapply(res23_top5, typeof) 
colnames(res23_top5) 
res23_top5$HPI <- as.numeric(as.character(res23_top5$HPI)) 
res23_top5$HPI <- as.factor(res23_top5$HPI)

# plot - geom point with log2foldchange and Alias as a label 
ggplot(data = res23_top5) + 
  geom_point(mapping = aes(x = HPI, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Hours post infection (hpi)") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  geom_text_repel(mapping = aes(x = HPI, y = log2FoldChange, group = Direction, label = Alias)) 
ggsave("MOK023_top5.jpeg")


# Top 5 genes up and down genes from each time point. ####
head(res124_24h_LF_sig) 
head(res124_24h_LF_sigDF) 
res124_24h_LF_sig_top10_up 
?top_n

library(tidyverse) 
library(RColorBrewer) 
library(ggrepel) 

res124_24h_LF_sig_top5_up <- res124_24h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK0124", 5), HPI = rep("24", 5), Direction = rep("Up", 5))
res124_24h_LF_sig_top5_down <- res124_24h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK0124", 5), HPI = rep("24", 5), Direction = rep("Down", 5))

res124_48h_LF_sig_top5_up <- res124_48h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK0124", 5), HPI = rep("48", 5), Direction = rep("Up", 5))
res124_48h_LF_sig_top5_down <- res124_48h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK0124", 5), HPI = rep("48", 5), Direction = rep("Down", 5)) 

res124_72h_LF_sig_top5_up <- res124_72h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK0124", 5), HPI = rep("72", 5), Direction = rep("Up", 5))
res124_72h_LF_sig_top5_down <- res124_72h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK0124", 5), HPI = rep("72", 5), Direction = rep("Down", 5)) 

res124_168h_LF_sig_top5_up <- res124_168h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange > 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK0124", 5), HPI = rep("168", 5), Direction = rep("Up", 5))
res124_168h_LF_sig_top5_down <- res124_168h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(log2FoldChange < 0) %>% top_n(-5, padj) %>% select(-baseMean, -lfcSE, -stat, -pvalue) %>% 
  add_column(Group = rep("MOK0124", 5), HPI = rep("168", 5), Direction = rep("Down", 5))

res124_top5 <- rbind(res124_24h_LF_sig_top5_up, res124_24h_LF_sig_top5_down, res124_48h_LF_sig_top5_up, res124_48h_LF_sig_top5_down, 
                    res124_72h_LF_sig_top5_up, res124_72h_LF_sig_top5_down, res124_168h_LF_sig_top5_up, res124_168h_LF_sig_top5_down)
View(res124_top5) 
lapply(res124_top5, typeof) 
colnames(res124_top5) 
res124_top5$HPI <- as.numeric(as.character(res124_top5$HPI)) 
res124_top5$HPI <- as.factor(res124_top5$HPI) 

write.csv(res124_top5, "MOK124_top5_LF_sig_genes.csv")

# plot - geom point with log2foldchange and Alias as a label 
ggplot(data = res124_top5) + 
  geom_point(mapping = aes(x = HPI, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Hours post infection (hpi)") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  geom_text_repel(mapping = aes(x = HPI, y = log2FoldChange, group = Direction, label = Alias)) 
ggsave("MOK124_top5.jpeg") 

# Make these one plot !!!! ####
# I make some cosmetic changes - italics for gene names 

plotA <- ggplot(data = res23_top5) + 
  geom_point(mapping = aes(x = HPI, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Hours post infection (hpi)") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  geom_text_repel(mapping = aes(x = HPI, y = log2FoldChange, group = Direction, label = Alias), fontface = "italic") 

plotB <- ggplot(data = res124_top5) + 
  geom_point(mapping = aes(x = HPI, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Hours post infection (hpi)") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  geom_text_repel(mapping = aes(x = HPI, y = log2FoldChange, group = Direction, label = Alias), fontface = "italic") 

library(ggpubr) 
plotAB <- ggarrange(plotA, plotB, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2) 
plotAB
ggsave("top5_both_strains.jpeg") 



# Milk protein genes - both strains ###### 
library(tidyverse) 

milk_genes <- c("LEPR", "LEP", "IGF1", "ABCG2", "OPN", "PPARGC1A", "CSN1", "CSN2", "CSN1S2", "CSN10", "DGAT1", "BDNF", "FTO", "GHR", "PRLR", "PRL",  "AGPAT6", "a-LACTA") 

res23_24h_LF_sig_milk <- res23_24h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(Alias %in% milk_genes) 
res23_48h_LF_sig_milk <- res23_48h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(Alias %in% milk_genes) 
res23_72h_LF_sig_milk <- res23_72h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(Alias %in% milk_genes) 
res23_168h_LF_sig_milk <- res23_168h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(Alias %in% milk_genes) 

res124_24h_LF_sig_milk <- res124_24h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(Alias %in% milk_genes) 
res124_48h_LF_sig_milk <- res124_48h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(Alias %in% milk_genes) 
res124_72h_LF_sig_milk <- res124_72h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(Alias %in% milk_genes) 
res124_168h_LF_sig_milk <- res124_168h_LF_sigDF %>% tibble::rownames_to_column() %>% filter(Alias %in% milk_genes) 

# For a heatmap I need a matrix, row names are alias, and then sample is a column - logs FC is a value inside 
# Join files with suffixes, and only select log2FC 
# Join all files into one 
milk_logs_23 <- full_join(res23_24h_LF_sig_milk, res23_48h_LF_sig_milk, by="Alias", suffix = c("_24h", "_48h"))
milk_logs_23 <- full_join(milk_logs_23, res23_72h_LF_sig_milk, by = "Alias", suffix = c("_48h", "_72h"))
milk_logs_23 <- full_join(milk_logs_23, res23_168h_LF_sig_milk, by = "Alias", suffix = c("_72h", "_168h")) %>% select(Alias, contains("log2FoldChange"))
View(milk_logs_23) 
write.csv(milk_logs_23, "Milk_DEgene_logs_MOK023.csv")

milk_logs_124 <- full_join(res124_24h_LF_sig_milk, res124_48h_LF_sig_milk, by="Alias", suffix = c("_24h", "_48h"))
milk_logs_124 <- full_join(milk_logs_124, res124_72h_LF_sig_milk, by = "Alias", suffix = c("_48h", "_72h"))
milk_logs_124 <- full_join(milk_logs_124, res124_168h_LF_sig_milk, by = "Alias", suffix = c("_72h", "_168h")) %>% select(Alias, contains("log2FoldChange"))
View(milk_logs_124) 
write.csv(milk_logs_124, "Milk_DEgene_logs_MOK124.csv")

milk_logs_both <- full_join(milk_logs_23, milk_logs_124, by = "Alias", suffix = c("_MOK023", "_MOK124"))
write.csv(milk_logs_both, "Milk_logs_both.csv")
# A table was then included in the thesis 
