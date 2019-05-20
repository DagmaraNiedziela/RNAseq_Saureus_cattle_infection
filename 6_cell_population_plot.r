# Import cell results from Paul ####


library(readr)
dagmara_cell_results <- read_csv("dagmara_cell.results.csv")
View(dagmara_cell_results) 
# The file is organised like countdata in DESeq - each sample is a column, and column names are sample names, rows are cell types 

# Create a file with long format 
library(reshape2)
cell_results_long <- melt(dagmara_cell_results) 
?melt
colnames(dagmara_cell_results)
View(cell_results_long)
colnames(cell_results_long) <- c("Cell_type", "Sample", "Proportion")  
nrow(cell_results_long)

# Import a file that describes samples - contains sample name, and then columns with cow, group, HPI, SCC 
sample_data <- read.csv("Sample_data.csv")
head(sample_data)
colnames(sample_data)

# Join with sample data ####
library(dplyr)
sample_data <- sample_data %>% select(Sample, Cow_ID, HPI, Group, SCC) 

cell_results_long_samples <- inner_join(cell_results_long, sample_data)
View(cell_results_long_samples)
nrow(cell_results_long_samples) 
colnames(cell_results_long_samples)
write.csv(cell_results_long_samples, "Cell_results_long_samples_SCC.csv")

library(ggplot2)
library(RColorBrewer)
brewer.pal.info 
# Pick a color palette 
display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
                   colorblindFriendly=FALSE) 
display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
                   colorblindFriendly=TRUE)

# Because here we don't want to see a time course, I set HPI to factor 
cell_results_long_samples$HPI <- as.factor(cell_results_long_samples$HPI)  
levels(cell_results_long_samples$HPI) 

levels(cell_results_long$Cell_type) 
# If I need a time course later, I can use + scale_x_continuous(breaks = c(0,24,48,72,168)) 

# There is 11 cell types - make Other gray, have a manual scale 
# From http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12 - this is the Paired palette: 
# '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928' 
# #ff7f00 is the orange for Other 
# grey #808080 
# '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#808080','#cab2d6','#6a3d9a','#ffff99','#b15928' 

my_cell_colours <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#808080','#cab2d6','#6a3d9a','#ffff99','#b15928') 
my_cell_colours2 <- c('#b2df8a','#33a02c','#e31a1c','#fdbf6f','#808080','#cab2d6','#6a3d9a','#ffff99','#b15928')

# Make cell averages
cell_averages <- cell_results_long_samples %>% group_by(HPI, Group, Cell_type) %>% 
  summarize(Proportion_mean = mean(Proportion, na.rm = TRUE)) 
View(cell_averages) 
write.csv(cell_averages, "cell_averages.csv")

# ** Average proportions ####
ggplot(data = cell_averages) + 
  geom_col(data = cell_averages, 
           mapping = aes(x = HPI, y = Proportion_mean, fill = Cell_type)) + 
  facet_wrap(cell_averages$Group) + 
  scale_fill_manual(values = my_cell_colours, 
                    name="Cell type") + 
  xlab("Hours post infection (hpi)") + 
  ylab("Proportion of cells") + 
  theme(text = element_text(size = 14, family = "Calibri")) 
?scale_color_manual

ggsave("Cell_proportions_group.jpeg") 

# RNAseq cell numbers - with SCC #### 
View(cell_results_long_samples)
nrow(cell_results_long_samples) 
colnames(cell_results_long_samples) 

library(dplyr) 

cell_results_long_abs <- cell_results_long_samples %>% mutate(Absolute = Proportion * SCC * 1000) 
write.csv(cell_results_long_abs, "cell_results_long_abs.csv")

# Make cell averages
cell_abs_averages <- cell_results_long_abs %>% group_by(HPI, Group, Cell_type) %>% 
  summarize(Absolute_mean = mean(Absolute, na.rm = TRUE)) 
View(cell_abs_averages) 
write.csv(cell_abs_averages, "cell_abs_averages.csv") 

library(ggplot2)

ggplot(data = cell_abs_averages) + 
  geom_col(data = cell_abs_averages, mapping = aes(x = HPI, y = Absolute_mean, fill = Cell_type), position ="dodge") + 
  facet_wrap(cell_averages$Group) + 
  scale_fill_manual(values = my_cell_colours, 
                    name="Cell type") + 
  xlab("Hours post infection (hpi)") + 
  ylab("Absolute cell count (cells/ml)") + 
  theme(text = element_text(size = 14, family = "Calibri")) 

ggsave("Cell_absolutes_group.jpeg") 

# Log scale - this is the plot used for thesis 
ggplot(data = cell_abs_averages) + 
  geom_col(data = cell_abs_averages, mapping = aes(x = HPI, y = Absolute_mean, fill = Cell_type), position ="dodge") + 
  facet_wrap(cell_averages$Group) + 
  scale_fill_manual(values = my_cell_colours, 
                    name="Cell type") + 
  xlab("Hours post infection (hpi)") + scale_y_log10() +
  ylab("Absolute cell count (log10 cells/ml)") + 
  theme(text = element_text(size = 14, family = "Calibri")) 

ggsave("Cell_absolutes_group_log10.jpeg") 

# Select only neutrophils, m1 and m2 macrophages - this not used for thesis but was used to have a better look 
cell_abs_averages3 <- cell_abs_averages %>% filter(Cell_type %in% c("Neutrophil", "Macrophage M1", "Macrophage M2")) 

ggplot(data = cell_abs_averages3) + 
  geom_col(data = cell_abs_averages3, mapping = aes(x = HPI, y = Absolute_mean, fill = Cell_type), position ="dodge") + 
  facet_wrap(cell_abs_averages3$Group) + 
  scale_fill_manual(values = my_cell_colours2, 
                    name="Cell type") + 
  xlab("Hours post infection (hpi)") + scale_y_log10() +
  ylab("Absolute cell count (log10 cells/ml)") + 
  theme(text = element_text(size = 14, family = "Calibri")) 

ggsave("Cell_absolutes_group_log10_neuts_and_monos.jpeg") 

# Plot together ####
library(ggplot2) 
library(ggpubr)
plotA <- ggplot(data = cell_averages) + 
  geom_col(data = cell_averages, 
           mapping = aes(x = HPI, y = Proportion_mean, fill = Cell_type)) + 
  facet_wrap(cell_averages$Group) + 
  scale_fill_manual(values = my_cell_colours, 
                    name="Cell type") + 
  xlab("Hours post infection (hpi)") + 
  ylab("Proportion of cells") + 
  theme(text = element_text(size = 14, family = "Calibri")) 

# A modification for the final plot - remove all absolute values below 1 (less than 1 cell/ml) by changing them to 0 - to make sure that they stay on the plot
library(dplyr)
colnames(cell_abs_averages) 
cell_abs_averages4 <- cell_abs_averages 
typeof(cell_abs_averages4$Absolute_mean)
cell_abs_averages4$Absolute_mean2 <- ifelse(cell_abs_averages4$Absolute_mean > 1, cell_abs_averages4$Absolute_mean, 0)
View(cell_abs_averages4) 
# Changing to 1 worked cause log 1 is 0 
cell_abs_averages4$Absolute_mean[cell_abs_averages4$Absolute_mean<1]=1

plotB <- ggplot(data = cell_abs_averages4) + 
  geom_col(data = cell_abs_averages4, mapping = aes(x = HPI, y = Absolute_mean, fill = Cell_type), position ="dodge") + 
  facet_wrap(cell_abs_averages4$Group) + 
  scale_fill_manual(values = my_cell_colours, 
                    name="Cell type") + 
  xlab("Hours post infection (hpi)") + scale_y_log10(breaks = c(10, 1000, 100000), labels = c(1,3,5)) +
  ylab("Absolute cell count (log10 cells/ml)") + 
  theme(text = element_text(size = 14, family = "Calibri")) 

plotAB <- ggarrange(plotA, plotB, 
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2) 
plotAB
ggsave("Cell_populations_AB.jpeg") 
