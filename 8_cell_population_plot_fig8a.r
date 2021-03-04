# Import cell results ####

library(tidyverse)
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
nrow(cell_results_long) #550 

# Import a file that describes samples - contains sample name, and then columns with cow, group, HPI, SCC 
sample_data <- read.csv("Sample_data.csv")
head(sample_data)
colnames(sample_data)

# Join with sample data ####
library(dplyr)
sample_data <- sample_data %>% select(Sample, Cow_ID, HPI, Group, SCC) 

cell_results_long_samples <- inner_join(cell_results_long, sample_data)
View(cell_results_long_samples)
nrow(cell_results_long_samples) #550
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

levels(cell_results_long$Cell_type) # NULL 
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
plotA2 <- ggplot(data = cell_averages) + 
  geom_col(data = cell_averages, 
           mapping = aes(x = HPI, y = Proportion_mean, fill = Cell_type)) + 
  facet_wrap(~cell_averages$Group) + 
  scale_fill_manual(values = my_cell_colours, 
                    name="Cell type") + 
  xlab("Hours post infection (hpi)") + 
  ylab("Proportion of cells") + 
  theme(text = element_text(size = 14, family = "Calibri")) 

ggsave("Cell_proportions_group.jpeg") 

# ** Add epithelial cell marker plot here #### 
# plot AC - plot2 comes from file 9 - epithelial cell markers 
# Add null space between the two plots 

plotAC <- ggarrange(plotA2, NULL, plot2, 
                    labels = c("A", "", "B"), heights = c(1,0.05,1.1),
                    ncol = 1, nrow = 3) 
plotAC
ggsave("Cell_populations_AC.jpeg", width = 24.2, height = 28, units = "cm", dpi = 1200) 

