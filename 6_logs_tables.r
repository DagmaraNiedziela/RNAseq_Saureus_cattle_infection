# Make LOGS tables of my significant genes ##### 
# LOGS is a data frame that has log2 fold changes only and ENSEMBL Gene IDs as row names - this will be used for pathview 
# Use tidyverse 
library(DESeq2)
library(dplyr) 
library(tidyverse)
detach(package:pathview) # pathview and dplyr have conflicting select functions, so I sometimes have to detach one or the other to progress with the code 
detach(package:dplyr)

head(res23_168h_LF_sigDF)
nrow(res23_24h_LF_sig) 
colnames(res23_24h_LF_sig) 

logs_23_24h <- res23_24h_LF_sigDF %>% select(log2FoldChange) 
head(logs_23_48h)
logs_23_48h <- res23_48h_LF_sigDF %>% select(log2FoldChange) 
logs_23_72h <- res23_72h_LF_sigDF %>% select(log2FoldChange) 
logs_23_168h <- res23_168h_LF_sigDF %>% select(log2FoldChange) 

# Convert row names to a new column called name? 
head(logs_23_48h)
logs_23_24h <- tibble::rownames_to_column(logs_23_24h)
logs_23_48h <- tibble::rownames_to_column(logs_23_48h)
logs_23_72h <- tibble::rownames_to_column(logs_23_72h)
logs_23_168h <- tibble::rownames_to_column(logs_23_168h) 

# Join time point per strain!!!! All single files still have rowname kept 

# 24hpi 
logs_24hpi <- full_join(logs_23_24h, logs_124_24h, by="rowname", suffix = c("_MOK023", "_MOK124"))
head(logs_24hpi) 
# Create rownames 
rownames(logs_24hpi) <- logs_24hpi$rowname 
logs_24hpi <- logs_24hpi[,-1] 

# 48 hpi 
logs_48hpi <- full_join(logs_23_48h, logs_124_48h, by="rowname", suffix = c("_MOK023", "_MOK124"))
head(logs_48hpi) 
# Create rownames
rownames(logs_48hpi) <- logs_48hpi$rowname 
logs_48hpi <- logs_48hpi[,-1]  

# 72 hpi 
logs_72hpi <- full_join(logs_23_72h, logs_124_72h, by="rowname", suffix = c("_MOK023", "_MOK124"))
head(logs_72hpi) 
# Create rownames
rownames(logs_72hpi) <- logs_72hpi$rowname 
logs_72hpi <- logs_72hpi[,-1] 

# 168 hpi 
logs_168hpi <- full_join(logs_23_168h, logs_124_168h, by="rowname", suffix = c("_MOK023", "_MOK124"))
head(logs_168hpi) 
# Create rownames 
rownames(logs_168hpi) <- logs_168hpi$rowname 
logs_168hpi <- logs_168hpi[,-1] 

# Join all log2fold changes 
logs_MOK023 <- full_join(logs_23_24h, logs_23_48h, by="rowname", suffix = c("_24h", "_48h"))
logs_MOK023 <- full_join(logs_MOK023, logs_23_72h, by = "rowname", suffix = c("_48h", "_72h"))
logs_MOK023 <- full_join(logs_MOK023, logs_23_168h, by = "rowname", suffix = c("_72h", "_168h"))
head(logs_MOK023) 
nrow(logs_MOK023) #7841 
colnames(logs_MOK023)
write.csv(logs_MOK023, "logs_MOK023.csv") 

# Create rownames 
rownames(logs_MOK023) <- logs_MOK023$rowname 
logs_MOK023 <- logs_MOK023[,-1] 

# Logs files for MOK124 #####
head(res124_24h_LF_sigDF)

logs_124_24h <- res124_24h_LF_sigDF %>% select(log2FoldChange) 
logs_124_48h <- res124_48h_LF_sigDF %>% select(log2FoldChange) 
logs_124_72h <- res124_72h_LF_sigDF %>% select(log2FoldChange) 
logs_124_168h <- res124_168h_LF_sigDF %>% select(log2FoldChange) 

# Convert row names to a new column called name? no need 
head(logs_124_48h)
logs_124_24h <- tibble::rownames_to_column(logs_124_24h)
logs_124_48h <- tibble::rownames_to_column(logs_124_48h)
logs_124_72h <- tibble::rownames_to_column(logs_124_72h)
logs_124_168h <- tibble::rownames_to_column(logs_124_168h) 

# Join all log2fold changes 
logs_MOK124 <- full_join(logs_124_24h, logs_124_48h, by="rowname", suffix = c("_24h", "_48h"))
logs_MOK124 <- full_join(logs_MOK124, logs_124_72h, by = "rowname", suffix = c("_48h", "_72h"))
logs_MOK124 <- full_join(logs_MOK124, logs_124_168h, by = "rowname", suffix = c("_72h", "_168h"))
head(logs_MOK124) 
nrow(logs_MOK124) #7841 
colnames(logs_MOK124)
write.csv(logs_MOK124, "logs_MOK124.csv") 

# Create rownames 
rownames(logs_MOK124) <- logs_MOK124$rowname 
logs_MOK124 <- logs_MOK124[,-1] 

# Make files for pathview - put rownames back, I decided to select them out of joined files 

colnames(logs_MOK023)
logs_23_24h_rows <- logs_MOK023 %>% select(log2FoldChange_24h) 
logs_23_48h_rows <- logs_MOK023 %>% select(log2FoldChange_48h)
logs_23_72h_rows <- logs_MOK023 %>% select(log2FoldChange_72h) 
logs_23_168h_rows <- logs_MOK023 %>% select(log2FoldChange_168h) 
head(logs_23_24h_rows) 

logs_124_24h_rows <- logs_MOK124 %>% select(log2FoldChange_24h) 
logs_124_48h_rows <- logs_MOK124 %>% select(log2FoldChange_48h)  
logs_124_72h_rows <- logs_MOK124 %>% select(log2FoldChange_72h) 
logs_124_168h_rows <- logs_MOK124 %>% select(log2FoldChange_168h) 
head(logs_124_168h_rows)


# Convert logs files to human for heatmaps ##### 

temp <- rownames_to_column(logs_MOK023) 
head(temp) 
logs_MOK023_human <- inner_join(species_1to1, temp, by = "rowname")
head(logs_MOK023_human) 
nrow(logs_MOK023_human) # without duplicates it is 6724 
typeof(logs_MOK023_human$Human.gene.stable.ID) # "integer" 
filter(logs_MOK023_human, Human.gene.stable.ID =="ENSG00000065054") 
# why is there duplicate rows for some genes???? and exactly the same? 
View(logs_MOK023_human)
logs_MOK023_human <- select(logs_MOK023_human, Human.gene.stable.ID, starts_with("log_")) 
# Remove duplicate rows 
duplicated(logs_MOK023_human) 
logs_MOK023_human <- logs_MOK023_human[!duplicated(logs_MOK023_human), ]
# Create rownames 
rownames(logs_MOK023_human) <- logs_MOK023_human$Human.gene.stable.ID 
logs_MOK023_human <- logs_MOK023_human[,-1] 

temp2 <- rownames_to_column(logs_MOK124) 
nrow(temp2) # 8889
head(temp2) 
logs_MOK124_human <- inner_join(species_1to1, temp2, by = "rowname")
head(logs_MOK124_human) 
nrow(logs_MOK124_human) # 8609, duplicated, later 7710 
typeof(logs_MOK124_human$Human.gene.stable.ID) # "integer" 
# why is there duplicate rows for some genes???? and exactly the same? 
View(logs_MOK124_human)
logs_MOK124_human <- select(logs_MOK124_human, Human.gene.stable.ID, starts_with("log_")) 
duplicated(logs_MOK124_human) 
rownames(logs_MOK124_human) <- logs_MOK124_human$Human.gene.stable.ID 
logs_MOK124_human <- logs_MOK124_human[!duplicated(logs_MOK124_human$rowname),]
# Create rownames 
rownames(logs_MOK124_human) <- logs_MOK124_human$Human.gene.stable.ID 
logs_MOK124_human <- logs_MOK124_human[,-1] 

colnames(logs_MOK023_human)
logs_23_24h_rows_human <- logs_MOK023_human %>% select(log2FoldChange_24h) 
logs_23_48h_rows_human <- logs_MOK023_human %>% select(log2FoldChange_48h)
logs_23_72h_rows_human <- logs_MOK023_human %>% select(log2FoldChange_72h) 
logs_23_168h_rows_human <- logs_MOK023_human %>% select(log2FoldChange_168h) 
head(logs_23_24h_rows_human) 

logs_124_24h_rows_human <- logs_MOK124_human %>% select(log2FoldChange_24h) 
logs_124_48h_rows_human <- logs_MOK124_human %>% select(log2FoldChange_48h)  
logs_124_72h_rows_human <- logs_MOK124_human %>% select(log2FoldChange_72h) 
logs_124_168h_rows_human <- logs_MOK124_human %>% select(log2FoldChange_168h) 
head(logs_124_168h_rows_human)
