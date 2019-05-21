# 20 05 2019
# This file is not finished - read through and remove reversing log2 fold changes because they are now correct 
# Pick only the pathway images that I used and maybe some examples? 

# Make LOGS tables of my significant genes ##### 
# LOGS is a data frame that has log2 fold changes only and ENSEMBL Gene IDs as row names - this will be used for pathview 
# Use tidyverse 
library(DESeq2)
library(dplyr) 
library(tidyverse)
detach(package:pathview) # pathview and dplyr have conflicting select functions, so I sometimes have to detach one or the other to progress with the code 
detach(package:dplyr)

head(res23_24h_LF_sigDF)
nrow(res23_24h_LF_sig) 
colnames(res23_24h_LF_sig) 

logs_23_24h <- res23_24h_LF_sig %>% select(log2FoldChange) 
logs_23_48h <- res23_48h_LF_sig %>% select(log2FoldChange) 
logs_23_72h <- res23_72h_LF_sig %>% select(log2FoldChange) 
logs_23_168h <- res23_168h_LF_sig %>% select(log2FoldChange) 

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
head(res124_24h_LF_sig)

logs_124_24h <- res124_24h_LF_sig %>% select(log2FoldChange) 
logs_124_48h <- res124_48h_LF_sig %>% select(log2FoldChange) 
logs_124_72h <- res124_72h_LF_sig %>% select(log2FoldChange) 
logs_124_168h <- res124_168h_LF_sig %>% select(log2FoldChange) 

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
logs_23_24h_rows <- logs_MOK023 %>% select(log_24h) 
logs_23_48h_rows <- logs_MOK023 %>% select(log_48h)
logs_23_72h_rows <- logs_MOK023 %>% select(log_72h) 
logs_23_168h_rows <- logs_MOK023 %>% select(log_168h) 
head(logs_23_24h_rows) 

logs_124_24h_rows <- logs_MOK124 %>% select(log_24h) 
logs_124_48h_rows <- logs_MOK124 %>% select(log_48h)  
logs_124_72h_rows <- logs_MOK124 %>% select(log_72h) 
logs_124_168h_rows <- logs_MOK124 %>% select(log_168h) 
head(logs_124_168h_rows)


# Convert logs files to human for pathview ##### 

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
logs_23_24h_rows_human <- logs_MOK023_human %>% select(log_24h) 
logs_23_48h_rows_human <- logs_MOK023_human %>% select(log_48h)
logs_23_72h_rows_human <- logs_MOK023_human %>% select(log_72h) 
logs_23_168h_rows_human <- logs_MOK023_human %>% select(log_168h) 
head(logs_23_24h_rows_human) 

logs_124_24h_rows_human <- logs_MOK124_human %>% select(log_24h) 
logs_124_48h_rows_human <- logs_MOK124_human %>% select(log_48h)  
logs_124_72h_rows_human <- logs_MOK124_human %>% select(log_72h) 
logs_124_168h_rows_human <- logs_MOK124_human %>% select(log_168h) 
head(logs_124_168h_rows_human)

# Plot various combinations of bovine KEGG pathways with pathview ##### 
library("pathview") 
detach(package:dplyr) 
# Pathview saves a PNG within the project home directory automatically 

# ** Strain - all HPI pathways ####
# NF-Kb in MOK023 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK023, pathway.id = "04064", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK023_NF-kB_2", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red2", cpd = "red"), 
                   na.col = "transparent") 
# I changed NAs to grey, because they were transparent/white, and it seems that genes that the genes that 
# Then I changed it back to default and zeros to white - clearer image 
# are NA on the list come up as grey (between 1 and -1) 
sum(is.na(logs_MOK124$log_24h)) # 1690 
sum(logs_MOK124$log_24h < 1, na.rm=TRUE) # 3646 
nrow(logs_MOK124) #8889 
sum(logs_MOK124$log_24h > 1, na.rm=TRUE) # 3553 
3553 + 3646 + 1690 # 8889 

# NF-Kb in MOK124 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK124, pathway.id = "04064", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK124_NF-kB_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red2", cpd = "red"), 
                   na.col = "transparent") 

# Hippo in MOK124 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK124, pathway.id = "04390", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK124_Hippo_2", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# Hippo in MOK023 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK023, pathway.id = "04390", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK023_Hippo_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# TNF signalling in MOK023 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK023, pathway.id = "04668", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK023_TNF_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# TNF in MOK124 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK124, pathway.id = "04668", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK124_TNF_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# Chemokine signalling in MOK023 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK023, pathway.id = "04062", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK023_chemokine_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# Chemokine in MOK124 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK124, pathway.id = "04062", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK124_chemokine_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# Cytokine receptor interaction in MOK023 (significant at 48 hpi) - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK023, pathway.id = "04060", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK023_cytokine_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# Cytokine receptor interaction in MOK124 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK124, pathway.id = "04060", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK124_cytokine_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 


# Rap1 in MOK023 (NS) - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK023, pathway.id = "04015", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK023_Rap1_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# Rap1 in MOK124 - log2 fold change from 8 to -8, PNG image - 
#this is a pathway only significant in MOK124, 
# at 24, 48 and 72 hpi 
pv.out <- pathview(gene.data = logs_MOK124, pathway.id = "04015", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK124_Rap1_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# ECM in MOK023 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK023, pathway.id = "04512", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK023_ECM_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# ECM in MOK124 - log2 fold change from 8 to -8, PNG image - 
# this is a pathway only significant in MOK124, at 48 hpi
pv.out <- pathview(gene.data = logs_MOK124, pathway.id = "04512", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK124_ECM_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# Tight junction in MOK023 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK023, pathway.id = "04530", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK023_tight_junction_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# Tight junction in MOK124 - log2 fold change from 8 to -8, PNG image - 
# this is a pathway only significant in MOK124, 
# at 48, 72 and 168 hpi 
pv.out <- pathview(gene.data = logs_MOK124, pathway.id = "04530", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK124_tight_junction_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# NON SIGNIFICANT - TLR signalling? 

# TLR signalling in MOK023 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK023, pathway.id = "04620", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK023_TLR_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# TLR in MOK124 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK124, pathway.id = "04620", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK124_TLR_allHPI", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red"))

# ** Check what are min and max log2foldchanges #####

# MOK124 
max(logs_MOK124$log_24h, na.rm = TRUE) # 9.793361 
min(logs_MOK124$log_24h, na.rm = TRUE) # -5.356991
sum(logs_MOK124$log_24h > 8, na.rm=TRUE) # 4 

max(logs_MOK124$log_48h, na.rm = TRUE) #  8.90111 
min(logs_MOK124$log_48h, na.rm = TRUE) # -9.897618
sum(logs_MOK124$log_48h > 8, na.rm=TRUE) # 2 
sum(logs_MOK124$log_48h < -8, na.rm=TRUE) # 20 

max(logs_MOK124$log_72h, na.rm = TRUE) # 8.817051
min(logs_MOK124$log_72h, na.rm = TRUE) # -8.821492
sum(logs_MOK124$log_72h > 8, na.rm=TRUE) # 3
sum(logs_MOK124$log_72h < -8, na.rm=TRUE) # 7 

max(logs_MOK124$log_168h, na.rm = TRUE) # 11.42208
min(logs_MOK124$log_168h, na.rm = TRUE) # -9.22833
sum(logs_MOK124$log_168h > 8, na.rm=TRUE) # 3 
sum(logs_MOK124$log_168h < -8, na.rm=TRUE) # 9 

# MOK023 
max(logs_MOK023$log_24h, na.rm = TRUE) # 9.038831
min(logs_MOK023$log_24h, na.rm = TRUE) # -6.530827
sum(logs_MOK023$log_24h > 8, na.rm=TRUE) # 4 

max(logs_MOK023$log_48h, na.rm = TRUE) #  8.509834
min(logs_MOK023$log_48h, na.rm = TRUE) # -6.692438
sum(logs_MOK023$log_48h > 8, na.rm=TRUE) # 3
sum(logs_MOK023$log_48h < -8, na.rm=TRUE) # 0 

max(logs_MOK023$log_72h, na.rm = TRUE) # 7.448415
min(logs_MOK023$log_72h, na.rm = TRUE) # -7.32069
sum(logs_MOK023$log_72h > 8, na.rm=TRUE) # 0
sum(logs_MOK023$log_72h < -8, na.rm=TRUE) # 0

max(logs_MOK023$log_168h, na.rm = TRUE) # 7.634158
min(logs_MOK023$log_168h, na.rm = TRUE) # -8.984974
sum(logs_MOK023$log_168h > 8, na.rm=TRUE) # 0 
sum(logs_MOK023$log_168h < -8, na.rm=TRUE) # 1  
# Way less of them in MOK023, but not too many over 8 or < -8, potentially it is ok to leave it as is in pathview. 



# ** Single time point - MOK023, get a separate vector with rownames  ####
pv.out <- pathview(gene.data = logs_23_24h_rows, pathway.id = "04064", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK023_24hpi_NF-kB", kegg.native = F, same.layer = F, 
                   sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 
# multi.state - T not necessary here, it only has one value - so this could be false

# Output a png image - pathway as above (kegg.native = T)
pv.out <- pathview(gene.data = logs_23_24h_rows, pathway.id = "04064", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK023_24hpi_NF-kB", kegg.native = T, same.layer = F, 
                   sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 
#  "Note: 684 of 7841 unique input IDs unmapped."

# Hippo in MOK124, 168 hpi - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_124_168h_rows, pathway.id = "04390", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "MOK124_168hpi_Hippo_2", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"),
                   high = list(gene = "red2", cpd = "red"), 
                   na.col = "transparent") 
# [1] "Note: 748 of 8889 unique input IDs unmapped." 
# Nicer to have all the NAs as white. 

# ** Both strains - HPI - Pathways at specific time point ####
# NF-kB at 24 hpi 
pv.out <- pathview(gene.data = logs_24hpi, pathway.id = "04064", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "24hpi_NF-kB", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red2", cpd = "red"), 
                   na.col = "transparent") 

# Hippo at 168 hpi 
pv.out <- pathview(gene.data = logs_168hpi, pathway.id = "04390", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "168hpi_Hippo", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"),
                   high = list(gene = "red2", cpd = "red"), 
                   na.col = "transparent") 

# ** NOTES for pathview: #####

# Error in UseMethod("select_") --> because plyr select is masking it! 
detach(package:dplyr) 
# Supposedly works better if select_ does not apply to objects of class... 
detach(package:dplyr, unload = TRUE) 

# Important changes: 
# change color scheme? - done 
# Gene range not from -1 to 1? - done (limit is log2foldchange and bins are divisions of color range)

# LAST: create separate vectors for use with time points (with row names) - done 
# These are taken from the compiled file, so they have NA rows - pathview treats them as zeros 
# so I changed all zeros and NAs to white 

# Pathview on human pathways ##### 

# NF-Kb in MOK023 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK023_human, pathway.id = "04064", gene.idtype = "ENSEMBL", 
                   species = "hsa", out.suffix = "MOK023_NF-kB_2", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red2", cpd = "red"), 
                   na.col = "transparent") 

# Hippo in MOK124 - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_MOK124_human, pathway.id = "04390", gene.idtype = "ENSEMBL", 
                   species = "hsa", out.suffix = "MOK124_Hippo_2", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red"))

# Single time point - MOK023, NF-kB 
pv.out <- pathview(gene.data = logs_23_24h_rows_human, pathway.id = "04064", gene.idtype = "ENSEMBL", 
                   species = "hsa", out.suffix = "MOK023_24hpi_NF-kB", kegg.native = T, same.layer = F, 
                   sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red3", cpd = "red")) 

# Hippo in MOK124, 168 hpi - log2 fold change from 8 to -8, PNG image  
pv.out <- pathview(gene.data = logs_124_168h_rows_human, pathway.id = "04390", gene.idtype = "ENSEMBL", 
                   species = "hsa", out.suffix = "MOK124_168hpi_Hippo_2", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"),
                   high = list(gene = "red2", cpd = "red"), 
                   na.col = "transparent") 

# For these I have no data frames for the moment ! 
# Pathways at specific time point 
# NF-kB at 24 hpi 
pv.out <- pathview(gene.data = logs_24hpi_human, pathway.id = "04064", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "24hpi_NF-kB", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"), 
                   high = list(gene = "red2", cpd = "red"), 
                   na.col = "transparent") 

# Hippo at 168 hpi 
pv.out <- pathview(gene.data = logs_168hpi_human, pathway.id = "04390", gene.idtype = "ENSEMBL", 
                   species = "bta", out.suffix = "168hpi_Hippo", keys.align = "y", 
                   kegg.native = T, match.data = F, multi.state = T, same.layer = F,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i], 
                   limit = list(gene = 8, cpd = 8), bins = list(gene = 16, cpd = 16), 
                   low = list(gene = "blue", cpd = "blue"), mid = list(gene = "white", cpd = "white"),
                   high = list(gene = "red2", cpd = "red"), 
                   na.col = "transparent") 