library(DESeq2)
library(ggplot2) 

# Plot single markers #####

# ** Plot cell markers MOK023 ####

# CK19, EPCAM - epithelial cell 
CK19Counts23 <- plotCounts(dds23, gene = "ENSBTAG00000004905", intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE)
CK19Counts23$Cow_ID <- as.factor(CK19Counts23$Cow_ID) 
plot <- ggplot(CK19Counts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() 
print(plot + ggtitle("CK19"))

EPCAMCounts23 <- plotCounts(dds23, gene = "ENSBTAG00000006474", intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE)
EPCAMCounts23$Cow_ID <- as.factor(EPCAMCounts23$Cow_ID) 
plot <- ggplot(EPCAMCounts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() 
print(plot + ggtitle("EPCAM")) 

# NF-kB subunit 1 - not for epithelial cells! 
NFKB1Counts23 <- plotCounts(dds23, gene = "ENSBTAG00000020270", intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE)
NFKB1Counts23$Cow_ID <- as.factor(NFKB1Counts23$Cow_ID) 
plot <- ggplot(NFKB1Counts23, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) + 
  scale_y_log10() + geom_point(size = 3) + geom_line() 
print(plot + ggtitle("NF-kB subunit 1")) 


# ** Plot cell markers MOK124 ####

# CK19, EPCAM - epithelial cell 
CK19Counts <- plotCounts(dds124, gene = "ENSBTAG00000004905", intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE)
CK19Counts$Cow_ID <- as.factor(CK19Counts$Cow_ID) 
plot <- ggplot(CK19Counts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() 
print(plot + ggtitle("CK19"))

EPCAMCounts <- plotCounts(dds124, gene = "ENSBTAG00000006474", intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE)
EPCAMCounts$Cow_ID <- as.factor(EPCAMCounts$Cow_ID) 
plot <- ggplot(EPCAMCounts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() 
print(plot + ggtitle("EPCAM"))

# NF-kB subunit 1
NFKB1Counts <- plotCounts(dds124, gene = "ENSBTAG00000020270", intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE)
NFKB1Counts$Cow_ID <- as.factor(NFKB1Counts$Cow_ID) 
plot <- ggplot(NFKB1Counts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) + 
  scale_y_log10() + geom_point(size = 3) + geom_line() 
print(plot + ggtitle("NF-kB subunit 1")) 

# After initial data exploration create files for other epithelial cell markers ###### 

# CK18 - epithelial cells ENSBTAG00000001517 
CK18Counts23 <- plotCounts(dds23, gene = "ENSBTAG00000001517", intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE)
CK18Counts23$Cow_ID <- as.factor(CK18Counts23$Cow_ID) 
plot <- ggplot(CK18Counts, aes(x = HPI, y = count, color = Cow_ID, group = Cow_ID)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() 
print(plot + ggtitle("CK18, MOK023")) 

CK18Counts <- plotCounts(dds124, gene = "ENSBTAG00000001517", intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE)
CK18 <- rbind(CK18Counts23, CK18Counts)

# Cytokeratin 14 (KRT14), ENSBTAG00000007583 
CK14Counts23 <- plotCounts(dds23, gene = "ENSBTAG00000007583", intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE)
CK14Counts <- plotCounts(dds124, gene = "ENSBTAG00000007583", intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE)
CK14 <- rbind(CK14Counts23, CK14Counts)

# E-cadherin (CDH1) ENSBTAG00000015991, PCDH1 ENSBTAG00000012518 
CDH1Counts23 <- plotCounts(dds23, gene = "ENSBTAG00000015991", intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE)
CDH1Counts <- plotCounts(dds124, gene = "ENSBTAG00000015991", intgroup = c("HPI", "Group", "Cow_ID"), returnData = TRUE)
CDH1 <- rbind(CDH1Counts23, CDH1Counts)

# Join the rest of gene files ##### 

CK19 <- rbind(CK19Counts23, CK19Counts) 
EPCAM <- rbind(EPCAMCounts23, EPCAMCounts)
NFKB1 <- rbind(NFKB1Counts, NFKB1Counts23) 

# Make plots with groups means ####### 

# First, define a data summary function (copied from R cookbook) - also this can be done using dplyr summarize function more easily but this was the code used here 

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
} 

library(RColorBrewer)

# Make means of my data - CK19 ####
tgcCK19 <- summarySE(CK19, measurevar="count", groupvars=c("Group","HPI"), na.rm = TRUE)
tgcCK19

# Plot with Standard error of the mean
ggplot(tgcCK19, aes(x=HPI, y=count, colour=Group)) + 
  geom_errorbar(aes(ymin=count-se, ymax=count+se), colour="black", width=0.1) +
  geom_line(aes(group=Group)) +
  geom_point() + scale_color_brewer(palette="Set1", direction = -1) + theme_bw() + 
  ggtitle("CK19 - epithelial cells")
## width - This is the width of the crossbars!
# in color_brewer - direction = -1 reverses the order of the colors :D 
# geom_line is dropping missing values 
ggsave("CK19.jpeg") 

tgc_EPCAM <- summarySE(EPCAM, measurevar="count", groupvars=c("Group","HPI"), na.rm = TRUE) 

ggplot(tgc_EPCAM, aes(x=HPI, y=count, colour=Group)) + 
  geom_errorbar(aes(ymin=count-se, ymax=count+se), colour="black", width=0.1) +
  geom_line(aes(group=Group)) +
  geom_point() + scale_color_brewer(palette="Set1", direction = -1) + theme_bw() + 
  ggtitle("EPCAM - epithelial cells") 
ggsave("EPCAM_group_means.jpeg") 

# CK18 
tgc_CK18 <- summarySE(CK18, measurevar="count", groupvars=c("Group","HPI"), na.rm = TRUE) 

ggplot(tgc_CK18, aes(x=HPI, y=count, colour=Group)) + 
  geom_errorbar(aes(ymin=count-se, ymax=count+se), colour="black", width=0.1) +
  geom_line(aes(group=Group)) +
  geom_point() + scale_color_brewer(palette="Set1", direction = -1) + theme_bw() + 
  ggtitle("CK18 - epithelial cells") 
ggsave("CK18_group_means.jpeg") 

# CK14 
tgc_CK14 <- summarySE(CK14, measurevar="count", groupvars=c("Group","HPI"), na.rm = TRUE) 

ggplot(tgc_CK14, aes(x=HPI, y=count, colour=Group)) + 
  geom_errorbar(aes(ymin=count-se, ymax=count+se), colour="black", width=0.1) +
  geom_line(aes(group=Group)) +
  geom_point() + scale_color_brewer(palette="Set1", direction = -1) + theme_bw() + 
  ggtitle("CK14 - epithelial cells") 
ggsave("CK14_group_means.jpeg") 

# CDH1 
tgc_CDH1 <- summarySE(CDH1, measurevar="count", groupvars=c("Group","HPI"), na.rm = TRUE) 

ggplot(tgc_CDH1, aes(x=HPI, y=count, colour=Group)) + 
  geom_errorbar(aes(ymin=count-se, ymax=count+se), colour="black", width=0.1) +
  geom_line(aes(group=Group)) +
  geom_point() + scale_color_brewer(palette="Set1", direction = -1) + theme_bw() + 
  ggtitle("CDH1 - epithelial cells") 
ggsave("CDH1_group_means.jpeg") 

# Make means of my data - NFKB1  
tgcNFKB1 <- summarySE(NFKB1, measurevar="count", groupvars=c("Group","HPI"), na.rm = TRUE) 

# Fix the hours post infection being a factor - as.character is crucial here! 
typeof(tgcNFKB1$HPI2) 
tgcNFKB1$HPI2 <- as.numeric(as.character(tgcNFKB1$HPI))

#** Plot NFkB for data exploration figure ####
NFkB_plot <- ggplot(tgcNFKB1, aes(x=HPI2, y=count, colour=Group)) + 
  geom_errorbar(aes(ymin=count-se, ymax=count+se), colour="black", width=2) +
  geom_line(aes(group=Group)) + scale_x_continuous(breaks = c(0,24,48,72,168)) +
  scale_y_continuous(labels = scales::comma) +
  geom_point() + scale_color_brewer(palette="Set1", direction = -1) + theme_gray() + 
  ggtitle("NFKB1") + 
  xlab("Hours post infection (hpi)") + ylab("Normalised gene count") +
  theme(text = element_text(size = 16, family = "Calibri"), plot.title = element_text(face = "italic")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
NFkB_plot
ggsave("NF-kB_group_mean_continuous_x2.jpeg") 

# # Join epithelial cell plots in a grid ##### 
install.packages("ggtext")
library(ggtext)

CK19plot2 <- ggplot(tgcCK19, aes(x=HPI, y=count, colour=Group)) + 
  geom_errorbar(aes(ymin=count-se, ymax=count+se), colour="black", width=0.1) +
  geom_line(aes(group=Group)) +
  geom_point() + scale_color_brewer(palette="Set1", direction = -1) + theme_bw() +
  ggtitle("KRT19") + labs(x="HPI",y="Normalised counts - *KRT19*") + 
    theme(plot.title = element_text(face = "italic"), 
          axis.title.y = ggtext::element_markdown()) 
CK19plot2
ggsave("CK19.jpeg") 

EPCAMplot2 <- ggplot(tgc_EPCAM, aes(x=HPI, y=count, colour=Group)) + 
  geom_errorbar(aes(ymin=count-se, ymax=count+se), colour="black", width=0.1) +
  geom_line(aes(group=Group)) +
  geom_point() + scale_color_brewer(palette="Set1", direction = -1) + theme_bw() +
  ggtitle("EPCAM") + labs(x="HPI",y="Normalised counts - *EPCAM*") + 
    theme(plot.title = element_text(face = "italic"), 
          axis.title.y = ggtext::element_markdown()) 

EPCAMplot2
ggsave("EPCAM_group_means.jpeg") 

# CK18 
CK18plot2 <- ggplot(tgc_CK18, aes(x=HPI, y=count, colour=Group)) + 
  geom_errorbar(aes(ymin=count-se, ymax=count+se), colour="black", width=0.1) +
  geom_line(aes(group=Group)) +
  geom_point() + scale_color_brewer(palette="Set1", direction = -1) + theme_bw() +
 ggtitle("KRT18") + labs(x="HPI",y="Normalised counts - *KRT18*") + 
  theme(plot.title = element_text(face = "italic"), 
        axis.title.y = ggtext::element_markdown()) 
CK18plot2
ggsave("CK18_group_means.jpeg") 

# CK14 
CK14plot2 <- ggplot(tgc_CK14, aes(x=HPI, y=count, colour=Group)) + 
  geom_errorbar(aes(ymin=count-se, ymax=count+se), colour="black", width=0.1) +
  geom_line(aes(group=Group)) +
  geom_point() + scale_color_brewer(palette="Set1", direction = -1) + theme_bw() +
 ggtitle("KRT14") + labs(x="HPI",y="Normalised counts - *KRT14*") + 
  theme(plot.title = element_text(face = "italic"), 
        axis.title.y = ggtext::element_markdown()) 
CK14plot2
ggsave("CK14_group_means.jpeg") 

# CDH1 
CDH1plot2 <- ggplot(tgc_CDH1, aes(x=HPI, y=count, colour=Group)) + 
  geom_errorbar(aes(ymin=count-se, ymax=count+se), colour="black", width=0.1) +
  geom_line(aes(group=Group)) +
  geom_point() + scale_color_brewer(palette="Set1", direction = -1) + theme_bw() +
  ggtitle("CDH1") + labs(x="HPI",y="Normalised counts - *CDH1*") + 
  theme(plot.title = element_text(face = "italic"), 
        axis.title.y = ggtext::element_markdown()) 
CDH1plot2
ggsave("CDH1_group_means.jpeg") 

plot2 <- ggarrange(CK14plot2, CK18plot2, CK19plot2, CDH1plot2, EPCAMplot2, 
                  #labels = c("A", "B", "C", "D", "E"),
                  ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom") 
plot2 
ggsave("Epithelial_cell_markers_notitles.jpeg")
# These are combined with cell composition in code file 8
