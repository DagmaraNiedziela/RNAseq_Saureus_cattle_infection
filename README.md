# RNAseq_Saureus_cattle_infection 
All code associated with my RNAseq analysis of two groups of cattle infected with different strains of Staphylococcus aureus. There is 5-6 cows in each group and transcriptomics data was generated from milk somatic cell RNA at 0, 24, 48, 72 and 168 hours post infection. The library preparation (polyA selected mRNA stranded libraries) and sequencing (paired end, sequencing depth 60M reads) was done by Macrogen. I used FastQC for transcriptome control, and STAR for alignment to Bos taurus genome (UMD3.1). 

This code was used as a part of a PhD project in University College Dublin and Teagasc, entitled "Strain-specific virulence of Staphylococcus aureus". 
This RNA-Seq analysis was performed with advice from Dr Paul Cormican in Teagasc. All code was created by me apart from the perl file which was used to align fastq files to the bovine genome using STAR. 
