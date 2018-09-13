#get genome 
# Download Bos taurus reference genome, version UMD3.1.1 from NCBI:
# nohup = "no hang up". A great way to send a task to the background and to free up terminal shell. Always include "&" after code.
nohup wget -o logfile -r -nd \
"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/055/GCF_000003055.6_Bos_taurus_UMD_3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna.gz" \
&
gunzip GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna.gz #cant seem to unzip that 

# Download annotation file for UMD3.1.1 NCBI Bos taurus Annotation Release 105:
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000003055.6_Bos_taurus_UMD_3.1.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff.gz
gunzip GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff.gz 

wget ftp://ftp.ensembl.org/pub/current_gtf/bos_taurus/Bos_taurus.UMD3.1.92.gtf.gz #checked today, Ensembl 

#New assembly? 
get new genome from NCBI 

##### Step 1 #### generate indexes - we are in STAR directory 
#trimmed fastq sequences are used for STAR ----- this isn't done because Paul already has indexes generated 
module load STAR 
STAR --runThreadN 20 #some nodes have 24, some have 48 threads, Paul uses 16 as a count here; if you use too many STAR will just give you less 
--runMode genomeGenerate
--genomeDir ./bovine_indexes 
--genomeFastaFiles ./bovine_genome/
--sjdbGTFfile ./bovine_annotation/Bos_taurus.UMD3.1.92.gtf 
--sjdbGTFtagExonParentTranscript Parent 
--sjdbOverhang 99 
--outFileNamePrefix ./bovine_indexes/Btau-UMD3.1.1

#### Step2 ##### map to genome 
#use trimmed fastq sequences! 
STAR --option1-name option1-value(s)--option2-name option2-value(s) 

STAR --runThreadN 10 
#--genomeLoad LoadAndKeep #this is if I run multiple STAR jobs, the genome will be in the environment for other jobs to use it 
--genomeDir ./bovine_indexes #if I manage to download sth there and index it 
--readFilesCommand zcat #if files are compressed as .gz - gunzip -c also works 
--readFilesIn ../1776L-11-7-AM_1.fastq.gz ../1776L-11-7-AM_2.fastq.gz  #comma for single reads, space for paired end reads 
#we probably want to use cleaned sequences 
--outFileNamePrefix ./samplename_mapped
--outSAMtype BAM Unsorted SortedByCoordinate #Paul only has sorted 
--quantMode GeneCounts #counts reads per gene - this gives me the nice output tables 

#Paul's code for one sample - we are in Run1 or Run2 folder 
module load STAR/2.5.2
STAR --readFilesCommand zcat --readFilesIn 1776L-11-7-AM_R1.clean.fastq.gz 1776L-11-7-AM_R2.clean.fastq.gz --runMode alignReads --genomeDir /data/tgsc1/shared/genomes/igenomes/bovine/ensembl/Bos_taurus/Ensembl/UMD3.1/Sequence/WholeGenomeFasta/ --runThreadN 16 --outFileNamePrefix ./align/s_1776L-11-7-AM_ --outSAMmode Full --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outFilterMismatchNoverLmax 0.02 --quantMode GeneCounts --sjdbGTFfile /data/tgsc1/shared/genomes/igenomes/bovine/ensembl/Bos_taurus/Ensembl/UMD3.1/Annotation/Genes/genes.gtf --limitBAMsortRAM 20000000000.0 --outSAMtype BAM SortedByCoordinate
#Paul's code ends 

####### EXPLANATIONS 
--runMode alignReads #this is default for STAR 
#look up all these 
--outSAMmode Full #Full is default 
--outSAMstrandField intronMotif #For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute, which STAR will generate with this ---- but my samples are stranded! 
--outReadsUnmapped Fastx #output of unmapped and partially mapped (i.e. mapped only one mate of a paired end read) reads in separate fastx files)
--outFilterMismatchNoverLmax 0.02 # default: 0.3, real: alignment will be output only if its ratio of mismatches to *mapped* length is less than or equal to this value --- what does it mean if its lower? 
--sjdbGTFfile /data/tgsc1/shared/genomes/igenomes/bovine/ensembl/Bos_taurus/Ensembl/UMD3.1/Annotation/Genes/genes.gtf #path to the GTF File with annotations
--limitBAMsortRAM 20000000000.0 #int>=0: maximum available RAM (bytes) for sorting BAM. If =0, it will be set to the genome index size. 0 value can only be used with {genomeLoad) NoSharedMemory option --- I totes dont get this 

#### GET MY DATA OUT 
grep "Uniquely mapped reads %" *_Log.final.out >> uniquely_mapped_Run1.txt #this gives me a file with all sample names and uniquely mapped read % for all! 
grep "Uniquely mapped reads number" *_Log.final.out >> uniquely_mapped_No_Run2.txt

#get rid of the start (Unmapped reads etc) from each Reads per Gene file - only reads left now. 
for i in `ls *ReadsPerGene.out.tab`;do grep "ENSBTAG" $i > ${i}.clean;done

mkdir STAR_ReadsPerGene #in RNAseq_Fastq
cp ../Run1/align/*10-7-AM_ReadsPerGene.out.tab.clean ./0_HPI
mkdir 0_HPI # and others, and I will only copy those files there 

#put in a header file for all files in them folder 
for i in `ls *clean`;do echo Processing ${i}; filename=${i%_ReadsPerGene.out.tab.clean}; echo -e "Genename\tall\tstrand1\t$filename" > tmp;cat tmp $i > ${i}.2;mv ${i}.2 $i;  done

####### THIS!!!!!! in each HPI folder (in STAR_ReadsPerGene) 
#put in a header file for all files in them folder 
for i in `ls *clean`;do echo Processing ${i}; filename=${i%_ReadsPerGene.out.tab.clean}; echo -e "Genename\tall\tstrand1\t$filename" > tmp;cat tmp $i > ${i}.2;mv ${i}.2 $i;  done ;
cut -f 1 s_568L-17-7_ReadsPerGene.out.tab.clean > Genename_ReadsPerGene.out.tab.clean.cut ;
for i in `ls *clean`;do cut -f 4 ${i} > ${i}.cut ; done ; 
paste *_ReadsPerGene.out.tab.clean.cut >> 168_HPI_GeneCounts.txt ; ####end of code - change date in 568 file and 168_HPI into folder name for each folder 

#### some notes for making GeneCounts tables, after the loop for columns 
ls #check if Genename is first as a file 
paste *_ReadsPerGene.out.tab.clean.cut >> 24_HPI_GeneCounts.txt

#this was not used 
paste echo *.ReadsPerGene.out.tab.clean >> STAR_ReadsPerGene_sample names.txt 
paste #Paul's code in a for loop 
cat STAR_ReadsPerGene_sample names.txt 0_HPI_sample_table.txt #not used 

#get my files in rows --- to do this first copy my reads per gene files from both runs into 1 folder! --- me trying 
paste <(cut -f 1,4 s_568L-17-7_ReadsPerGene.out.tab.clean) <(cut -f 1,4 s_616L-17-7_ReadsPerGene.out.tab.clean) | head #Paul's sample code 
paste <(cut -f 1,4 s_1776L-10-7-AM_ReadsPerGene.out.tab.clean) > 0_HPI_sample_table.txt 
paste 0_HPI_sample_table.txt *<(cut -f 4 *10-7-AM_ReadsPerGene.out.tab.clean) >> 0_HPI_GeneCounts.txt #####this works!!!! I just cant see past the 3 columns, and 1776 is in twice 
for i in `ls *clean`;do paste 0_HPI_sample_table.txt <(cut -f 4 ${i}) >> 0_HPI_sample_tableX.txt; done ###none of these work 