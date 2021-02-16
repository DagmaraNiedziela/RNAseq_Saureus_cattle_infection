# This file contains notes - please don't run everything! 
# The for loops are the only parts that need to be ran 

pwd 
ls 
mkdir RNAseq_fastq 
cd RNAseq_fastq 
mkdir Run1 
cd Run1 
#files copied by Filezilla 

#Check md5sums to see if files were downloaded properly 
md5sum -c md5sum.txt > md5sum_check.txt #should return file name and OK 

#checked one file just in case 
echo d3a8a9ec4d5cf9cad0277516d4b6c5a1  1776L-11-7-AM_1.fastq.gz | md5sum -c
#answer: 1776L-11-7-AM_1.fastq.gz: OK 

#FAST QC
lsload #check node usage 
ssh compute01 #get onto node 
module avail #check available software 
module load fastqc #load software 
module unload fastqc #when I finish with software 

#MAIN CODE - FOR LOOP 
for i in `ls *.gz`;do echo -e "module load fastqc\nfastqc $i" > ${i}.fastqc.sh;bsub -e ${i}.fastqc.sh.e -o ${i}.fastqc.sh.o < ${i}.fastqc.sh;done 

#this is the loop that does fastqc in both my folders (I cd into each folder to do it) 
it will write (echo)
module load fastqc 
fastqc "sample" 
into a script ${i}.fastqc.sh 
#then you queue it onto server in separate jobs (bsub), files .e and .o are reports (.o - output), and it runs all scripts onto bsub 

#Run fastQC on trimmed sequences (in Run1 and Run2)
for i in `ls *.clean.fastq.gz`;do echo -e "module load fastqc\nfastqc $i" > ${i}.fastqc.sh;bsub -e ${i}.fastqc.sh.e -o ${i}.fastqc.sh.o < ${i}.fastqc.sh;done

bsub #is submitting jobs to a queue 
bqueues #check how many jobs queued 
bjobs #check how many of my jobs are running or pending 

nohup fastqc -o fastqc_cleaned *.clean.fastq.gz &
