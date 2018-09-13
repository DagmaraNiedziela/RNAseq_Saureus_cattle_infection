#make a folder with summaries, unzip all files and put text summaries in one folder with a name corresponding to sample name 
#one sample example 
mkdir fastqc_summary 
unzip 1776L-11-7-AM_1_fastqc.zip 
cp 1776L-11-7-AM_1_fastqc/summary.txt fastqc_summary/1776L-11-7-AM_1.txt #for one sample

#WORKING SCRIPT 
#Location of code: dniedziela/RNAseq_fastq/Run1/fastqc 
mkdir fastqc_summary 

for file in `ls *.zip`; \
do outfile=`basename $file | perl -p -e 's/\.zip//'`; \
echo -e "unzip $file; \
cp $outfile/summary.txt ./fastqc_summary/${outfile}_summary.txt" >> summary.sh; \
done 

sh summary.sh #run shell script 
#after I do it I can use cat to combine all summaries into one file 
cd fastqc_summary 
ls 
echo * #to make sure it will list all files 
cat * > merged_Fastqc_summaries_Run1.txt 
#end of working script 
#the same done for Run2 (merged file name ends with Run2) 

##########NOTES
#try for all samples? 
for i in `ls *.zip`;do echo -e "unzip ${i};cp ${i}/summary.txt fastqc_summary/${i}.txt" > ${i}.summary.sh;done 
#the summary file is: 
unzip 1922L-13-7_2_fastqc.zip.zip;cp 1922L-13-7_2_fastqc.zip/summary.txt fastqc_summary/1922L-13-7_2_fastqc.zip.txt 

#after the unzip had wrong directories I fixed without the unzip bit it still works 
for file in `ls *.zip`; \
do outfile=`basename $file | perl -p -e 's/\.zip//'`; \
echo -e "cp $outfile/summary.txt ./fastqc_summary/${outfile}_summary.txt" >> summary2.sh; \
done 

#unzip all first - this is wrong code 
unzip *.zip 
for i in `*fastqc`;do echo -e "cp ${i}fastqc/summary.txt fastqc_summary/${i}.txt" > ${i}.summary.sh;done 
for i in `ls *_fastqc/summary.txt`;do echo -e "cp ${i} fastqc_summary/${i}" > ../${i}.summary.sh;done #this does it for all files inside that folder not just the summary file 

sh *.summary.sh #run all shell scripts 
#after I do it I can use cat to combine all summaries into one file 
cd fastqc_summary 
ls 
echo * #to make sure it will list all files 
cat * > merged_Fastqc_summaries_Run1.txt 
paste * > columned_fastqc_summary.txt #will merge files column by column 