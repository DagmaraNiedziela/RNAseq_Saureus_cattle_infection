for i in $(ls *_1.*);do  sn=`echo $i | cut -d '_' -f 1`; r2=`echo $i | sed 's/_1./_2./'`; echo  -e "module load fastp\nfastp --in1 $i --in2 $r2 --out1=${sn}_R1.clean.fastq.gz --out2=${sn}_R2.clean.fastq.gz -j ${sn}.json -h ${sn}.html -p -c\nmodule unload fastp" > ${sn}.fastp.sh && bsub -e ${sn}.fastp.sh.e -o ${sn}.fastp.sh.o < ${sn}.fastp.sh;done 
#this is done in folder Run1 and Run2, multiple sh scripts are generated and then ran by bsub 
sn=`echo $i | cut -d '_' -f 1`; r2=`echo $i | sed 's/_1./_2./'` #what is this? define variables for input files, sn got rid of underscore ("cut") in my file name to tidy it in output, and r2 defined read 2 files by "sed" command replacing _1 with _2 

#one sample 
module load fastp
fastp --in1 616L-13-7_1.fastq.gz --in2 616L-13-7_2.fastq.gz --out1=616L-13-7_R1.clean.fastq.gz --out2=616L-13-7_R2.clean.fastq.gz -j 616L-13-7.json -h 616L-13-7.html -p -c
module unload fastp 

#check options -j -h -p -c 
# -c base correction for PE data enabled 
# -p enable overrepresented sequence analysis 
#important -g - can specify for polyG correction - it is a default for NovaSeq sequencing, because no signal in those systems means G 
#no default correction options were disabled 
# -j and -h specify reporting 
# -j, --json the json format report file name (string [=fastp.json])
# -h, --html the html format report file name (string [=fastp.html]) 

######### COMBINE .JSON FILES INTO TABLE ###### 
#time to grep json files! 
#also, have a look at why I have no graphs 

grep "total_reads" 1776L-10-7-AM.json >> json_1776.txt;
grep "q20_rate" 1776L-10-7-AM.json >> json_1776.txt ;
grep "q30_rate" 1776L-10-7-AM.json >> json_1776.txt;
grep "gc_content" 1776L-10-7-AM.json >> json_1776.txt ;

#do these in Run2 and then in Run1 folders 
#Run1 has json files in a fastp_results directory - and some out of it, account for that 
for i in $(ls *.json);do echo  -e "echo ${i} >> ${i}_reads.txt; grep "total_reads" $i >> ${i}_reads.txt;grep "q20_rate" $i >> ${i}_reads.txt;grep "q30_rate" $i >> ${i}_reads.txt;grep "gc_content" $i  >> ${i}_reads.txt" >> json_summary.sh;done #it was giving me only the last one because the for loop was replacing and not appending script 
mkdir json_summary ;
mv *json_reads.txt ./json_summary ;
cd json_summary 
ls 
echo * #to make sure it will list all files 
paste * > columned_json_summary.txt #will merge files column by column 
#optionally 
cat * > merged_json_summaries_Run1.txt 

