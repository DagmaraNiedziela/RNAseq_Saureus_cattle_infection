#!/usr/bin/perl



#open(my $fh1, ">>load_genome.sh") or die "Could not open file 'loadgenome' $!";

#print $fh1 "STAR --genomeLoad LoadAndExit --genomeDir /home/cormicp/people/genomes/bovine_ensembl3.1_feb2015\n";


open (READ, "samples") || die "Cannot open listdirs\n";
while (<READ>) {
    $contents = $_;
    chomp $contents;

    push (@list, $contents);


}


$count_job = '0';

foreach $sam (@list){
$r1=$sam;
$search= "R1";
$replace= "R2";

$sam =~ s/$search/$replace/;
$r2=$sam;
 

print "$r1\t";

@name2 = split /_/, $sam;
$sample = "s" . "_$name2[0]" ;

$runscript1 = "runscript_1"  . ".sh";
#$runscript2 = "runscript_2"  . ".sh";
#$runscript3 = "runscript_3"  . ".sh";
#$runscript4 = "runscript_4"  . ".sh";

print "$sample\t$runscript\n";

$staroutput = "$sample" . "_";

$samfile = "$sample" . "_Aligned.out.sam";

$bamfile = "$sample" . "_Aligned.out.bam";

$flagstat = "$bamfile" . ".flagstat";
 

	$count_job++;







open(my $fh, ">>$runscript1") or die "Could not open file '$runscript1' $!";
print $fh "module load bbmap/36.67\n";
print $fh "module load STAR/2.5.2\n";
print $fh "module load samtools/1.3.1\n";
print $fh "module load fastqc/11.5\n";


print $fh "STAR --readFilesCommand zcat --readFilesIn $r1 $r2 --runMode alignReads --genomeDir /data/tgsc1/shared/genomes/igenomes/bovine/ensembl/Bos_taurus/Ensembl/UMD3.1/Sequence/WholeGenomeFasta/ --runThreadN 16 --outFileNamePrefix ./align/$staroutput --outSAMmode Full --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outFilterMismatchNoverLmax 0.02 --quantMode GeneCounts --sjdbGTFfile /data/tgsc1/shared/genomes/igenomes/bovine/ensembl/Bos_taurus/Ensembl/UMD3.1/Annotation/Genes/genes.gtf --limitBAMsortRAM 20000000000.0 --outSAMtype BAM SortedByCoordinate\n";




#print $fh "samtools view -b -o ./align/$bamfile ./align/$samfile\n";
#print $fh "rm ./align/$samfile\n";
#print $fh "samtools flagstat  ./align/$bamfile > ./align/$flagstat\n";
#print $fh "gzip ./fastq/$clean1\n";
#print $fh "fastqc -t 2 ./fastq/$clean1 ./fastq/$r1\n";

print $fh "module unload bbmap/36.67\n";
print $fh "module unload STAR/2.5.2\n";
print $fh "module unload samtools/1.3.1\n";
print $fh "module unload fastqc/11.5\n";

close($fh);







#system("bsub -o $runscript.o -e $runscript.e < $runscript");

#sleep (600);
    }

