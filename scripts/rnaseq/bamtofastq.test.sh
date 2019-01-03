#!/bin/bash

#-----------------------------------------------------------------
## Converting bam to fastq and checking fastqc
## Running this on GBM file from SU2C so that I can compare the fastq
#-----------------------------------------------------------------

bam="/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/RNAseq/GBM.PFA.HF/A61519_Dirks_G523_LineAligned.sortedByCoord.out.bam "

#-----------------------------
## Bedtools 
#-----------------------------

cmd1=" module load bedtools/2.23.0;
module load samtools/1.8;

samtools sort -n \"$bam\" -o /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/RNAseq/GBM.PFA.HF/A61519_Dirks_G523_LineAligned.sorted.bam;
bedtools bamtofastq -i /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/RNAseq/GBM.PFA.HF/A61519_Dirks_G523_LineAligned.sorted.bam -fq /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/RNAseq/GBM.PFA.HF/A61519_Dirks_G523_Line_R1.bedtools.fq -fq2 /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/RNAseq/GBM.PFA.HF/A61519_Dirks_G523_Line_R2.bedtools.fq ; 
"
echo $cmd1 > scripts/rnaseq/bedtools.G523.bam2fastq.sh
sed -i 's/"//g' scripts/rnaseq/bedtools.G523.bam2fastq.sh

qsub -q all.q -cwd -N bedtools.G523.bam2fastq -e stdout/rnaseq -o stdout/rnaseq scripts/rnaseq/bedtools.G523.bam2fastq.sh

#-----------------------------
## Picard
#-----------------------------

name="A61519_Dirks_G523" ;  

cmd2=" 
module load picard/2.6.0 ; 
java -jar -Xmx4g /mnt/work1/software/picard/2.6.0/picard.jar SamToFastq I=$bam FASTQ=/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/RNAseq/GBM.PFA.HF/\"$name\"_R1.picard.fastq.gz  SECOND_END_FASTQ=/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/RNAseq/GBM.PFA.HF/\"$name\"_R2.picard.fastq.gz VALIDATION_STRINGENCY=LENIENT ; 
java -jar -Xmx4g /mnt/work1/software/picard/2.6.0/picard.jar SamToFastq I=/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/RNAseq/GBM.PFA.HF/A61519_Dirks_G523_LineAligned.sorted.bam FASTQ=/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/RNAseq/GBM.PFA.HF/\"$name\"_R1.sorted.picard.fastq.gz SECOND_END_FASTQ=/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/RNAseq/GBM.PFA.HF/\"$name\"_R2.sorted.picard.fastq.gz VALIDATION_STRINGENCY=LENIENT ; 
"

echo $cmd2 > scripts/rnaseq/picard.G523.bam2fastq.sh
sed -i 's/"//g' scripts/rnaseq/picard.G523.bam2fastq.sh

qsub -q all.q -cwd -N picard.G523.bam2fastq -e stdout/rnaseq -o stdout/rnaseq scripts/rnaseq/picard.G523.bam2fastq.sh




