#!/bin/bash
# Convert bams to fastqs

for f in /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/RNAseq/LSCp/*.bam ;
do 
    name=$(echo $f | awk -F"/" '{ print $NF}' | sed -e 's/_Aligned.sortedByCoord.out.bam//g' )
    cmd="module load picard/2.6.0 ; java -jar -Xmx4g $picard_dir/picard.jar SamToFastq I=$f FASTQ=\"$name\"_R1.fastq.gz  SECOND_END_FASTQ=\"$name\"_R2.fastq.gz VALIDATION_STRINGENCY=LENIENT" ;
    echo $cmd > scripts/rnaseq/"$name".bam2fastq.sh
    sed -i 's/"//g' scripts/rnaseq/"$name".bam2fastq.sh
    qsub -q all.q -N LSC.bam2fastq.$name -e stdout/rnaseq/ -o stdout/rnaseq/ -cwd scripts/rnaseq/"$name".bam2fastq.sh
done

