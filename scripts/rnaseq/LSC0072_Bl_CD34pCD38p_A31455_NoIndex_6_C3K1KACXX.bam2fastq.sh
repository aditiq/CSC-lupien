module load picard/2.6.0 ; java -jar -Xmx4g /mnt/work1/software/picard/2.6.0/picard.jar SamToFastq I=/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/RNAseq/LSCp/LSC0072_Bl_CD34pCD38p_A31455_NoIndex_6_C3K1KACXX_Aligned.sortedByCoord.out.bam FASTQ=LSC0072_Bl_CD34pCD38p_A31455_NoIndex_6_C3K1KACXX_R1.fastq.gz SECOND_END_FASTQ=LSC0072_Bl_CD34pCD38p_A31455_NoIndex_6_C3K1KACXX_R2.fastq.gz VALIDATION_STRINGENCY=LENIENT