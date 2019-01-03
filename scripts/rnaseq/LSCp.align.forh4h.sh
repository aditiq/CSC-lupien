#!/bin/bash


module load igenome-human/hg38 
module load  STAR/2.4.2a

## Generating genome indexes files -- Needs to be done only once for each genome/annotation combination. But is readlength dependant. 
## So we need to make sure that LSCp is 75bo in length -- They all are. So we ddont need to generate new genome index files

## Map reads to genome

STAR --runMode alignReads \
--runThreadN 30 \
--genomeDir /cluster/tools/data/genomes/human/hg38/STARIndex_hg38_genomeSAsparseD2 \
--readFilesIn $1 $2 \
--readFilesCommand zcat \
--outFileNamePrefix path/to/output/dir/prefix \
--outSAMtype BAM SortedByCoordinate \
--outSAMprimaryFlag AllBestScore \
--outSAMunmapped Within \
--outSAMstrandField intronMotif \
--chimSegmentMin 20 \
--twopassMode Basic \
--outFilterIntronMotifs RemoveNoncanonical \
--quantMode TranscriptomeSAM GeneCounts \
--genomeSAsparseD 2
