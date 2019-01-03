#!/bin/bash
#$ -cwd
#$ -N C3D
#$ -o logs/c3d.log
#$ -S /bin/bash

module load bedtools/2.23.0
module load R/3.4.1

# Cross Cell-type Correlation in DNaseI hypersensitivity (C3D)
# Written by: Tahmid Mehdi
# Princess Margaret Cancer Centre - University Health Network, July 25, 2016

# Cross Cell-type Correlation in DNaseI hypersensitivity (C3D)
# Takes correlations between open regions of chromatin based on DNaseI hypersensitivity signals
# Regions with high correlations are candidates for 3D interactions
# Performs association tests on each candidate & adjusts p-values
# Identifies transcription factor motifs which overlap open regions
# Produces interaction landscapes and motif tracks in PDF format


## Changing the script to accomodate different input variables

reference=$1
db=$2
anchor=$3
outDirectory=$4
matrix=$5

# # make output directory
# mkdir -p $4
# timestamp function
timestamp() {
  date +"%Y-%m-%d_%H-%M-%S"
}

if the anchor file does not have 5 fields, create a new formatted one
anchorCols=$(awk '{print NF}' $3 | sort -nu | tail -n 1)

if [ "$anchorCols" -ne "5" ]; then
    awk '{print $1"\t"$2"\t"$3"\t.\t"$1":"$2"-"$3}' $3 > $4/anchors.bed
else
	cat $3 > $4/anchors.bed
fi

#
# # if matrix is missing, map bedgraphs to reference
# if [ "$5" = "" ]; then
#         # Map all background files to one reference sample/peak catalogue
#         echo "$(timestamp): Mapping peak files"
#         cut -f 1-3 $1 > $4/ref.bed
#
#         counter=1
#         cat $2 | while read i; do
#                 echo "Mapping ${i} to file $counter.map.bed"
#                 mapBed -a $4/ref.bed -b ${i} -c 4 -o max -null 0 | awk 'BEGIN{ OFS="\t" }{ print $1, $2, $3, $4 }' > $4/$counter.map.bed
#                 counter=$((counter + 1))
#         done
# else
# 	tail -n +2 $5 | awk '{print $1}' | awk -F'[:-]' '{print $1"\t"$2"\t"$3}' > $4/ref1.bed
# 	sort -k1,1 -k2,2n $4/ref1.bed > $4/ref.bed
# 	rm $4/ref1.bed
# fi

# echo "----------------------"
# echo "Starting R script"
# echo "----------------------"
# echo "                      "
# echo "                      "
#
# if [ "$5" = "" ]; then
#     matrixfilecheck="n"
# fi

#Rscript /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/scripts/C3D/c3d.R $4 $4 $4/anchors.bed $2 $(timestamp) "$matrixfilecheck";
