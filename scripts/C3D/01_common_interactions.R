#============================================================================================================================================================
# Objective: Analyse C3D results
#============================================================================================================================================================

    # 1. What is the distribution of the co-accessible regions ?
        # A. E-E , E-P and P-P

    # 2. Co-accessibility regions of LSC17 genes

    # 3. Do LSC, do PFA share co-accessibility hubs ?
        # A. How many of them map to the same gene -- shared chromosomal structure
        # B. How many map to different genes -- different model of regulation
        # C. Identify no. of E-P interactions mapping to common genes and different genes

#============================================================================================================================================================

## Identify no. of E-P interactions mapping to common genes and different genes
#
# for f in results/PCSC1/C3D/finalresults/*combined.results.txt ;
# do
#   name=$(echo $f | awk -F"/" '{ print $NF}'  | sed -e 's/.combined.results.txt//g')
#   awk 'BEGIN  {FS=OFS="\t"} { if($NF <0.01 && $7 >= 0.7 && $1!="COORD_1") print $1,$4,$7,$8,$9,$10}' $f > results/PCSC1/C3D/finalresults/$name.sig.results.txt
# done

library(data.table)
library(ggplot2)

gbm <- fread("results/PCSC1/C3D/finalresults/GBM.sig.results.txt", header=F, sep="\t", stringsAsFactors=F, data.table=F)
lsc <- fread("results/PCSC1/C3D/finalresults/LSCp.sig.results.txt", header=F, sep="\t", stringsAsFactors=F, data.table=F)
pfa <- fread("results/PCSC1/C3D/finalresults/PFA.sig.results.txt", header=F, sep="\t", stringsAsFactors=F, data.table=F)
