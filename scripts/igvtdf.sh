#!/bin/bash

## Generate tdf files for iGV browser for the FE tracks

files="100091_+_+_TAAGGCGA_L008_R1
8001_TAAGGCGA_L001_R1
8011_TAAGGC_L005_R1
8021_TAAGGCGA_L004_R1
8031_TAAGGCGA_L004_R1
8032_CGTACTAG_L004_R1
8041_TCCTGAGC_L002_R1
8051_TCCTGAGC_L006_R1
8056_TAAGGCGA_L008_R1
8057_CGTACTAG_L008_R1
8076_TCCTGAGC_L006_R1
8077_TAAGGCGA_L007_R1
8081_TAAGGCGA_L008_R1
8082_CGTACTAG_L003_R1
8083_TCCTGAGC_L008_R1
8096_TCCTGAGC_L004_R1
8101_TAAGGCGA_L007_R1
90239_-_-_TAAGGCGA_L007_R1
G361
G411
G489
G498
G523
G549
G551
G566
G567
G571
G577
G583
G594
G613
G620
G648
G702
G705
G706
G719
G729
G744
G797
PFA2
PFA3
PFA4
PFA5
PFA6
PFA7
PFA8"

for f in $files ;
do

  cmd1="module load igvtools/2.3.32 ; igvtools toTDF data/peaks/\"$f\"/\"$f\"_FE.bedgraph data/peaks/\"$f\"/\"$f\"_FE.tdf hg38" ; 
  echo $cmd1 > scripts/igvtdf/gentdf."$f".sh
  sed -i 's/"//g' scripts/igvtdf/gentdf."$f".sh
  qsub -cwd -N igvtdf."$f" -q light.q -e stdout/igvtdf/ -o stdout/igvtdf/ scripts/igvtdf/gentdf."$f".sh
  
done