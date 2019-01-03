#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# R-3.5.0
# mordor
# Objective : Differential test using VGAM binomial from package monocle. Similar to cistopic but also provides a qvalue
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

#!/bin/bash

module load R/3.5.0

files="CSC.ESC
GBM.BRAIN
GBM.ESC
GBM.HF
GBM.PFA.HF
LSCp.ESC
LSCp.HEMATDIFF
LSCp.HEMATHSC
LSCp.LSCb
LSCposandneg
PFA.BRAIN
PFA.ESC"

for f in $files ;
do

  qsub -cwd -N difftest."$f" -q lupiengroup \
    -e stdout/monocle.difftest -o stdout/monocle.difftest \
    -b y /mnt/work1/software/R/3.5.0/Rscript scripts/20_monocledifftest.R data/ConsensusSet/PCSC1/"$f".Consensus.Catalogue.Binarymat.txt results/PCSC1/monocle.difftest/"$f".anno.txt results/PCSC1/monocle.difftest/"$f".difftest.txt
done



## KitchenSink2
qsub -cwd -N difftest.KitchenSink2 -q lupiengroup -e stdout/monocle.difftest/ -o stdout/monocle.difftest/ -b y /mnt/work1/software/R/3.5.0/Rscript scripts/20_monocledifftest.R data/ConsensusSet/KitchenSink2/KitchenSink2.Consensus.Catalogue.Binarymat.txt results/KitchenSink2/monocle.difftest/KitchenSink2.anno.txt results/KitchenSink2/monocle.difftest/KitchenSink2.difftest.txt

## KitchenSink2 Tissue Adjusted
qsub -cwd -N difftest.KitchenSink2.TissueAdj -q lupiengroup -e stdout/monocle.difftest/ -o stdout/monocle.difftest/ -b y /mnt/work1/software/R/3.5.0/Rscript scripts/20_monocledifftest.KitchenSink2.R data/ConsensusSet/KitchenSink2/KitchenSink2.Consensus.Catalogue.Binarymat.txt results/KitchenSink2/monocle.difftest/KitchenSink2.anno.txt results/KitchenSink2/monocle.difftest/KitchenSink2.TissueAdj.difftest.txt

## Beta patch
files="CSC.ESC
GBM.BRAIN
GBM.ESC
GBM.HF
GBM.PFA.HF
LSCp.ESC
LSCp.HEMATDIFF
LSCp.HEMATHSC
LSCp.LSCb
LSCposandneg
PFA.BRAIN
PFA.ESC"

for f in $files ;
do

  qsub -cwd -N difftest.betapatch."$f" -q lupiengroup \
    -e stdout/monocle.difftest -o stdout/monocle.difftest \
    -b y /mnt/work1/software/R/3.5.0/Rscript scripts/20_monocledifftest.betapatch.R data/ConsensusSet/PCSC1/"$f".Consensus.Catalogue.Binarymat.txt results/PCSC1/monocle.difftest/"$f".anno.txt results/PCSC1/monocle.difftest/"$f".difftest.betapatch.txt
done

## KitchenSink2
qsub -cwd -N difftest.betapatch.KitchenSink2 -q lupiengroup -e stdout/monocle.difftest/ -o stdout/monocle.difftest/ -b y /mnt/work1/software/R/3.5.0/Rscript scripts/20_monocledifftest.betapatch.R data/ConsensusSet/KitchenSink2/KitchenSink2.Consensus.Catalogue.Binarymat.txt results/KitchenSink2/monocle.difftest/KitchenSink2.anno.txt results/KitchenSink2/monocle.difftest/KitchenSink2.difftest.betapatch.txt

## KitchenSink2 Tissue Adjusted
qsub -cwd -N difftest.betapatch.KitchenSink2.TissueAdj -q lupiengroup -e stdout/monocle.difftest/ -o stdout/monocle.difftest/ -b y /mnt/work1/software/R/3.5.0/Rscript scripts/20_monocledifftest.KitchenSink2.betapatch.R data/ConsensusSet/KitchenSink2/KitchenSink2.Consensus.Catalogue.Binarymat.txt results/KitchenSink2/monocle.difftest/KitchenSink2.anno.txt results/KitchenSink2/monocle.difftest/KitchenSink2.TissueAdj.difftest.betapatch.txt
