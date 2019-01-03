#!/bin/bash


#GBM
qsub -N GBM.repeat.fisher.jaccard -e stdout/repeat/ -o stdout/repeat/ -cwd \
scripts/repeat/03_fishertestfilegen.sh \
data/ConsensusSet/PCSC1/GBM.Consensus.Catalogue.narrowPeak.sorted \
data/ConsensusSet/Repeatanalysis/Background.merged.bed \
results/PCSC1/Repeatenrich/fishertest/ "GBM"

# PFA
qsub -N PFA.repeat.fisher.jaccard -e stdout/repeat/ -o stdout/repeat/ -cwd \
scripts/repeat/03_fishertestfilegen.sh \
data/ConsensusSet/PCSC1/PFA.Consensus.Catalogue.narrowPeak.sorted \
data/ConsensusSet/Repeatanalysis/Background.merged.bed \
results/PCSC1/Repeatenrich/fishertest/ "PFA"

#HF
qsub -N HF.repeat.fisher.jaccard -e stdout/repeat/ -o stdout/repeat/ -cwd \
scripts/repeat/03_fishertestfilegen.sh \
data/ConsensusSet/PCSC1/HF.Consensus.Catalogue.narrowPeak.sorted \
data/ConsensusSet/Repeatanalysis/Background.merged.bed \
results/PCSC1/Repeatenrich/fishertest/ "HF"

#LSCp
qsub -N LSCp.repeat.fisher.jaccard -e stdout/repeat/ -o stdout/repeat/ -cwd \
scripts/repeat/03_fishertestfilegen.sh \
data/ConsensusSet/PCSC1/LSCp.Consensus.Catalogue.narrowPeak.sorted \
data/ConsensusSet/Repeatanalysis/Background.merged.bed \
results/PCSC1/Repeatenrich/fishertest/ "LSCp"

#LSC-neg
qsub -N LSCneg.repeat.fisher.jaccard -e stdout/repeat/ -o stdout/repeat/ -cwd \
scripts/repeat/03_fishertestfilegen.sh \
data/ConsensusSet/PCSC1/LSC.neg.Consensus.Catalogue.narrowPeak.sorted \
data/ConsensusSet/Repeatanalysis/Background.merged.bed \
results/PCSC1/Repeatenrich/fishertest/ "LSCneg"

#LSC-bulk
qsub -N LSCbulkrepeat.fisher.jaccard -e stdout/repeat/ -o stdout/repeat/ -cwd \
scripts/repeat/03_fishertestfilegen.sh \
data/ConsensusSet/PCSC1/LSC.bulk.Consensus.Catalogue.narrowPeak.sorted \
data/ConsensusSet/Repeatanalysis/Background.merged.bed \
results/PCSC1/Repeatenrich/fishertest/ "LSCbulk"

#Hemat - HSCs
qsub -N HSC.repeat.fisher.jaccard -e stdout/repeat/ -o stdout/repeat/ -cwd \
scripts/repeat/03_fishertestfilegen.sh \
/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/RyanCorces.Hemat/hg38/peakCalls/ConsensusSet/HSC.ConsensusSet.bed.sorted \
data/ConsensusSet/Repeatanalysis/Background.merged.bed \
results/PCSC1/Repeatenrich/fishertest/ "HSC"

#Hemat - Progenitor
qsub -N Hemat.Progen.repeat.fisher.jaccard -e stdout/repeat/ -o stdout/repeat/ -cwd \
scripts/repeat/03_fishertestfilegen.sh \
/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/RyanCorces.Hemat/hg38/peakCalls/ConsensusSet/Hemat.Progenitor.ConsensusSet.bed.sorted \
data/ConsensusSet/Repeatanalysis/Background.merged.bed \
results/PCSC1/Repeatenrich/fishertest/ "Hemat.Progen"


#Hemat - Differentiated
qsub -N Hemat.differentiated.repeat.fisher.jaccard -e stdout/repeat/ -o stdout/repeat/ -cwd \
scripts/repeat/03_fishertestfilegen.sh \
/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/RyanCorces.Hemat/hg38/peakCalls/ConsensusSet/Hemat.differentiated.ConsensusSet.bed.sorted \
data/ConsensusSet/Repeatanalysis/Background.merged.bed \
results/PCSC1/Repeatenrich/fishertest/ "Hemat.differentiated"

# hESCs
qsub -N ESC.repeat.fisher.jaccard -e stdout/repeat/ -o stdout/repeat/ -cwd \
scripts/repeat/03_fishertestfilegen.sh \
/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ESC/hg38/Consensus.ESC.bed.sorted \
data/ConsensusSet/Repeatanalysis/Background.merged.bed \
results/PCSC1/Repeatenrich/fishertest/ "ESC"

# Brain Atlas + ENCODE brain
qsub -N NormalBrain.repeat.fisher.jaccard -e stdout/repeat/ -o stdout/repeat/ -cwd \
scripts/repeat/03_fishertestfilegen.sh \
/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/Brain.Consensus/Brain.Consensus.sortBed.sorted \
data/ConsensusSet/Repeatanalysis/Background.merged.bed \
results/PCSC1/Repeatenrich/fishertest/ "NormalBrain"
