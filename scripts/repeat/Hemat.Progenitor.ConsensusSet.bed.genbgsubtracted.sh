module load bedtools/2.23.0 ; intersectBed -a data/ConsensusSet/Repeatanalysis/Background.merged.bed -b /mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/RyanCorces.Hemat/hg38/peakCalls/ConsensusSet/Hemat.Progenitor.ConsensusSet.bed -v | sortBed -i stdin > data/ConsensusSet/Repeatanalysis/Bg.subtracted.Hemat.Progenitor.ConsensusSet.bed.bed
