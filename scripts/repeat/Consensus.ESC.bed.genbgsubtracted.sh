module load bedtools/2.23.0 ; intersectBed -a data/ConsensusSet/Repeatanalysis/Background.merged.bed -b /mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ESC/hg38/Consensus.ESC.bed -v | sortBed -i stdin > data/ConsensusSet/Repeatanalysis/Bg.subtracted.Consensus.ESC.bed.bed
