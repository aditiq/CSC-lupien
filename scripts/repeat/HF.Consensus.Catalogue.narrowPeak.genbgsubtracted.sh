module load bedtools/2.23.0 ; intersectBed -a data/ConsensusSet/Repeatanalysis/Background.merged.bed -b data/ConsensusSet/PCSC1/HF.Consensus.Catalogue.narrowPeak -v | sortBed -i stdin > data/ConsensusSet/Repeatanalysis/Bg.subtracted.HF.Consensus.Catalogue.narrowPeak.bed