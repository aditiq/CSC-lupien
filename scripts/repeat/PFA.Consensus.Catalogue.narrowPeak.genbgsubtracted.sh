module load bedtools/2.23.0 ; intersectBed -a data/ConsensusSet/Repeatanalysis/Background.merged.bed -b data/ConsensusSet/PCSC1/PFA.Consensus.Catalogue.narrowPeak -v | sortBed -i stdin > data/ConsensusSet/Repeatanalysis/Bg.subtracted.PFA.Consensus.Catalogue.narrowPeak.bed