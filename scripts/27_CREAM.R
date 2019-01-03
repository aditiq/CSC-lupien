#!/bin/bash

## Run CREAM to identify blocks of chromatin accessibility for CSCs
## This is to identify regions that have shared enhancers/promoters in one block but not at the exact same location
## Maybe CSCs cluster well on CORES identified from CREAM apart from HKG regions
## Thus CREAM will also be run on normal differentiated samples --- Kitchensink1 and 2


library(CREAM)

## CREAM on PCSC1
pcsc1.cream <- CREAM("data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.narrowPeak", WScutoff = 1.5, MinLength = 1000, peakNumMin = 2)
write.table(pcsc1.cream, file="results/PCSC1/CREAM/PCSC1.CREAM.bed", sep="\t", quote = F,row.names=F, col.names=F)

## CREAM on KitchenSink1
ks1.cream <- CREAM("data/ConsensusSet/KitchenSink1/KitchenSink1.Consensus.Catalogue.narrowPeak", WScutoff = 1.5, MinLength = 1000, peakNumMin = 2)
write.table(ks1.cream, file="results/PCSC1/CREAM/KitchenSink1.CREAM.bed", sep="\t",quote = F, row.names=F, col.names=F)

## CREAM on KitchenSink2
ks2.cream <- CREAM("data/ConsensusSet/KitchenSink2/KitchenSink2.Consensus.Catalogue.narrowPeak", WScutoff = 1.5, MinLength = 1000, peakNumMin = 2)
write.table(ks2.cream, file="results/PCSC1/CREAM/KitchenSink2.CREAM.bed", sep="\t",quote = F, row.names=F, col.names=F)

## CREAM on LSCp
ks2.cream <- CREAM("data/ConsensusSet/PCSC1/LSCp.Consensus.Catalogue.narrowPeak", WScutoff = 1.5, MinLength = 1000, peakNumMin = 2)
write.table(ks2.cream, file="results/PCSC1/CREAM/LSCp.CREAM.bed", sep="\t",quote = F, row.names=F, col.names=F)

## CREAM on GBM
ks2.cream <- CREAM("data/ConsensusSet/PCSC1/GBM.Consensus.Catalogue.narrowPeak", WScutoff = 1.5, MinLength = 1000, peakNumMin = 2)
write.table(ks2.cream, file="results/PCSC1/CREAM/GBM.CREAM.bed", sep="\t",quote = F, row.names=F, col.names=F)

## CREAM on PFA
ks2.cream <- CREAM("data/ConsensusSet/PCSC1/PFA.Consensus.Catalogue.narrowPeak", WScutoff = 1.5, MinLength = 1000, peakNumMin = 2)
write.table(ks2.cream, file="results/PCSC1/CREAM/PFA.CREAM.bed", sep="\t",quote = F, row.names=F, col.names=F)
