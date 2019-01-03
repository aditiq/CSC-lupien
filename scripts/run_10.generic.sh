#!/bin/bash

module load R/3.4.1

Rscript scripts/10_flexmix.generic.R $1 $2 $3 $4

#for f in {2..100} ; do qsub -N flexmix.binomial.enhancer."$f" -cwd -e stdout/flexmix/ -o stdout/flexmix/ -q all.q scripts/run_10.generic.sh data/ConsensusSet/PCSC1/PCSC1.Consensus.Catalogue.Enhancer.Binarymat.txt "Enhancer" $f $f ; done


