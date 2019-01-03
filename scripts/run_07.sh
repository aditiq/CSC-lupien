#!/bin/bash

module load R/3.4.1

Rscript scripts/07_mcc.R $1 $2

#for f in data/ConsensusSet/PCSC1/enh.split/PCSC1.Consensus.Catalogue.Enhancer.Binarymat* ; do name=$(echo $f| awk -F"/" '{ print $NF}' | sed -e 's/PCSC1.Consensus.Catalogue.Enhancer.Binarymat//g' ) ; qsub -cwd -N MCC.Enh."$name" -e stdout/mcc -o stdout/mcc/ -q all.q scripts/run_07.sh $f $name ; done  
