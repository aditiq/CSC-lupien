#!/bin/bash

module load R/3.4.1

Rscript scripts/08_jaccard.R $1 $2
#for f in data/ConsensusSet/PCSC1/enh.split/PCSC1.Consensus.Catalogue.Enhancer.Binarymat* ; do name=$(echo $f| awk -F"/" '{ print $NF}' | sed -e 's/PCSC1.Consensus.Catalogue.Enhancer.Binarymat//g' ) ; qsub -cwd -N Jaccard.Enh."$name" -e stdout/jaccard/ -o stdout/jaccard/ -q light.q scripts/run_08.sh $f $name ; done

