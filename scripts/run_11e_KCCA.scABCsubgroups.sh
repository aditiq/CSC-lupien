#!/bin/bash

module load R/3.4.1

Rscript scripts/11e_KCCA.scABCsubgroups.R

#qsub -cwd -N scABC.KCCA -e stdout/ -o stdout/ scripts/run_11e_KCCA.scABCsubgroups.sh 
