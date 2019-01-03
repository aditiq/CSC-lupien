#!/bin/bash

module load R/3.4.1
# Arg 1 - No. of clusters

Rscript scripts/06_KCCA.Enhancer.R $1
