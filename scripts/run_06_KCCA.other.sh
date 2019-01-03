#!/bin/bash

module load R/3.4.1

# $1= binarydsfile
# $2= opname
# $3= numK
Rscript scripts/06_KCCA.other.R $1 $2 $3
