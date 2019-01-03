#!/bin/bash

module load deeptools/2.4.2     

computeMatrix reference-point -R $1 -S $2 \
-o $3.matrix.gz --upstream $4 --downstream $5 --binSize $6