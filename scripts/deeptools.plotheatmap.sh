#!/bin/bash

module load deeptools/2.4.2     

plotHeatmap -m $1 --outFileName $2.pdf --outFileNameMatrix $2.heatmap.mat 