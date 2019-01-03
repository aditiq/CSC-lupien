#!/bin/bash

## Generate table of alignment stats 
## Script Run in directory /mnt/work1/users/lupiengroup/People/qamraa99/HG38.HG38.Pancancer.CSC


stats <- dir("data/bams/", pattern="stats.txt", all.files= T,    full.names= T, recursive= T,   include.dirs= T)
peaks <- dir("data/peaks/", pattern=".narrowPeak", all.files= T,    full.names= T, recursive= T,   include.dirs= T)
meta <- read.table("data/metadata.txt", header=T, sep="\t", stringsAsFactors=F)