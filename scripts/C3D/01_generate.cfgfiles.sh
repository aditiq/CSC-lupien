#!/bin/bash

## Run script to set up directory structure and anchor and ref files


## Generate signal matrix for ref.bed against db files


## Use that signal matrix for running the R script




#
#
# files="LSCp
# GBM
# PFA" ;
#
# for f in $files;
# do
#     mkdir -p /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/C3D/\"$f\"/
#     for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/c3danchors/c3d.anchors.gencodev24.500bparoundTSS.chr*;
#     do
#       name=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/c3d.anchors.gencodev24.500bparoundTSS.//g' | sed -e 's/.bed//g' )
#       cmd1="reference=/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/PCSC1/\"$f\".Consensus.Catalogue.narrowPeak.sorted
#             db=/mnt/work1/users/lupiengroup/People/qamraa99/common.data/ATAC-Catalogue/ATAC-ENCODE/hg38/atac.encode.hg38.bdg.files.txt
#             anchor=\"$k\"
#             outDirectory=/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/C3D/\"$f\"/
#             window=500000
#             correlationThreshold=0.7
#             correlationMethod=pearson
#             figures=n
#             figureWidth=1000000
#             zoom=100000
#             colours=red,blue,green,purple
#             tracks=y
#             assembly=hg38" ;
#       echo "${cmd1}" > scripts/C3D/"$f"."$name".cfg  ;
#       sed -i 's/"//g' scripts/C3D/"$f"."$name".cfg    ;
#     done
# done
#
# files="LSCp.dnase
# GBM.dnase
# PFA.dnase" ;
#
# for f in $files;
# do
#     mkdir -p /mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/C3D/\"$f\"/
#     for k in /mnt/work1/users/lupiengroup/People/qamraa99/common.data/c3danchors/c3d.anchors.gencodev24.500bparoundTSS.chr*;
#     do
#       name=$(echo $k | awk -F"/" '{ print $NF}' | sed -e 's/c3d.anchors.gencodev24.500bparoundTSS.//g' | sed -e 's/.bed//g' )
#       cmd1="reference=/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/data/ConsensusSet/PCSC1/\"$f\".Consensus.Catalogue.narrowPeak.sorted
#             db=/mnt/work1/users/lupiengroup/Projects/CommonData/ENCODE/hg38/dnase.hg38.bedgraph.txt
#             anchor=\"$k\"
#             outDirectory=/mnt/work1/users/lupiengroup/People/qamraa99/HG38.Pancancer.CSC/results/PCSC1/C3D/\"$f\"/
#             window=500000
#             correlationThreshold=0.7
#             correlationMethod=pearson
#             figures=n
#             figureWidth=1000000
#             zoom=100000
#             colours=red,blue,green,purple
#             tracks=y
#             assembly=hg38" ;
#       echo "${cmd1}" > scripts/C3D/"$f"."$name".cfg  ;
#       sed -i 's/"//g' scripts/C3D/"$f"."$name".cfg    ;
#     done
# done
