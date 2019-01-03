awk 'BEGIN {FS=OFS="\t"} { if ($1 < 0.05) print $3}' scABC.Combined.LSCnegvsLSCpos.txt  | awk -F"_" '{ print $1"\t"$2"\t"$3}' > scABC.Combined.LSCnegvsLSCpos.LSCneg.p0.05.bed
awk 'BEGIN {FS=OFS="\t"} { if ($2 < 0.05) print $3}' scABC.Combined.LSCnegvsLSCpos.txt  | awk -F"_" '{ print $1"\t"$2"\t"$3}' > scABC.Combined.LSCnegvsLSCpos.LSCpos.p0.05.bed

