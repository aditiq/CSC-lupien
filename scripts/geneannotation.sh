for f in results/PCSC1/Cluster/scABC/*.p0.05.bed ;
do
        sortBed -i $f | closestBed -a stdin -b ../common.data/gencode.v24.annotation.TSS.wdgeneinfo.sorted.bed -d |  \
        awk 'BEGIN {FS=OFS="\t"} { if ($NF<=100) print $12}' | sort -k1,1 -u  > "$f".Genes.GencodeTSS500bp.bed
done


