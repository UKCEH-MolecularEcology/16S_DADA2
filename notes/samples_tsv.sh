#!/usr/bin/env bash

## initial table
## header
#echo -e "Sample\tR1\tR2" > config/samples.tsv
## samples
#fpath="/mnt/isilon/projects/ecosystem_biology/circles/raw/fastq/"
#find "${fpath}" -type f -name "*_R1_001.fastq.gz" -exec basename -s '.fastq.gz' {} \; | sed 's/_16S.*//' | sort | while read sid; do
#    r1=$(find "${fpath}" -type f -name "${sid}_*_R1_001.fastq.gz")
#    r2=$(find "${fpath}" -type f -name "${sid}_*_R2_001.fastq.gz")
#    echo -e "${sid}\t${r1}\t${r2}"
#done >> config/samples.tsv


for file in /mnt/data/sbusi/vip/data/fastq/*.gz
do 
    echo $file
done | paste - - | awk '{print $0=$1"\t"$1"\t"$2}' | \
    awk -v "OFS=\t" '{$1=$1;sub(/^.*MB/,"MB",$1); print}' | \
    awk -v OFS="\t" '{ split($1, arr, "_"); $1=arr[1]"_"arr[2] }1' | \
    awk -v "OFS=\t" '{$1=$1;sub("_","-",$1); print}' | \
    sed $'1 i\\\nSample\tR1\tR2' > config/samples.tsv

# 2021.02.19: update (issue #2)
python notes/samples_tsv_02.py
