#!/bin/bash

module load extenv/fg
module load platypus

APPDIR=/ccc/cont007/dsku/lautrec/home/app/fg/fg/products/platypus-0.8.1/bin

NORMALBAM=/ccc/scratch/cont007/fg0094/soudadel/MESO/BAM_post_al/normal/MYBAMID.pa.bam
OUTPUTVCF=/ccc/store/cont007/fg0094/soudadel/MESO/calling_germline/MYBAMID.GERMLINE.vcf
REFFILES=/ccc/work/cont007/fg0094/soudadel/MESO/calling_somatic/ref_files

$APPDIR/Platypus.py callVariants --bamFiles=$NORMALBAM --output=$OUTPUTVCF --refFile=$FG_BIOBANK/by-name/Homo_sapiens/hs38dh/hs38dh_all_chr.fasta --nCPU=12 --regions=$REFFILES/hs38dh.bed --badReadsThreshold=0 --qdThreshold=0 --rmsmqThreshold=20 --hapScoreThreshold=10 --scThreshold=0.99
