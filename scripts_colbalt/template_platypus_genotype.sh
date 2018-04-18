#!/bin/bash

module load extenv/fg
module load platypus

APPDIR=/ccc/cont007/dsku/lautrec/home/app/fg/fg/products/platypus-0.8.1/bin

TUMORBAM=/ccc/scratch/cont007/fg0094/soudadel/MESO/BAM_post_al/tumor/MYVCFID.pa.bam
INPUTVCF=/ccc/store/cont007/fg0094/soudadel/MESO/calling_germline_without_genotype/MYVCFID.GERMLINE.vcf.gz
OUTPUTVCF=/ccc/store/cont007/fg0094/soudadel/MESO/calling_germline_with_genotype/MYVCFID.GERMLINE.genotype.vcf
REFFILES=/ccc/work/cont007/fg0094/soudadel/MESO/calling_somatic/ref_files

$APPDIR/Platypus.py callVariants --bamFiles=tumor1.bam,tumor2.bam --refFile=$FG_BIOBANK/by-name/Homo_sapiens/hs38dh/hs38dh_all_chr.fasta --regions=$REFFILES/hs38dh.bed --nCPU=12 --output=$OUTPUTVCF --source=$INPUTVCF --minPosterior=0 --getVariantsFromBAMs=0
