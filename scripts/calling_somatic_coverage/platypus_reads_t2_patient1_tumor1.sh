#!/bin/bash

module load extenv/fg
module load platypus

APPDIR=/ccc/cont007/dsku/lautrec/home/app/fg/fg/products/platypus-0.8.1/bin

TUMORBAM=/path/to/BAM_post_al/tumor/patient1_tumor1.pa.bam
INPUTVCF=/path/to/somatic_read_count_t2/input/patient1_tumor2.PASS.vcf.gz
OUTPUTVCF=/path/to/somatic_read_count_t2/output/patient1_tumor1.other_tumor_positions.vcf
REFFILES=/path/to/ref_files

$APPDIR/Platypus.py callVariants --bamFiles=$TUMORBAM --refFile=$FG_BIOBANK/by-name/Homo_sapiens/hs38dh/hs38dh_all_chr.fasta --regions=$REFFILES/hs38dh.bed --nCPU=12 --output=$OUTPUTVCF --source=$INPUTVCF --minPosterior=0 --getVariantsFromBAMs=0
