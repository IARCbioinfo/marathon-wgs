#!/bin/bash

module load extenv/fg
module load platypus

APPDIR=/ccc/cont007/dsku/lautrec/home/app/fg/fg/products/platypus-0.8.1/bin

NORMALBAM=/path/to/BAM_post_al/normal
TUMORBAM=/path/to/BAM_post_al/tumor
INPUTVCF=/path/to/calling_germline_without_genotype/patient1_normal.GERMLINE.vcf.gz
OUTPUTVCF=/path/to/calling_germline_with_genotype/patient1_normal.GERMLINE.genotype.vcf
REFFILES=/path/to/calling_somatic/ref_files

$APPDIR/Platypus.py callVariants --bamFiles=$NORMALBAM/patient1_normal.pa.bam,$TUMORBAM/patient1_tumor1.pa.bam,$TUMORBAM/patient1_tumor2.pa.bam --refFile=$FG_BIOBANK/by-name/Homo_sapiens/hs38dh/hs38dh_all_chr.fasta --regions=$REFFILES/hs38dh.bed --nCPU=12 --output=$OUTPUTVCF --source=$INPUTVCF --minPosterior=0 --getVariantsFromBAMs=0
