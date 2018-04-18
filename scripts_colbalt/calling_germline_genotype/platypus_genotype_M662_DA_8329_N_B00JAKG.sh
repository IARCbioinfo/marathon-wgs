#!/bin/bash

module load extenv/fg
module load platypus

APPDIR=/ccc/cont007/dsku/lautrec/home/app/fg/fg/products/platypus-0.8.1/bin

NORMALBAM=/ccc/scratch/cont007/fg0094/soudadel/MESO/BAM_post_al/normal
TUMORBAM=/ccc/scratch/cont007/fg0094/soudadel/MESO/BAM_post_al/tumor
INPUTVCF=/ccc/store/cont007/fg0094/soudadel/MESO/calling_germline_without_genotype/M662_DA_8329_N_B00JAKG.GERMLINE.vcf.gz
OUTPUTVCF=/ccc/store/cont007/fg0094/soudadel/MESO/calling_germline_with_genotype/M662_DA_8329_N_B00JAKG.GERMLINE.genotype.vcf
REFFILES=/ccc/work/cont007/fg0094/soudadel/MESO/calling_somatic/ref_files

$APPDIR/Platypus.py callVariants --bamFiles=$NORMALBAM/M662_DA_8329_N_B00JAKG.pa.bam,$TUMORBAM/M662_DA_8329_T_B00JAKF.pa.bam,$TUMORBAM/M662_DA_8329_T_B00JAKH.pa.bam --refFile=$FG_BIOBANK/by-name/Homo_sapiens/hs38dh/hs38dh_all_chr.fasta --regions=$REFFILES/hs38dh.bed --nCPU=12 --output=$OUTPUTVCF --source=$INPUTVCF --minPosterior=0 --getVariantsFromBAMs=0
