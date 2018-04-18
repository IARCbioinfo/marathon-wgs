#!/bin/bash

module load extenv/fg
module load platypus

APPDIR=/ccc/cont007/dsku/lautrec/home/app/fg/fg/products/platypus-0.8.1/bin

TUMORBAM=/ccc/scratch/cont007/fg0094/soudadel/MESO/BAM_post_al/tumor/M662_DA_19392_T_B00JAKQ.pa.bam
INPUTVCF=/ccc/work/cont007/fg0094/soudadel/MESO/somatic_read_count_t2/input/M662_DA_19392_T_B00JAKP.PASS.vcf.gz
OUTPUTVCF=/ccc/work/cont007/fg0094/soudadel/MESO/somatic_read_count_t2/output/M662_DA_19392_T_B00JAKQ.other_tumor_positions.vcf
REFFILES=/ccc/work/cont007/fg0094/soudadel/MESO/calling_somatic/ref_files

$APPDIR/Platypus.py callVariants --bamFiles=$TUMORBAM --refFile=$FG_BIOBANK/by-name/Homo_sapiens/hs38dh/hs38dh_all_chr.fasta --regions=$REFFILES/hs38dh.bed --nCPU=12 --output=$OUTPUTVCF --source=$INPUTVCF --minPosterior=0 --getVariantsFromBAMs=0
