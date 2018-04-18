#!/usr/bin/env bash

for vcf in `ls /home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/*.vcf`
do

vcf_iden_id=`basename ${vcf} | cut -d. -f1 `
echo $vcf_iden_id
/home/pgm/Workspace/MPM/R_controle_qualite/germline/germline_fractions_alleliques.R $vcf_iden_id

done
