#!/usr/bin/env bash

cd ~/Workspace/MPM/calling_somatic_result/STRELKA/

for calling_dir in `ls ~/Workspace/MPM/calling_somatic_result/STRELKA`

do

echo `ls $calling_dir/results/variants/somatic.snvs.vcf.gz`
cp $calling_dir/results/variants/somatic.snvs.vcf.gz ~/Workspace/MPM/calling_somatic_result/$calling_dir.SOMATIC.vcf.gz
gunzip -d ~/Workspace/MPM/calling_somatic_result/$calling_dir.SOMATIC.vcf.gz

done
