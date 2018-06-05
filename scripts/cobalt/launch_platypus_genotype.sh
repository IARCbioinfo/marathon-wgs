#!/bin/bash

#Run Post Alignment
#the ccc_msub command does not accept argument , use a template script

cd /ccc/work/cont007/fg0094/soudadel/MESO/calling_germline/add_genotype

COUNTER=1
for vcf_file in `ls /ccc/store/cont007/fg0094/soudadel/MESO/calling_germline_without_genotype/*.vcf.gz`

do

vcf_file_list[$COUNTER]=$vcf_file
file_id=`basename ${vcf_file_list[$COUNTER]} | cut -d. -f1 `
echo $file_id
sed "s/MYVCFID/${file_id}/g"  template_platypus_genotype.sh > platypus_genotype_${file_id}.sh
#ccc_msub -A fg0094 -o platypus_genotype_${file_id}.log -e platypus_genotype_${file_id}.err -T 40000 -c 28 -N 2 ./platypus_genotype_${file_id}.sh
let COUNTER+=1

done
