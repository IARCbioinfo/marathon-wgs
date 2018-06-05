#!/bin/bash

#Run Post Alignment
#the ccc_msub command does not accept argument , use a template script

cd /ccc/work/cont007/fg0094/soudadel/MESO/calling_germline

COUNTER=1
for bam_file in `ls /ccc/scratch/cont007/fg0094/soudadel/MESO/BAM_post_al/normal/*.bam`

do

bam_file_list[$COUNTER]=$bam_file
bam_id=`basename ${bam_file_list[$COUNTER]} | cut -d. -f1 `
echo $bam_id
sed "s/MYBAMID/${bam_id}/g"  template_platypus.sh> platypus_${bam_id}.sh
ccc_msub -A fg0094 -o platypus_${bam_id}.log -e platypus_${bam_id}.err -T 40000 -c 28 -N 2 ./platypus_${bam_id}.sh
let COUNTER+=1

done
