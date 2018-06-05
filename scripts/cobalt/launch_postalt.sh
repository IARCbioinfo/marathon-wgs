#!/bin/bash

#Run Post Alignment
#the ccc_msub command does not accept argument , use a template script

cd /ccc/work/cont007/fg0094/soudadel/MESO/post_align

COUNTER=1
for bam_file in `ls /ccc/scratch/cont007/fg0094/soudadel/MESO/BAM_raw/*.bam`

do

bam_file_list[$COUNTER]=$bam_file
bam_id=`basename ${bam_file_list[$COUNTER]} | cut -d. -f1 `
echo $bam_id
sed "s/MYBAMID/${bam_id}/g"  template_postalt.sh> postalt_${bam_id}.sh
ccc_msub -A fg0094 -o postalt_${bam_id}.log -e postalt_${bam_id}.err -T 40000 -c 28 -N 2 ./postalt_${bam_id}.sh
let COUNTER+=1

done
