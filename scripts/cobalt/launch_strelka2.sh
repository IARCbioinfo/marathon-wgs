#!/usr/bin/env bash
# Run Strelka2
# the ccc_msub command does not accept argument , use a template script

cd /ccc/work/cont007/fg0094/soudadel/MESO/calling_somatic

COUNTER=1
for bam_file in `ls /ccc/scratch/cont007/fg0094/soudadel/MESO/BAM_post_al/tumor/*.bam`

do

bam_file_list[$COUNTER]=$bam_file
bam_id=`basename ${bam_file_list[$COUNTER]} | cut -d"." -f1 `
echo $bam_id
sed "s/TUMORBAMID/${bam_id}/g" template_strelka2.sh > strelka2_${bam_id}.sh
#ccc_msub -A fg0094 -o strelka2_${bam_id}.log -e strelka2_${bam_id}.err -T 80000 -n 28 ./strelka2_${bam_id}.sh
let COUNTER+=1

done
