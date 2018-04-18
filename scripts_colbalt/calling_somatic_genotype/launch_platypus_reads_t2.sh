#!/usr/bin/env bash
# Run platypus
# the ccc_msub command does not accept argument , use a template script

cd /ccc/work/cont007/fg0094/soudadel/MESO/somatic_read_count_t2

ccc_msub -A fg0094 -o platypus_reads_t2_M662_DA_12323_T_B00JAKC.log -e platypus_reads_t2_M662_DA_12323_T_B00JAKC.err -T 80000 -n 28 ./platypus_reads_t2_M662_DA_12323_T_B00JAKC.sh

ccc_msub -A fg0094 -o platypus_reads_t2_M662_DA_12323_T_B00JAKD.log -e platypus_reads_t2_M662_DA_12323_T_B00JAKD.err -T 80000 -n 28 ./platypus_reads_t2_M662_DA_12323_T_B00JAKD.sh

ccc_msub -A fg0094 -o platypus_reads_t2_M662_DA_19392_T_B00JAKP.log -e platypus_reads_t2_M662_DA_19392_T_B00JAKP.err -T 80000 -n 28 ./platypus_reads_t2_M662_DA_19392_T_B00JAKP.sh

ccc_msub -A fg0094 -o platypus_reads_t2_M662_DA_19392_T_B00JAKQ.log -e platypus_reads_t2_M662_DA_19392_T_B00JAKQ.err -T 80000 -n 28 ./platypus_reads_t2_M662_DA_19392_T_B00JAKQ.sh

ccc_msub -A fg0094 -o platypus_reads_t2_M662_DA_5009_T_B00JAJB.log -e platypus_reads_t2_M662_DA_5009_T_B00JAJB.err -T 80000 -n 28 ./platypus_reads_t2_M662_DA_5009_T_B00JAJB.sh

ccc_msub -A fg0094 -o platypus_reads_t2_M662_DA_5009_T_B00JAJC.log -e platypus_reads_t2_M662_DA_5009_T_B00JAJC.err -T 80000 -n 28 ./platypus_reads_t2_M662_DA_5009_T_B00JAJC.sh

ccc_msub -A fg0094 -o platypus_reads_t2_M662_DA_6063_T_B00JAJE.log -e platypus_reads_t2_M662_DA_6063_T_B00JAJE.err -T 80000 -n 28 ./platypus_reads_t2_M662_DA_6063_T_B00JAJE.sh

ccc_msub -A fg0094 -o platypus_reads_t2_M662_DA_6063_T_B00JAJF.log -e platypus_reads_t2_M662_DA_6063_T_B00JAJF.err -T 80000 -n 28 ./platypus_reads_t2_M662_DA_6063_T_B00JAJF.sh

ccc_msub -A fg0094 -o platypus_reads_t2_M662_DA_8329_T_B00JAKF.log -e platypus_reads_t2_M662_DA_8329_T_B00JAKF.err -T 80000 -n 28 ./platypus_reads_t2_M662_DA_8329_T_B00JAKF.sh

ccc_msub -A fg0094 -o platypus_reads_t2_M662_DA_8329_T_B00JAKH.log -e platypus_reads_t2_M662_DA_8329_T_B00JAKH.err -T 80000 -n 28 ./platypus_reads_t2_M662_DA_8329_T_B00JAKH.sh

ccc_msub -A fg0094 -o platypus_reads_t2_M662_DA_LB110287_T_B00JAK0.log -e platypus_reads_t2_M662_DA_LB110287_T_B00JAK0.err -T 80000 -n 28 ./platypus_reads_t2_M662_DA_LB110287_T_B00JAK0.sh

ccc_msub -A fg0094 -o platypus_reads_t2_M662_DA_LB110287_T_B00JAK2.log -e platypus_reads_t2_M662_DA_LB110287_T_B00JAK2.err -T 80000 -n 28 ./platypus_reads_t2_M662_DA_LB110287_T_B00JAK2.sh



#cd /ccc/work/cont007/fg0094/soudadel/MESO/somatic_read_count_t2

#COUNTER=1
#for bam_file in `ls /ccc/scratch/cont007/fg0094/soudadel/MESO/BAM_post_al/tumor/*.bam`

#do

#bam_file_list[$COUNTER]=$bam_file
#bam_id=`basename ${bam_file_list[$COUNTER]} | cut -d"." -f1 `
#echo $bam_id
#sed "s/TUMORBAMID/${bam_id}/g" template_platypus_reads_t2.sh > platypus_reads_t2_${bam_id}.sh
#ccc_msub -A fg0094 -o platypus_reads_t2_${bam_id}.log -e platypus_reads_t2_${bam_id}.err -T 80000 -n 28 ./platypus_reads_t2_#${bam_id}.sh
#let COUNTER+=1

#done
