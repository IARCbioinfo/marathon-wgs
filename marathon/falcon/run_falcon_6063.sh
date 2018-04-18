#!/bin/bash

log_folder="/data/soudadel/MPM/falcon/output/patient_6063/"

for chr in {1..22}
do

mkdir -p ${log_folder}chr${chr}
bsub -q "normal" -J "Falcon_6063_${chr}" -n 1 -R "rusage[mem=4000]" -M 4000 -oo ${log_folder}patient_6063.chr${chr}.log -eo ${log_folder}patient_6063.chr${chr}.error.log "/home/soudadel/MPM/falcon/falcon.R /data/soudadel/MPM/falcon/input/M662_DA_6063_N_B00JAJG_chromosomes/M662_DA_6063_N_B00JAJG.GERMLINE.chr${chr}.vcf B00JAJE B00JAJF ${chr} ${log_folder}chr${chr}/"

done

for chr in "X" "Y"
do

mkdir -p ${log_folder}chr${chr}
bsub -q "normal" -J "Falcon_6063_${chr}" -n 1 -R "rusage[mem=4000]" -M 4000 -oo ${log_folder}patient_6063.chr${chr}.log -eo ${log_folder}patient_6063.chr${chr}.error.log "/home/soudadel/MPM/falcon/falcon.R /data/soudadel/MPM/falcon/input/M662_DA_6063_N_B00JAJG_chromosomes/M662_DA_6063_N_B00JAJG.GERMLINE.chr${chr}.vcf B00JAJE B00JAJF ${chr} ${log_folder}chr${chr}/"

done
