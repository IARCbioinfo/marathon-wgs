#!/bin/bash

log_folder="/data/soudadel/MPM/falcon/output/patient_LB110287/"

for chr in {1..22}
do

mkdir -p ${log_folder}chr${chr}
bsub -q "normal" -J "Falcon_LB110287_${chr}" -n 1 -R "rusage[mem=4000]" -M 4000 -oo ${log_folder}patient_LB110287.chr${chr}.log -eo ${log_folder}patient_LB110287.chr${chr}.error.log "/home/soudadel/MPM/falcon/falcon.R /data/soudadel/MPM/falcon/input/M662_DA_LB110287_N_B00JAK1_chromosomes/M662_DA_LB110287_N_B00JAK1.GERMLINE.chr${chr}.vcf B00JAK0 B00JAK2 ${chr} ${log_folder}chr${chr}/"

done

for chr in "X" "Y"
do

mkdir -p ${log_folder}chr${chr}
bsub -q "normal" -J "Falcon_LB110287_${chr}" -n 1 -R "rusage[mem=4000]" -M 4000 -oo ${log_folder}patient_LB110287.chr${chr}.log -eo ${log_folder}patient_LB110287.chr${chr}.error.log "/home/soudadel/MPM/falcon/falcon.R /data/soudadel/MPM/falcon/input/M662_DA_LB110287_N_B00JAK1_chromosomes/M662_DA_LB110287_N_B00JAK1.GERMLINE.chr${chr}.vcf B00JAK0 B00JAK2 ${chr} ${log_folder}chr${chr}/"

done
