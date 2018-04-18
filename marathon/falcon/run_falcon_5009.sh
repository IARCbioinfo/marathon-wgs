#!/bin/bash

log_folder="/data/soudadel/MPM/falcon/output/patient_5009/"

for chr in {1..22}
do

mkdir -p ${log_folder}chr${chr}
bsub -q "normal" -J "Falcon_5009_${chr}" -n 1 -R "rusage[mem=4000]" -M 4000 -oo ${log_folder}patient_5009.chr${chr}.log -eo ${log_folder}patient_5009.chr${chr}.error.log "/home/soudadel/MPM/falcon/falcon.R /data/soudadel/MPM/falcon/input/M662_DA_5009_N_B00JAJD_chromosomes/M662_DA_5009_N_B00JAJD.GERMLINE.chr${chr}.vcf B00JAJB B00JAJC ${chr} ${log_folder}chr${chr}/"

done

for chr in "X" "Y"
do

mkdir -p ${log_folder}chr${chr}
bsub -q "normal" -J "Falcon_5009_${chr}" -n 1 -R "rusage[mem=4000]" -M 4000 -oo ${log_folder}patient_5009.chr${chr}.log -eo ${log_folder}patient_5009.chr${chr}.error.log "/home/soudadel/MPM/falcon/falcon.R /data/soudadel/MPM/falcon/input/M662_DA_5009_N_B00JAJD_chromosomes/M662_DA_5009_N_B00JAJD.GERMLINE.chr${chr}.vcf B00JAJB B00JAJC ${chr} ${log_folder}chr${chr}/"

done
