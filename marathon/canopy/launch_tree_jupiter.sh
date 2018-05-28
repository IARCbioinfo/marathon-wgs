#!/bin/bash

log_folder="/data/soudadel/MPM/canopy/output/generated_clustering/"

#bsub -q "normal" -J "CanopyTree_5009" -n 1 -R "rusage[mem=8000]" -M 8000 -oo ${log_folder}tree_patient_5009.log -eo ${log_folder}tree_patient_5009.error.log "/data/soudadel/MPM/canopy/scripts/canopy_tree.R 5009 /data/soudadel/MPM/canopy/output/generated_clustering/5009/"

#bsub -q "normal" -J "CanopyTree_6063" -n 1 -R "rusage[mem=8000]" -M 8000 -oo ${log_folder}tree_patient_6063.log -eo ${log_folder}tree_patient_6063.error.log "/data/soudadel/MPM/canopy/scripts/canopy_tree.R 6063 /data/soudadel/MPM/canopy/output/generated_clustering/6063/"

#bsub -q "normal" -J "CanopyTree_8329" -n 1 -R "rusage[mem=8000]" -M 8000 -oo ${log_folder}tree_patient_8329.log -eo ${log_folder}tree_patient_8329.error.log "/data/soudadel/MPM/canopy/scripts/canopy_tree.R 8329 /data/soudadel/MPM/canopy/output/generated_clustering/8329/"

bsub -q "normal" -J "CanopyTree_12323" -n 1 -R "rusage[mem=8000]" -M 8000 -oo ${log_folder}tree_patient_12323.log -eo ${log_folder}tree_patient_12323.error.log "/data/soudadel/MPM/canopy/scripts/canopy_tree.R 12323 /data/soudadel/MPM/canopy/output/generated_clustering/12323/"

#bsub -q "normal" -J "CanopyTree_19392" -n 1 -R "rusage[mem=8000]" -M 8000 -oo ${log_folder}tree_patient_19392.log -eo ${log_folder}tree_patient_19392.error.log "/data/soudadel/MPM/canopy/scripts/canopy_tree.R 19392 /data/soudadel/MPM/canopy/output/generated_clustering/19392/"

#bsub -q "normal" -J "CanopyTree_LB110287" -n 1 -R "rusage[mem=8000]" -M 8000 -oo ${log_folder}tree_patient_LB110287.log -eo ${log_folder}tree_patient_LB110287.error.log "/data/soudadel/MPM/canopy/scripts/canopy_tree.R LB110287 /data/soudadel/MPM/canopy/output/generated_clustering/LB110287/"
