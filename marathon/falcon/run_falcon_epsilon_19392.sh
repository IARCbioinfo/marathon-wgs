#!/bin/bash

for chr in {1..22}
do

Rscript falcon_epsilon.R "/home/pgm/Workspace/MPM/marathon/falcon/output/patient_19392/chr$chr/falcon.patient_19392.tumor_B00JAKP.chr_$chr.Falcon.rda" "/home/pgm/Workspace/MPM/marathon/falcon/output/patient_19392/chr$chr/falcon.patient_19392.tumor_placeholder.chr_$chr.output.txt" "/home/pgm/Workspace/MPM/marathon/falcon/output/"

done

for chr in "X" "Y"
do

Rscript falcon_epsilon.R "/home/pgm/Workspace/MPM/marathon/falcon/output/patient_19392/chr$chr/falcon.patient_19392.tumor_B00JAKP.chr_$chr.Falcon.rda" "/home/pgm/Workspace/MPM/marathon/falcon/output/patient_19392/chr$chr/falcon.patient_19392.tumor_placeholder.chr_$chr.output.txt" "/home/pgm/Workspace/MPM/marathon/falcon/output/"

done
