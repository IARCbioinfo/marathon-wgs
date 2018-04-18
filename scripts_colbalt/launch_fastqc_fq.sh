#!/bin/bash

#Run fastqc on fq
#the ccc_msub command does not accept argument , use a template script


COUNTER=1

cd /ccc/work/cont007/fg0094/soudadel/MESO/fastqc

for fq_file in `find /ccc/store/cont007/fg0094/fg0094/rawdata/projet_MESOTHEWG_662 -type f -name '*.fastq.gz'`

do

fq_file_list[$COUNTER]=$fq_file
fq_identif=`basename -s .fastq.gz  ${fq_file_list[$COUNTER]}`
fq_id=`echo ${fq_file_list[$COUNTER]} | cut -d/ -f2 `
echo $fq_file
echo $fq_identif
sed "s+FQFILE+${fq_file}+g" template_fastqc_fq.sh > fastqc_fq_${fq_identif}.sh
#ccc_msub -A fg0094 -o fastqc_fq_${fq_identif}.log -e fastqc_fq_${fq_identif}.err -T 10000 -n 15 ./fastqc_fq_${fq_identif}.sh
let COUNTER+=1

done
