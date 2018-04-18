#!/bin/sh

annovar="/mnt/beegfs/old_appli57/annovar/annovar_2017jul16/";

##### ANNOTATION with Annovar

COUNTER=1

for strelka_file in `ls /data/soudadel/MPM/calling_somatic/*.vcf`
do
  strelka_file_list[$COUNTER]=$strelka_file
  individual_id=`basename ${strelka_file_list[$COUNTER]}`
  echo $strelka_file
  echo $individual_id

  #grep "#CHROM" $strelka_file > ${strelka_file}_header.txt
  #grep "PASS" $strelka_file > ${strelka_file}.txt
  #awk -F "\t" '{print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11}' ${strelka_file}.txt > ${strelka_file}_OK.txt
  #cat ${strelka_file}_header.txt ${strelka_file}_OK.txt > ${strelka_file}_final.txt

  #bsub -q "normal" -m "cn01 cn02 cn03 cn04 cn05 cn06 cn07 cn08 cn09 cn10 cn11 cn12" -J ${individual_id} -n 8 -R "rusage[mem=8000]" -M 8000 -oo ${log_folder}/joblog_convert_snvs_${individual_id}.txt -eo ${log_folder}/joberror_convert_snvs_${individual_id}.txt "${annovar}convert2annovar.pl ${strelka_file}_final.txt -format vcf4old -includeinfo -outfile ${strelka_file}_converted_snvs"
  bsub -q "normal" -m "cn01 cn02 cn03 cn04 cn05 cn06 cn07 cn08 cn09 cn10 cn11 cn12" -J ${individual_id} -n 8 -R "rusage[mem=8000]" -M 8000 -oo ${log_folder}/joblog_annot_snvs_${individual_id}.txt -eo ${log_folder}/joberror_annot_snvs_${individual_id}.txt "${annovar}table_annovar.pl ${strelka_file}_converted_snvs /mnt/beegfs/old_appli57/annovar/Annovar_DB/hg38db -buildver hg38 -out ${strelka_file}_multianno -remove -protocol refGene -operation g -otherinfo -nastring NA -otherinfo"

  #rm ${strelka_file}_header.txt
  #rm ${strelka_file}.txt
  #rm ${strelka_file}_OK.txt

  let COUNTER+=1

done
