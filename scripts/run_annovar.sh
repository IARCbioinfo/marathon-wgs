#!/bin/bash

#annovar="/mnt/beegfs/old_appli57/annovar/";

if [ "$#" -eq  "0" ]
then
  echo "###############################################################"
  echo "## Run Annovar (variant annotation)                          ##"
  echo "##                                                           ##"
  echo "## Missing argument \$1 : path to Annovar folder              ##"
  echo "###############################################################"
else
  annovar = $1

  COUNTER=1

  for strelka_file in `ls *.vcf`
  do
    strelka_file_list[$COUNTER]=$strelka_file
    individual_id=`basename ${strelka_file_list[$COUNTER]}`
    echo $strelka_file
    echo $individual_id

    grep "#CHROM" $strelka_file > ${strelka_file}_header.txt
    grep "PASS" $strelka_file > ${strelka_file}.txt
    awk -F "\t" '{print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11}' ${strelka_file}.txt > ${strelka_file}_OK.txt
    cat ${strelka_file}_header.txt ${strelka_file}_OK.txt > ${strelka_file}_final.txt

    ${annovar}annovar_2017jul16/convert2annovar.pl ${strelka_file}_final.txt -format vcf4old -includeinfo -outfile ${strelka_file}_converted_snvs
    ${annovar}annovar_2017jul16/table_annovar.pl ${strelka_file}_converted_snvs ${annovar}Annovar_DB/hg38db -buildver hg38 -out ${strelka_file}_multianno -remove -protocol refGene -operation g -otherinfo -nastring NA -otherinfo

    rm ${strelka_file}_header.txt
    rm ${strelka_file}.txt
    rm ${strelka_file}_OK.txt

    let COUNTER+=1

  done
fi
