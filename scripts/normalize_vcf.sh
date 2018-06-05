#!/bin/bash
if [ "$#" -eq  "0" ]
then
  echo "###############################################################"
  echo "## Normalize VCF                                             ##"
  echo "##                                                           ##"
  echo "## Missing argument \$1 : path to reference file              ##"
  echo "###############################################################"
else
  for vcffile in `ls *.vcf.gz`
  do
    vcfiden=`basename $vcffile | cut -d. -f1 `
    echo $vcfiden
    vt decompose -s $vcffile | vt decompose_blocksub -a - | vt normalize -r $1 -q - | vt uniq - > ./$vcfiden.normalized.vcf
  done
fi
