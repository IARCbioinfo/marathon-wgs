for vcffile in `ls ./*.vcf`

do

vcfiden=`basename $vcffile | cut -d. -f1 `
echo $vcfiden
cat $vcffile | grep "PASS" > ./${vcfiden}.PASS.vcf

done
