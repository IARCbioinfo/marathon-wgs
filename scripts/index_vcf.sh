for vcffile in `ls ./*.vcf`

do

echo $vcffile

sed 's/^##FORMAT=<ID=NV,Number=.,/##FORMAT=<ID=NV,Number=A,/1g' $vcffile | sed 's/^##FORMAT=<ID=NR,Number=.,/##FORMAT=<ID=NR,Number=A,/1g' > ${vcffile}_sedded.vcf #correct header format
rm $vcffile
mv ${vcffile}_sedded.vcf $vcffile
cp $vcffile ${vcffile}_bck

bgzip ${vcffile}
tabix -p vcf $vcffile.gz
mv ${vcffile}_bck $vcffile

done
