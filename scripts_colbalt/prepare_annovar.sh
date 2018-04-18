rm -f M662_DA_12323_T_B00JAKC.normalized.final.converted.txt
rm -f M662_DA_12323_T_B00JAKC.normalized.final.converted.txt.avinput
rm -f M662_DA_12323_T_B00JAKC.normalized.final.txt
rm -f M662_DA_12323_T_B00JAKC.normalized.header.txt
rm -f M662_DA_12323_T_B00JAKC.normalized.OK.txt
rm -f M662_DA_12323_T_B00JAKC.normalized.txt

grep "#CHROM" M662_DA_12323_T_B00JAKC.normalized.vcf > M662_DA_12323_T_B00JAKC.normalized.header.txt
grep "PASS" M662_DA_12323_T_B00JAKC.normalized.vcf > M662_DA_12323_T_B00JAKC.normalized.txt
awk -F "\t" '{print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11}' M662_DA_12323_T_B00JAKC.normalized.txt > M662_DA_12323_T_B00JAKC.normalized.OK.txt
cat M662_DA_12323_T_B00JAKC.normalized.header.txt M662_DA_12323_T_B00JAKC.normalized.OK.txt > M662_DA_12323_T_B00JAKC.normalized.final.txt

rm -f M662_DA_12323_T_B00JAKC.normalized.header.txt
rm -f M662_DA_12323_T_B00JAKC.normalized.OK.txt
rm -f M662_DA_12323_T_B00JAKC.normalized.txt
