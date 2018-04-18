declare -a arr=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for i in "${arr[@]}"
do
   #tabix -h /home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/M662_DA_5009_N_B00JAJD.GERMLINE.vcf.gz $i > /home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/M662_DA_5009_N_B00JAJD_chromosomes/M662_DA_5009_N_B00JAJD.GERMLINE.$i.vcf
   tabix -h /home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/M662_DA_12323_N_B00JAKE.GERMLINE.vcf.gz $i > /home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/M662_DA_12323_N_B00JAKE_chromosomes/M662_DA_12323_N_B00JAKE.GERMLINE.$i.vcf
   tabix -h /home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/M662_DA_19392_N_B00JAKR.GERMLINE.vcf.gz $i > /home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/M662_DA_19392_N_B00JAKR_chromosomes/M662_DA_19392_N_B00JAKR.GERMLINE.$i.vcf
   tabix -h /home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/M662_DA_6063_N_B00JAJG.GERMLINE.vcf.gz $i > /home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/M662_DA_6063_N_B00JAJG_chromosomes/M662_DA_6063_N_B00JAJG.GERMLINE.$i.vcf
   tabix -h /home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/M662_DA_8329_N_B00JAKG.GERMLINE.vcf.gz $i > /home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/M662_DA_8329_N_B00JAKG_chromosomes/M662_DA_8329_N_B00JAKG.GERMLINE.$i.vcf
   tabix -h /home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/M662_DA_LB110287_N_B00JAK1.GERMLINE.vcf.gz $i > /home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/M662_DA_LB110287_N_B00JAK1_chromosomes/M662_DA_LB110287_N_B00JAK1.GERMLINE.$i.vcf
done
