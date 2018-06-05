for vcffile in `ls /home/pgm/Workspace/MPM/calling_germline_with_genotype/indexed/*.vcf.gz`

do

vcfiden=`basename $vcffile | cut -d. -f1 `
echo $vcfiden
~/Programs/vt/vt decompose -s $vcffile | ~/Programs/vt/vt decompose_blocksub -a - | ~/Programs/vt/vt normalize -r /home/pgm/Workspace/MPM/ref_files/hs38dh_all_chr.fasta -q - | ~/Programs/vt/vt uniq - > /home/pgm/Workspace/MPM/calling_germline_with_genotype/normalized/$vcfiden.normalized.vcf

done
