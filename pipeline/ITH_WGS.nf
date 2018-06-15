#!/usr/bin/env nextflow

params.input_BAM = null
params.input_ref = null
params.input_regions = null

// BAM correspondance file
bams = file(params.input_BAM)

bams_repo = Channel.fromPath(bams)
  .splitCsv(header: ['patient_id', 'normal', 'tumor1', 'tumor2'], skip: 1 )
  .map {  row ->
  [row.patient_id ,
    [file(row.normal), file(row.normal.replace(".bam",".bam.bai"))],
    [file(row.tumor1), file(row.tumor1.replace(".bam",".bam.bai"))],
    [file(row.tumor2), file(row.tumor2.replace(".bam",".bam.bai"))]
  ] }

// Print repo
//repo = bams_repo.map {path -> path }
//repo.subscribe {println it}

/* Data access
patient_id : path[0]
normal array : path[1]
normal bam : path[1][0]
normal bai : path[1][1]
tumor1 array : path[2]
tumor2 array : path[3]
*/

bams_repo.into { bams_repo ; bam ; tumor1 ; tumor2 }
normal = bam.map { path -> path[1][0] }
tumor1 = tumor1.map { path -> path[2][0] }
tumor2 = tumor2.map { path -> path[3][0] }
all_bams = normal.mix(tumor1, tumor2)
patient_bams = bam.map { path -> [ path[1][0], path[2][0], path[3][0] ] }

/** PROCESSES **/

process reads_QC {
  input:
  // fastQC

  output:
  // some HTML files

  shell:
  // TODO
}

// Etape de post-alignement à finir
process post_alignment {
  input:
  file i from all_bams

  output:
  file("${i}.HEAD") into BAM_head

  shell:
  '''
  samtools view -h ${i} | k8 bwa-postalt.js $APPDIR/bwakit-0.7.15/resource-GRCh38/hs38DH.fa.alt | \
  sambamba view -S -f bam -l 0 /dev/stdin | \
  sambamba sort -t 50 -m 6G --tmpdir=/ccc/scratch/cont007/fg0094/soudadel/MESO_tmp -o /ccc/store/cont007/fg0094/soudadel/MESO/post_align/MYBAMID_pa.bam /dev/stdin

  //head -n1 !{i} > !{i}.HEAD
  '''
}

process germline_calling {
  input:
  file i from normal
  file input_ref
  file input_regions

  output:
  file("${i}.germline.vcf") into VCF_germline

  shell:
  '''
  Platypus.py callVariants --bamFiles=!{i} --output=!{i}.germline.vcf --refFile=!{input_ref} --regions=!{input_regions} --badReadsThreshold=0 --qdThreshold=0 --rmsmqThreshold=20 --hapScoreThreshold=10 --scThreshold=0.99
  '''
}

process somatic_calling {
  input:
  file i from patient_bams
  file input_ref
  file input_regions

  output:
  file("${i}.somatic.vcf") into VCF_somatic

  shell:
  '''
  strelka2.sh .......
  '''
}

process germline_tumor_coverage {
  // Platypus
}

process somatic_tumor_coverage {
  // Platypus
}

process quality_control {
  // 4 scripts R
}

process split_VCF {
  // split VCF by chromosome
}

process copy_numbers {
  // run falcon
}

process copy_numbers_epsilon {
  // run falcon_epsilon (to compute errors)
}

process compile_copy_numbers {
  // compile falcon results in a TSV file
}

process MCMC {
  // canopy pre-clustering & MCMC sampling
}

process tree {
  // canopy tree
}

process filter {
  // filter informative events
}
