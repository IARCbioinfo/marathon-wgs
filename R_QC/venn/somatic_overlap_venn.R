#!/usr/bin/env Rscript

library(VennDiagram)
source("https://gist.githubusercontent.com/mfoll/a4dfbb92068dc559f130/raw/714dc8c2e97987fd4385dcef2722b3ef986d38d6/get_vcf_data.r")
setwd("/home/pgm/Workspace/MPM/VCF_finaux/somatic_sandbox")

args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) { stop("Input name missing!\n", call.=FALSE) }

tumor1_name = args[1]
tumor2_name = args[2]

patient_id = unlist(strsplit(tumor1_name, "_"))[3]
sample1_id = unlist(strsplit(tumor1_name, "_"))[5]
sample2_id = unlist(strsplit(tumor2_name, "_"))[5]

tumor1_VCFfile = paste(tumor1_name, '.normalized.vcf_multianno.hg38_multianno.txt', sep='')
tumor2_VCFfile = paste(tumor2_name, '.normalized.vcf_multianno.hg38_multianno.txt', sep='')

tumor1_VCFcontent = read.table(tumor1_VCFfile, sep="\t", header=FALSE, stringsAsFactors=FALSE)
tumor2_VCFcontent = read.table(tumor2_VCFfile, sep="\t", header=FALSE, stringsAsFactors=FALSE)

generated_files = '/home/pgm/Workspace/MPM/R_controle_qualite/venn/generated_files_somatic/'


############################################
## Generate unique ID and intersect
############################################
tumor1_VCFcontent$uniqid <- paste(tumor1_VCFcontent$V1, tumor1_VCFcontent$V2, tumor1_VCFcontent$V4, tumor1_VCFcontent$V5, sep='_')
tumor2_VCFcontent$uniqid <- paste(tumor2_VCFcontent$V1, tumor2_VCFcontent$V2, tumor2_VCFcontent$V4, tumor2_VCFcontent$V5, sep='_')

inter = intersect(tumor1_VCFcontent$uniqid, tumor2_VCFcontent$uniqid)

nb_tumor1 = length(tumor1_VCFcontent$uniqid)
nb_tumor2 = length(tumor2_VCFcontent$uniqid)
nb_inter  = length(inter)


############################################
## Stacked bars with tumor 1 & 2
############################################
svg(filename=paste(generated_files, patient_id, "_", sample2_id, "_on_", sample1_id, ".svg", sep=''))
grid.newpage()
draw.pairwise.venn(nb_tumor1, nb_tumor2, nb_inter, category = c(sample1_id, sample2_id), lty = rep("blank", 2), fill = c(rgb(1,0,0,0.4), rgb(0,0,1,0.3)), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
dev.off() # close the device
############################################


# STRELKA : An additional feature of the somatic caller is that it uses two calling tiers to reduce false positives.
# The first tier (tier1) is a set of input data filtration and model parameter settings with relatively
# stringent values, whereas the second tier (tier2) uses more permissive settings. All calls are
# initially  made  using  tier1  settings,  after  which  the  variant  is  called  again using tier2.

# refCounts = Value of FORMAT column $REF + “U” (e.g. if REF="A" then use the value in FOMRAT/AU)
# altCounts = Value of FORMAT column $ALT + “U” (e.g. if ALT="T" then use the value in FOMRAT/TU)
# tier1RefCounts = First comma-delimited value from $refCounts
# tier1AltCounts = First comma-delimited value from $altCounts
# Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
