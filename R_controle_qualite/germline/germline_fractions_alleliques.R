#!/usr/bin/env Rscript

source("https://gist.githubusercontent.com/mfoll/a4dfbb92068dc559f130/raw/714dc8c2e97987fd4385dcef2722b3ef986d38d6/get_vcf_data.r")
generated_files = '/home/pgm/Workspace/MPM/R_controle_qualite/germline/generated_files/'
setwd("/home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox")

args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) { stop("Input name missing!\n", call.=FALSE) }

# sample_name = 'M662_DA_12323_N_B00JAKE'
sample_name = args[1]

cat("Loading VCF file...")
VCFfile = paste(sample_name, '.GERMLINE.vcf', sep='')
VCFcontent = read.table(VCFfile, sep="\t", header=FALSE, stringsAsFactors=FALSE)
cat("Loaded.")

############################################
## Calculate allelic fractions
############################################
DPCounts <- unlist(lapply(1:nrow(VCFcontent), function(i) get_genotype(VCFcontent[i,10], VCFcontent[i, 9], "NR")))
altCounts <- unlist(lapply(1:nrow(VCFcontent), function(i) get_genotype(VCFcontent[i,10], VCFcontent[i, 9], "NV")))

allelic_fractions = altCounts / DPCounts

svg(filename=paste(generated_files, sample_name, ".svg", sep=''))
hist(allelic_fractions, breaks=50, main="Allelic fractions", xlim=c(0,1))
dev.off() # close the device
############################################

# GT       # Un-phased genotype calls
# GL       # Genotype log-likelihoods (natural log) for AA,AB and BB genotypes, where A = ref and B = variant. Only applicable for bi-allelic sites
# GOF      # Phred-scaled goodness-of-fit score for the genotype call
# GQ       # Phred-scaled quality score for the genotype call
# NR  = DP # Number of reads covering variant position in this sample
# NV  = AO (alt ) # Number of reads at variant position which support the called variant in this sample