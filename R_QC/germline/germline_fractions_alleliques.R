#!/usr/bin/env Rscript
source("https://gist.githubusercontent.com/mfoll/a4dfbb92068dc559f130/raw/714dc8c2e97987fd4385dcef2722b3ef986d38d6/get_vcf_data.r")

args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) { stop("Input missing!\n", call.=FALSE) }

path_to_germline_VCF     = args[1]
path_to_output_directory = args[2]

path_to_germline_VCF_spl = unlist(strsplit(path_to_germline_VCF, "/"))
sample_name = gsub(".vcf", "", path_to_germline_VCF_spl[length(path_to_germline_VCF_spl)])

cat("Loading VCF file...")
VCFcontent = read.table(path_to_germline_VCF, sep="\t", header=FALSE, stringsAsFactors=FALSE)
cat("Loaded.")

############################################
## Calculate allelic fractions
############################################
DPCounts <- unlist(lapply(1:nrow(VCFcontent), function(i) get_genotype(VCFcontent[i,10], VCFcontent[i, 9], "NR")))
altCounts <- unlist(lapply(1:nrow(VCFcontent), function(i) get_genotype(VCFcontent[i,10], VCFcontent[i, 9], "NV")))

allelic_fractions = altCounts / DPCounts

svg(filename=paste(path_to_output_directory, "/", sample_name, ".svg", sep=''))
hist(allelic_fractions, breaks=50, main="Allelic fractions", xlim=c(0,1))
dev.off() # close the device
############################################

# GT       # Un-phased genotype calls
# GL       # Genotype log-likelihoods (natural log) for AA,AB and BB genotypes, where A = ref and B = variant. Only applicable for bi-allelic sites
# GOF      # Phred-scaled goodness-of-fit score for the genotype call
# GQ       # Phred-scaled quality score for the genotype call
# NR  = DP # Number of reads covering variant position in this sample
# NV  = AO (alt ) # Number of reads at variant position which support the called variant in this sample
