#!/usr/bin/env Rscript

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

generated_files = '/home/pgm/Workspace/MPM/R_controle_qualite/intersect_2_tumors/generated_files/'
print(paste("Generating chart for patient ", patient_id, "...", sep=""))

############################################
## Calculate allelic fractions for tumor 1 & 2
############################################
calculate_allelic_fractions = function(VCFcontent) {
  DPCounts <- unlist(lapply(1:nrow(VCFcontent), function(i) get_genotype(VCFcontent[i,21], VCFcontent[i, 19], "DP")))
  refCounts <- unlist(lapply(1:nrow(VCFcontent), function(i) get_genotype(VCFcontent[i,21], VCFcontent[i, 19], paste(VCFcontent[i,4] ,"U", sep=''), FALSE)))
  altCounts <- unlist(lapply(1:nrow(VCFcontent), function(i) get_genotype(VCFcontent[i,21], VCFcontent[i, 19], paste(VCFcontent[i,5] ,"U", sep=''), FALSE)))
  
  allelic_fractions = vector()
  for (i in 1:length(refCounts)) { # for each mutation
    ref_allele_count = as.numeric(unlist(strsplit(refCounts[i], ','))[1])
    alt_allele_count = as.numeric(unlist(strsplit(altCounts[i], ','))[1])
    allelic_fractions[i] = alt_allele_count / (alt_allele_count + ref_allele_count)
  }
  return(allelic_fractions)
}
tumor1_allfrac = calculate_allelic_fractions(tumor1_VCFcontent)
tumor2_allfrac = calculate_allelic_fractions(tumor2_VCFcontent)


############################################
## Generate unique ID and intersect
############################################
tumor1_VCFcontent$uniqid = paste(tumor1_VCFcontent$V1, tumor1_VCFcontent$V2, tumor1_VCFcontent$V4, tumor1_VCFcontent$V5, sep='_')
tumor2_VCFcontent$uniqid = paste(tumor2_VCFcontent$V1, tumor2_VCFcontent$V2, tumor2_VCFcontent$V4, tumor2_VCFcontent$V5, sep='_')
inter = intersect(tumor1_VCFcontent$uniqid, tumor2_VCFcontent$uniqid)
inter_t1 = tumor1_VCFcontent[tumor1_VCFcontent$uniqid %in% inter, ]
inter_t2 = tumor2_VCFcontent[tumor2_VCFcontent$uniqid %in% inter, ]

tumor1_inter_allfrac = calculate_allelic_fractions(inter_t1)
tumor2_inter_allfrac = calculate_allelic_fractions(inter_t2)


############################################
## Stacked bars with tumor 1 & 2
############################################
breaks = seq(0,1,0.01)
svg(filename=paste(generated_files, patient_id, "_", sample2_id, "_on_", sample1_id, ".svg", sep=''))

par(family="Times",las=1)

h1a = hist(tumor1_allfrac, breaks=breaks, plot=FALSE)
ylim1 = ceiling(max(h1a$counts)/50) * 50 # j'arrondis à la cinquantaine du dessus sinon c'est moche

h2a = hist(tumor2_allfrac, breaks=breaks, plot=FALSE)
ylim2 = ceiling(max(h2a$counts)/50) * 50 # j'arrondis à la cinquantaine du dessus sinon c'est moche
h2a$counts = -h2a$counts # on inverse les valeurs

# ylimv = max(c(ylim1, ylim2)) # je garde le max des max
ylimv = 300 # je force l'axe des Y pour avoir la même échelle pour tous
plot(h1a, ylim=c(-ylimv, ylimv), xlim=c(0,1), xlab="Allelic fractions", border=rgb(1,0,0,0.4), main=paste("Patient", patient_id, "(tumor VS tumor)"))
plot(h2a, ylim=c(-ylimv, ylimv), xlim=c(0,1), border=rgb(0,0,1,0.3), add=TRUE)

hist(tumor1_inter_allfrac, breaks=breaks, xlim=c(0,1), ylim=c(-ylimv, ylimv), col=rgb(0.5,0,0.5,0.6), add=TRUE)

h2b = hist(tumor2_inter_allfrac, breaks=breaks, plot=FALSE)
h2b$counts = -h2b$counts
plot(h2b, add=TRUE, ylim=c(-ylimv, ylimv), col=rgb(0.5,0,0.5,0.6))

legend("topright", c(sample1_id, sample2_id, "overlap"), col=c(rgb(1,0,0,0.4), rgb(0,0,1,0.3), rgb(0.5,0,0.5,0.6)), lwd=10)
box()
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


# pour mettre mes deux somatiques ensemble
# hist(rnorm(100))
# hist(rnorm(100),ylim=c(-20,20))
# h = hist(rnorm(100),add=T,plot=F)
# h$counts = -h$counts
# plot(h,add=T)
