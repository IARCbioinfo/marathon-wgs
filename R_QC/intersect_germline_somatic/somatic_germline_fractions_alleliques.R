#!/usr/bin/env Rscript

source("https://gist.githubusercontent.com/mfoll/a4dfbb92068dc559f130/raw/714dc8c2e97987fd4385dcef2722b3ef986d38d6/get_vcf_data.r")

args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) { stop("Input name missing!\n", call.=FALSE) }

path_to_normal_VCF       = args[1]
path_to_output_directory = args[2]
patient_id               = args[3]
sample1_id               = args[4]
sample2_id               = args[5]


print(paste("> processing ", patient_id, sep=''))
print("> loading germline VCF... ")
normal_VCFcontent = read.table(path_to_normal_VCF, sep="\t", header=FALSE, stringsAsFactors=FALSE)
print("DONE.")


############################################
## Calculate allelic fractions for tumor 1 & 2
############################################
calculate_germline_allelic_fractions = function(VCFcontent) {
  # 10 is normal
  DPCounts <- unlist(lapply(1:nrow(VCFcontent), function(i) get_genotype(VCFcontent[i,10], VCFcontent[i, 9], "NR")))
  altCounts <- unlist(lapply(1:nrow(VCFcontent), function(i) get_genotype(VCFcontent[i,10], VCFcontent[i, 9], "NV")))
  normal = altCounts / DPCounts

  # 11 is tumor1
  DPCounts <- unlist(lapply(1:nrow(VCFcontent), function(i) get_genotype(VCFcontent[i,11], VCFcontent[i, 9], "NR")))
  altCounts <- unlist(lapply(1:nrow(VCFcontent), function(i) get_genotype(VCFcontent[i,11], VCFcontent[i, 9], "NV")))
  t1 = altCounts / DPCounts

  # 12 is tumor2
  DPCounts <- unlist(lapply(1:nrow(VCFcontent), function(i) get_genotype(VCFcontent[i,12], VCFcontent[i, 9], "NR")))
  altCounts <- unlist(lapply(1:nrow(VCFcontent), function(i) get_genotype(VCFcontent[i,12], VCFcontent[i, 9], "NV")))
  t2 = altCounts / DPCounts

  return(list(normal, t1, t2))
}
print("> calculating allelic fractions... ")
allFracs = calculate_germline_allelic_fractions(normal_VCFcontent)
allFracs_normal = allFracs[[1]]
allFracs_tumor1 = allFracs[[2]]
allFracs_tumor2 = allFracs[[3]]
print("DONE.")


############################################
## Save computed data
############################################
save.image(file=paste(path_to_output_directory, "/", patient_id, ".germline_VS_tumor_in_germline.rda", sep=''))
# load(file=paste(generated_files, patient_id, "germline_VS_tumor_in_germline.rda", sep=''))


############################################
## Stacked bars with normal & tumor
############################################
print("> generating chart...")
breaks = seq(0,1,0.02)
svg(filename=paste(path_to_output_directory, "/", patient_id, ".germline_VS_tumors.svg", sep=''))

par(family="Times",las=1)

h1 = hist(allFracs_normal, breaks=breaks, plot=FALSE)
ylim1 = ceiling(max(h1$counts)/50) * 50
plot(h1, ylim=c(-ylim1, ylim1), xlab="Allelic fractions", col=rgb(0,0,0,0.4), border=rgb(1,1,1,1), main=paste("Patient", patient_id, "(germline VS tumors)"))

h2 = hist(allFracs_tumor1, breaks=breaks, plot=FALSE)
h2$counts = -h2$counts # on inverse les valeurs
plot(h2, ylim=c(-ylimv, ylimv), xlim=c(0,1), col=rgb(1,0,0,0.4), border=rgb(1,1,1,1), add=TRUE)

h3 = hist(allFracs_tumor2, breaks=breaks, plot=FALSE)
h3$counts = -h3$counts # on inverse les valeurs
plot(h3, ylim=c(-ylimv, ylimv), xlim=c(0,1), col=rgb(0,0,1,0.3), border=rgb(1,1,1,1), add=TRUE)

legend("topleft", c("germline", sample1_id, sample2_id, "overlap"), col=c(rgb(0,0,0,0.4), rgb(1,0,0,0.4), rgb(0,0,1,0.3), rgb(0.5,0,0.5,0.6)), lwd=10)

dev.off() # close the device
print("DONE.")
############################################
