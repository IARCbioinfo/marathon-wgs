#!/usr/bin/env Rscript

library("Canopy")
source("https://gist.githubusercontent.com/mfoll/a4dfbb92068dc559f130/raw/714dc8c2e97987fd4385dcef2722b3ef986d38d6/get_vcf_data.r")


##########################################
## Retrieve arguments
##########################################
args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) { stop("Input name missing!\n", call.=FALSE) }

##
input_folder                      = args[1]
patient_id                        = args[2]
tumor1_id                         = args[3]
tumor2_id                         = args[4]
input_somatic_VCF_t1              = args[5]
input_somatic_VCF_t2              = args[6]
input_somatic_VCF_t1_positions_t2 = args[7]
input_somatic_VCF_t2_positions_t1 = args[8]
output_path                       = args[9]
K                                 = args[10]
param_writeskip                   = as.numeric(args[11])
param_thin                        = as.numeric(args[12])

##
# input_folder                      = "/home/pgm/Workspace/MPM/VCF_finaux/falcon/patient_5009/"
# patient_id                        = "5009"
# tumor1_id                         = "B00JAJB"
# tumor2_id                         = "B00JAJC"
# input_somatic_VCF_t1              = "/home/pgm/Workspace/MPM/VCF_finaux/somatic_sandbox/M662_DA_5009_T_B00JAJB.normalized.vcf_multianno.hg38_multianno.txt"
# input_somatic_VCF_t2              = "/home/pgm/Workspace/MPM/VCF_finaux/somatic_sandbox/M662_DA_5009_T_B00JAJC.normalized.vcf_multianno.hg38_multianno.txt"
# input_somatic_VCF_t1_positions_t2 = "/home/pgm/Workspace/MPM/scripts_colbalt/calling_somatic_genotype/output/M662_DA_5009_T_B00JAJB.other_tumor_positions.vcf"
# input_somatic_VCF_t2_positions_t1 = "/home/pgm/Workspace/MPM/scripts_colbalt/calling_somatic_genotype/output/M662_DA_5009_T_B00JAJC.other_tumor_positions.vcf"
# output_path                       = paste("/home/pgm/Workspace/MPM/marathon/canopy/generated_clustering/5009/", patient_id, ".K", 3, sep='')
# K                                 = 2

##
cat("####### ARGUMENTS #######\n")
cat(paste("input_folder: ", input_folder, "\n", sep=''))
cat(paste("patient_id: ", patient_id, "\n", sep=''))
cat(paste("tumor1_id: ", tumor1_id, "\n", sep=''))
cat(paste("tumor2_id: ", tumor2_id, "\n", sep=''))
cat(paste("K (nb subclones): ", K, "\n", sep=''))


##########################################
## Debug function
##########################################
# write_matrices = function() {
#   write.table(C, quote = F, file = "/home/pgm/Workspace/MPM/marathon/canopy/5009_C.tsv")
#   write.table(Y, quote = F, file = "/home/pgm/Workspace/MPM/marathon/canopy/5009_Y.tsv")
#   write.table(R, quote = F, file = "/home/pgm/Workspace/MPM/marathon/canopy/5009_R.tsv")
#   write.table(X, quote = F, file = "/home/pgm/Workspace/MPM/marathon/canopy/5009_X.tsv")
#   write.table(WM, quote = F, file = "/home/pgm/Workspace/MPM/marathon/canopy/5009_WM.tsv")
#   write.table(Wm, quote = F, file = "/home/pgm/Workspace/MPM/marathon/canopy/5009_Wm.tsv")
#   write.table(epsM, quote = F, file = "/home/pgm/Workspace/MPM/marathon/canopy/5009_epsM.tsv")
#   write.table(epsm, quote = F, file = "/home/pgm/Workspace/MPM/marathon/canopy/5009_epsm.tsv")
# }


##########################################
## Prepare CNAs input for canopy
##########################################
compile_ascn_tumor = function(patient_id, tumor_id) {
  chromosomes = c(seq(1:22), "X", "Y")
  header = c("chr", "st_snp", "end_snp", "st_bp", "end_bp", "Minor_copy", "Major_copy", "Minor.sd", "Major.sd")
  ascn = data.frame(matrix(, nrow=0, ncol=length(header)))
  colnames(ascn) = header
  for (chr in chromosomes) {
    input_folder_chr = paste(input_folder, "chr", chr, "/", sep="")
    input_file_name = paste("falcon.patient_", patient_id, ".tumor_", tumor_id, ".chr_", chr, ".output.txt", sep="") # falcon.patient_5009.tumor_B00JAJB.chr_Y.output.txt
    ascn_chr = data.frame(read.table(paste(input_folder_chr, input_file_name, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE))
    ascn = rbind(ascn, ascn_chr)
  }
  ## Remove rows with Minor.sd = Major.sd = NA
  ascn = ascn[(!is.na(ascn$Minor.sd) & !is.na(ascn$Major.sd)),]
  return (ascn)
}

ascn_primary = compile_ascn_tumor(patient_id, tumor1_id)
ascn_relapse = compile_ascn_tumor(patient_id, tumor2_id)


##########################################
## Build C matrix (CNA over regions)
## See https://groups.google.com/forum/#!searchin/canopy_phylogeny/overlap|sort:date/canopy_phylogeny/lT2S6t4QNT4/2XjzeCUcBQAJ
##########################################
ascn_total = rbind(ascn_primary, ascn_relapse)
ascn_total$chr[ascn_total$chr == "X"] = 23 # to convert chr as numeric and sort it
ascn_total$chr[ascn_total$chr == "Y"] = 24
ascn_total$chr = as.numeric(ascn_total$chr)
ascn_total = ascn_total[order( ascn_total[,1], ascn_total[,2] ), ] # sort by chromosome
rownames(ascn_total) = seq(1, nrow(ascn_total))
chromosomes_impacted = unique(ascn_total$chr)
all_regions = c()
for (chr in chromosomes_impacted) {
  sub_ascn = ascn_total[ascn_total$chr == chr, ] # all ascn of a chromosome
  starts   = sub_ascn$st_bp
  ends     = sub_ascn$end_bp
  points   = sort(unique(c(starts, ends)))
  new_starts = points[1:(length(points)-1)] # all points except last one are starts
  new_ends   = points[2:(length(points))] # all points except first one are ends
  chr_regions = paste("chr", chr, ":", new_starts, "-", new_ends, sep="")
  all_regions = c(all_regions, chr_regions)
}

# prepare empty C matrix
C = as.data.frame(matrix(0, nrow=length(all_regions), ncol=nrow(ascn_total)))
rownames(C) = all_regions
colnames(C) = paste("CNA", seq(1:nrow(ascn_total)), sep="_")

# process all copy numbers and fill C matrix with 1 or 0 for each copy number
for (i in seq(1:nrow(ascn_total))) {
  copy_number = ascn_total[i, ]
  
  for (j in seq(1:nrow(C))) {
    C_chr   = unlist(strsplit(rownames(C[j, ]), ":"))[1]
    C_start = as.numeric(unlist(strsplit(unlist(strsplit(rownames(C[j, ]), ":"))[2], "-"))[1])
    C_end   = as.numeric(unlist(strsplit(unlist(strsplit(rownames(C[j, ]), ":"))[2], "-"))[2])
    
    if (paste("chr", copy_number$chr, sep="") == C_chr) {
      if (C_start >= as.numeric(copy_number$st_bp) && C_start < as.numeric(copy_number$end_bp)) {
        C[j, i] = 1
      }
    }
  }
}
C = C[rowSums(C) > 0, ] # Remove regions with no CNA event


##########################################
## Build W & eps matrices
##########################################
WM   = as.data.frame(matrix(NA, nrow=nrow(C), ncol=2))
Wm   = as.data.frame(matrix(NA, nrow=nrow(C), ncol=2))
epsM = as.data.frame(matrix(NA, nrow=nrow(C), ncol=2))
epsm = as.data.frame(matrix(NA, nrow=nrow(C), ncol=2))
rownames(WM) = rownames(Wm) = rownames(epsM) = rownames(epsm) = rownames(C)
colnames(WM) = colnames(Wm) = colnames(epsM) = colnames(epsm) = c(tumor1_id, tumor2_id)

build_matrices = function(ascn_tumor, ascn_tumor_index, WM, Wm, epsM, epsm) {
  for (i in seq(1:nrow(WM))) {
    C_chr   = unlist(strsplit(rownames(WM[i, ]), ":"))[1]
    C_start = as.numeric(unlist(strsplit(unlist(strsplit(rownames(WM[i, ]), ":"))[2], "-"))[1])
    C_end   = as.numeric(unlist(strsplit(unlist(strsplit(rownames(WM[i, ]), ":"))[2], "-"))[2])
    
    for (j in seq(1:nrow(ascn_tumor))) {
      copy_number = ascn_tumor[j, ]
      copy_number$chr = gsub("X", "23", copy_number$chr)
      copy_number$chr = gsub("Y", "24", copy_number$chr)
      if (paste("chr", copy_number$chr, sep="") == C_chr && C_start >= as.numeric(copy_number$st_bp) && C_start < as.numeric(copy_number$end_bp) && is.na(WM[i, ascn_tumor_index])) {
        WM[i, ascn_tumor_index] = copy_number$Major_copy
        Wm[i, ascn_tumor_index] = copy_number$Minor_copy
        epsM[i, ascn_tumor_index] = copy_number$Major.sd
        epsm[i, ascn_tumor_index] = copy_number$Minor.sd
      }
    }
  }
  return(list(WM, Wm, epsM, epsm))
}

matrices = build_matrices(ascn_primary, 1, WM, Wm, epsM, epsm)
WM = matrices[[1]]
Wm = matrices[[2]]
epsM = matrices[[3]]
epsm = matrices[[4]]
matrices = build_matrices(ascn_relapse, 2, WM, Wm, epsM, epsm)
WM = matrices[[1]]
Wm = matrices[[2]]
epsM = matrices[[3]]
epsm = matrices[[4]]

# Complete missing W & eps with other tumor data recomputed
WM[is.na(WM)] = 1 # not found copy numbers are 1
Wm[is.na(Wm)] = 1

for (i in 1:nrow(epsM)) {
  
  ascn_coord = rownames(epsM[i,])[1]
  ascn_coord_spl = unlist(strsplit(ascn_coord, ':'))
  ascn_coord_chr = sub("chr", "", ascn_coord_spl[1])
  ascn_coord_spl = unlist(strsplit(ascn_coord_spl[2], '-'))
  ascn_coord_start_bp = ascn_coord_spl[1]
  ascn_coord_end_bp = ascn_coord_spl[2]
  
  if (ascn_coord_chr == "23") { ascn_coord_chr = "X" }
  if (ascn_coord_chr == "24") { ascn_coord_chr = "Y" }
  
  if (is.na(epsM[i,1])) { # primary is NA
    input_folder_chr = paste(input_folder, "chr", ascn_coord_chr, "/", sep="")
    ascn_complementary_file = paste("falcon.patient_", patient_id, ".tumor_", tumor1_id, ".chr_", ascn_coord_chr, ".output_epsilon.txt", sep="")
    ascn_complementary_chr = data.frame(read.table(paste(input_folder_chr, ascn_complementary_file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE))
    ascn_complementary_row = ascn_complementary_chr[which(ascn_complementary_chr$chr == ascn_coord_chr & ascn_complementary_chr$st_bp == ascn_coord_start_bp),]
    if (nrow(ascn_complementary_row) == 0) { ascn_complementary_row = ascn_complementary_chr[which(ascn_complementary_chr$chr == ascn_coord_chr & ascn_coord_start_bp >= ascn_complementary_chr$st_bp & ascn_coord_start_bp <= ascn_complementary_chr$end_bp),] }
    if (nrow(ascn_complementary_row) == 1) {
      epsM[i,1] = ascn_complementary_row$Major.sd
      epsm[i,1] = ascn_complementary_row$Minor.sd
    }
  }
  
  # this code duplication is a shame...
  if (is.na(epsM[i,2])) { # relapse is NA
    input_folder_chr = paste(input_folder, "chr", ascn_coord_chr, "/", sep="")
    ascn_complementary_file = paste("falcon.patient_", patient_id, ".tumor_", tumor2_id, ".chr_", ascn_coord_chr, ".output_epsilon.txt", sep="")
    ascn_complementary_chr = data.frame(read.table(paste(input_folder_chr, ascn_complementary_file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE))
    ascn_complementary_row = ascn_complementary_chr[which(ascn_complementary_chr$chr == ascn_coord_chr & ascn_complementary_chr$st_bp == ascn_coord_start_bp),]
    if (nrow(ascn_complementary_row) == 0) { ascn_complementary_row = ascn_complementary_chr[which(ascn_complementary_chr$chr == ascn_coord_chr & ascn_coord_start_bp >= ascn_complementary_chr$st_bp & ascn_coord_start_bp <= ascn_complementary_chr$end_bp),] }
    if (nrow(ascn_complementary_row) == 1) {
      epsM[i,2] = ascn_complementary_row$Major.sd
      epsm[i,2] = ascn_complementary_row$Minor.sd
    }
  }
}

# Just for security, but normally no NA at this step
epsM[is.na(epsM)] = 0.001 # not found standard deviations are 0.001
epsm[is.na(epsm)] = 0.001


##########################################
## SNA input : R & X for VAF
##########################################
VCF_somatic_tumor1 = read.table(input_somatic_VCF_t1, sep = "\t", header=FALSE, stringsAsFactors=FALSE)
VCF_somatic_tumor1$uniqID <- do.call(paste, c(VCF_somatic_tumor1[c("V1", "V2", "V3", "V4", "V5")], sep = "_")) 
VCF_somatic_tumor2 = read.table(input_somatic_VCF_t2, sep = "\t", header=FALSE, stringsAsFactors=FALSE)
VCF_somatic_tumor2$uniqID <- do.call(paste, c(VCF_somatic_tumor2[c("V1", "V2", "V3", "V4", "V5")], sep = "_")) 

# Retrieve all DP & AO counts for tumor1
DPCounts_tumor1 = unlist(lapply(1:nrow(VCF_somatic_tumor1), function(i) get_genotype(VCF_somatic_tumor1[i,21], VCF_somatic_tumor1[i, 19], "DP")))
AOCounts_tumor1 = unlist(lapply(1:nrow(VCF_somatic_tumor1), function(i) get_genotype(VCF_somatic_tumor1[i,21], VCF_somatic_tumor1[i, 19], paste(VCF_somatic_tumor1[i,5] ,"U", sep=''), FALSE)))
AOCounts_tumor1_spl = strsplit(AOCounts_tumor1, ',')
AOCounts_tumor1 = as.numeric(unlist(AOCounts_tumor1_spl)[2*(1:length(AOCounts_tumor1))-1])

# Retrieve all DP & AO counts for tumor2
DPCounts_tumor2 = unlist(lapply(1:nrow(VCF_somatic_tumor2), function(i) get_genotype(VCF_somatic_tumor2[i,21], VCF_somatic_tumor2[i, 19], "DP")))
AOCounts_tumor2 = unlist(lapply(1:nrow(VCF_somatic_tumor2), function(i) get_genotype(VCF_somatic_tumor2[i,21], VCF_somatic_tumor2[i, 19], paste(VCF_somatic_tumor2[i,5] ,"U", sep=''), FALSE)))
AOCounts_tumor2_spl = strsplit(AOCounts_tumor2, ',')
AOCounts_tumor2 = as.numeric(unlist(AOCounts_tumor2_spl)[2*(1:length(AOCounts_tumor2))-1])

# Load tumor1 counts with positions of tumor2
VCF_somatic_tumor1_positions_t2 = as.data.frame(read.table(input_somatic_VCF_t1_positions_t2, sep = "\t", header=FALSE, stringsAsFactors=FALSE))
VCF_somatic_tumor1_positions_t2_uniqID <- do.call(paste, c(VCF_somatic_tumor1_positions_t2[c("V1", "V2", "V2", "V4", "V5")], sep = "_")) 
VCF_somatic_tumor1_positions_t2 = cbind(VCF_somatic_tumor1_positions_t2, VCF_somatic_tumor1_positions_t2_uniqID)

# Retrieve AO and DP counts in tumor1, with tumor2 positions
DPCounts_positions2_in_tumor1 = unlist(lapply(1:nrow(VCF_somatic_tumor1_positions_t2), function(i) get_genotype(VCF_somatic_tumor1_positions_t2[i,10], VCF_somatic_tumor1_positions_t2[i, 9], "NR")))
AOCounts_positions2_in_tumor1 = unlist(lapply(1:nrow(VCF_somatic_tumor1_positions_t2), function(i) get_genotype(VCF_somatic_tumor1_positions_t2[i,10], VCF_somatic_tumor1_positions_t2[i, 9], "NV")))

# Load tumor2 counts with positions of tumor1
VCF_somatic_tumor2_positions_t1 = as.data.frame(read.table(input_somatic_VCF_t2_positions_t1, sep = "\t", header=FALSE, stringsAsFactors=FALSE))
VCF_somatic_tumor2_positions_t1_uniqID <- do.call(paste, c(VCF_somatic_tumor2_positions_t1[c("V1", "V2", "V2", "V4", "V5")], sep = "_")) 
VCF_somatic_tumor2_positions_t1 = cbind(VCF_somatic_tumor2_positions_t1, VCF_somatic_tumor2_positions_t1_uniqID)

# Retrieve AO and DP counts in tumor2, with tumor1 positions
DPCounts_positions1_in_tumor2 = unlist(lapply(1:nrow(VCF_somatic_tumor2_positions_t1), function(i) get_genotype(VCF_somatic_tumor2_positions_t1[i,10], VCF_somatic_tumor2_positions_t1[i, 9], "NR")))
AOCounts_positions1_in_tumor2 = unlist(lapply(1:nrow(VCF_somatic_tumor2_positions_t1), function(i) get_genotype(VCF_somatic_tumor2_positions_t1[i,10], VCF_somatic_tumor2_positions_t1[i, 9], "NV")))

# Generate R & X inputs
Rt1 = as.data.frame(cbind(AOCounts_tumor1, AOCounts_positions1_in_tumor2, VCF_somatic_tumor1$uniqID))
rownames(Rt1) = VCF_somatic_tumor1$uniqID
Rt2 = as.data.frame(cbind(AOCounts_positions2_in_tumor1, AOCounts_tumor2, VCF_somatic_tumor2$uniqID))
rownames(Rt2) = VCF_somatic_tumor2$uniqID
Rt2_dedup = Rt2[!(Rt2$V3 %in% Rt1$V3),] # keep only tumor2 rows which are not in t1
colnames(Rt1) = colnames(Rt2_dedup) = c(tumor1_id, tumor2_id)
R = rbind(Rt1, Rt2_dedup)
R = R[,1:ncol(R)-1]
for (i in 1:ncol(R)) {
  R[, i] = as.numeric(as.character(R[, i]))
}

Xt1 = as.data.frame(cbind(DPCounts_tumor1, DPCounts_positions1_in_tumor2, VCF_somatic_tumor1$uniqID))
rownames(Xt1) = VCF_somatic_tumor1$uniqID
Xt2 = as.data.frame(cbind(DPCounts_positions2_in_tumor1, DPCounts_tumor2, VCF_somatic_tumor2$uniqID))
rownames(Xt2) = VCF_somatic_tumor2$uniqID
Xt2_dedup = Xt2[!(Xt2$V3 %in% Xt1$V3),]
colnames(Xt1) = colnames(Xt2_dedup) = c(tumor1_id, tumor2_id)
X = rbind(Xt1, Xt2_dedup)
X = X[,1:ncol(X)-1]
for (i in 1:ncol(X)) {
  X[, i] = as.numeric(as.character(X[, i]))
}

# Generate Y input
Y_column_names = rownames(C)
Y_row_names    = rownames(X)
Y = as.data.frame(matrix(0, nrow=length(Y_row_names), ncol=length(Y_column_names)))
rownames(Y) = Y_row_names
colnames(Y) = Y_column_names

# process all positions and fill Y matrix with 1 when this position is in a CNA region
for (i in seq(1:nrow(Y))) {
  Ysna_chr      = unlist(strsplit(rownames(Y)[i], '_'))[1]
  if (Ysna_chr == "chrX") { Ysna_chr = "chr23" }
  if (Ysna_chr == "chrY") { Ysna_chr = "chr24" }
  Ysna_position = as.numeric(unlist(strsplit(rownames(Y)[i], '_'))[2])
  
  for (j in seq(1:ncol(Y))) {
    Ycna = unlist(strsplit(colnames(Y)[j], ':'))
    Ycna_chr   = Ycna[1]
    Ycna_start = as.numeric(unlist(strsplit(Ycna[2], "-"))[1])
    Ycna_end   = as.numeric(unlist(strsplit(Ycna[2], "-"))[2])
    
    if (Ysna_chr == Ycna_chr && Ysna_position >= Ycna_start && Ysna_position <= Ycna_end) {
      Y[i, j] = 1 # [row, col]
    }
  }
}

# Finally update non cna regions
Y = cbind(rep(0, nrow(Y)), Y)
colnames(Y)[1] = "non-cna_region"
Y[rowSums(Y) == 0, 1] = 1




##########################################
##########################################
##                                      ##
##          Process Canopy              ##
##                                      ##
##########################################
##########################################



##########################################
## Binomial pre-clustering of SNAs
##########################################

# A multivariate binomial mixture clustering step can be applied to the SNAs before MCMC 
# sampling. We show in our paper via simulations that this pre-clustering method helps 
# the Markov chain converge faster with smaller estimation error (especially when mutations 
# show clear cluster patterns by visualization). This clustering step can also remove likely 
# false positives before feeding the mutations to the MCMC algorithm.

if(K==2) Kclusters = 2:3 # Range of number of clusters to run
if(K==3) Kclusters = 2:5 # Range of number of clusters to run
if(K==4) Kclusters = 2:7 # Range of number of clusters to run
if(K>=5) Kclusters = 2:10 # Range of number of clusters to run
num_run = 10 # How many EM runs per clustering step for each mutation cluster wave
canopy.cluster = canopy.cluster(R = data.matrix(R),
                                X = data.matrix(X),
                                num_cluster = Kclusters,
                                num_run = num_run)

bic_output  = canopy.cluster$bic_output # BIC for model selection (# of clusters)
Mu          = canopy.cluster$Mu # VAF centroid for each cluster
Tau         = canopy.cluster$Tau  # Prior for mutation cluster, with a K+1 component
sna_cluster = canopy.cluster$sna_cluster # cluster identity for each mutation

# Plot BIC values
BIC = data.frame(as.numeric(bic_output))
BIC = cbind(Kclusters, BIC)
rownames(BIC) = Kclusters
colnames(BIC) = c("Kclust", "BIC")
optimalBIC = BIC[which.max(BIC$BIC), ]

svg(filename=paste(output_path, ".optimal_cluster_number.svg", sep=''))
plot(BIC, xlab="Kclust (nb clusters)", xaxt="n", main="Optimal number of SNAs clusters")
lines(BIC)
abline(v = optimalBIC$Kclust, col="red")
axis(side=1, at=Kclusters, labels=Kclusters)
legend("bottomleft", c(paste("optimal Kclust = ", optimalBIC$Kclust, sep='')), col=c(rgb(1,0,0,1)), lty=1:2, cex=0.8, box.lty=0)
dev.off() # close the device


##########################################
## MCMC sampling
##########################################

# Canopy samples in subtree space with varying number of subclones(denoted as K) by a 
# Markov chain Monte Carlo (MCMC) method. A plot of posterior likelihood (pdf format) 
# will be generated for each subtree space and we recommend users to refer to the plot 
# as a sanity check for sampling convergence and to choose the number of burn-ins and 
# thinning accordingly. Note that this step can be time-consuming, especially with 
# larger number of chains numchain specifies the number of chains with random initiations, 
# a larger value of which is in favor of not getting stuck in local optima) and longer 
# chains (simrun specifies number of iterations per chain). MCMC sampling is the most 
# computationally heavy step in Canopy. It is recommended that jobs are run in parallel 
# on high-performance cluster.


# This function is for cases where SNAs are pre-clustered by the Binomial mixture EM algorithm
numchain     = 10 # number of chains with random initiations
max.simrun   = 50000 # to increase later
min.simrun   = 10000 # to increase later
writeskip    = param_writeskip
K = K:K

source("/data/soudadel/MPM/canopy/scripts/custom_canopy.sample.cluster.R")
# source("/home/pgm/Workspace/MPM/marathon/libs")
sampchain = custom_canopy.sample.cluster(as.matrix(R), as.matrix(X), 
                                         sna_cluster, 
                                         as.matrix(WM), as.matrix(Wm), 
                                         as.matrix(epsM), as.matrix(epsm), 
                                         C=as.matrix(C),
                                         Y=as.matrix(Y), 
                                         K, 
                                         max.simrun = max.simrun, min.simrun = min.simrun, numchain = numchain,
                                         writeskip = writeskip, 
                                         projectname = output_path, 
                                         cell.line = TRUE, 
                                         plot.likelihood = TRUE)
save.image(file = paste(output_path, '.postmcmc_image.rda', sep=''), compress = 'xz')


##########################################
## BIC for model selection
##########################################
# load(file = "/home/pgm/Workspace/MPM/marathon/canopy/generated_clustering/5009/5009.K5.postmcmc_image.rda")
# Canopy uses BIC as a model selection criterion to determine to optimal number of subclones.

# burnin: initial fraction of the MCMC samples to be discarded as burn-in. Must be a value in [0, 1).
# thin positive integer. If thin = n, only every n-th realizations of the Markov chain is kept.

burnin = 90
thin = param_thin # If there is error in the bic and canopy.post step below, make sure
# burnin and thinning parameters are wisely selected so that there are
# posterior trees left.

bic = canopy.BIC(sampchain = sampchain, projectname = output_path, K = K,
                 numchain = numchain, burnin = burnin, thin = thin, pdf = TRUE)
write.table(bic, file = paste(output_path, '.BIC', sep=""))




