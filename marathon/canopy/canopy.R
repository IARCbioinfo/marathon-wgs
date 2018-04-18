#!/usr/bin/env Rscript

library("Canopy")
source("https://gist.githubusercontent.com/mfoll/a4dfbb92068dc559f130/raw/714dc8c2e97987fd4385dcef2722b3ef986d38d6/get_vcf_data.r")


load('/home/pgm/Workspace/MPM/marathon/canopy/demo/preprocessed.rda')


##########################################
## Retrieve arguments
##########################################
args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) { stop("Input name missing!\n", call.=FALSE) }

##
input_folder          = args[1]
patient_id            = args[2]
tumor1_id             = args[3]
tumor2_id             = args[4]
input_somatic_VCF_t1  = args[5]
input_somatic_VCF_t2  = args[6]
input_mpileup_t1      = args[7]
input_mpileup_t2      = args[8]

##
input_folder          = "/home/pgm/Workspace/MPM/VCF_finaux/falcon/patient_5009/"
patient_id            = "5009"
tumor1_id             = "B00JAJB"
tumor2_id             = "B00JAJC"
input_somatic_VCF_t1  = "/home/pgm/Workspace/MPM/VCF_finaux/somatic_sandbox/M662_DA_5009_T_B00JAJB.normalized.vcf_multianno.hg38_multianno.txt"
input_somatic_VCF_t2  = "/home/pgm/Workspace/MPM/VCF_finaux/somatic_sandbox/M662_DA_5009_T_B00JAJC.normalized.vcf_multianno.hg38_multianno.txt"
input_mpileup_t1      = "/home/pgm/Workspace/MPM/VCF_finaux/other_tumor_coverage/M662_DA_5009_T_B00JAJB.mpileup.tsv"
input_mpileup_t2      = "/home/pgm/Workspace/MPM/VCF_finaux/other_tumor_coverage/M662_DA_5009_T_B00JAJC.mpileup.tsv"

##
cat("####### ARGUMENTS #######\n")
cat(paste("input_folder: ", input_folder, "\n", sep=''))
cat(paste("patient_id: ", patient_id, "\n", sep=''))
cat(paste("tumor1_id: ", tumor1_id, "\n", sep=''))
cat(paste("tumor2_id: ", tumor2_id, "\n", sep=''))


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
  starts   = sub_ascn$st_snp
  ends     = sub_ascn$end_snp
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
      if (C_start >= as.numeric(copy_number$st_snp) && C_start < as.numeric(copy_number$end_snp)) {
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
      if (paste("chr", copy_number$chr, sep="") == C_chr && C_start >= as.numeric(copy_number$st_snp) && C_start < as.numeric(copy_number$end_snp) && is.na(WM[i, ascn_tumor_index])) {
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

# All NAs copy numbers are 1 & standard deviations are 0
WM[is.na(WM)] = 1
Wm[is.na(Wm)] = 1
epsM[is.na(epsM)] = 0
epsm[is.na(epsm)] = 0


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

# Init R & X with all tumor1 position, and complete tumor2 counts
R = data.frame(cbind(AOCounts_tumor1, rep(NA, length(AOCounts_tumor1))))
X = data.frame(cbind(DPCounts_tumor1, rep(NA, length(DPCounts_tumor1))))
rownames(R) = rownames(X) = VCF_somatic_tumor1$uniqID
mpileup_tumor2 = read.table(input_mpileup_t2, sep = "\t", header=FALSE, stringsAsFactors=FALSE)
mpileup_tumor2_uniqID <- do.call(paste, c(mpileup_tumor2[c("V1", "V2", "V3", "V4", "V5")], sep = "_")) 
mpileup_tumor2_AO = data.frame(mpileup_tumor2$V6)
mpileup_tumor2_DP = data.frame(mpileup_tumor2$V7)
rownames(mpileup_tumor2_AO) = rownames(mpileup_tumor2_DP) = mpileup_tumor2_uniqID
for (i in seq(1:nrow(mpileup_tumor2_AO))) { # for each mpileup, put count with its tumor1 equivalent
  position_name = rownames(mpileup_tumor2_AO)[i]
  R[position_name, 2] = mpileup_tumor2_AO[i, 1]
  X[position_name, 2] = mpileup_tumor2_DP[i, 1]
}

# Keep only tumor2 positions which are not already in R (= all except the intersected ones)
AOCounts_tumor2 = as.data.frame(AOCounts_tumor2) # all AO counts in tumor 2
DPCounts_tumor2 = as.data.frame(DPCounts_tumor2) # all DP counts in tumor 2
rownames(AOCounts_tumor2) = rownames(DPCounts_tumor2) = VCF_somatic_tumor2$uniqID
colnames(DPCounts_tumor2) = colnames(DPCounts_tumor2) = c("t2")
AOCounts_tumor2_not_in_t1 = subset(AOCounts_tumor2, ! rownames(AOCounts_tumor2) %in% rownames(R)) # all AO counts in tumor2 and not in tumor1
DPCounts_tumor2_not_in_t1 = subset(DPCounts_tumor2, ! rownames(DPCounts_tumor2) %in% rownames(X)) # all DP counts in tumor2 and not in tumor1
Rt2 = data.frame(cbind(rep(NA, length(AOCounts_tumor2_not_in_t1)), AOCounts_tumor2_not_in_t1))
Xt2 = data.frame(cbind(rep(NA, length(DPCounts_tumor2_not_in_t1)), DPCounts_tumor2_not_in_t1))
mpileup_tumor1 = read.table(input_mpileup_t1, sep = "\t", header=FALSE, stringsAsFactors=FALSE)
mpileup_tumor1_uniqID <- do.call(paste, c(mpileup_tumor1[c("V1", "V2", "V3", "V4", "V5")], sep = "_")) 
mpileup_tumor1_AO = data.frame(mpileup_tumor1$V6)
##
mpileup_tumor1_AO = subset(AOCounts_tumor2, ! rownames(AOCounts_tumor2) %in% rownames(R))
##
mpileup_tumor1_DP = data.frame(mpileup_tumor1$V7)
rownames(mpileup_tumor1_AO) = mpileup_tumor1_uniqID
for (i in seq(1:nrow(mpileup_tumor1_AO))) { # for each mpileup, put count with its tumor1 equivalent
  position_name = rownames(mpileup_tumor1_AO)[i]
  Rt2[position_name, 1] = mpileup_tumor1_AO[i, 1]
  Xt2[position_name, 1] = mpileup_tumor1_DP[i, 1]
}

colnames(R) = colnames(Rt2) = colnames(X) = colnames(Xt2) = c(tumor1_id, tumor2_id)
R = rbind(R, Rt2)
X = rbind(X, Xt2)




## Remove NAs
# R_t1 = na.omit(SNA_in$R) # nb reads supporting mutant allele
# X_t1 = na.omit(SNA_in$X) # nb reads total


##########################################
## Build Y matrix (SNA over regions)
##########################################









# ##########################################
# ## SNA input : R & X for VAF
# ##########################################
# VCF_somatic_content_t1 = read.table(input_somatic_VCF_t1, sep = "\t", header=FALSE, stringsAsFactors=FALSE)
# VCF_somatic_header = c("Chr", "Start", "End", "Ref", "Alt" ,"Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S00716_N", "S00716_T1")
# colnames(VCF_somatic_content_t1) = VCF_somatic_header
# 
# VCF_somatic_content_t1_extract_sna =  VCF_somatic_content_t1[tail(seq_along(VCF_somatic_content_t1), 11)]
# 
# source("/home/pgm/Workspace/MPM/marathon/libs/readVCFforCanopy.R")
# SNA_t1 = readVCFforCanopy(VCF_somatic_content_t1_extract_sna)


## Remove NAs
# R_t1 = na.omit(SNA_in$R) # nb reads supporting mutant allele
# X_t1 = na.omit(SNA_in$X) # nb reads total


# data("MDA231")
# projectname = MDA231$projectname ## name of project
# R = MDA231$R; R ## mutant allele read depth (for SNAs)
# X = MDA231$X; X ## total depth (for SNAs)
# WM = MDA231$WM; WM ## observed major copy number (for CNA regions)
# Wm = MDA231$Wm; Wm ## observed minor copy number (for CNA regions)
# epsilonM = MDA231$epsilonM ## standard deviation of WM, pre-fixed here
# epsilonm = MDA231$epsilonm ## standard deviation of Wm, pre-fixed here
# ## whether CNA regions harbor specific CNAs (only needed for overlapping CNAs)
# C = MDA231$C; C
# Y = MDA231$Y; Y ## whether SNAs are affected by CNAs



