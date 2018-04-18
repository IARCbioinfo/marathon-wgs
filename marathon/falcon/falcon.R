#!/usr/bin/env Rscript

library("falcon")
source("https://gist.githubusercontent.com/mfoll/a4dfbb92068dc559f130/raw/714dc8c2e97987fd4385dcef2722b3ef986d38d6/get_vcf_data.r")
# source("/home/pgm/Workspace/MPM/marathon/libs/falcon.output.R")
# source("/home/pgm/Workspace/MPM/marathon/libs/falcon.qc.R")
source("/home/soudadel/MPM/falcon/libs/falcon.output.R")
source("/home/soudadel/MPM/falcon/libs/falcon.qc.R")


##########################################
## Retrieve arguments
##########################################
args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) { stop("Input name missing!\n", call.=FALSE) }

input_file = args[1]
tumor1_sample_id = args[2]
tumor2_sample_id = args[3]
chr = args[4]
generated_files = args[5]

cat("####### ARGUMENTS #######\n")
cat(paste("input_file: ", input_file, "\n", sep=''))
cat(paste("tumor1_sample_id: ", tumor1_sample_id, "\n", sep=''))
cat(paste("tumor2_sample_id: ", tumor2_sample_id, "\n", sep=''))
cat(paste("chr: ", chr, "\n", sep=''))
cat(paste("generated_files: ", generated_files, "\n\n", sep=''))

# input_file = "/home/pgm/Workspace/MPM/VCF_finaux/germline_sandbox/M662_DA_5009_N_B00JAJD_chromosomes/M662_DA_5009_N_B00JAJD.GERMLINE.chrY.vcf"
# tumor1_sample_id = "B00JAJB"
# tumor2_sample_id = "B00JAJC"
# chr = "Y"
# generated_files = "/home/pgm/Workspace/MPM/marathon/falcon/generated_files/"


##########################################
## Set parameters
##########################################
file_name = unlist(strsplit(input_file, "_chromosomes/"))[2]

input_name_data  = unlist(strsplit(file_name, "_"))
input_name_data2 = unlist(strsplit(input_name_data[5], ".GERMLINE"))
patient_id = input_name_data[3]
normal_sample_id = input_name_data2[1]


##########################################
## Load germline VCF from 1 patient (2 tumors + 1 normal)
## Column order in VCF : TUMOR1  TUMOR2  NORMAL
##########################################
header = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", tumor1_sample_id, tumor2_sample_id, normal_sample_id)
cat(paste("\nLoading germline VCF with genotypes: ", file_name, "... "))
VCF_germline_content = read.table(input_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(VCF_germline_content) = header
cat(paste("Loaded.\n"))


##########################################
## Format data from the 2 tumor columns
##########################################
buildInput = function(VCF_germline_content, normal_sample_id, tumor_sample_id) {
  input_header = c("Tumor_ReadCount_Alt", "Tumor_ReadCount_Ref", "Tumor_ReadCount_Total", "Normal_ReadCount_Alt", "Normal_ReadCount_Ref", "Normal_ReadCount_Total", "Reference_Allele", "TumorSeq_Allele1", "TumorSeq_Allele2", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2", "Chromosome", "Start_position", "End_position")
  ## Create empy data frame
  input_tumor = data.frame(matrix(, nrow=nrow(VCF_germline_content), ncol=length(input_header)))
  colnames(input_tumor) = input_header
  
  cat("\nRetrieving read counts & genotype info... ")
  ## FORMAT:NR = DP
  ## FORMAT:NV = AO
  tumor_total_counts  = unlist(lapply(1:nrow(VCF_germline_content), function(i) get_genotype(VCF_germline_content[i, tumor_sample_id], VCF_germline_content[i, "FORMAT"], "NR")))
  tumor_alt_counts    = unlist(lapply(1:nrow(VCF_germline_content), function(i) get_genotype(VCF_germline_content[i, tumor_sample_id], VCF_germline_content[i, "FORMAT"], "NV")))
  tumor_genotype      = unlist(lapply(1:nrow(VCF_germline_content), function(i) get_genotype(VCF_germline_content[i, tumor_sample_id], VCF_germline_content[i, "FORMAT"], "GT", FALSE)))
  normal_total_counts = unlist(lapply(1:nrow(VCF_germline_content), function(i) get_genotype(VCF_germline_content[i, normal_sample_id], VCF_germline_content[i, "FORMAT"], "NR")))
  normal_alt_counts   = unlist(lapply(1:nrow(VCF_germline_content), function(i) get_genotype(VCF_germline_content[i, normal_sample_id], VCF_germline_content[i, "FORMAT"], "NV")))
  normal_genotype     = unlist(lapply(1:nrow(VCF_germline_content), function(i) get_genotype(VCF_germline_content[i, normal_sample_id], VCF_germline_content[i, "FORMAT"], "GT", FALSE)))
  cat("DONE.")
  
  cat("\nFilling input data frame... ")
  for (i in 1:nrow(VCF_germline_content)) { # for each mutation
    REF = VCF_germline_content[i, "REF"]
    ALT = VCF_germline_content[i, "ALT"]
    
    input_tumor[i, "Tumor_ReadCount_Alt"]     = tumor_alt_counts[i]
    input_tumor[i, "Tumor_ReadCount_Ref"]     = tumor_total_counts[i] - tumor_alt_counts[i]
    input_tumor[i, "Tumor_ReadCount_Total"]   = tumor_total_counts[i]
    input_tumor[i, "Normal_ReadCount_Alt"]    = normal_alt_counts[i]
    input_tumor[i, "Normal_ReadCount_Ref"]    = normal_total_counts[i] - normal_alt_counts[i]
    input_tumor[i, "Normal_ReadCount_Total"]  = normal_total_counts[i]
    
    tumor_allele1 = unlist(strsplit(tumor_genotype[i], "/"))[1]
    tumor_allele2 = unlist(strsplit(tumor_genotype[i], "/"))[2]
    input_tumor[i, "TumorSeq_Allele1"]        = if (tumor_allele1==0) REF else ALT
    input_tumor[i, "TumorSeq_Allele2"]        = if (tumor_allele2==0) REF else ALT
    
    normal_allele1 = unlist(strsplit(normal_genotype[i], "/"))[1]
    normal_allele2 = unlist(strsplit(normal_genotype[i], "/"))[2]
    input_tumor[i, "Match_Norm_Seq_Allele1"]  = if (normal_allele1==0) REF else ALT
    input_tumor[i, "Match_Norm_Seq_Allele2"]  = if (normal_allele2==0) REF else ALT
    
    input_tumor[i, "Reference_Allele"]        = REF
    input_tumor[i, "Chromosome"]              = gsub("chr", "", VCF_germline_content[i, "CHROM"])
    input_tumor[i, "Start_position"]          = VCF_germline_content[i, "POS"]
    input_tumor[i, "End_position"]            = VCF_germline_content[i, "POS"]
  }
  cat("DONE.")
  return(input_tumor)
}

primary = buildInput(VCF_germline_content, normal_sample_id, tumor1_sample_id)
relapse = buildInput(VCF_germline_content, normal_sample_id, tumor2_sample_id)

# save image file.
save.image(file=paste(generated_files, file_name, '.input_falcon.rda', sep=''))  


##########################################
## Run falcon
##########################################

# calculate depth ratio (total read counts of tumor versus normal)
rdep_relapse=sum(relapse$Tumor_ReadCount_Total)/sum(relapse$Normal_ReadCount_Total)
rdep_primary=sum(primary$Tumor_ReadCount_Total)/sum(primary$Normal_ReadCount_Total)

# Factorize Falcon code... never seen such a big code duplication...
process_chromosome = function(tumor_content, chr, patient_id, sample_id, rdep_tumor) {
  ###########################################
  # Focus on germline heterozygous variants.
  ###########################################
  
  # remove variants with missing genotype
  tumor_content=tumor_content[tumor_content[,'Match_Norm_Seq_Allele1']!=' ',]
  tumor_content=tumor_content[tumor_content[,'Match_Norm_Seq_Allele2']!=' ',]
  tumor_content=tumor_content[tumor_content[,'Reference_Allele']!=' ',]
  tumor_content=tumor_content[tumor_content[,'TumorSeq_Allele1']!=' ',]
  tumor_content=tumor_content[tumor_content[,'TumorSeq_Allele2']!=' ',]
  
  # get germline heterozygous loci (normal allele1 != normal allele2)
  tumor_content=tumor_content[(as.matrix(tumor_content[,'Match_Norm_Seq_Allele1'])!=as.matrix(tumor_content[,'Match_Norm_Seq_Allele2'])),]
  
  
  ############################################################
  # QC procedures to remove false neg and false pos variants.
  # The thresholds can be adjusted.
  ############################################################
  
  # remove indels (this can be relaxed but we think indels are harder to call than SNPs)
  indel.filter1=nchar(as.matrix(tumor_content[,'Reference_Allele']))<=1
  indel.filter2=nchar(as.matrix(tumor_content[,'Match_Norm_Seq_Allele1']))<=1
  indel.filter3=nchar(as.matrix(tumor_content[,'Match_Norm_Seq_Allele2']))<=1
  indel.filter4=nchar(as.matrix(tumor_content[,'TumorSeq_Allele1']))<=1
  indel.filter5=nchar(as.matrix(tumor_content[,'TumorSeq_Allele2']))<=1
  tumor_content=tumor_content[indel.filter1 & indel.filter2 & indel.filter3 & indel.filter4 & indel.filter5,]
  
  # total number of reads greater than 30 in both tumor and normal
  depth.filter1=(tumor_content[,"Normal_ReadCount_Ref"]+tumor_content[,"Normal_ReadCount_Alt"])>=30
  depth.filter2=(tumor_content[,"Tumor_ReadCount_Ref"]+tumor_content[,"Tumor_ReadCount_Alt"])>=30
  tumor_content=tumor_content[depth.filter1 & depth.filter2,]
  
  
  #########################
  # Generate FALCON input.
  #########################
  
  # Data frame with four columns: tumor ref, tumor alt, normal ref, normal alt.
  readMatrix.tumor_content=as.data.frame(tumor_content[,c('Tumor_ReadCount_Ref',
                                                  'Tumor_ReadCount_Alt',
                                                  'Normal_ReadCount_Ref',
                                                  'Normal_ReadCount_Alt')])
  colnames(readMatrix.tumor_content)=c('AT','BT','AN','BN')
  dim(readMatrix.tumor_content); dim(tumor_content)
  
  
  ###############################
  # Run FALCON and view results.
  ###############################
  
  if (nrow(tumor_content) > 0) {
    tauhat.tumor_content=getChangepoints(readMatrix.tumor_content)
    cn.tumor_content = getASCN(readMatrix.tumor_content, tauhat=tauhat.tumor_content, rdep = rdep_tumor, threshold = 0.3)
    
    # Chromosomal view of segmentation results.
    filename = paste(generated_files, 'falcon.patient_', patient_id, '.tumor_', sample_id, '.chr_',chr,sep='')
    pdf(file=paste(filename,'.pdf',sep=''), width=5, height=8)
    view(cn.tumor_content,pos=tumor_content[,'Start_position'], rdep = rdep_tumor)
    dev.off()
    
    # save image file.
    save.image(file=paste(filename, '.Falcon.rda', sep=''))  
    
    
    ########################################
    # Further curate FALCON's segmentation.
    ########################################
    
    # From the pdf above, we see that:
    # (1) There are small segments that need to be removed;
    # (2) Consecutive segments with similar allelic cooy number states need to be combined.
    if(length(tauhat.tumor_content)>0){
      length.thres=10^6  # Threshold for length of segments, in base pair.
      delta.cn.thres=0.3  # Threshold of absolute copy number difference between consecutive segments.
      falcon.qc.list = falcon.qc(readMatrix = readMatrix.tumor_content,
                                 tauhat = tauhat.tumor_content,
                                 cn = cn.tumor_content,
                                 st_bp = tumor_content[,"Start_position"],
                                 end_bp = tumor_content[,"End_position"],
                                 rdep = rdep_tumor,
                                 length.thres = length.thres,
                                 delta.cn.thres = delta.cn.thres)
      
      tauhat.tumor_content=falcon.qc.list$tauhat
      cn.tumor_content=falcon.qc.list$cn
    }
    
    # Chromosomal view of QC'ed segmentation results.
    pdf(file=paste(filename,'.QC.pdf',sep=''),width=5,height=8)
    view(cn.tumor_content,pos=tumor_content[,'Start_position'], rdep = rdep_tumor)
    dev.off()
    
    
    #################################################
    # Generate Canopy's input with s.d. measurement.
    #################################################
    
    # This is to generate table output including genomic locations for 
    # segment boudaries.
    # For Canopy's input, we use Bootstrap-based method to estimate the
    # standard deviations for the allele-specific copy numbers.
    falcon.output=falcon.output(readMatrix = readMatrix.tumor_content,
                                tauhat = tauhat.tumor_content,
                                cn = cn.tumor_content,
                                st_bp = tumor_content[,"Start_position"],
                                end_bp = tumor_content[,"End_position"],
                                nboot = 5000)
    falcon.output = cbind(chr=rep(chr,nrow(falcon.output)), falcon.output)
    write.table(falcon.output, file=paste(filename,'.output.txt',sep=''), col.names =T, row.names = F, sep='\t', quote = F)  
  }
}


#################################################
## Falcon processes each chromosome separately
#################################################
cat(paste("\nProcessing chromosome ", chr, "...", sep=""))
primary_chr=primary[which(primary[,'Chromosome']==chr),]
relapse_chr=relapse[which(relapse[,'Chromosome']==chr),]

cat(paste("\nProcessing primary tumor ", tumor1_sample_id, " of patient ", patient_id, "...\n", sep=""))
process_chromosome(primary_chr, chr, patient_id, tumor1_sample_id, rdep_primary)
cat(paste("\nProcessing relapse tumor ", tumor2_sample_id, " of patient ", patient_id, "...\n", sep=""))
process_chromosome(relapse_chr, chr, patient_id, tumor2_sample_id, rdep_relapse)
cat(paste("\nChromosome", chr, " DONE.\n\n", sep=""))
cat("########################################\n\n")
