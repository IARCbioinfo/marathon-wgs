#!/usr/bin/env Rscript

library("falcon")


##########################################
## Retrieve arguments
##########################################
args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) { stop("Input name missing!\n", call.=FALSE) }

data_file_missing_epsilon = args[1]
data_file_coordinates     = args[2]
output_files              = args[3]

## test data
# data_file_missing_epsilon = "/home/pgm/Workspace/MPM/marathon/falcon/output/patient_5009/chr6/falcon.patient_5009.tumor_B00JAJB.chr_6.Falcon.rda"
# data_file_coordinates     = "/home/pgm/Workspace/MPM/marathon/falcon/output/patient_5009/chr6/falcon.patient_5009.tumor_placeholder.chr_6.output.txt"
# output_files              = "/home/pgm/Workspace/MPM/marathon/falcon/output/"
##


##########################################
## Load falcon data
##########################################
load(data_file_missing_epsilon)


cat("####### ARGUMENTS #######\n")
cat(paste("data_file_missing_epsilon: ", data_file_missing_epsilon, "\n", sep=''))
cat(paste("data_file_coordinates: ", data_file_coordinates, "\n", sep=''))
cat(paste("patient_id: ", patient_id, "\n", sep=''))
cat(paste("tumor1_sample_id_primay: ", tumor1_sample_id, "\n", sep=''))
cat(paste("tumor2_sample_id_relapse: ", tumor2_sample_id, "\n", sep=''))
cat(paste("chr: ", chr, "\n", sep=''))
cat(paste("output_files: ", output_files, "\n\n", sep=''))


##########################################
## Run falcon
##########################################

# calculate depth ratio (total read counts of tumor versus normal)
rdep_relapse=sum(relapse$Tumor_ReadCount_Total)/sum(relapse$Normal_ReadCount_Total)
rdep_primary=sum(primary$Tumor_ReadCount_Total)/sum(primary$Normal_ReadCount_Total)

process_chromosome = function(tumor_content, tOri_to_filter_content, chr, patient_id, sample_id, rdep_tumor, coordinates) {

  ###########################################
  # Focus on germline heterozygous variants.
  ###########################################
  
  tumor_ori_filtered = tOri_to_filter_content
  
  # remove variants with missing genotype
  tumor_content=tumor_content[tumor_ori_filtered[,'Match_Norm_Seq_Allele1']!=' ',]
  tumor_ori_filtered=tumor_ori_filtered[tumor_ori_filtered[,'Match_Norm_Seq_Allele1']!=' ',]
  tumor_content=tumor_content[tumor_ori_filtered[,'Match_Norm_Seq_Allele2']!=' ',]
  tumor_ori_filtered=tumor_ori_filtered[tumor_ori_filtered[,'Match_Norm_Seq_Allele2']!=' ',]
  tumor_content=tumor_content[tumor_ori_filtered[,'Reference_Allele']!=' ',]
  tumor_ori_filtered=tumor_ori_filtered[tumor_ori_filtered[,'Reference_Allele']!=' ',]
  tumor_content=tumor_content[tumor_ori_filtered[,'TumorSeq_Allele1']!=' ',]
  tumor_ori_filtered=tumor_ori_filtered[tumor_ori_filtered[,'TumorSeq_Allele1']!=' ',]
  tumor_content=tumor_content[tumor_ori_filtered[,'TumorSeq_Allele2']!=' ',]
  tumor_ori_filtered=tumor_ori_filtered[tumor_ori_filtered[,'TumorSeq_Allele2']!=' ',]
  
  # get germline heterozygous loci (normal allele1 != normal allele2)
  tumor_content=tumor_content[(as.matrix(tumor_ori_filtered[,'Match_Norm_Seq_Allele1'])!=as.matrix(tumor_ori_filtered[,'Match_Norm_Seq_Allele2'])),]
  tumor_ori_filtered=tumor_ori_filtered[(as.matrix(tumor_ori_filtered[,'Match_Norm_Seq_Allele1'])!=as.matrix(tumor_ori_filtered[,'Match_Norm_Seq_Allele2'])),]
  
  
  ############################################################
  # QC procedures to remove false neg and false pos variants.
  # The thresholds can be adjusted.
  ############################################################
  
  # remove indels (this can be relaxed but we think indels are harder to call than SNPs)
  indel.filter1=nchar(as.matrix(tumor_ori_filtered[,'Reference_Allele']))<=1
  indel.filter2=nchar(as.matrix(tumor_ori_filtered[,'Match_Norm_Seq_Allele1']))<=1
  indel.filter3=nchar(as.matrix(tumor_ori_filtered[,'Match_Norm_Seq_Allele2']))<=1
  indel.filter4=nchar(as.matrix(tumor_ori_filtered[,'TumorSeq_Allele1']))<=1
  indel.filter5=nchar(as.matrix(tumor_ori_filtered[,'TumorSeq_Allele2']))<=1
  tumor_content=tumor_content[indel.filter1 & indel.filter2 & indel.filter3 & indel.filter4 & indel.filter5,]
  tumor_ori_filtered=tumor_ori_filtered[indel.filter1 & indel.filter2 & indel.filter3 & indel.filter4 & indel.filter5,]
  
  # total number of reads greater than 30 in both tumor and normal
  depth.filter1=(tumor_ori_filtered[,"Normal_ReadCount_Ref"]+tumor_ori_filtered[,"Normal_ReadCount_Alt"])>=30
  depth.filter2=(tumor_ori_filtered[,"Tumor_ReadCount_Ref"]+tumor_ori_filtered[,"Tumor_ReadCount_Alt"])>=30
  tumor_content=tumor_content[depth.filter1 & depth.filter2,]
  tumor_ori_filtered=tumor_ori_filtered[depth.filter1 & depth.filter2,]
  
  
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
    
    for (i in seq(nrow(coordinates))) {
      
      if (!is.na(coordinates[i,]$Major.sd) || !is.na(coordinates[i,]$Minor.sd)) {
        tauhat.tumor_content = c(coordinates[i,]$st_snp, coordinates[i,]$end_snp)
        source("/home/pgm/Workspace/MPM/marathon/libs/falcon.getASCN.epsilon.R")
        cn.tumor_content = falcon.getASCN.epsilon(readMatrix.tumor_content, tauhat=tauhat.tumor_content, rdep = rdep_tumor, threshold = 0.3)
        
        # falcon bugfix
        if (ncol(cn.tumor_content$ascn) <= 1 ) { # if only 1 ascn found, ther is a bug, so we duplicate
          cn.tumor_content$ascn = cbind(cn.tumor_content$ascn, cn.tumor_content$ascn[,1])
          cn.tumor_content$Haplotype[[2]] = cn.tumor_content$Haplotype[[1]]
        }
        
        #################################################
        # Generate Canopy's input with s.d. measurement.
        #################################################
        
        # This is to generate table output including genomic locations for 
        # segment boudaries.
        # For Canopy's input, we use Bootstrap-based method to estimate the
        # standard deviations for the allele-specific copy numbers.
        source("/home/pgm/Workspace/MPM/marathon/libs/falcon.output.R")
        falcon.output=falcon.output(readMatrix = readMatrix.tumor_content,
                                    tauhat = tauhat.tumor_content,
                                    cn = cn.tumor_content,
                                    st_bp = tumor_content[,"Start_position"],
                                    end_bp = tumor_content[,"End_position"],
                                    nboot = 5000,
                                    ascn_to_1 = 1)
        
        falcon.output = cbind(chr=rep(chr,nrow(falcon.output)), falcon.output)
        falcon.output[which(falcon.output == 0)] = NA
        filename = paste(output_files, "patient_", patient_id, "/chr", chr, "/falcon.patient_", patient_id, ".tumor_", sample_id, ".chr_", chr, ".output_epsilon.txt", sep='')
        
        if (!file.exists(filename)) {
          write.table(na.omit(falcon.output), file=filename, col.names =T, row.names = F, sep='\t', quote = F, append=F)  
        } else {
          write.table(na.omit(falcon.output), file=filename, col.names =F, row.names = F, sep='\t', quote = F, append=T)  
        }
      } 
    }
  }
}


#################################################
## Falcon processes each chromosome separately
#################################################
cat(paste("\nProcessing chromosome ", chr, "...", sep=""))
primary_chr=primary[which(primary[,'Chromosome']==chr),]
relapse_chr=relapse[which(relapse[,'Chromosome']==chr),]

cat("\nProcessing tumor1 with tumor2 tauhat")
data_file_coordinates_t2 = sub("placeholder", tumor2_sample_id, data_file_coordinates)
coordinates_t2 = na.omit(read.table(data_file_coordinates_t2, header=T))
if (nrow(coordinates_t2)) {
  process_chromosome(primary_chr, relapse_chr, chr, patient_id, tumor1_sample_id, rdep_primary, coordinates_t2)
}

cat("\nProcessing tumor2 with tumor1 tauhat")
data_file_coordinates_t1 = sub("placeholder", tumor1_sample_id, data_file_coordinates)
coordinates_t1 = na.omit(read.table(data_file_coordinates_t1, header=T))
if (nrow(coordinates_t1)) {
  process_chromosome(relapse_chr, primary_chr, chr, patient_id, tumor2_sample_id, rdep_relapse, coordinates_t1)
}

cat(paste("\nChromosome", chr, " DONE.\n\n", sep=""))
cat("########################################\n\n")

