#!/usr/bin/env Rscript

library("falcon")
source("https://gist.githubusercontent.com/mfoll/a4dfbb92068dc559f130/raw/714dc8c2e97987fd4385dcef2722b3ef986d38d6/get_vcf_data.r")
source("/home/pgm/Workspace/MPM/marathon/libs/falcon.output.R")
# source("/home/soudadel/MPM/falcon/libs/falcon.output.R")
# source("/home/soudadel/MPM/falcon/libs/falcon.qc.R")


##########################################
## Retrieve arguments
##########################################
args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) { stop("Input name missing!\n", call.=FALSE) }

data_file_missing_epsilon = args[1]
data_file_coordinates     = args[2]
output_files              = args[3]

## test data
data_file_missing_epsilon = "/home/pgm/Workspace/MPM/marathon/falcon/output/patient_5009/chr6/falcon.patient_5009.tumor_B00JAJB.chr_6.Falcon.rda"
data_file_coordinates     = "/home/pgm/Workspace/MPM/marathon/falcon/output/patient_5009/chr6/falcon.patient_5009.tumor_placeholder.chr_6.output.txt"
output_files              = "/home/pgm/Workspace/MPM/marathon/falcon/output/"
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

# process_chromosome = function(tumor_content, tOri_to_filter_content, chr, patient_id, sample_id, rdep_tumor, coordinates) {
  
  tumor_content = relapse_chr
  tOri_to_filter_content = primary_chr
  sample_id = tumor2_sample_id
  rdep_tumor = rdep_primary
  coordinates = coordinates_t1
  
  
  ###########################################
  # Focus on germline heterozygous variants.
  ###########################################
  
  tumor_ori_filtered = tOri_to_filter_content
  
  # remove variants with missing genotype
  tumor_ori_filtered=tumor_ori_filtered[tumor_ori_filtered[,'Match_Norm_Seq_Allele1']!=' ',]
  tumor_ori_filtered=tumor_ori_filtered[tumor_ori_filtered[,'Match_Norm_Seq_Allele2']!=' ',]
  tumor_ori_filtered=tumor_ori_filtered[tumor_ori_filtered[,'Reference_Allele']!=' ',]
  tumor_ori_filtered=tumor_ori_filtered[tumor_ori_filtered[,'TumorSeq_Allele1']!=' ',]
  tumor_ori_filtered=tumor_ori_filtered[tumor_ori_filtered[,'TumorSeq_Allele2']!=' ',]
  
  # get germline heterozygous loci (normal allele1 != normal allele2)
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
  tumor_ori_filtered=tumor_ori_filtered[indel.filter1 & indel.filter2 & indel.filter3 & indel.filter4 & indel.filter5,]
  
  # total number of reads greater than 30 in both tumor and normal
  depth.filter1=(tumor_ori_filtered[,"Normal_ReadCount_Ref"]+tumor_ori_filtered[,"Normal_ReadCount_Alt"])>=30
  depth.filter2=(tumor_ori_filtered[,"Tumor_ReadCount_Ref"]+tumor_ori_filtered[,"Tumor_ReadCount_Alt"])>=30
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
        cn.tumor_content = getASCN(readMatrix.tumor_content, tauhat=tauhat.tumor_content, rdep = rdep_tumor, threshold = 0.3)
        
        
        ##################
        # function (readMatrix, rdep = NULL, tauhat = NULL, threshold = 0.15, 
        #           pOri = c(0.49, 0.51), error = 1e-05, maxIter = 1000) 
        # {
        
        
        readMatrix = readMatrix.tumor_content
        rdep = rdep_tumor
        tauhat = tauhat.tumor_content
        threshold = 0.3
        pOri = c(0.49, 0.51)
        error = 1e-05
        maxIter = 1000
        
        
          AN = readMatrix$AN
          BN = readMatrix$BN
          AT = readMatrix$AT
          BT = readMatrix$BT
          if (is.null(rdep)) 
            rdep = median(AT + BT)/median(AN + BN)
          if (is.null(tauhat)) 
            tauhat = getChangepoints(readMatrix, pOri = pOri, error = error, 
                                     maxIter = maxIter)
          N = length(AT)
          tau = sort(unique(c(1, tauhat, N)))
          K = length(tau) - 1
          pa = pb = rep(0, K)
          for (i in 1:K) {
            ids = tau[i]:(tau[i + 1] - 1)
            if (i == K) 
              ids = tau[i]:tau[i + 1]
            p = as.numeric(.Call("GetP", as.numeric(AT[ids]), as.numeric(BT[ids]), 
                                 as.numeric(AN[ids]), as.numeric(BN[ids]), as.numeric(error), 
                                 as.numeric(maxIter), as.numeric(pOri), PACKAGE = "falcon"))
            pa[i] = p[1]
            pb[i] = p[2]
            if (diff(p) < 0.1) {
              temp = as.numeric(.Call("LikH", as.numeric(AT[ids]), 
                                      as.numeric(BT[ids]), as.numeric(AN[ids]), as.numeric(BN[ids]), 
                                      as.numeric(p), PACKAGE = "falcon"))
              p2 = sum(AT[ids] + BT[ids])/sum(AT[ids] + BT[ids] + 
                                                AN[ids] + BN[ids])
              temp2 = as.numeric(.Call("Lik", as.numeric(AT[ids]), 
                                       as.numeric(BT[ids]), as.numeric(AN[ids]), as.numeric(BN[ids]), 
                                       as.numeric(rep(p2, 2)), PACKAGE = "falcon"))
              if (!is.na(temp)[1] && !is.na(temp[2]) && !is.na(temp2)) {
                bic = temp[1] - temp2 - temp[2]/2 + log(p2 * 
                                                          (1 - p2) * sum(AT[ids] + BT[ids] + AN[ids] + 
                                                                           BN[ids]))/2 + log(2 * pi)/2
              }
              else if (!is.na(temp)[1] && !is.na(temp2)) {
                bic = temp[1] - temp2 + log(p2 * (1 - p2) * 
                                              sum(AT[ids] + BT[ids] + AN[ids] + BN[ids]))/2 + 
                  log(2 * pi)/2
              }
              if (bic < 0) {
                pa[i] = p2
                pb[i] = p2
              }
            }
          }
          rawcns1 = pa/(1 - pa)/rdep
          rawcns2 = pb/(1 - pb)/rdep
          cns1 = hardthres(rawcns1, low = 1 - threshold, high = 1 + 
                             threshold)
          cns2 = hardthres(rawcns2, low = 1 - threshold, high = 1 + 
                             threshold)
          Haplotype = list()
          for (i in 1:K) {
            # if (cns1[i] != cns2[i]) {
              ids = tau[i]:(tau[i + 1] - 1)
              if (i == K) 
                ids = tau[i]:tau[i + 1]
              gt = 1/(1 + (pb[i]/pa[i])^(AT[ids] - BT[ids]) * 
                        ((1 - pb[i])/(1 - pa[i]))^(AN[ids] - BN[ids]))
              temp3 = rep("A", length(ids))
              temp3[which(gt > 0.5)] = "B"
              Haplotype[[i]] = temp3
            # }
          }
          cn.tumor_content = list(tauhat = tauhat, ascn = rbind(cns1, cns2), Haplotype = Haplotype, 
                      readMatrix = readMatrix)
        ##################
        hardthres = function(v, low=0.9, high=1.1){
          n = length(v)
          for (i in 1:n){ if (v[i]>low && v[i]<high) v[i] = 1 }
          v
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
        filename = paste(output_files, "patient_", patient_id, "/chr", chr, "/falcon.patient_", patient_id, ".tumor_", sample_id, ".chr_", chr, ".output_epsilon.txt", sep='')
        
        if (!file.exists(filename)) {
          write.table(na.omit(falcon.output), file=filename, col.names =T, row.names = F, sep='\t', quote = F, append=F)  
        } else {
          write.table(na.omit(falcon.output), file=filename, col.names =F, row.names = F, sep='\t', quote = F, append=T)  
        }
      } 
    }
  }
# }


#################################################
## Falcon processes each chromosome separately
#################################################
cat(paste("\nProcessing chromosome ", chr, "...", sep=""))
primary_chr=primary[which(primary[,'Chromosome']==chr),]
relapse_chr=relapse[which(relapse[,'Chromosome']==chr),]

cat("\nProcessing tumor1 with tumor2 tauhat")
data_file_coordinates_t2 = sub("placeholder", tumor2_sample_id, data_file_coordinates)
coordinates_t2 = na.omit(read.table(data_file_coordinates_t2, header=T))
process_chromosome(primary_chr, relapse_chr, chr, patient_id, tumor1_sample_id, rdep_primary, coordinates_t2)
cat("\nProcessing tumor2 with tumor1 tauhat")
data_file_coordinates_t1 = sub("placeholder", tumor1_sample_id, data_file_coordinates)
coordinates_t1 = na.omit(read.table(data_file_coordinates_t1, header=T))
process_chromosome(relapse_chr, primary_chr, chr, patient_id, tumor2_sample_id, rdep_relapse, coordinates_t1)
cat(paste("\nChromosome", chr, " DONE.\n\n", sep=""))
cat("########################################\n\n")

