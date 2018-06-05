#!/usr/bin/env Rscript


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
output_path                       = args[5]

##
cat("####### ARGUMENTS #######\n")
cat(paste("input_folder: ", input_folder, "\n", sep=''))
cat(paste("patient_id: ", patient_id, "\n", sep=''))
cat(paste("tumor1_id: ", tumor1_id, "\n", sep=''))
cat(paste("tumor2_id: ", tumor2_id, "\n", sep=''))


##########################################
## Prepare CNAs input
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
## Write data
##########################################
write.table(ascn_primary, file=paste(output_path, patient_id, '.', tumor1_id, '.CNA.tsv', sep=''), quote=F)
write.table(ascn_relapse, file=paste(output_path, patient_id, '.', tumor2_id, '.CNA.tsv', sep=''), quote=F)
