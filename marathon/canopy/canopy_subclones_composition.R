#!/usr/bin/env Rscript


##########################################
## Retrieve arguments
##########################################
args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) { stop("Input name missing!\n", call.=FALSE) }

##
patient_id            = args[1]
data_path             = args[2]
file_name             = args[3]
input_somatic_VCF_t1  = args[4]
input_somatic_VCF_t2  = args[5]
only_exonic           = args[6]

##
cat("####### ARGUMENTS #######\n")
cat(paste("patient_id: ", patient_id, "\n", sep=''))
cat(paste("data_path: ", data_path, "\n", sep=''))
cat(paste("file_name: ", file_name, "\n", sep=''))


##########################################
## Load somatic VCF
##########################################
VCF_somatic_tumor1 = read.table(input_somatic_VCF_t1, sep = "\t", header=FALSE, stringsAsFactors=FALSE)
VCF_somatic_tumor1$uniqID <- do.call(paste, c(VCF_somatic_tumor1[c("V1", "V2", "V3", "V4", "V5")], sep = "_"))
VCF_somatic_tumor2 = read.table(input_somatic_VCF_t2, sep = "\t", header=FALSE, stringsAsFactors=FALSE)
VCF_somatic_tumor2$uniqID <- do.call(paste, c(VCF_somatic_tumor2[c("V1", "V2", "V3", "V4", "V5")], sep = "_"))
cols = c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene", "AAChange.refGene", "Otherinfo")
colnames(VCF_somatic_tumor1) = cols
colnames(VCF_somatic_tumor1)[ncol(VCF_somatic_tumor1)] = "uniqID"
colnames(VCF_somatic_tumor2) = cols
colnames(VCF_somatic_tumor2)[ncol(VCF_somatic_tumor2)] = "uniqID"


##########################################
## Load subclones
##########################################
output = data.frame(matrix(ncol = 7, nrow = 0))
col_names = c("subclone", "data", "chr", "start", "end", "alt", "ref")
colnames(output) = col_names

file_lines = readLines(paste(data_path, patient_id, "/", file_name, sep=''))
for(i in 1: length(file_lines)) {
  raw_subclone = unlist(strsplit(file_lines[i], ': '))
  raw_mutations = sub(" ", "", unlist(strsplit(raw_subclone[2], ', ')))

  for (j in 1:length(raw_mutations)) {
    found_mutation = VCF_somatic_tumor1[which(VCF_somatic_tumor1$uniqID == raw_mutations[j]),]
    if (nrow(found_mutation) == 0) found_mutation = VCF_somatic_tumor2[which(VCF_somatic_tumor2$uniqID == raw_mutations[j]),]
    if (nrow(found_mutation) == 0) {
      chr_tmp = unlist(strsplit(raw_mutations[j], ':'))
      chr = chr_tmp[1]
      start_tmp = unlist(strsplit(chr_tmp[2], '-'))
      start = start_tmp[1]
      end_tmp = unlist(strsplit(start_tmp[2], '_'))
      end = end_tmp[1]
      new_row = unname(data.frame(raw_subclone[1],
                                  paste('copy number:', raw_mutations[j]),
                                  chr,
                                  start,
                                  end,
                                  '',
                                  ''))
      colnames(new_row) = col_names
      output = rbind(output, new_row)
    } else if (only_exonic == 1 & found_mutation$Func.refGene == "exonic" & found_mutation$ExonicFunc.refGene != "synonymous SNV") {
      new_row = unname(data.frame(raw_subclone[1],
                                  paste(found_mutation$Gene.refGene, " (", found_mutation$ExonicFunc.refGene, ")", sep=""),
                                  found_mutation$Chr,
                                  as.character(found_mutation$Start),
                                  as.character(found_mutation$End),
                                  found_mutation$Ref,
                                  found_mutation$Alt))
      colnames(new_row) = col_names
      output = rbind(output, new_row)
    } else if (only_exonic == 0) {
      new_row = unname(data.frame(raw_subclone[1],
                                  paste(found_mutation$Gene.refGene, " (", found_mutation$Func.refGene, ")", sep=""),
                                  found_mutation$Chr,
                                  as.character(found_mutation$Start),
                                  as.character(found_mutation$End),
                                  found_mutation$Ref,
                                  found_mutation$Alt))
      colnames(new_row) = col_names
      output = rbind(output, new_row)
    }
  }
}

if (only_exonic == 0) output_file_name = "_TREE_most_likely.tsv" else output_file_name = "_TREE_most_likely.exonic.tsv"
write.table(output, file = paste(data_path, patient_id, "/", patient_id, output_file_name, sep=''), row.names=FALSE, sep="\t")
