#Script to calculate averages of biological replicates with appropriate controls

rm(list=ls())

setwd("/home/ubuntu/Drive0/Plastid_browser/plageb_6_0_arabidopsis_with_psoralen")

load_data <- function(file, strand, name) #function to load data from input files into a data frame
{
  input_data <- read.table(paste("./datasets/", file, sep = ""), sep='\t', stringsAsFactors = FALSE, colClasses = c("character", "integer", "integer", "integer", "double"))
  if(strand == '-') {    strand_value <- -1       
  } else{                strand_value <- 1 }
  coverage <- data.frame(input_data [,4], input_data[,5], name, strand_value, (input_data[,5]+1) / sum(input_data[,5]) * 1000000)
  colnames(coverage) <- c('position', 'data', 'name', 'strand', 'normalized')
  #cat(paste("Load_data ", file, "\n"))
  return(coverage)
}

data1 <- load_data("Psoralen_Cipro_sig6_biotin_1_chloroplast_mapping_cov.bed", 1, "Cipro_sig6_1_b")
data2 <- load_data("Psoralen_Cipro_sig7_biotin_2_chloroplast_mapping_cov.bed", 1, "Cipro_sig6_2_b")
data3 <- load_data("Psoralen_Cipro_sig8_biotin_3_chloroplast_mapping_cov.bed", 1, "Cipro_sig6_3_b")

input1 <- load_data("Psoralen_Cipro_sig6_input_1_chloroplast_mapping_cov.bed", 1, "Cipro_sig6_1_i")
input2 <- load_data("Psoralen_Cipro_sig7_input_2_chloroplast_mapping_cov.bed", 1, "Cipro_sig6_2_i")
input3 <- load_data("Psoralen_Cipro_sig8_input_3_chloroplast_mapping_cov.bed", 1, "Cipro_sig6_3_i")



pctinput <- data.frame(position = data1$position)
pctinput$rep1 <- data1$normalized / input1$normalized
pctinput$rep2 <- data2$normalized / input2$normalized
pctinput$rep3 <- data3$normalized / input3$normalized

pctinput$average <- apply(pctinput[, 2:4], 1, mean)
pctinput$sd <- apply(pctinput[, 2:4], 1, sd)

output <- cbind("Pt", "0", "154477", pctinput$position, pctinput$average, pctinput$sd)

write.table(output, file = "/home/ubuntu/Drive0/Plastid_browser/plageb_6_0_arabidopsis_with_psoralen/datasets/2nd_Psoralen_Cipro_sig6_avg.bed", quote = FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

#test_col0_avg.bed
#test_sig2_avg.bed
#test_sig6_avg.bed