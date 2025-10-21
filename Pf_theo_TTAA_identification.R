library(tidyverse)
library(IRanges)
library(Biostrings) 
library(ShortRead)
library(openxlsx)
library(dplyr)
library(data.table)
library(rtracklayer)
library(parallel)

###input
###output
###functions

#######################################################recognize TTAA by R instead of using homer
#######################################################recognize TTAA by R instead of using homer
#######################################################recognize TTAA by R instead of using homer
#######################################################recognize TTAA by R instead of using homer
#######################################################recognize TTAA by R instead of using homer
Genome.fasta_all <- readDNAStringSet("./Input/Genome/PlasmoDB-63_Pfalciparum3D7_Genome.fasta") 
names(Genome.fasta_all)
#No contigs in PF3D7 genome v3
TTAA <- "TTAA"
Genome_theo_TTAA <- vmatchPattern(TTAA, Genome.fasta_all, fixed = TRUE)
#transform Granges object into bed files
export.bed(Genome_theo_TTAA,con='./Output/PF3D7_theo_TTAA_R.bed')
Genome_theo_TTAA_bed <- read.table('./Output/PF3D7_theo_TTAA_R.bed')
Genome_theo_TTAA_bed_modified <- Genome_theo_TTAA_bed %>% dplyr::select(V1, V10, V11, V7, V13, V14)
Genome_theo_TTAA_bed_modified <- Genome_theo_TTAA_bed_modified %>% dplyr::rename(V1=V1, V2=V10, V3=V11, V4=V7, V5=V13, V6=V14)
nrow(Genome_theo_TTAA_bed_modified) ##330136

Genome_theo_TTAA_bed_modified_final <- data.frame(Genome_theo_TTAA_bed_modified, freq = 2)
#repeat each row for adding strandness column later
Genome_theo_TTAA_bed_modified_final<- Genome_theo_TTAA_bed_modified_final[rep(row.names(Genome_theo_TTAA_bed_modified_final), Genome_theo_TTAA_bed_modified_final$freq), ]
Genome_theo_TTAA_bed_modified_final <- Genome_theo_TTAA_bed_modified_final[, c(1:6)]
#change row names into order of number
row.names(Genome_theo_TTAA_bed_modified_final) <- c(1: nrow(Genome_theo_TTAA_bed_modified_final))
nrow(Genome_theo_TTAA_bed_modified_final)
#add strandness '+' and '-'
Genome_theo_TTAA_bed_modified_final$V6 <- rep(c('+', '-'), nrow(Genome_theo_TTAA_bed_modified))
Genome_theo_TTAA_bed_modified_final$V4 <- "TTAA"

write.table(Genome_theo_TTAA_bed_modified_final, './Output/PF3D7_theo_TTAA_R_modified_final.bed', col.names = F, quote = F, row.names = F, sep = '\t')


