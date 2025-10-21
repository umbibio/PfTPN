library(tidyverse)
library(IRanges)
library(Biostrings) 
library(ShortRead)
library(openxlsx)
library(dplyr)
library(data.table)

cm_converted <- read.xlsx("./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon202412_202502_202504.xlsx")
cm_converted2 <- cm_converted%>%dplyr::select(Chrom,Site,pTN9_Day20_C,Day17,Day20,Day17_C_0504,Day18_C_0504,Day20_C_0504_2,Day20_C_0504) ###7 samples after Day15

cm_converted2$Total <- rowSums(cm_converted2[,3:9])
bed_data <- data.frame(V1=cm_converted2$Chrom,
                       V2=cm_converted2$Site-1,
                       V3=cm_converted2$Site+4,
                       V4=cm_converted2$Total,
                       V5=NaN,
                       V6=rep(c("+", "-"), times = nrow(cm_converted2)/2))

####To filter out 0 counts####
bed_data2 <- bed_data%>%filter(V4!=0)

write.table(
  bed_data2,
  file = "./Output/PlasmoDB/SP_allWTafterDay15.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

TTAA_hit <- read.table("./Output/TTAAhitsPf/TTAAhits_R_gtf_include_contigs.bed")

