library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(bedtoolsr)
library(scales)
library(ggpmisc)
library(edgeR)
library(venneuler)
library(grid)
library(rtracklayer)


cm <- read.xlsx("./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon_all.xlsx")
##################Adding Annotation from UTR regions###################
TTAA_sites <- read.table("/Users/sidaye/Documents/R/API_TnSeq/Output/TTAAhitsPf/PF3D7_theo_TTAA_R_modified_final.bed")
gff_file <- "/Users/sidaye/Documents/R/API_TnSeq/Input/Genome/PlasmoDB-63_Pfalciparum3D7.gff"
gff <- import(gff_file)

function_gff_to_utr_bed <- function(x){
  # Filter the GFF data to extract UTR information
  utr <- subset(gff, type %in% c("five_prime_UTR", "three_prime_UTR"))
  head(utr)
  
  #####Turn utr region from Granges object into bed file format#####
  seqnames <- as.character(seqnames(utr))
  ranges <- as.data.frame(ranges(utr))
  strand <- as.character(strand(utr))
  type <- mcols(utr)$type
  ID <- unlist(mcols(utr)$Parent)
  
  utr_bed <- data.frame(
    seqnames = seqnames,
    start = ranges$start,
    end = ranges$end,
    ID = ID,
    type = type,
    strand = strand,
    width = ranges$width)
  ###Use \\. to represent the literal dot
  utr_bed$ID <- unlist(lapply(strsplit(utr_bed$ID,"\\."),'[[',1))
  return(utr_bed)
}

utr_bed <- function_gff_to_utr_bed(gff)
utr5_bed1 <- utr_bed%>%dplyr::filter(type=="five_prime_UTR")
utr3_bed1 <- utr_bed%>%dplyr::filter(type=="three_prime_UTR")
utr5_bed <- data.frame(V1=utr5_bed1$seqnames, 
                       V2=utr5_bed1$start,
                       V3=utr5_bed1$end,
                       V4=utr5_bed1$ID,
                       V5=utr5_bed1$width,
                       V6=utr5_bed1$strand)

utr3_bed <- data.frame(V1=utr3_bed1$seqnames, 
                       V2=utr3_bed1$start,
                       V3=utr3_bed1$end,
                       V4=utr3_bed1$ID,
                       V5=utr3_bed1$width,
                       V6=utr3_bed1$strand)

options(bedtools.path = "/opt/homebrew/bin/")
intersect_bed_utr5 <- bedtoolsr::bt.intersect(a = TTAA_sites, b = utr5_bed, u = T)
intersect_bed_utr3 <- bedtoolsr::bt.intersect(a = TTAA_sites, b = utr3_bed, u = T)

intersect_bed_utr5$TTAA_ID <- paste(intersect_bed_utr5$V1,":",intersect_bed_utr5$V2)
intersect_bed_utr3$TTAA_ID <- paste(intersect_bed_utr3$V1,":",intersect_bed_utr3$V2)


TTAA_sites_index <- cm%>%dplyr::select(Chrom, Site)
TTAA_sites_index$TTAA_ID <- paste(TTAA_sites_index$Chrom,":",TTAA_sites_index$Site)

cm$Location[TTAA_sites_index$TTAA_ID%in%intersect_bed_utr5$TTAA_ID] <- "UTR5"
cm$Location[TTAA_sites_index$TTAA_ID%in%intersect_bed_utr3$TTAA_ID] <- "UTR3"
table(cm$Location)
dim(cm)
write.xlsx(cm, "./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon_UTR_annotation_all.xlsx")



