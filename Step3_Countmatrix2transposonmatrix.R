library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(doParallel)
library(rtracklayer)

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

#######################This script is just for turning original count matrix into transposon matrix#######################
###The transposon matrix will tell the overlapped exons
###!!!important: Do not use exon converted count matrix to be turned into transposon matrix, use the original one###############

###Input overlapped exons rows
overlaps <- read.xlsx('./Output/count_matrix/overlaps/count_matrix_all_overlaps.xlsx')

#!!!!! No gene description column

###Input the original count matrix with gene.description without total column and exon conversion
#countmatrix<- read.xlsx("./Output/count_matrix/Pf_count_matrix_all_Bg_removed_siteslevel_final.xlsx")
countmatrix<- read.xlsx("./Output/count_matrix/Pf_count_matrix_all_Bg_removed_siteslevel_cpm_final.xlsx")

######make sure to remove the duplicated rows for overlaps########
countmatrix <- countmatrix[1:660272,]

###Input sample names(colnames), only include WT samples
sample_names <- as.data.frame(colnames(countmatrix)[c(13:ncol(countmatrix))])

colnames(sample_names) <- "V1"

countmatrix2transposonmatrix <- function(tmp, samples){
  tmp_sense <- tmp[seq(1,nrow(tmp),2),]
  tmp_antisense <- tmp[seq(2,nrow(tmp),2),]
  UII_no_strandness <- data.frame(Chrom = tmp_sense$Chrom,
                                  Site =  tmp_sense$Site,
                                  Assigned_location = rep(NA, nrow(tmp_sense)),
                                  Location_sense = tmp_sense$Location,
                                  Location_antisense = tmp_antisense$Location
  )
  ####!different from counts_theor_prepare.R 
  
  UII_no_strandness$Location_sense[which(UII_no_strandness$Location_sense == 'intergenic')] <- 0
  UII_no_strandness$Location_sense[which(UII_no_strandness$Location_sense == 'intron')] <- 1
  UII_no_strandness$Location_sense[which(UII_no_strandness$Location_sense =='exon')] <- 3
  
  UII_no_strandness$Location_antisense[which(UII_no_strandness$Location_antisense == 'intergenic')] <- 0
  UII_no_strandness$Location_antisense[which(UII_no_strandness$Location_antisense == 'intron')] <- 1
  UII_no_strandness$Location_antisense[which(UII_no_strandness$Location_antisense =='exon')] <- 3
  
  UII_no_strandness$Location_sense <- as.numeric(UII_no_strandness$Location_sense)
  UII_no_strandness$Location_antisense <- as.numeric(UII_no_strandness$Location_antisense)
  
  UII_no_strandness$Assigned_location <- UII_no_strandness %>% transmute(Assigned_location = (Location_sense + Location_antisense))
  
  UII_count_matrix <- UII_no_strandness
  table(UII_count_matrix$Assigned_location) ####confirm the number of exon-exon overlaps
  UII_no_strandness$Assigned_location[which(UII_no_strandness$Assigned_location == 0),] <- 'intergenic'
  UII_no_strandness$Assigned_location[which(UII_no_strandness$Assigned_location == 1),] <- 'intron'
  UII_no_strandness$Assigned_location[which(UII_no_strandness$Assigned_location == 2),] <- 'intron'
  UII_no_strandness$Assigned_location[which(UII_no_strandness$Assigned_location == 3),] <- 'exon'
  UII_no_strandness$Assigned_location[which(UII_no_strandness$Assigned_location == 6),] <- 'exon'
  UII_no_strandness$Assigned_location[which(UII_no_strandness$Assigned_location == 4),] <- 'exon&intron'
  
  geneID <- data.frame(sense_geneID = tmp_sense$GeneID,
                       antisen_geneID = tmp_antisense$GeneID,
                       geneID = rep(NA, nrow(tmp_sense))
                       
  )
  
  #gtf_genome_58 <- read.table('/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/genome/PlasmoDB-58_PknowlesiH.gtf', sep = '\t')
  #table(gtf_genome_58$V3)
  
  #sense_geneID_ls <- geneID %>% dplyr::filter(sense_geneID !='NA')
  #sense_geneID_ls <- sense_geneID_ls$sense_geneID
  #antisen_geneID_ls <-geneID %>% dplyr::filter(antisen_geneID !='NA') 
  #antisen_geneID_ls <- antisen_geneID_ls$antisen_geneID
  
  #just to append all the geneID in count matrix and to see the number of unique geneIDs
  #gene_ID_ls <- append(sense_geneID_ls,antisen_geneID_ls)
  #TotalNo_inserted_theo_genes <- length(unique(gene_ID_ls))
  
  
  #merge sense matrix and antisense matrix
  ##be careful here, make sure if it is 13 since we change the format of matrix for several times 
  count_merge_sense <- tmp_sense %>% dplyr::select(contains('TPN'))
  count_merge_antisense <- tmp_antisense%>% dplyr::select(contains('TPN'))
  
  #make an empty dataframe
  #if it is only 1 sample:
  #count_merge <- as.data.frame(matrix(nrow = length(count_merge_sense), ncol =(length(as.list(samples$V1))) ))
  
  ###!!!!! if it is more than 1 sample
  count_merge <- as.data.frame(matrix(nrow = nrow(count_merge_sense), ncol =(length(as.list(samples$V1))) ))
  
  colnames(count_merge) <-as.list(samples$V1)
  
  #for(i in 1:nrow(count_merge_sense)){
  #  for(j in 1:nrow(samples)){
  #    count_merge[i,j] <- (count_merge_sense[i,j]+count_merge_antisense[i,j])
  #  }
  #}
  for(i in 1:nrow(samples)){
    count_merge[,i] <- (count_merge_sense[,i]+count_merge_antisense[,i])
  }
  
  transposon_count_matrix <- data.frame(Chrom = UII_no_strandness$Chrom,
                                        Site = UII_no_strandness$Site,
                                        Location = UII_no_strandness$Assigned_location,
                                        sense_geneID = geneID$sense_geneID,
                                        antisen_geneID = geneID$antisen_geneID,
                                        Present_in_any_samples = rep(NA, nrow(UII_no_strandness)))
  
  
  count_merge_total <- as.data.frame(rowSums(count_merge))
  colnames(count_merge_total) <- 'Total'
  transposon_count_matrix <- cbind(transposon_count_matrix,count_merge, count_merge_total)
  transposon_count_matrix$Present_in_any_samples = transposon_count_matrix$Total
  
  transposon_count_matrix$Present_in_any_samples[which(transposon_count_matrix$Present_in_any_samples!=0)] <- 'yes'
  transposon_count_matrix$Present_in_any_samples[which(transposon_count_matrix$Present_in_any_samples==0)] <- 'no'
  table(transposon_count_matrix$Present_in_any_samples)
  return(transposon_count_matrix)
}

transposon_count_matrix_all <- countmatrix2transposonmatrix(tmp=countmatrix, samples=sample_names)

###double check overlapped exons
overlapped_tm <- transposon_count_matrix_all%>% dplyr::filter((!is.na(sense_geneID) & !is.na(antisen_geneID)))


#################To label UTR regions###############
#####Need to excluding the sites at UTR regions
TTAA_sites <- read.table("./Output/TTAAhitsPf/PF3D7_theo_TTAA_R_modified_final.bed")


gff_file <- "./Input/Genome/PlasmoDB-63_Pfalciparum3D7.gff"
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

TTAA_sites_index <- transposon_count_matrix_all%>%dplyr::select(Chrom, Site)
TTAA_sites_index$TTAA_ID <- paste(TTAA_sites_index$Chrom,":",TTAA_sites_index$Site)

transposon_count_matrix_all$Assigned_location[TTAA_sites_index$TTAA_ID%in%intersect_bed_utr5$TTAA_ID] <- "UTR5"
transposon_count_matrix_all$Assigned_location[TTAA_sites_index$TTAA_ID%in%intersect_bed_utr3$TTAA_ID] <- "UTR3"
table(transposon_count_matrix_all$Assigned_location)


#write.xlsx(transposon_count_matrix_all, "./Output/transposon_matrix/all/transposon_count_matrix_all_Bg_removed_siteslevel.xlsx", na.string='NA', keepNA=F)
#write.xlsx(transposon_count_matrix_all, "./Output/transposon_matrix/all/transposon_count_matrix_all_Bg_removed_siteslevel_final.xlsx", na.string='NA', keepNA=F)
write.xlsx(transposon_count_matrix_all, "./Output/transposon_matrix/all/transposon_count_matrix_all_Bg_removed_siteslevel_cpm_final.xlsx", na.string='NA', keepNA=F)
#write.xlsx(overlapped_tm, "./Output/transposon_matrix/overlaps/transposon_count_matrix_all_overlapped_exons_Bg_removed_siteslevel.xlsx", na.string='NA', keepNA=F)


