library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(scales)

###note:
#1 solve the problem of overlapped exons
#2 repeat the GeneID and gene.description in exons only as well for the convenience of downstream analysis
#3 add the antisense strand of overlapped exons to the tail of dataframe
getwd()
setwd('/Users/sidaye/Documents/R/API_TnSeq')


cm1 <- read.xlsx("./Input/count_matrix202412/counts.distri.Pf.all.xlsx")
cm2 <- read.xlsx("./Input/count_matrix202502/counts.distri.Pf.all_202502.xlsx")
cm3 <- read.xlsx("./Input/count_matrix202504/counts.distri.Pf.all_202504.xlsx")
cm4 <- read.xlsx("./Input/count_matrix202505/counts.distri.Pf.all_202505.xlsx")
cm5 <- read.xlsx("./Input/count_matrix202507/counts.distri.Pf.all_202507.xlsx")
cm6 <- read.xlsx("./Input/count_matrix202508/counts.distri.Pf.all_202508.xlsx")

colnames(cm4)[c(15, 17, 22, 23, 25, 27)] <- c("Day14_C_0505_2","Day14_M_0505_2","Day18_C_0505","Day18_M_0505","Day20_C_0505_2","Day20_M_0505_2")
colnames(cm1)[c(12:ncol(cm1))] <- paste0(colnames(cm1[,c(12:ncol(cm1))]), "_TPN")
colnames(cm2)[c(12:ncol(cm2))] <- paste0(colnames(cm2[,c(12:ncol(cm2))]), "_TPN")
colnames(cm3)[c(12:ncol(cm3))] <- paste0(colnames(cm3[,c(12:ncol(cm3))]), "_TPN")
colnames(cm4)[c(12:ncol(cm4))] <- paste0(colnames(cm4[,c(12:ncol(cm4))]), "_TPN")

cm <- cbind(cm1[,c(12:ncol(cm1))],cm2[,c(12:ncol(cm2))])
cm <- cbind(cm,cm3[,12:ncol(cm3)])
cm <- cbind(cm,cm4[,12:ncol(cm4)])
cm <- cbind(cm,cm5[,12:ncol(cm5)])
cm <- cbind(cm,cm6[,12:ncol(cm6)])
colnames(cm)

cm_tmp <- read.xlsx("./Output/count_matrix/Pf_count_matrix_202412_202502_202504_202505.xlsx")

cm_tmp <- cm_tmp[,1:11]
cm <- cbind(cm_tmp, cm)
dim(cm)
########To add gene.description column from PlasmoDB##########
########gene.description column should be at column6 after GeneID column###########
########There are multiple transcripts corresponds to one geneID, need to remove duplication######
gene.description <- read.csv("./Input/5720_total_Pf_product_description.csv")
dim(gene.description)
gene.description <-gene.description[,c(1,3)]
colnames(gene.description)[1] <- "GeneID"
colnames(gene.description)[2] <- "gene.description"
gene.description <- gene.description[!duplicated(gene.description$GeneID), ]
dim(gene.description) #####In total, 5720 unique geneID
cm_i <- left_join(cm[,1:5], gene.description, by="GeneID")
cm_final <- cbind(cm_i, cm[,6:ncol(cm)])
colnames(cm_final)

write.xlsx(cm_final, "./Output/count_matrix/Pf_count_matrix_all.xlsx")

####Colnames for cm4 needs to be changed and every sample name should add "_TPN" at tail in order to be filtered####
####Need to add gene description column first

#cm_all <- read.xlsx("./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged.xlsx")
#cm_bgremoved <- read.xlsx("./Output/count_matrix/all/cm_75Pk_essentialomeOnly_Bg_removed_siteslevel_per_sample.xlsx")
########there are overlapped TTAA both locates on two strands' exons

#mode1: check exon-exon overlaps
#mode2: check exon-intron or intron exon overlaps
#mode3: check intron-intron or intron exon overlaps
cm <- read.xlsx("./Output/count_matrix/Pf_count_matrix_all.xlsx") ####No UTR annotation####
#colnames(cm)[6] <- "gene.description"
#write.xlsx(cm, "./Output/count_matrix/Pf_count_matrix_all.xlsx")
table(cm$Location)
check_overlaps <- function(cm, mode){
  cm$Index <- seq(from = 1, to = nrow(cm), by=1)
  ###Check where the overlaps are
  ###Create an intermediate checking table
  Intermediate <- cm%>%dplyr::select(Chrom, Site, GeneID, gene.description, Location, Index)
  sense_Intermediate <- Intermediate[seq(from = 1, to = nrow(cm), by=2),]
  names(sense_Intermediate)[grep('Location',names(sense_Intermediate))] <- 'Location_sense'
  names(sense_Intermediate)[grep('GeneID',names(sense_Intermediate))] <- 'GeneID_sense'
  names(sense_Intermediate)[grep('gene.description',names(sense_Intermediate))] <- 'gene.description_sense'
  
  antisense_Intermediate <- Intermediate[seq(from = 0, to = nrow(cm), by=2),]
  names(antisense_Intermediate)[grep('Location',names(antisense_Intermediate))] <- 'Location_antisense'
  names(antisense_Intermediate)[grep('GeneID',names(antisense_Intermediate))] <- 'GeneID_antisense'
  names(antisense_Intermediate)[grep('gene.description',names(antisense_Intermediate))] <- 'gene.description_antisense'
  names(antisense_Intermediate)[grep('Index',names(antisense_Intermediate))] <- 'Index_antisense'
  antisense_Intermediate <- antisense_Intermediate%>%dplyr::select(GeneID_antisense, gene.description_antisense, Location_antisense, Index_antisense)
  checking_df <- cbind(sense_Intermediate,antisense_Intermediate)
  if (mode==1){
    overlaps <- checking_df %>% dplyr::filter(Location_sense == 'exon' & Location_antisense == 'exon')
  }else if (mode==2){
    overlaps <- checking_df %>% dplyr::filter((Location_sense == 'exon' & Location_antisense == 'intron')|(Location_sense == 'intron' & Location_antisense == 'exon'))
  }else{
    overlaps <- checking_df %>% dplyr::filter((Location_sense == 'intron' & Location_antisense == 'intron'))
  }
  
  return(overlaps)
}

overlaps <- check_overlaps(cm, mode=1)
nrow(overlaps) #####5061
write.xlsx(overlaps, "./Output/count_matrix/overlaps/count_matrix_all_overlaps.xlsx", na.string='NA', keepNA=F)
overlaps <- read.xlsx('./Output/count_matrix/overlaps/count_matrix_all_overlaps.xlsx')

####Attached the duplicated overlaps sites
conversion_cm <- function(cm, overlaps){
  cm$strand <- rep(c("+", "-"), length.out = nrow(cm))
  overlapped_antisense <- cm[overlaps$Index_antisense,]
  overlapped_sense <- cm[overlaps$Index,]
  ###convert exons
  index_exon <- grep("exon", cm$Location)
  ####For unoverlapped exons, if the exon locates at sense strand, +1: change the next row, if antisense strand, -1: change the previous row)
  strand_vector <- cm$strand[index_exon]
  converted_vector <- ifelse(strand_vector == "+", 1, -1)
  index_converted_exon <- index_exon+converted_vector
  table(cm$Location[index_converted_exon]) #### If there are overlapped exons existed, can be showed and No_exon/2 = No_overlaps
  cm$Location[index_converted_exon] <- "exon"
  
  ###convert geneID and gene description as well
  length(index_exon)
  ###remove overlapped antisense exon and will add it in the tails of the dataframe
  index_exon2 <- index_exon[!(index_exon %in% overlaps$Index_antisense)]
  strand_vector2 <- strand_vector[!(index_exon %in% overlaps$Index_antisense)]
  length(index_exon2)
  length(strand_vector2)
  converted_vector2 <- ifelse(strand_vector2 == "+", 1, -1)
  index_converted_exon2 <- index_exon2+converted_vector2
  table(cm$Location[index_converted_exon2]) ###check the converted  index are all non-exons and with No_overlaps of exons
  cm$GeneID[index_converted_exon2] <- cm$GeneID[index_exon2]
  cm$gene.description[index_converted_exon2] <- cm$gene.description[index_exon2]
  
  ###only change geneID, gene.description of the overlaps_sense into overlaps_antisense, and rbind the dataframe
  overlapped_sense$GeneID <- overlapped_antisense$GeneID
  overlapped_sense$gene.description <- overlapped_antisense$gene.description
  cm_new <- rbind(cm,overlapped_sense)
  cm_new <- rbind(cm_new,overlapped_antisense)
  cm_new <- cm_new[,-ncol(cm_new)]
  return(cm_new)
}

cm_converted <- conversion_cm(cm, overlaps)
#cm_converted_all <- conversion_cm(cm_all, overlaps)
#cm_converted_bgremoved <- conversion_cm(cm_bgremoved, overlaps)

write.xlsx(cm_converted, "./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon_all.xlsx")


########To check how many final usable reads in the count matrix for each sample
cm<- read.xlsx("./Output/count_matrix/Pf_count_matrix_all.xlsx")
#cm_all <- cm_all[,1:(ncol(cm_all)-1)]
#colnames(cm_all)[12:ncol(cm_all)] <- paste0(colnames(cm_all)[12:ncol(cm_all)], "_TPN")
#cm5 <- read.xlsx("./Input/count_matrix202507/counts.distri.Pf.all_202507.xlsx")
#cm5 <- cm5%>%dplyr::select(contains("TPN"))
#cm6 <- read.xlsx("./Input/count_matrix202508/counts.distri.Pf.all_202508.xlsx")
#cm6 <- cm6%>%dplyr::select(contains("TPN"))
#cm56 <- cbind(cm5,cm6)
#colnames(cm56)
#cm_all2 <- cbind(cm_all,cm56)

Total_count <- function(countmatrix){
  countmatrix2 <- countmatrix%>%dplyr::select(contains("TPN"))
  total <- colSums(countmatrix2)
  df_total <- as.data.frame(total)
  df_total <- rownames_to_column(df_total, var = "Sample")
  return(df_total)
}

PfTPN_totalcount <- Total_count(cm)

PfTPN_totalcount$Sample <- factor(
  PfTPN_totalcount$Sample,
  levels = PfTPN_totalcount$Sample
)

ggplot(PfTPN_totalcount, aes(x = Sample, y = total)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "",
       x = "Sample ID",
       y = "Final usable reads") +
  scale_y_continuous(labels = label_comma()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
