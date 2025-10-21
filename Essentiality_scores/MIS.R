library(tidyverse)
library(IRanges)
library(Biostrings) 
library(ShortRead)
library(openxlsx)
library(dplyr)
library(data.table)
library(UpSetR)
library(ggVennDiagram)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(bedtoolsr)
library(mixtools)
library(scales)
library(colorspace)
library(cowplot)
library(rtracklayer)
library(ggrepel)
library(edgeR)
library(cowplot)

#This script is to pipe the table of normalized reads after bg noise removed in MIS model construction, remain CDS length for benmark analysis
#The input can be count matrix(bidirectional) or transposon matrix(unidirectional)
##########################################Start from here###################
##########################################Start from here###################
##########################################Start from here###################

####CPM normalization does not work##########
transposon_count_matrix <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix_all_Bg_removed_siteslevel_final.xlsx")
#transposon_count_matrix <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix_all_Bg_removed_siteslevel_cpm_final.xlsx")
colnames(transposon_count_matrix) ####Check if there is total column


sites_ID <- transposon_count_matrix[,c(1:6)]
cm_matrix <- transposon_count_matrix[,c(7:(ncol(transposon_count_matrix)-1))]

#####optional:CPM normalization########
####cpm normalization can be performed on count matrix or on transposon matrix###############
#cpm_cm_matrix <- cpm(cm_matrix)
#transposon_count_matrix <- cbind(sites_ID, cpm_cm_matrix)
#colnames(transposon_count_matrix)
#table(transposon_count_matrix$Assigned_location)
#####optional:normalization by JA method###############
scale_factor <- 25000000
count_matrix_norm <- sweep(cm_matrix, 2, colSums(cm_matrix), FUN = "/") * scale_factor
transposon_count_matrix <- cbind(sites_ID, count_matrix_norm)
colnames(transposon_count_matrix)
#####optional:normalization by JA method###############

####Extract sample
#sample_names <- as.data.frame(colnames(transposon_count_matrix))
#sample_names <- as.data.frame(colnames(transposon_count_matrix)[c(7,8,10,12,14,16,18,20,21)])
#colnames(sample_names) <- "V1"
#transposon_count_matrix_WT <- transposon_count_matrix[,c(1:6,7,9,11,13,15,17,19,20,22,23,24,25,28,30,32,33,36,38)] ###11 WT samples 
#transposon_count_matrix_WT <- transposon_count_matrix[,c(1:6,7,9,15,17,19,24,25,28,30,32,33)] ###11 WT samples after Day9
#transposon_count_matrix_WT <- transposon_count_matrix[,c(1:6,9,19,32,33)] ###4 WT Day20 samples

####WT after day 14############
#transposon_count_matrix_WT <- transposon_count_matrix[,c(1:6,7,8,10,12,14,16,18,20,21,23,24,25,26,29,31,33,34,37,39,41,42,44,45,48,51,54,55)] ### All WT samples
transposon_count_matrix_WT <- transposon_count_matrix[,c(1:6,10,16,18,20,25,26,29,31,33,34,44,45,48,51,54,55,60,62,64,66,
                                                         68,72,74,76,78,80,84,86,88,90,94,96,98,101,107,110,114,115,123,
                                                         126,129,136,138,140)] ### WT samples after Day14


transposon_count_matrix_WT <- transposon_count_matrix[,c(1:6,10,16,18,20,25,26,29,31,33,34,44,45,48,51,54,55,60,62,64,66,
                                                         68,72,74,76,78,80,84,86,88,90,94,96,98,101,107,110,114,115,123,
                                                         126,129,136,138,140)] ### WT samples after Day14

colnames(transposon_count_matrix_WT)

#####Optional:MSD checking(similarity between samples)######
mat <- as.matrix(transposon_count_matrix_WT[,c(7:ncol(transposon_count_matrix_WT))])
#d <- dist(t(mat)) ###Since dist() calculate the distance between rows, need to transpose the matrix, dist(t(mat)) is for Eucldean distanc
d <- as.dist(1 - cor(mat, method = "pearson")) ###1 - Pearson correlation Focus more on patterns not absolute values

mds <- cmdscale(d, k = 2) #K=2 means 2D MDS plot

mds_df <- data.frame(
  x = mds[,1],
  y = mds[,2],
  sample = colnames(mat)
)

ggplot(mds_df, aes(x = x, y = y, label = sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2) +  # <- non-overlapping labels
  theme_minimal() +
  labs(x = "MDS1", y = "MDS2", title = "MDS Plot")
#####Optional:MSD checking(similarity between samples)######



#####Optional:column sum checking######
colSums(transposon_count_matrix_WT[,c(7:ncol(transposon_count_matrix_WT))])
check <- transposon_count_matrix_WT
check$Total <- rowSums(check%>%dplyr::select(contains('TPN')))
check <- check[,c(1:6,ncol(check))]
#####Optional:column sum checking######


####Treated group########Day20_M_0504_TPN(36), Day20_M_0505_TPN(57)#############
transposon_count_matrix_WT <- transposon_count_matrix[,c(1:6,11,17,19,30,32,35,49,52,56,57,59,61,63,65,67,69,73,
                                                         77,79,81,87,89,91,97,98,109,112,118,119,128,131,137,139,141)] ### MEV samples after Day15


#####To check if UTR regions has been labelled, if not, need to exclude the TTAA sites locating at UTR regions#####
table(transposon_count_matrix_WT$Assigned_location)
#####Need to excluding the sites at UTR regions
#TTAA_sites <- read.table("./Output/TTAAhitsPf/PF3D7_theo_TTAA_R_modified_final.bed")


#gff_file <- "./Input/Genome/PlasmoDB-63_Pfalciparum3D7.gff"
#gff <- import(gff_file)

#function_gff_to_utr_bed <- function(x){
#  # Filter the GFF data to extract UTR information
#  utr <- subset(gff, type %in% c("five_prime_UTR", "three_prime_UTR"))
#  head(utr)
  
  #####Turn utr region from Granges object into bed file format#####
#  seqnames <- as.character(seqnames(utr))
#  ranges <- as.data.frame(ranges(utr))
#  strand <- as.character(strand(utr))
#  type <- mcols(utr)$type
#  ID <- unlist(mcols(utr)$Parent)
  
#  utr_bed <- data.frame(
#    seqnames = seqnames,
#    start = ranges$start,
#    end = ranges$end,
#    ID = ID,
#    type = type,
#    strand = strand,
#    width = ranges$width)
  ###Use \\. to represent the literal dot
#  utr_bed$ID <- unlist(lapply(strsplit(utr_bed$ID,"\\."),'[[',1))
#  return(utr_bed)
#}

#utr_bed <- function_gff_to_utr_bed(gff)
#utr5_bed1 <- utr_bed%>%dplyr::filter(type=="five_prime_UTR")
#utr3_bed1 <- utr_bed%>%dplyr::filter(type=="three_prime_UTR")
#utr5_bed <- data.frame(V1=utr5_bed1$seqnames, 
#                       V2=utr5_bed1$start,
#                       V3=utr5_bed1$end,
#                       V4=utr5_bed1$ID,
#                       V5=utr5_bed1$width,
#                       V6=utr5_bed1$strand)

#utr3_bed <- data.frame(V1=utr3_bed1$seqnames, 
#                       V2=utr3_bed1$start,
#                       V3=utr3_bed1$end,
#                       V4=utr3_bed1$ID,
#                       V5=utr3_bed1$width,
#                       V6=utr3_bed1$strand)

#options(bedtools.path = "/opt/homebrew/bin/")
#intersect_bed_utr5 <- bedtoolsr::bt.intersect(a = TTAA_sites, b = utr5_bed, u = T)
#intersect_bed_utr3 <- bedtoolsr::bt.intersect(a = TTAA_sites, b = utr3_bed, u = T)

#intersect_bed_utr5$TTAA_ID <- paste(intersect_bed_utr5$V1,":",intersect_bed_utr5$V2)
#intersect_bed_utr3$TTAA_ID <- paste(intersect_bed_utr3$V1,":",intersect_bed_utr3$V2)

#TTAA_sites_index <- transposon_count_matrix%>%dplyr::select(Chrom, Site)
#TTAA_sites_index$TTAA_ID <- paste(TTAA_sites_index$Chrom,":",TTAA_sites_index$Site)

#transposon_count_matrix_WT$Assigned_location[TTAA_sites_index$TTAA_ID%in%intersect_bed_utr5$TTAA_ID] <- "UTR5"
#transposon_count_matrix_WT$Assigned_location[TTAA_sites_index$TTAA_ID%in%intersect_bed_utr3$TTAA_ID] <- "UTR3"
#table(transposon_count_matrix_WT$Assigned_location)

###gtf file contains no isoforms###
#input gtf file in order to calculate transcript length and CDS length
gtf <- read.table('./Input/Genome/PlasmoDB-63_Pfalciparum3D7.gtf', sep = '\t')
Total_gene_ls <-  gsub(' ', '', gsub(';', '', lapply(strsplit(gtf$V9, 'gene_id'), '[[', 2)))
length(unique(na.omit(Total_gene_ls))) # In total, there are 5720 unique geneID in Pf3D7 v63

######Including all transcript or the first transcript length to represent#####
###mode==1: only include XXXX.1 transcript
###mode==2: include all transcript length

trans_len <- function(gtf, mode){
  #extract transcript to get transcript length
  transcipt <- gtf %>% dplyr::filter(V3 == 'transcript') 
  #transcipt$V9 <- gsub(' ', '', gsub(';', '', lapply(strsplit(transcipt$V9, 'gene_id'), '[[', 2)))
  transcipt$V9 <- gsub('transcript_id ', '', gsub(';', '', lapply(strsplit(transcipt$V9, 'gene_id'), '[[', 1)))
  transcipt <- transcipt %>% dplyr::mutate(V10=(V5-V4+1))
  Total_transciptlength <- data.frame(geneID = transcipt$V9,
                                      Total.transcipt.length = transcipt$V10)
  Total_transciptlength$geneID <- gsub(" ","",Total_transciptlength$geneID)
  if (mode==1){
    Total_transciptlength2 <- Total_transciptlength %>% dplyr::filter(grepl("\\.1", geneID))
    Total_transciptlength2$geneID <- gsub('\\.1', '',Total_transciptlength2$geneID)
    return(Total_transciptlength2)
  }else{
    return(Total_transciptlength)
  }
}

Total_transciptlength <- trans_len(gtf, mode=1) 
dim(Total_transciptlength)
###In total, 5720 transcript in Pf genome###


###mode==1: only include XXXX.1 CDS length
###mode==2: include all transcript CDS length
###!!! Some genes such as pseudogene has no CDS annotation
####!!!To define exon or CDS in very tricky in Pf3D7, exon has 5720 unique rows, but CDS has 5318 unique rows
exon_len <- function(gtf, mode){
  #####For those genes have >1 exons 
  exon <- gtf %>% dplyr::filter(V3 == 'CDS')
  ###! those genes have only 1 exon has no CDS annotation
  #exons <- gtf %>% dplyr::filter(V3 == 'exon')
  #exon_counts <- exons %>%group_by(V9) %>%summarise(exon_count = n())
  exon$V9 <- gsub('transcript_id ', '', gsub(';', '', lapply(strsplit(exon$V9, 'gene_id'), '[[', 1)))
  #length of CDS
  exon <- exon %>% mutate(V10 = (V5-V4 + 1))
  exon$V10 <- as.numeric(exon$V10)
  #convert data frame to data table 
  df <- setDT(exon)
  #find sum of observed insertions with respect to each gene by the data.table package
  Total_exonlength <- df[ ,list(sum=sum(V10)), by=V9]
  colnames(Total_exonlength)[1] <- 'geneID'
  colnames(Total_exonlength)[2] <- 'Total.CDS.length'
  Total_exonlength$geneID <- gsub(' ','',Total_exonlength$geneID)
  if (mode==1){
    Total_exonlength2 <- Total_exonlength%>% dplyr::filter(grepl("\\.1", geneID))
    Total_exonlength2$geneID <- gsub('\\.1', '',Total_exonlength2$geneID)
    return(Total_exonlength2) 
  }else{
    return(Total_exonlength) 
  }
  #####For those genes have only 1 exon will have no CDS annotation
  #exon_tb <- gtf %>% dplyr::filter(V3 == 'exon')
  #exon_tb$V9 <- gsub('transcript_id ', '', gsub(';', '', lapply(strsplit(exon_tb$V9, 'gene_id'), '[[', 1)))
  #exon_tb2 <- exon_tb %>% dplyr::filter(grepl("\\.1", V9))
  #exon_tb2$V9 <- gsub('\\.1', '',exon_tb2$V9)
  #length of CDS
  #exon_tb2 <- exon_tb2 %>% mutate(V10 = (V5-V4 + 1))
  #exon_tb2$V10 <- as.numeric(exon_tb2$V10)
  #####To filter out genes has only 1 exon
  #exon_tb3 <- exon_tb2%>%group_by(V9)%>%summarise(V9=V9,V10=V10,V11=n())
  #exon_tb4 <-distinct(exon_tb3)
  #exon_tb4 <- exon_tb4%>%dplyr::filter(V11==1)
  
}

Total_exonlength <- exon_len(gtf, mode=1)
dim(Total_exonlength)
#######5318 genes#########
Total_exonlength_pcgenes <- Total_exonlength %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
dim(Total_exonlength_pcgenes) 
#######5285 protein-coding genes without API/MIT genes#########

Pf_MIS_MFS <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/MIS_MFS_Pf.xlsx", startRow=2)
setdiff(Pf_MIS_MFS$Gene_ID, Total_exonlength$geneID)

get_insertedgeneFreq <- function(tmp){
  present_yes_matrix <- tmp
  sense_geneID_ls_yes <- present_yes_matrix %>% dplyr::filter(sense_geneID !='NA')
  sense_geneID_ls_yes <- sense_geneID_ls_yes$sense_geneID
  antisen_geneID_ls_yes <-present_yes_matrix %>% dplyr::filter(antisen_geneID !='NA') 
  antisen_geneID_ls_yes <- antisen_geneID_ls_yes$antisen_geneID
  #just to append all the geneID in count matrix and to see the number of unique geneIDs
  gene_ID_ls_yes <- append(sense_geneID_ls_yes,antisen_geneID_ls_yes) 
  return(gene_ID_ls_yes)
}

####PF3D7_1031400 has wrong annotation on UTR and exons#######
theo_TTAA_CDS <- function(transposon_count_matrix){
  #extract exon.transposon.matrix from transposon_count_matrix
  exon.transposon.matrix <-transposon_count_matrix %>% dplyr::filter(Assigned_location == 'exon')
  #5343 genes covered with in exon.transposon.matrix
  gene_ID_ls <- get_insertedgeneFreq(exon.transposon.matrix)
  length(unique(gene_ID_ls))
  #Theoretically, genes with TTAA 
  gene_ID_ls.total <- get_insertedgeneFreq(transposon_count_matrix)
  length(unique(gene_ID_ls.total))
  
  #42 genes has only TTAA in their introns
  genelist.has.only.intron.TTAA <- setdiff(unique(gene_ID_ls.total),unique(gene_ID_ls))
  length(genelist.has.only.intron.TTAA)
  #############################################
  Theo.insertion.each.gene <- as.data.frame(table(gene_ID_ls))
  colnames(Theo.insertion.each.gene)[1] <- 'geneID'
  colnames(Theo.insertion.each.gene)[grep("Freq",colnames(Theo.insertion.each.gene))] <- 'Theo.num.unique.insertions'
  Theo.insertion.each.gene$Theo.num.unique.insertions[is.na(Theo.insertion.each.gene$Theo.num.unique.insertions)] <- 0 ########only in exons(no introns)
  return(Theo.insertion.each.gene)
}

####Need to double check theoretical TTAA numbers within CDS in IGV####
Theo.insertion.each.gene <- theo_TTAA_CDS(transposon_count_matrix_WT)
nrow(Theo.insertion.each.gene)

Total_samples <- 44 #all WT samples after Day14:202510

#Total_samples <- 10 #all treated samples after Day15

essential_geneslist <- read.xlsx("./Output/background_genelist_54outof65curated.xlsx")
essential_geneslist <- data.frame(V1=essential_geneslist$GeneID)
nrow(essential_geneslist)
total.product.Pk <- read.csv("./Input/5720_total_Pf_product_description.csv")
#TTAAhits.greater.than.transcript99.TTAA.ID <- read.xlsx("./Output/TTAA_greater_than_99transcript/Pf.TTAAhits.greater.than.transcript99.TTAA.ID.xlsx")
#nrow(TTAAhits.greater.than.transcript99.TTAA.ID)
TTAAhits.greater.than.transcript99.TTAA.ID <- read.xlsx("./Output/TTAA_greater_than_99transcript/Pf.TTAAhits.greater.than.cds99.TTAA.ID.xlsx")
nrow(TTAAhits.greater.than.transcript99.TTAA.ID)

MIS <- function(transposon_count_matrix, Total_samples, Total_exonlength, Total_transciptlength, Theo.insertion.each.gene, essential_geneslist, total.product.Pk, TTAAhits.greater.than.transcript99.TTAA.ID){
  transposon_count_matrix$Total <- rowSums(transposon_count_matrix%>%dplyr::select(contains('TPN')))
  ####optional: turn counts into integers to further remove noise####
  #transposon_count_matrix$Total <- round(transposon_count_matrix$Total)
  
  ####This step is important, need to recalculate Present_in_any_samples based on Total column since we may extracted different samples for calculation
  transposon_count_matrix <- transposon_count_matrix%>% mutate (Present_in_any_samples=ifelse(Total>0, "yes","no"))
  present_yes_matrix <- transposon_count_matrix %>% dplyr::filter(Present_in_any_samples == 'yes')
  exon_matrix <- present_yes_matrix %>% dplyr::filter(Assigned_location == 'exon'|Assigned_location == 'exon&intron')
  ####Rule out 99% transcript TTAA sites
  exon_matrix$TTAA.ID <- paste(exon_matrix$Chrom, exon_matrix$Site, sep = ":")
  exon_matrix.modified <- left_join(exon_matrix, TTAAhits.greater.than.transcript99.TTAA.ID, by = 'TTAA.ID', relationship = "many-to-many")
  table(exon_matrix.modified$greater.than.transcript99) #2350
  #nrow(exon_matrix) - 1005 
  nrow(exon_matrix.modified) - 2350 
  #filter out '1' labelled rows on greater.than.transcript99 column
  exon_matrix <- exon_matrix.modified[is.na(exon_matrix.modified$greater.than.transcript99), ]
  nrow(exon_matrix)
  
  length(unique(append(unique(exon_matrix$sense_geneID), unique(exon_matrix$antisense_geneID))))
  present.exon.no99.TTAA.ID <- data.frame(TTAA.ID = exon_matrix$TTAA.ID,
                                          TTAA.ID.yes = rep(1, length(exon_matrix$TTAA.ID)))
  
  df2 <- setDT(exon_matrix)
  ##To extract sense_geneID, antisen_geneID, Total only
  df2 <- df2%>%dplyr::select(sense_geneID, antisen_geneID, Total)
  
  #find sum of observed insertions for each gene
  ob.gene.insertions.sense <- df2[ ,list(sum=sum(Total)), by=sense_geneID]
  ob.gene.insertions.sense <- ob.gene.insertions.sense %>% dplyr::filter(sense_geneID != 'NA')
  colnames(ob.gene.insertions.sense)[1] <- 'geneID'
  ob.gene.insertions.antisense <- df2[ ,list(sum=sum(Total)), by=antisen_geneID]
  ob.gene.insertions.antisense <- ob.gene.insertions.antisense %>% dplyr::filter(antisen_geneID != 'NA')
  colnames(ob.gene.insertions.antisense)[1] <- 'geneID'
  
  #merge two table
  ob.gene.insertions <- rbind(ob.gene.insertions.sense, ob.gene.insertions.antisense)
  
  ###########################Then, merge the CDS(exon)length and theoretical insertions sites within exons (no introns)
  Total.df <- left_join(Total_exonlength, Theo.insertion.each.gene, by = 'geneID')
  #calculate the theo TTAA density by number of theo TTAA sites of specific genes * 1000/ corresponding length of CDS of specific genes
  #the TTAA density of the gene g, which is calculated as the number of theoretical TTAA per kb of the CDS
  Total.df$Theo.TTAA.density <- (Total.df$Theo.num.unique.insertions* 1000)/Total.df$Total.CDS.length
  Total.df <- left_join(Total.df, Total_transciptlength, by = 'geneID')
  
  
  Total.df <- left_join(Total.df, ob.gene.insertions, by = "geneID")
  #colnames(Total.df)[13] <- 'sum.observed.insertions'
  colnames(Total.df)[ncol(Total.df)] <- 'sum.observed.insertions'
  Total.df$sum.observed.insertions[is.na(Total.df$sum.observed.insertions)] <- 0
  ###!!!should remove CDS length, this is for modified MIS
  ###!!!should remove CDS length, this is for modified MIS
  ###!!!should remove CDS length, this is for modified MIS
  Total.df$Mg <- log10((Total.df$sum.observed.insertions + 1)/(Total.df$Theo.TTAA.density))
  #Total.df$Mg <- log10((Total.df$sum.observed.insertions + 1)/(Total.df$Theo.num.unique.insertions))
  ################################input validated essential genes list from Brendan
  #nrow(unique(essential_geneslist))
  
  essential_geneslist <- as.data.frame(unique(essential_geneslist$V1))
  colnames(essential_geneslist)[1] <- 'geneID'
  
  #essential_geneslistMg  extraction
  essential_geneslistMg <- left_join(essential_geneslist, Total.df, by = "geneID")
  ##if has already done the QC for the gold plus essential genes##
  cutoff <- max(essential_geneslistMg$sum.observed.insertions)+1
  
  #optional:
  #cutoff <- quantile(essential_geneslistMg$sum.observed.insertions)[[4]] ####75% percentile
  essential_geneslistMg.filtered <- essential_geneslistMg %>% dplyr::filter(sum.observed.insertions < cutoff)
  Outliers <- setdiff(essential_geneslistMg$geneID, essential_geneslistMg.filtered$geneID)
  
  #label background validated essential genes as 3, the  outliers of validated essential genes are 2
  Total.df$background <- ifelse(Total.df$geneID %in% essential_geneslist$geneID, 2, 1)
  Total.df$background <- ifelse(Total.df$geneID %in% essential_geneslistMg.filtered$geneID, 3, Total.df$background)
  
  #normalization
  u=mean(essential_geneslistMg.filtered$Mg)
  s=sd(essential_geneslistMg.filtered$Mg)
  
  #plot(Total.df$Mg)
  Total.df$new_score <- (Total.df$Mg-u)/s - max((essential_geneslistMg.filtered$Mg-u)/s)
  
  total.product.Pk <- total.product.Pk[,c(1,3)]
  colnames(total.product.Pk)[1] <- 'geneID'
  Total.df <- left_join(Total.df, total.product.Pk, by = "geneID",relationship = "many-to-many")%>%distinct()
  Total.df <- Total.df[order(Total.df$new_score), ] #no problem here, the sum.observed.insertions=0 is just due to the order of genes changed
  
  #should remove Missing values on new score, since normalmixEM function can not have missing value
  Total.df <-Total.df[grep('FALSE', is.na(Total.df$new_score)), ]
  gm <- normalmixEM(Total.df$new_score, k=2, lambda=c(0.5,0.5))
  
  ## take a look at the recovered values
  mu1_hat<- gm$mu[1]
  mu2_hat <- gm$mu[2]
  sigma1_hat <- gm$sigma[1]
  sigma2_hat <- gm$sigma[2]
  
  ## Now recover probability of each point in the vector x, comming from the first distibution
  post.probs <- gm$posterior[,1]
  #plot(sort(post.probs))
  #################################################################################plot
  # filter out genes without TTAA
  idx=sort(Total.df$new_score, index.return = TRUE)
  Total.df$MIS <-sort(post.probs)
  Total.df$geneIndex <- idx$ix
  print(mu1_hat)
  print(mu2_hat)
  print(sigma1_hat)
  print(sigma2_hat)
  print(gm$lambda[1])
  print(gm$lambda[2])
  return(Total.df)
}

set.seed(000001)
Total.df <- MIS(transposon_count_matrix_WT, Total_samples, Total_exonlength, Total_transciptlength, Theo.insertion.each.gene, essential_geneslist, total.product.Pk, TTAAhits.greater.than.transcript99.TTAA.ID)
#Total.df <- Total.df %>% dplyr::filter(Theo.num.unique.insertions>0)
Total.df2 <-Total.df
#write.xlsx(Total.df, "./Output/MIS/MIS_readyforplot_bgremoved_WT_afterDay14_final.xlsx", na.string='NA', keepNA=F)
write.xlsx(Total.df, "./Output/MIS/MIS_readyforplot_bgremoved_WT_afterDay14_final_removeCDS.xlsx", na.string='NA', keepNA=F)
write.xlsx(Total.df, "./Output/MIS/MIS_readyforplot_bgremoved_WT_afterDay14_final_cpm.xlsx", na.string='NA', keepNA=F)
#write.xlsx(Total.df, "./Output/MIS/MIS_readyforplot_bgremoved_MEV_afterDay15_final.xlsx", na.string='NA', keepNA=F)

#Total.df2 <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_WT_afterDay14_final.xlsx")
Total.df2 <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_WT_afterDay14_final_removeCDS.xlsx")
Total.df2 <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_WT_afterDay14_final_cpm.xlsx")
#Total.df2 <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_MEV_afterDay15_final.xlsx")
######################################
p_MSg_dis <- ggplot(Total.df2, aes(x = new_score)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "grey", binwidth = .2) + 
  geom_density(lwd = 1.2,
               linetype = 1,
               colour = 4,
               fill =4,
               alpha = 0.25,
               bw = .4) + theme_bw() + labs(x = "MSg", y="Density") +
  ggtitle('Distribution of MSg after normalization')

p_MSg_dis + theme(
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(colour = NA, fill = "transparent"), legend.text = element_text(size=12))+xlim(-10,10)+ylim(0,0.3)
####MIS(no cpm normalization)
mu1_hat=-1.737362
mu2_hat=1.00647
sigma1_hat=1.479701
sigma2_hat=1.128931
p1=0.7186227
p2=0.2813773
####MIS remove CDS(no cpm normalization)
mu1_hat=-3.172892
mu2_hat=-0.5111282
sigma1_hat=0.9614906
sigma2_hat=1.037281
p1=0.5085244
p2=0.4914756
####MIS cpm without removing CDS()
mu1_hat=-6.730373
mu2_hat=-0.9417126
sigma1_hat=0.6884221
sigma2_hat=1.587059
p1=0.04249076
p2=0.9575092


####Raw countmatrix WT 35 samples after Day15:202508
mu1_hat=-1.331092
mu2_hat=3.043834
sigma1_hat=2.124663
sigma2_hat=1.940507
p1=0.6473502
p2=0.3526498
####Raw countmatrix Mev 34 samples after Day15:202508
mu1_hat=-2.429206
mu2_hat= 0.1401968
sigma1_hat=1.163729
sigma2_hat=1.387114
p1=0.413863
p2=0.586137

####Raw countmatrix 11 samples after Day9
mu1_hat=-2.083014
mu2_hat=1.744716
sigma1_hat=1.49224
sigma2_hat=2.325121
p1=0.4190593
p2=0.5809407
####Raw countmatrix 4 Day20 samples after Day9
mu1_hat=-2.404878
mu2_hat=1.072718
sigma1_hat=1.347732
sigma2_hat=3.116712
p1=0.5099755
p2=0.4900245
####Raw countmatrix 4 Day20 samples after Day9
mu1_hat=-2.533779
mu2_hat=0.3353212
sigma1_hat=0.7148273
sigma2_hat=2.101922
p1=0.2811103
p2=0.7188897
####Raw countmatrix 33 samples excluding UTR sites
mu1_hat=-2.572495
mu2_hat=0.5305749
sigma1_hat=1.590832
sigma2_hat=2.153502
p1=0.3783098
p2=0.6216902

###loess normalization

####Excluding samples before Day8(inclusive)
mu1_hat=-3.370847
mu2_hat=0.1489451
sigma1_hat=0.9521778
sigma2_hat=1.934921
p1=0.1169401
p2=0.8830599

df.EM1 <- data.frame(y=dnorm(x=seq(-10,10, 0.01), mu1_hat, sigma1_hat)* p1, 
                     x=seq(-10,10, 0.01))
df.EM2 <- data.frame(y=dnorm(x=seq(-10,10, 0.01), mu2_hat, sigma2_hat)* p2, 
                     x=seq(-10,10, 0.01))

#corrected MSg distribution
p2 <- ggplot(Total.df2, aes(x = new_score)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "grey", binwidth = .2) + 
  geom_density(lwd = 1.2,
               linetype = 1,
               colour = 4,
               fill =4,
               alpha = 0.25,
               bw = .4) + theme_bw() + labs(x = "MSg", y="Density") +
  geom_line(data = df.EM1, aes(x = df.EM1$x, y=df.EM1$y,color = "Em.guassian1"), lty=4, lwd = 1.2, show.legend = FALSE) +
  geom_line(data = df.EM2, aes(x = df.EM2$x, y=df.EM2$y,color = "Em.guassian2"), lty=4, lwd = 1.2, show.legend = FALSE) +
  scale_colour_manual("", 
                      breaks = c("Original distribution", "Em.guassian1", "Em.guassian2"),
                      values = c("Original distribution"=4, "Em.guassian1"="#8ECFC9", 
                                 "Em.guassian2"="#FA7F6F")) +
  ggtitle('')

p2 + theme(
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 16),
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.25, 0.8), 
  legend.key = element_rect(colour = NA, fill = "transparent"), legend.text = element_text(size=12))+xlim(-10,10)+ylim(0,0.3)


####Portrait, 4 X 4 inches
#################################################


###########################################ready for plot###########################################
###########################################ready for plot###########################################
###########################################ready for plot###########################################
###########################################ready for plot###########################################
###########################################ready for plot###########################################
###########################################ready for plot###########################################
#Total.df2 <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_WT_afterDay14_final.xlsx")
Total.df2 <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_WT_afterDay14_final_removeCDS.xlsx")
#Total.df2 <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_WT_afterDay14_final_cpm.xlsx")

Pf.MIS.JA <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/MIS_MFS_Pf.xlsx", startRow=2)
Pf.MIS.JA <- Pf.MIS.JA%>%select(Gene_ID, MIS)
Pf.MIS.JA <- Pf.MIS.JA[order(Pf.MIS.JA$MIS), ]
Pf.MIS.JA$geneIndex <- 1:nrow(Pf.MIS.JA)
colnames(Pf.MIS.JA)[1] <- "geneID"
Total.df2 <- Pf.MIS.JA
###############################################
p.MIS <- ggplot(Total.df2, aes(x=geneIndex, y= MIS)) +
  geom_point(aes(colour = MIS)) +
  labs(x = "Rank-ordered genes", y="MIS")+
  scale_colour_gradient2(low = muted("blue"), mid = "white",
                         limits = c(0, 1),
                         breaks = c(0.00, 0.50, 1.00),
                         high = "red" , midpoint = 0.5,  name = "MIS")+
  ggtitle('Mutagenesis index score (MIS)') + scale_x_continuous(breaks=seq(0, 5000, 2500))+ylim(c(0,1))

p.MIS+theme(
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
  axis.text = element_text(size = 12),  axis.title=element_text(size=14), legend.background = element_blank())+theme_cowplot() 
#############################################ggplot background


p.MIS.background <- ggplot(Total.df2, aes(x=geneIndex, y= MIS)) +
  geom_point(color='grey') +
  labs(x = "Rank-ordered genes", y="MIS.SP")+
  ggtitle('Mutagenesis index score (MIS.SP)') + scale_x_continuous(breaks=seq(0, 5000, 2500)) + theme(
    plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+ylim(c(0,1))

p.MIS.background <- ggplot(Total.df2, aes(x=geneIndex, y= MIS)) +
  geom_point(color='grey') +
  labs(x = "Rank-ordered genes", y="MIS.JA")+
  ggtitle('Mutagenesis index score (MIS.JA)') + scale_x_continuous(breaks=seq(0, 5000, 2500)) + theme(
    plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+ylim(c(0,1))

#############################################ggplot, validated essential genes and outliers
colors <- c("validated essential genes" = "#bf0000", "outliers of essential genes list" = "#2e659f")

p.MIS.background + geom_point(data=Total.df2[Total.df2$background == 3, ], aes(x=geneIndex, y=MIS, color='validated essential genes'),
                              shape = 21, colour = "black", fill = "#bf0000", size = 5, stroke = 0.5) + theme(
                                plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
                                legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
                                axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+
  geom_point(data=Total.df2[Total.df2$background == 2, ], aes(x=geneIndex, y=MIS, color='outliers of essential genes list'), shape = 21, colour = "black", fill = "#2e659f", size = 5, stroke = 0.5)

goldplus <- Total.df2[Total.df$background == 2, ]
goldplus2 <- Total.df2[Total.df$background == 3, ]
goldplusAll <- rbind(goldplus, goldplus2)

#############################################Plot discrepant genes for Pf.JA.MIS and 
wronggenes <- read.xlsx("./Input/WrongGenes_PfJAMIS.xlsx")
wronggenes.df <- Total.df2[Total.df2$geneID %in% wronggenes$Gene_ID, ]
geneName.df <- wronggenes%>%select(Gene_ID,GeneName, wrongtype)
colnames(geneName.df)[1] <- "geneID"
wronggenes.df <- left_join(wronggenes.df, geneName.df, by="geneID")
wronggenes1 <- wronggenes.df%>% dplyr::filter(wrongtype=="KO.essential")
wronggenes2 <- wronggenes.df%>% dplyr::filter(wrongtype=="KO.dispensable")

Pf.MIS.SP <- p.MIS.background + 
  geom_point(data=Total.df2[Total.df2$geneID %in% wronggenes2$geneID, ], aes(x=geneIndex, y=MIS),
             shape = 21, colour = "black", fill = "#2e832c", size = 5, stroke = 0.5)+
  geom_point(data=Total.df2[Total.df2$geneID %in% wronggenes1$geneID, ], aes(x=geneIndex, y=MIS),
             shape = 21, colour = "black", fill = "red", size = 5, stroke = 0.5) +
  geom_text_repel(data=wronggenes.df, aes(x=geneIndex, y=MIS, label=GeneName, color=wrongtype), size=4, box.padding = unit(0.6, "lines"),
                  segment.linetype=2,
                  max.overlaps = Inf,
                  show.legend=F,
                  nudge_x=100,
                  nudge_y=-0.03,
                  #0 indicates left alignment, 0.5 indicates center alignment, and 1 indicates right alignment. Adjusting hjust allows you to control how the labels are positioned horizontally relative to the data points.
                  hjust = 0,
                  min.segment.length = 0,
                  force=1,
                  fontface="italic",
                  family="sans")+
  scale_color_manual(values = c("KO.essential" = "red", "KO.dispensable" = "#2e832c")) 
Pf.MIS.SP 
ggsave(filename = "./Output/Figures/13Wronggenes_PfSP_final.pdf", plot=Pf.MIS.SP, width = 4,height = 4, dpi = 300)



Pf.MIS.JA <- p.MIS.background + geom_point(data=Total.df2[Total.df2$geneID %in% wronggenes1$geneID, ], aes(x=geneIndex, y=MIS),
                                           shape = 21, colour = "black", fill = "red", size = 5, stroke = 0.5) +
  geom_point(data=Total.df2[Total.df2$geneID %in% wronggenes2$geneID, ], aes(x=geneIndex, y=MIS),
             shape = 21, colour = "black", fill = "#2e832c", size = 5, stroke = 0.5)+
  geom_text_repel(data=wronggenes.df, aes(x=geneIndex, y=MIS, label=GeneName, color=wrongtype), size=6, box.padding = unit(0.6, "lines"),
                  segment.linetype=5,
                  max.overlaps = Inf,
                  show.legend=F,
                  nudge_x=300,
                  nudge_y=0.007,
                  #0 indicates left alignment, 0.5 indicates center alignment, and 1 indicates right alignment. Adjusting hjust allows you to control how the labels are positioned horizontally relative to the data points.
                  hjust = 0,
                  min.segment.length = 0,
                  force=1,
                  fontface="italic",
                  family="sans")+
  scale_color_manual(values = c("KO.essential" = "red", "KO.dispensable" = "#2e832c"))

ggsave(filename = "./Output/Figures/13Wronggenes_PfJA_v2.pdf", plot=Pf.MIS.JA, width = 4,height = 4, dpi = 300)


#write.xlsx(Total.df2[Total.df2$geneID %in% golden_list_nonessential_Pk$V1, ], "/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Math_model_backgroundgenelist2/Pk/Pk_nonessential_golden_list.xlsx", na.string='NA', keepNA=F)
#write.xlsx(Total.df2[Total.df2$geneID %in% golden_list_nonessential_Pk2$V1, ], "/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Output/Math_model_backgroundgenelist2/Pk/Pk_nonessential_golden_list2_curated.xlsx", na.string='NA', keepNA=F)
#############################################check any gene in Pk only
Gene.name= 'PKNH_1032300'
Gene.Index <- Total.df2[grep(Gene.name, Total.df2$geneID),]$geneIndex
Gene.MIS <- Total.df2[grep(Gene.name, Total.df2$geneID),]$MIS
#Total.df2[grep(Gene.name, Total.df2$geneID),]
p.MIS.background + geom_point(aes(x=Gene.Index, y=Gene.MIS), shape = 21, colour = "black", fill = "#2e832c", size = 5, stroke = 1) +ggtitle(paste('MIS:', Gene.name))+
  geom_point(aes(x=300, y=1), shape = 21, colour = "black", fill = "#2e832c", size = 3, stroke = 1)+
  geom_text(aes(x = 1200, y = 1, label = Gene.name))


