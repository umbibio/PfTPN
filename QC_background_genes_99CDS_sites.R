library(tidyverse)
library(IRanges)
library(Biostrings) 
library(ShortRead)
library(openxlsx)
library(dplyr)
library(data.table)
library(rtracklayer)
library(parallel)
library(stringr)
######Calculate the sites locating at 99% of the transcript or 99% of the CDS

####Input should be exon conversion with overlaps followed at tails and raw counts
cm_converted <- read.xlsx("./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon_UTR_annotation_all.xlsx") ###Double check if includes total columns

######
gff <- read.table('./Input/Genome/PlasmoDB-63_Pfalciparum3D7.gff', sep = '\t')
gtf <- read.table('./Input/Genome/PlasmoDB-63_Pfalciparum3D7.gtf', sep = '\t')
Total_gene_ls <-  gsub(' ', '', gsub(';', '', lapply(strsplit(gtf$V9, 'gene_id'), '[[', 2)))
length(unique(Total_gene_ls)) # In total, there are unique 5270 gene ID in Pf v63 genome

#####There are mutiple genes isoforms are annotated in Pf v63 genome, to make it simple, only extract gene.1 transcript info
trans_len <- function(gtf){
  #extract transcript to get transcript length
  transcipt <- gtf %>% dplyr::filter(V3 == 'transcript')
  ####Only extract gene1.1 transcript
  transcipt <- transcipt %>% filter(str_detect(V9, "\\.1"))
  
  transcipt$V9 <- gsub(' ', '', gsub(';', '', lapply(strsplit(transcipt$V9, 'gene_id'), '[[', 2)))
  transcipt <- transcipt %>% dplyr::mutate(V10=(V5-V4+1))
  Total_transciptlength <- data.frame(geneID = transcipt$V9,
                                      Total.transcipt.length = transcipt$V10)
  return(Total_transciptlength)
}

Total_transciptlength <- trans_len(gtf) ###5270 unique genes

#exon_len <- function(gtf){
#  ####Only extract gene1.1 exons
#  gtf <- gtf %>% filter(str_detect(V9, "\\.1"))
  
#  exon <- gtf %>% dplyr::filter(V3 == 'exon')
#  exon$V9 <- gsub(' ', '', gsub(';', '', lapply(strsplit(exon$V9, 'gene_id'), '[[', 2)))
  #length of CDS
#  exon <- exon %>% mutate(V10 = (V5-V4 + 1))
#  exon$V10 <- as.numeric(exon$V10)
  
  #convert data frame to data table 
#  df <- setDT(exon)
  #To calculate the sum of all exons with respect to each gene by the data.table package
#  Total_exonlength <- df[ ,list(sum=sum(V10)), by=V9]
#  colnames(Total_exonlength)[1] <- 'geneID'
#  colnames(Total_exonlength)[2] <- 'Total.CDS.length'
#  dim(Total_exonlength) #5502, contains all the genes
#  return(Total_exonlength)
#}
#Total_exonlength <- exon_len(gtf)

##############To get CDS length for each gene#############################
###mode==1: only include XXXX.1 CDS length
###mode==2: include all transcript CDS length
###!!! Some genes such as pseudogene has no CDS annotation
exon_len <- function(gtf, mode){
  exon <- gtf %>% dplyr::filter(V3 == 'exon')
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
}
Total_exonlength <- exon_len(gtf, mode=1)
######rule out those theretical TTAA sites located at >99% CDS
######rule out those theretical TTAA sites located at >99% CDS
######rule out those theretical TTAA sites located at >99% CDS

#######input bed file of TTAA bed file overlapped with gtf
TTAAhits_R_gtf_include_contigs <- read.table('./Output/TTAAhitsPf/TTAAhits_R_gtf_include_contigs.bed', sep = '\t')
######modified TTAA loci only
TTAA_R <- read.table('./Output/TTAAhitsPf/PF3D7_theo_TTAA_R_modified_final.bed', sep = '\t')

#TTAA_sites_trans_tails <- function(TTAAhits_R_gtf_include_contigs, TTAA_R, trans_tails, Total_transciptlength, Total_exonlength, Total_theo_TTAA){
#  TTAA_R$V7 <- paste(TTAA_R$V1 , paste(TTAA_R$V2, TTAA_R$V3, sep= "-"), sep = ":")
  
#  TTAAhits_R_gtf_include_contigs$V16 <- paste(TTAAhits_R_gtf_include_contigs$V1 , paste(TTAAhits_R_gtf_include_contigs$V2, TTAAhits_R_gtf_include_contigs$V3, sep= "-"), sep = ":")
#  TTAAhits_R_gtf_include_contigs_trans <- TTAAhits_R_gtf_include_contigs %>% dplyr::filter(V9 == 'transcript')
#  TTAAhits_R_gtf_include_contigs_trans$V15 <- gsub(';', '', unlist(lapply(strsplit(TTAAhits_R_gtf_include_contigs_trans$V15, '; gene_id '), '[[', 2)))
#  TTAAhits_R_gtf_include_contigs_trans_filtered <- TTAAhits_R_gtf_include_contigs_trans %>% dplyr::filter(V6 == V13)
#  TTAAhits_R_gtf_include_contigs_trans_filtered <- TTAAhits_R_gtf_include_contigs_trans_filtered %>% dplyr::mutate(V17=ifelse(V6 == '+', V2-V10, V11-V3))
  #!!!!!!!!if one genes intron is too long, then TTAA sites located at >99% CDS is not the case we want to rule out. so the criteria should be >99% of transcript length
#  Total_transciptlength$cutoff.very.end.distance <- Total_transciptlength$Total.transcipt.length * trans_tails
  
#  colnames(TTAAhits_R_gtf_include_contigs_trans_filtered)[15] <- 'geneID'
#  TTAAhits_R_gtf_include_contigs_trans_filtered.merged <- left_join(TTAAhits_R_gtf_include_contigs_trans_filtered, Total_exonlength, by = 'geneID')
#  TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total <- left_join(TTAAhits_R_gtf_include_contigs_trans_filtered.merged, Total_transciptlength, by = 'geneID')
  #label those therectical TTAA sites located at >99% transcript
#  TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total$greater.than.transcript99 <- ifelse(TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total$V17 > TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total$cutoff.very.end.distance, 1, 0)
#  table(TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total$greater.than.transcript99)
#  TTAAhits.greater.than.transcript99 <- TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total %>% dplyr::filter(greater.than.transcript99 == '1')
#  print(nrow(TTAAhits.greater.than.transcript99)) #1072
#  print(nrow(TTAAhits.greater.than.transcript99)/Total_theo_TTAA) # 0.669% theoretical TTAA sites located at > 99% of transcript length
  
#  colnames(TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total)[17] <- "Distance.to.TSS"
#  colnames(TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total)[16] <- "TTAA.ID"
  
#  colnames(TTAAhits.greater.than.transcript99)[16] <- 'TTAA.ID'
#  colnames(TTAAhits.greater.than.transcript99)[17] <- "Distance.to.TSS"
#  TTAAhits.greater.than.transcript99$TTAA.ID <- unlist(lapply(strsplit(TTAAhits.greater.than.transcript99$TTAA.ID, '-'), '[[', 1))
#  list_of_dfs <- list(TTAAhits_R_gtf_include_contigs_trans_total = TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total, 
#                      TTAAhits.greater.than.transcript99 = TTAAhits.greater.than.transcript99)
#  return(list_of_dfs)
#}

#################To remove sites at 99% of CDS rather than transcripts#################
TTAA_sites_CDS_tails <- function(TTAAhits_R_gtf_include_contigs, TTAA_R, CDS_tails, Total_transciptlength, Total_exonlength, Total_theo_TTAA){
  TTAA_R$V7 <- paste(TTAA_R$V1 , paste(TTAA_R$V2, TTAA_R$V3, sep= "-"), sep = ":")
  
  TTAAhits_R_gtf_include_contigs$V16 <- paste(TTAAhits_R_gtf_include_contigs$V1 , paste(TTAAhits_R_gtf_include_contigs$V2, TTAAhits_R_gtf_include_contigs$V3, sep= "-"), sep = ":")
  TTAAhits_R_gtf_include_contigs_trans <- TTAAhits_R_gtf_include_contigs %>% dplyr::filter(V9 == 'transcript')
  TTAAhits_R_gtf_include_contigs_trans$V15 <- gsub(';', '', unlist(lapply(strsplit(TTAAhits_R_gtf_include_contigs_trans$V15, '; gene_id '), '[[', 2)))
  TTAAhits_R_gtf_include_contigs_trans_filtered <- TTAAhits_R_gtf_include_contigs_trans %>% dplyr::filter(V6 == V13)
  TTAAhits_R_gtf_include_contigs_trans_filtered <- TTAAhits_R_gtf_include_contigs_trans_filtered %>% dplyr::mutate(V17=ifelse(V6 == '+', V2-V10, V11-V3))
  #!!!!!!!!if one genes intron is too long, then TTAA sites located at >99% CDS is not the case we want to rule out. so the criteria should be >99% of transcript length
  Total_transciptlength$cutoff.very.end.distance <- Total_transciptlength$Total.transcipt.length * CDS_tails
  
  colnames(TTAAhits_R_gtf_include_contigs_trans_filtered)[15] <- 'geneID'
  TTAAhits_R_gtf_include_contigs_trans_filtered.merged <- left_join(TTAAhits_R_gtf_include_contigs_trans_filtered, Total_exonlength, by = 'geneID')
  TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total <- left_join(TTAAhits_R_gtf_include_contigs_trans_filtered.merged, Total_transciptlength, by = 'geneID')
  #label those therectical TTAA sites located at >99% transcript
  TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total$greater.than.transcript99 <- ifelse(TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total$V17 > TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total$cutoff.very.end.distance, 1, 0)
  table(TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total$greater.than.transcript99)
  TTAAhits.greater.than.transcript99 <- TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total %>% dplyr::filter(greater.than.transcript99 == '1')
  print(nrow(TTAAhits.greater.than.transcript99)) #1072
  print(nrow(TTAAhits.greater.than.transcript99)/Total_theo_TTAA) # 0.669% theoretical TTAA sites located at > 99% of transcript length
  
  colnames(TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total)[17] <- "Distance.to.TSS"
  colnames(TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total)[16] <- "TTAA.ID"
  
  colnames(TTAAhits.greater.than.transcript99)[16] <- 'TTAA.ID'
  colnames(TTAAhits.greater.than.transcript99)[17] <- "Distance.to.TSS"
  TTAAhits.greater.than.transcript99$TTAA.ID <- unlist(lapply(strsplit(TTAAhits.greater.than.transcript99$TTAA.ID, '-'), '[[', 1))
  list_of_dfs <- list(TTAAhits_R_gtf_include_contigs_trans_total = TTAAhits_R_gtf_include_contigs_trans_filtered.merged.total, 
                      TTAAhits.greater.than.transcript99 = TTAAhits.greater.than.transcript99)
  return(list_of_dfs)
}


list_of_dfs <- TTAA_sites_CDS_tails(TTAAhits_R_gtf_include_contigs, TTAA_R, CDS_tails=0.99,Total_transciptlength, Total_exonlength,Total_theo_TTAA=330136)
list_of_dfs$TTAAhits.greater.than.transcript99
write.xlsx(list_of_dfs$TTAAhits.greater.than.transcript99, "./Output/TTAA_greater_than_99transcript/Pf.TTAAhits.greater.than.cds99.TTAA.ID.xlsx", na.string='NA', keepNA=F)


TTAAhits.greater.than.transcript99 <- read.xlsx("./Output/TTAA_greater_than_99transcript/Pf.TTAAhits.greater.than.cds99.TTAA.ID.xlsx")

######Original essential gene list with high confidence by published data######
gene.table <- read.xlsx("./Input/ref_gene_list.xlsx", sheet = 1) #sheet1 is WT Pf api-associated essential genes
table(gene.table$Category)
#####Only extract essential gene list
gene.table0 <- gene.table[gene.table$Category=="ESSENTIAL ",]
essential_geneslist <- as.data.frame(gene.table0$Gene.ID)
colnames(essential_geneslist) <- "V1"
#make sure the gene in the list is unique
nrow(unique(essential_geneslist))

gene.table <- read.xlsx("./Input/ref_gene_list.xlsx", sheet = 3) #sheet3 is WT Pf essential genes
essential_geneslist <- as.data.frame(gene.table$Gene.ID)
colnames(essential_geneslist) <- "V1"
nrow(essential_geneslist)
#make sure the gene in the list is unique
nrow(unique(essential_geneslist)) ###68 curated essential genes in total

######recurated by Dr. Krithika, extended essential gene list########
###sheet2 is essential genes(no Apicomplast-associated genes)
gene.table <- read.xlsx("./Input/Pf_high_confidence_essential_dispensable_genelist.xlsx", sheet = 2) 
essential_geneslist <- as.data.frame(gene.table$Gene)
colnames(essential_geneslist) <- "V1"
nrow(essential_geneslist)
#make sure the gene in the list is unique
nrow(unique(essential_geneslist)) ###68 curated essential genes in total

########remove sites locates at the tail of transcript(99% of transcript)
###PKNH_0424500: have >1000 insertion at the vert end of the gene, but because the piggypac transposon insert seamlessly, if the end
#Sequence is TTAA, it will not affect the transcription of the gene
remove_TTAA_99transcript <- function(cm_converted, TTAAhits.greater.than.transcript99){
  #Notice that can only use exon_cm_converted to filter out 99% of transcript, since cm_exon_converted only modify the exon_rows, those exon-intron rows will mess up if use the total cm_converted to remove TTAA sites locating at 99% of transcript
  cm_converted_exon <- cm_converted%>%dplyr::filter(Location=='exon')
  cm_converted_exon$TTAA.ID <- paste(cm_converted_exon$Chrom, cm_converted_exon$Site,cm_converted_exon$GeneID, sep = ":")
  TTAAhits.greater.than.transcript99.TTAA.ID <- paste(TTAAhits.greater.than.transcript99$TTAA.ID, TTAAhits.greater.than.transcript99$geneID, sep = ":")
  cm_converted_exon_filtered <- cm_converted_exon %>% dplyr::filter(!(TTAA.ID %in% TTAAhits.greater.than.transcript99.TTAA.ID))
  #nrow(cm_converted_exon)-nrow(cm_converted_exon_filtered) #2122
  return(cm_converted_exon_filtered)
}

cm_converted_exon_filtered <- remove_TTAA_99transcript(cm_converted, TTAAhits.greater.than.transcript99)
cm_converted_essential_gl <- cm_converted_exon_filtered %>% dplyr::filter(GeneID %in% essential_geneslist$V1)
colnames(cm_converted_essential_gl)

#######Only extract WT columns after Day15
cm_converted_essential_gl <- cm_converted_essential_gl[,c(1:12,14,16,22,24,26,31,32,35,37,39,40,50,51,54,
                                                          57,60,61,66,68,70,72,74,78,80,82,84,86,90,92,94,96,100,
                                                          102,104,107,113,116,120,121,126,129,132,135,139,142,144,146)]
#######Only extract Mev treated columns
#cm_converted_essential_gl <- cm_converted_essential_gl%>% select(-matches("_C"),-matches("D0"))
###double check, no Mev treated samples
colnames(cm_converted_essential_gl)
cm_converted_essential_gl$Total <- rowSums(cm_converted_essential_gl[,c(13:20)]) ###13 includes gene.description column
df <- cm_converted_essential_gl%>%dplyr::select(Chrom, Site, GeneID, Total)

df_genelist <- df %>% group_by(GeneID) %>% summarise(sum.observed.counts=sum(Total))
dim(df_genelist)
quantile(df_genelist$sum.observed.counts,probs = seq(0, 1, 0.1))
#cutoff <- 55.25 ###75% quantile for WT samples only for API associated essential genes
#cutoff <- 79.25 ###75% quantile for Mev samples only for API associated essential genes

cutoff <- 13 ###75% quantile for WT samples only for essential genes


df_genelist_filtered <- df_genelist[df_genelist$sum.observed.counts<cutoff,]
nrow(df_genelist_filtered)#54
plot(sort(df_genelist$sum.observed.counts), xlab= '54 essential genes (putative)', ylab= 'Sum of observed counts', xlim=c(1, 70), pch = 16,
     cex.lab = 1.4,
     cex.axis = 1.2)
abline(h = 13, col = "red", lty=2,lwd = 3)
abline(v = nrow(df_genelist_filtered), col = "blue", lty=2,lwd = 3)

df_genelist$outlier <- ifelse(df_genelist$sum.observed.counts<cutoff,"NO","Yes")
Product.description <- read.csv("./Input/5720_total_Pf_product_description.csv") #PF3D7_1031400 repeated
Product.description <- Product.description[,c(1,3,4)]
colnames(Product.description)[1] <- "GeneID"
df_genelist <- left_join(df_genelist,Product.description, by="GeneID")
df_genelist <-df_genelist%>%distinct()
write.xlsx(df_genelist_filtered, "./Output/background_genelist_54outof65curated.xlsx", na.string='NA', keepNA=F)
write.xlsx(df_genelist, "./Output/background_genelist65_WTonly202510curated.xlsx", na.string='NA', keepNA=F)

table <- read.xlsx("./Output/background_genelist_54outof65curated.xlsx")



