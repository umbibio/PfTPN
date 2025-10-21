library(tidyverse)
library(openxlsx)
library(flexmix)
library(countreg)
library(doParallel)
library(mixtools)
library(data.table)
library(cowplot)

getwd()

####For MMIS, we choose to use weight function with a sigmoid drop at 99% of the transcript

####This script is for hybrid BM and MIS model###########test on local computer
####This script is for hybrid BM and MIS model###########test on local computer
####This script is for hybrid BM and MIS model###########test on local computer

####note: This script is just to fit a logistic regression for Gi
####The s.state showed in the final table is just for last sample of Gibbs sampling
#input exon converted count matrix with bg noise removed at site level and CPM normalization for all 75 samples for essentialome only
####This script is for hybrid BM and MIS model
###Did not need to remove 99% of the transcript
###No need to inout transposon matrix v2 with gene description

######Protein-coding transcript##########
transposon_count_matrix <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix_all_Bg_removed_siteslevel_final.xlsx")
#transposon_count_matrix <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202307_Novaseq/Tnseq202307/Output/transposon_matrix/Pk/transposon_count_matrix75essentialomeonly_run13.xlsx")


transposon_count_matrix_WT <- transposon_count_matrix[,c(1:6,10,16,18,20,25,26,29,31,33,34,44,45,48,51,54,55,60,62,64,66,
                                                         68,72,74,76,78,80,84,86,88,90,94,96,98,101,107,110,114,115,123,
                                                         126,129,136,138,140)] ### WT samples after Day14


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

#Pf_MIS_MFS <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/MIS_MFS_Pf.xlsx", startRow=2)

#setdiff(Pf_MIS_MFS$Gene_ID, Total_exonlength$geneID)

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
#########This part is to calculate TTAA metrics and apply the weight for each site based on its relative distance to TSS###################################
#########This part is to calculate TTAA metrics and apply the weight for each site based on its relative distance to TSS###################################
#########This part is to calculate TTAA metrics and apply the weight for each site based on its relative distance to TSS###################################
####We need to remove the introns and stitch all the exons together as well as the TTAA within exons
##Input modified df, which removed all the introns and concatenate all the exons for each gene
############V2, V3 are modified TTAA loci, V10 and V11 are modified transcipt loci
TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info <- read.xlsx("./Output/TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info.xlsx")
##V4 is original TTAA_ID, the unique TTAA identifier
##V2 is modified TTAA_ID after removing introns
##V10 is modified loci of start of transcript after removing introns
##V11 is modified loci of end of transcript after removing introns
##V13 is the strandness of genes
TTAA_metrics_ri <- function(TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info){
  TTAAhits_R_gtf <-TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info[,c(1:15)] 
  TTAAhits_R_gtf$V16 <- paste(TTAAhits_R_gtf$V1, TTAAhits_R_gtf$V2, sep=":")
  ############V12 is CDS length
  TTAAhits_R_gtf$V12 <- TTAAhits_R_gtf$V11-TTAAhits_R_gtf$V10+1
  TTAAhits_R_gtf_exons <- TTAAhits_R_gtf
  
  ###################metric4: To calculate relative distance to 5' end for every TTAA and ###################
  #####V17 is distance to TSS for every TTAA sites within exons
  TTAAhits_R_gtf_exons <- as.data.frame(TTAAhits_R_gtf_exons)
  #### double check with IGV, distance for sense strandness should +1, no change for antisense strandness
  TTAAhits_R_gtf_exons <- TTAAhits_R_gtf_exons%>%mutate(V17=ifelse(V13=="+", (V2-V10+1), (V11-V3)))
  #### For both sense and antisense strand, overlap of TTAA at TSS will not disrupt gene's expression theoretically
  #### For those sites distance to 5' <= 0, the w should be 0, we can turn those sites' distances into modified trans length, which will let the NW=0 in the downstream analysis
  TTAAhits_R_gtf_exons_mod <- TTAAhits_R_gtf_exons %>% mutate(V17=ifelse(V17<=0, V12, V17))
  TTAAhits_R_gtf_exons_mod$R_i <- abs(TTAAhits_R_gtf_exons_mod$V17-TTAAhits_R_gtf_exons_mod$V12)/(TTAAhits_R_gtf_exons_mod$V12)-0.5
  return(TTAAhits_R_gtf_exons_mod)
}

TTAAhits_R_gtf_exons <- TTAA_metrics_ri(TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info)
#write.xlsx(TTAAhits_R_gtf_exons,"./Output/HM/TTAAhits_R_gtf_exons_relative_loci_normalized_CDS.xlsx")

# Sigmoid function for the drop at the 99% tail
sigmoid <- function(x, midpoint, slope) {
  y <- 1 / (1 + exp(-slope * (x - midpoint)))
  return(y)
}

# Quadratic function
quadratic_function <- function(x) {
  return(-8/5 * x^2 + 2/5 * x + 1)
}

weight_function <- function(Ri){
  W_quadratic <- quadratic_function(Ri)
  # Define sigmoid parameters for the drop at the beginning
  midpoint_drop <- -0.48  # Adjust as needed
  slope_drop <- -1000     # Adjust as needed for a drop at the beginning
  # Apply sigmoid drop at the beginning
  sigmoid_drop <- sigmoid(Ri, midpoint_drop, slope_drop)
  W_with_drop <- W_quadratic * (1 - sigmoid_drop)
  return(list(W_with_drop=W_with_drop, W_quadratic=W_quadratic))
}
#########Need to make sure the Present_in_any_samples == 'yes' has been correctly performed
Total_samples <- 44
essential_geneslist <- read.xlsx("./Output/background_genelist_54outof65curated.xlsx")
essential_geneslist <- data.frame(V1=essential_geneslist$GeneID)
nrow(essential_geneslist)
total.product.Pk <- read.csv("./Input/5720_total_Pf_product_description.csv")
modified_MIS <- function(transposon_count_matrix, Total_samples, Total_exonlength, Total_transciptlength, TTAAhits_R_gtf_exons, Theo.insertion.each.gene, essential_geneslist, total.product.Pk){
  transposon_count_matrix$Total <- rowSums(transposon_count_matrix%>%dplyr::select(contains('TPN')))
  transposon_count_matrix <- transposon_count_matrix%>%mutate(Present_in_any_samples=ifelse(Total>0, "yes","no"))
  present_yes_matrix <- transposon_count_matrix %>% dplyr::filter(Present_in_any_samples == 'yes')
  exon_matrix <- present_yes_matrix  %>% dplyr::filter(Assigned_location == 'exon')
  ###overlapped exons which have insertions in the observed dataset
  overlapped_tm <- exon_matrix%>% dplyr::filter((!is.na(sense_geneID) & !is.na(antisen_geneID)))
  print(nrow(exon_matrix))
  ####58991
  
  exon_TTAA_ID <-exon_matrix[,c(1,2,4,5,ncol(exon_matrix))]
  nrow(exon_TTAA_ID)
  
  length(unique(exon_TTAA_ID$sense_geneID))+length(unique(exon_TTAA_ID$antisen_geneID))-2 #2 NAs
  ###!is.na
  exon_TTAA_ID_sense <- exon_TTAA_ID[!is.na(exon_TTAA_ID$sense_geneID),]
  exon_TTAA_ID_sense <- exon_TTAA_ID_sense%>% dplyr::mutate(V4 = paste(Chrom,Site,sep = ":"))
  ###is.na
  #exon_TTAA_ID_antisense <- exon_TTAA_ID[is.na(exon_TTAA_ID$sense_geneID),]
  exon_TTAA_ID_antisense <- exon_TTAA_ID[!is.na(exon_TTAA_ID$antisen_geneID),]
  exon_TTAA_ID_antisense <- exon_TTAA_ID_antisense%>% dplyr::mutate(V4 = paste(Chrom,Site,sep = ":"))
  
  TTAAhits_R_gtf_exons_sense <- TTAAhits_R_gtf_exons%>% dplyr::filter(V13=="+")
  TTAAhits_R_gtf_exons_sense <- TTAAhits_R_gtf_exons_sense[seq(from=1, to=nrow(TTAAhits_R_gtf_exons_sense), by=2),]
  TTAAhits_R_gtf_exons_sense$V4 <- unlist(lapply(strsplit(TTAAhits_R_gtf_exons_sense$V4,"-"),"[[",1))
  TTAAhits_R_gtf_exons_sense_d5 <- TTAAhits_R_gtf_exons_sense[,c(1:4,grep("R_i",colnames(TTAAhits_R_gtf_exons)))]
  
  TTAAhits_R_gtf_exons_antisense <- TTAAhits_R_gtf_exons%>% dplyr::filter(V13=="-")
  TTAAhits_R_gtf_exons_antisense <- TTAAhits_R_gtf_exons_antisense[seq(from=2, to=nrow(TTAAhits_R_gtf_exons_antisense), by=2),]
  TTAAhits_R_gtf_exons_antisense$V4 <- unlist(lapply(strsplit(TTAAhits_R_gtf_exons_antisense$V4,"-"),"[[",1))
  TTAAhits_R_gtf_exons_antisense_d5 <- TTAAhits_R_gtf_exons_antisense[,c(1:4,grep("R_i",colnames(TTAAhits_R_gtf_exons)))]
  ####V4 is TTAA ID/original TTAA identifier
  exon_TTAA_ID_sense <- left_join(exon_TTAA_ID_sense,TTAAhits_R_gtf_exons_sense_d5, by="V4")
  nrow(exon_TTAA_ID_sense)
  exon_TTAA_ID_antisense <- left_join(exon_TTAA_ID_antisense,TTAAhits_R_gtf_exons_antisense_d5, by="V4")
  nrow(exon_TTAA_ID_antisense)
  #28645+30347
  
  #####merge all
  #####V17 is relative distance to TSS for every site, xi
  exon_TTAA_ID_WX <- rbind(exon_TTAA_ID_sense,exon_TTAA_ID_antisense)
  nrow(exon_TTAA_ID_WX)
  #####weight column based on W(Ri)=-3Ri^2+1/2Ri+1
  #exon_TTAA_ID_WX$W <- -3 * (exon_TTAA_ID_WX$R_i)^2+(1/2)*(exon_TTAA_ID_WX$R_i)+1
  
  #####weight column based on W(Ri)=-5/8Ri^2+2/5Ri+1 passing(0.5,0.4),(0,1),(0.5,0.8)
  #########Applied weight function###########
  #########Applied weight function###########
  #########Applied weight function###########
  ##Mode1:
  #exon_TTAA_ID_WX$W <- (-8/5) * (exon_TTAA_ID_WX$R_i)^2+(2/5)*(exon_TTAA_ID_WX$R_i)+1
  ##Mode2:
  W <- weight_function(exon_TTAA_ID_WX$R_i)
  exon_TTAA_ID_WX$W <- W$W_with_drop
  #########Applied weight function###########
  #########Applied weight function###########
  #########Applied weight function###########
  #####NW column=W(Ri)*Total nomalized reads
  exon_TTAA_ID_WX$NW <- exon_TTAA_ID_WX$Total * exon_TTAA_ID_WX$W
  nrow(exon_TTAA_ID_WX)
  head(exon_TTAA_ID_WX)
  tail(exon_TTAA_ID_WX)
  #######change NW columns in order to pipe in the algorithm easily
  colnames(exon_TTAA_ID_WX)[grep("Total",colnames(exon_TTAA_ID_WX))] <- "Total_reads"
  colnames(exon_TTAA_ID_WX)[grep("NW",colnames(exon_TTAA_ID_WX))] <- "Total"
  
  #convert data frame to data table 
  df2 <- setDT(exon_TTAA_ID_WX)
  #df2 <- df2[, c(4,6,79)]
  #df2 <- df2[, c(4,6,84)]
  #df2 <- df2[, c(3,4,12)]
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
  Total.df$Mg <- log10((Total.df$sum.observed.insertions + 1)/(Total.df$Theo.num.unique.insertions))
  
  ################################input validated essential genes list from Brendan
  #nrow(unique(essential_geneslist))
  
  essential_geneslist <- as.data.frame(unique(essential_geneslist$V1))
  colnames(essential_geneslist)[1] <- 'geneID'
  
  #essential_geneslistMg  extraction
  essential_geneslistMg <- left_join(essential_geneslist, Total.df, by = "geneID")
  cutoff <- quantile(essential_geneslistMg$sum.observed.insertions)[[4]]
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
  Total.df <- left_join(Total.df, total.product.Pk, by = "geneID")
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
  Total.df$MMIS <-sort(post.probs)
  Total.df$geneIndex <- idx$ix
  print(mu1_hat)
  print(mu2_hat)
  print(sigma1_hat)
  print(sigma2_hat)
  print(gm$lambda[1])
  print(gm$lambda[2])
  return(Total.df)
}

Total.df <- modified_MIS(transposon_count_matrix, Total_samples, Total_exonlength, Total_transciptlength, TTAAhits_R_gtf_exons, Theo.insertion.each.gene, essential_geneslist, total.product.Pk)
#Total.df <- Total.df %>% dplyr::filter(Theo.num.unique.insertions>0)
write.xlsx(Total.df, "./Output/MMIS/Total_df_modified_MIS.xlsx")
Total.df2 <- read.xlsx("./Output/MMIS/Total_df_modified_MIS.xlsx")

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
  legend.key = element_rect(colour = NA, fill = "transparent"), legend.text = element_text(size=12))+xlim(-10,10)+ylim(0,0.5)
####MIS(no cpm normalization)
mu1_hat=-2.152826
mu2_hat=0.04334797
sigma1_hat=0.4623344
sigma2_hat=0.416218
p1=0.5559969
p2=0.4440031

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
  legend.key = element_rect(colour = NA, fill = "transparent"), legend.text = element_text(size=12))+xlim(-10,10)+ylim(0,0.5)

###############################################
p.MIS <- ggplot(Total.df2, aes(x=geneIndex, y= MMIS)) +
  geom_point(aes(colour = MMIS)) +
  labs(x = "Rank-ordered genes", y="MMIS")+
  scale_colour_gradient2(low = muted("blue"), mid = "white",
                         limits = c(0, 1),
                         breaks = c(0.00, 0.50, 1.00),
                         high = "red" , midpoint = 0.5,  name = "MMIS")+
  ggtitle('Mutagenesis index score (MIS)') + scale_x_continuous(breaks=seq(0, 5000, 2500))+ylim(c(0,1))

p.MIS+theme(
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
  axis.text = element_text(size = 12),  axis.title=element_text(size=14), legend.background = element_blank())+theme_cowplot()
#############################################ggplot background
p.MIS.background <- ggplot(Total.df2, aes(x=geneIndex, y= MMIS)) +
  geom_point(color='grey') +
  labs(x = "Rank-ordered genes", y="MIS.SP")+
  ggtitle('Mutagenesis index score (MIS.SP)') + scale_x_continuous(breaks=seq(0, 5000, 2500)) + theme(
    plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+ylim(c(0,1))


wronggenes <- read.xlsx("./Input/WrongGenes_PfJAMIS.xlsx")
wronggenes.df <- Total.df2[Total.df2$geneID %in% wronggenes$Gene_ID, ]
geneName.df <- wronggenes%>%dplyr::select(Gene_ID,GeneName, wrongtype)
colnames(geneName.df)[1] <- "geneID"
wronggenes.df <- left_join(wronggenes.df, geneName.df, by="geneID")
wronggenes1 <- wronggenes.df%>% dplyr::filter(wrongtype=="KO.essential")
wronggenes2 <- wronggenes.df%>% dplyr::filter(wrongtype=="KO.dispensable")

Pf.MIS.SP <- p.MIS.background + 
  geom_point(data=Total.df2[Total.df2$geneID %in% wronggenes2$geneID, ], aes(x=geneIndex, y=MMIS),
             shape = 21, colour = "black", fill = "#2e832c", size = 5, stroke = 0.5)+
  geom_point(data=Total.df2[Total.df2$geneID %in% wronggenes1$geneID, ], aes(x=geneIndex, y=MMIS),
             shape = 21, colour = "black", fill = "red", size = 5, stroke = 0.5) +
  geom_text_repel(data=wronggenes.df, aes(x=geneIndex, y=MMIS, label=GeneName, color=wrongtype), size=4, box.padding = unit(0.6, "lines"),
                  segment.linetype=2,
                  max.overlaps = Inf,
                  show.legend=F,
                  nudge_x=100,
                  nudge_y=0.00,
                  #0 indicates left alignment, 0.5 indicates center alignment, and 1 indicates right alignment. Adjusting hjust allows you to control how the labels are positioned horizontally relative to the data points.
                  hjust = 0,
                  min.segment.length = 0,
                  force=1,
                  fontface="italic",
                  family="sans")+
  scale_color_manual(values = c("KO.essential" = "red", "KO.dispensable" = "#2e832c")) 
Pf.MIS.SP 
ggsave(filename = "./Output/Figures/13Wronggenes_PfSP_final.pdf", plot=Pf.MIS.SP, width = 4,height = 4, dpi = 300)

