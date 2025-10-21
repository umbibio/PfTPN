library(tidyverse)
library(IRanges)
library(Biostrings) 
library(ShortRead)
library(openxlsx)
library(dplyr)
library(data.table)
library(rtracklayer)
library(parallel)
setwd("/Users/sidaye/Documents/R/API_TnSeq")
####This script is for generating countmatrix template and mapping reads to TTAA sites
#######################################################bedtoos intersect -a -b -wa -wb > .bed !!!must use gtf file instead of gff file
#check the total number of protein_coding_genes of gtf files are whether the same as the DB gennome Info


#######################################################create empty matrix
samples2 <- read.table('/mathspace/data01/mathbio_lab/Tnseq/202410_API_Pf3D7_NextSeq/statistics/Input/Sample_list.txt', header = F)
#input.dir.bed<- "../Input/mapped_pairedRead1_bed/"

Pf_theo <- read.table('/mathspace/data01/mathbio_lab/Tnseq/202410_API_Pf3D7_NextSeq/statistics/Input/PF3D7_theo_TTAA_R_modified_final.bed', sep = '\t')
Pf_theo_anti <- Pf_theo %>% dplyr::filter(V6 == '-')
Pf_theo_sense <- Pf_theo %>% dplyr::filter(V6 == '+')
index_Pf_anti <- which(Pf_theo$V6 == '-')
index_Pf_sense <- which(Pf_theo$V6 == '+')

Total_Pf_theo <- nrow(Pf_theo)
Total_Pf_theo_R_gtf <- Total_Pf_theo

counts_distribution_Pf_samples <- as.data.frame(matrix(nrow = Total_Pf_theo, ncol =(length(as.list(samples2$V1))) ))
colnames(counts_distribution_Pf_samples) <-as.list(samples2$V1)
#mutate add a column, transmute create a new column
counts_distribution_Pf <- data.frame(Chrom = rep(NA, Total_Pf_theo),
                                     Site = rep(NA, Total_Pf_theo),
                                     R1_pointing_downsteam = rep(NA, Total_Pf_theo),
                                     R1_pointing_upstream = rep(NA, Total_Pf_theo),
                                     GeneID = rep(NA, Total_Pf_theo),
                                     Location = rep(0, Total_Pf_theo),
                                     Location1 = rep(0, Total_Pf_theo),
                                     Location2 = rep(0, Total_Pf_theo),
                                     Present_in_any_samples = rep(NA, Total_Pf_theo),
                                     Run_repeats = rep(NA, Total_Pf_theo),
                                     Theoretical = rep(1, Total_Pf_theo)
)
counts.distri.Pf.all <- cbind(counts_distribution_Pf, counts_distribution_Pf_samples)

counts.distri.Pf.all$Chrom <- Pf_theo$V1
counts.distri.Pf.all$Site <- Pf_theo$V2


counts.distri.Pf.all$R1_pointing_downsteam[index_Pf_anti] <- Pf_theo_anti$V6
counts.distri.Pf.all$R1_pointing_upstream[index_Pf_sense] <- Pf_theo_sense$V6


#############################################################################according to TTAA bed files to input the count matrix
###########!!!!must use gtf file to intersect 
###########!!!!must use gtf file to intersect
###########!!!!must use gtf file to intersect
###########!!!!must use gtf file to intersect

#Exons = transcript - introns; CDS = transcript- introns - UTRs; CDS = Exons - UTRs; intergenic = !transcript
#！！！bedtools intersect only shows intersect rows, do not show all the genes' exons.
TTAAhits_Pf <- read.table('/mathspace/data01/mathbio_lab/Tnseq/202410_API_Pf3D7_NextSeq/statistics/Input/TTAAhits_R_gtf_include_contigs.bed', sep = '\t' )
#extract transcript and Exons
TTAAhits_Pf_subset <- TTAAhits_Pf %>% dplyr::filter(V9 == "transcript"| V9 == "exon") 
nrow(TTAAhits_Pf_subset )#1004660
Pf_GeneID <- lapply(strsplit(TTAAhits_Pf_subset$V15, ' '), '[[',4)
Pf_GeneID <- gsub( ';', '', Pf_GeneID)
TTAAhits_Pf_subset$V15 <- Pf_GeneID
## Only remain those column V6=V13
TTAAhits_Pf_subset <- TTAAhits_Pf_subset %>% dplyr::filter(V6 == V13) #make it half and reduce redundancy
TTAAhits_Pf_transcript <- TTAAhits_Pf_subset %>% dplyr::filter(V9 == "transcript")
nrow(TTAAhits_Pf_subset) # 502330

TTAAhits_Pf_transcript<- TTAAhits_Pf_transcript[, c(1,2,6,15)]
nrow(TTAAhits_Pf_transcript) #259364

TTAAhits_Pf_exon <- TTAAhits_Pf_subset %>% dplyr::filter(V9 == "exon")
TTAAhits_Pf_exon<- TTAAhits_Pf_exon[, c(1,2,6,15)]
nrow(TTAAhits_Pf_exon) #242966

counts.distri.Pf.all$R1_pointing_downsteam[which(is.na(counts.distri.Pf.all$R1_pointing_downsteam) == 'TRUE')] <- 'NA'
counts.distri.Pf.all$R1_pointing_upstream[which(is.na(counts.distri.Pf.all$R1_pointing_upstream) == 'TRUE')] <- 'NA'

#for(i in 1:nrow(counts.distri.Pf.all)){
#  TTAAhits_Pf_transcript
#  TTAAhits_Pf_exon
#  #transverse counts.distri.Pf.all
#  subset_exon <- TTAAhits_Pf_exon %>% dplyr::filter(V2 == counts.distri.Pf.all[i ,]$Site)
#  #considering the strandness
#  subset_exon <- subset_exon %>% dplyr::filter(V1 == counts.distri.Pf.all[i ,]$Chrom & (V6 == counts.distri.Pf.all[i ,]$R1_pointing_downsteam | V6 == counts.distri.Pf.all[i ,]$R1_pointing_upstream))
#  subset_transcript <- TTAAhits_Pf_transcript %>% dplyr::filter(V2 == counts.distri.Pf.all[i ,]$Site)
#  subset_transcript <- subset_transcript %>% dplyr::filter(V1 == counts.distri.Pf.all[i ,]$Chrom & (V6 == counts.distri.Pf.all[i ,]$R1_pointing_downsteam | V6 == counts.distri.Pf.all[i ,]$R1_pointing_upstream))
#  if(nrow(subset_exon) == 1){
#    counts.distri.Pf.all[i, ]$Location1 = 1 # 1 means exon
#    counts.distri.Pf.all[i, ]$Location2 = 2 # 2 means transcript, and if it is exon, it must be in transcript
#    counts.distri.Pf.all[i, ]$GeneID <- subset_exon$V15
#  } 
#  if(nrow(subset_transcript) == 1){
#    counts.distri.Pf.all[i, ]$Location2 = 2 # 2 means transcript
#    counts.distri.Pf.all[i, ]$GeneID <- subset_transcript$V15
#  } 
#  cat(paste('processing column', i))
#  cat('\n')
#}
##wait too long time,about up to 10 hours

###Use parallel::mclapply() instead of lapply() to run tasks across multiple cores
# Check the number of cores
detectCores()
# Define the function to process each row
process_row <- function(i) {
  subset_exon <- TTAAhits_Pf_exon %>%
    filter(V2 == counts.distri.Pf.all[i, ]$Site) %>%
    filter(V1 == counts.distri.Pf.all[i, ]$Chrom & 
             (V6 == counts.distri.Pf.all[i, ]$R1_pointing_downsteam | 
                V6 == counts.distri.Pf.all[i, ]$R1_pointing_upstream))
  
  subset_transcript <- TTAAhits_Pf_transcript %>%
    filter(V2 == counts.distri.Pf.all[i, ]$Site) %>%
    filter(V1 == counts.distri.Pf.all[i, ]$Chrom & 
             (V6 == counts.distri.Pf.all[i, ]$R1_pointing_downsteam | 
                V6 == counts.distri.Pf.all[i, ]$R1_pointing_upstream))
  
  result <- counts.distri.Pf.all[i, ]  # Copy the row to modify
  
  # Update result based on subset conditions
  if (nrow(subset_exon) == 1) {
    result$Location1 <- 1  # exon
    result$Location2 <- 2  # transcript (if in exon, it's in transcript)
    result$GeneID <- subset_exon$V15
  } 
  if (nrow(subset_transcript) == 1) {
    result$Location2 <- 2  # transcript
    result$GeneID <- subset_transcript$V15
  }
  
  cat(paste('processing row', i, '\n'))
  return(result)
}

# Set up parallel processing using 16 cores
num_cores <- 4
results <- mclapply(1:nrow(counts.distri.Pf.all), process_row, mc.cores = num_cores)

# Combine the results back into a single data frame
counts.distri.Pf.all <- bind_rows(results)


# loci in transcript but not in exon is in intron; introns = transcript - exons
#if location1+location2=0 ： intergenic 0+0
#if location1+location2=2 ： intron 0+2
#if location1+location2=3 ： exon 1+2

counts.distri.Pf.all$Location <- counts.distri.Pf.all %>% transmute(Location = (Location1 + Location2))
saveRDS(counts.distri.Pf.all, "./Input/Rdata_new/counts_distri_Pf_theo_only.RData")
table(counts.distri.Pf.all$Location)

##################################################################Input bed file to generate count matrix
##################################################################Input bed file to generate count matrix
##################################################################Input bed file to generate count matrix
##################################################################Input bed file to generate count matrix
##################################################################Input bed file to generate count matrix
##################################################################Input bed file to generate count matrix



counts.distri.Pf.all <- readRDS("./Input/Rdata_new/counts_distri_Pf_theo_only.RData")




input.dir.bed<- "/Users/sidaye/Documents/R/Tnseq/202207/Tnseq/Input/sam_markdup_before_extract/bed/"
samples2 <- read.table('./Input/Sample_list_Pf.txt', header = F)
table(counts.distri.Pf.all$Location)
Total_Pf_theo <- nrow(counts.distri.Pf.all)
Total_Pf_theo
#make sure number of Total_Pf_theo is right
for(i in 1:length(as.list(samples2$V1))){
  counts.distri.Pf.all[, i+11] <- rep(0, Total_Pf_theo)
  input.dir.bed
  sample_bed <- paste(input.dir.bed, samples2[i, 1], ".mapped.sorted.bed", sep = "")
  mapped <- read.table(sample_bed, sep = '\t')
  #discard V4
  mapped<- mapped[,c(1,2,3,6)]
  #sense strand need +0 and antisense strand need -4,which is different from Bd samples
  mapped1 <- mapped %>% mutate(V2 = ifelse(V6 == '+', V2, V3-4))
  mapped1 <- mapped1[,c(1,2,4)]
  #count duplicated rows
  mapped1 <- setDT(mapped1)[,list(Count=.N),names(mapped1)]
  names(mapped1) <- c('Chrom', 'Site', 'Strand', 'Count')
  
  #create a dataframe to record rows of mapped1 can not be assigned to a UII
  unmapped <- data.frame(V1 = rep(0, nrow(mapped1)),
                         V2 = rep(0, nrow(mapped1)),
                         V3 = rep(0, nrow(mapped1)),
                         Count = rep(0, nrow(mapped1))
  )
  
  for (j in 1:nrow(mapped1)){
    #in this step will get rid of mapped reads which is not at theoretical sites
    ind_count <- which(counts.distri.Pf.all$Site==mapped1[j, ]$Site & counts.distri.Pf.all$Chrom == mapped1[j, ]$Chrom)
    if(length(ind_count)==0){
      unmapped[j, ] <- mapped1[j, ]
    }else{
      if(mapped1[j, ]$Strand == '+'){ind_count_all <- ind_count[1]}else{ind_count_all <- ind_count[2]}
      counts.distri.Pf.all[ind_count_all, (i+11)] <- mapped1[j, ]$Count
    }
  }
  output_for_unmapped <- './Output/unmapped_goodReads_remove_duplicates_markdup/'
  write.xlsx(unmapped, paste0(output_for_unmapped, samples2[i, 1], '_unmapped_goodReads.xlsx'), na.string='NA', keepNA=F)
  cat(paste('processing file', samples2[i, 1]))
  cat('\n')
}

saveRDS(counts.distri.Pf.all, "./Input/Rdata_new/remove_duplicates_markdup/counts_distri_Pf_all.RData")

#optional: use function rowSums
counts.distri.Pf.all$Present_in_any_samples <- counts.distri.Pf.all %>% transmute(Present_in_any_samples  = Sample5+Sample7+Sample8+Sample16+Sample24+Sample25+Sample26+Sample27+Sample28+Sample29+Sample30+Sample31+Sample32+Sample33+Sample54+Sample55+Sample56+Sample57+Sample58+Sample59)
counts.distri.Pf.all$Present_in_any_samples <- counts.distri.Pf.all$Present_in_any_samples[[1]]
counts.distri.Pf.all$Present_in_any_samples[which(counts.distri.Pf.all$Present_in_any_samples!=0)] <- 'yes'
counts.distri.Pf.all$Present_in_any_samples[which(counts.distri.Pf.all$Present_in_any_samples==0)] <- 'no'
counts.distri.Pf.all$Location <- counts.distri.Pf.all %>% transmute(Location = (Location1 + Location2))
counts.distri.Pf.all$Location <- counts.distri.Pf.all$Location[[1]]
table(counts.distri.Pf.all$Location)
counts.distri.Pf.all$Location[which(counts.distri.Pf.all$Location == 0)] <- 'intergenic'
counts.distri.Pf.all$Location[which(counts.distri.Pf.all$Location == 2)] <- 'intron'
counts.distri.Pf.all$Location[which(counts.distri.Pf.all$Location == 3)] <- 'exon'
#error in is.nan(tmp) : default method not implemented for type 'list', because transmute create list
table(counts.distri.Pf.all$Location)#this is the only the statistics for reads locations 
write.xlsx(counts.distri.Pf.all, "./Output/reads_count_matrix_remove_duplicates_markdup/counts.distri.Pf.all.xlsx", na.string='NA', keepNA=F)
#write.xlsx(counts.distri.Pf.all, "./Input/Rdata_new/counts.distri.Pf.all.xlsx", na.string='NA', keepNA=F)



###################################Check with PlasmoDB,total theo TTAA insertion sites

PlasmoDB_TTAA <- read.csv('./Input/PlasmoDB_Pf.csv')
PlasmoDB_TTAA_contigs <- unlist(lapply(strsplit(PlasmoDB_TTAA$Genomic.Location, ':'), '[[',1))
names(Pf.fasta_all)
PlasmoDB_TTAA.df <- data.frame(Chrom = PlasmoDB_TTAA_contigs)
num_PlasmoDB_TTAA <- as.data.frame(table(PlasmoDB_TTAA.df))
#rule out contigs
num_PlasmoDB_TTAA_no_contigs <- num_PlasmoDB_TTAA[122:137, ]
total_num_PlasmoDB_TTAA <- sum(num_PlasmoDB_TTAA_no_contigs$Freq)
total_num_TTAA_R <- nrow(counts.distri.Pf.all)
total_num_df <- data.frame(class=c('PlasmoDB_TTAA', 'R_TTAA'), num=c(total_num_PlasmoDB_TTAA))

#plot
total_num.histo <- ggplot(total_num_df, aes(class, num)) + 
  geom_bar(stat='identity') +
  ylab('Count')+
  xlab('Comparison between PlasmoDB and R code results for TTAA recognition')+
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) 


print(total_num.histo)

#######################################

#######################Sample_list
samples2 <- read.table('./Input/Sample_list_Pf.txt', header = F)

input.dir.bed<- "./Input/sam_markdup_before_extract/bed/"

counts.distri.Pf.all <- readRDS("./Input/counts_distri_Pf_theo_only.RData")
table(counts.distri.Pf.all$Location)
Total_Pf_theo <- nrow(counts.distri.Pf.all)
####!!!!
####!!!!
####!!!!
####!!!!need to change number of samples and sample name according to different batch's samplesheet!!!
n <- 7
counts.distri.Pf.all <- counts.distri.Pf.all[, 1:(11+n)]
names(counts.distri.Pf.all)[(11+1):(11+n)] <- samples2$V1



#make sure number of Total_Pf_theo is right
for(i in 1:length(as.list(samples2$V1))){
  counts.distri.Pf.all[, i+11] <- rep(0, Total_Pf_theo)
  input.dir.bed
  sample_bed <- paste(input.dir.bed, samples2[i, 1], ".mapped.sorted.bed", sep = "")
  mapped <- read.table(sample_bed, sep = '\t')
  #discard V4
  mapped<- mapped[,c(1,2,3,6)]
  #sense strand need +0 and antisense strand need -4,which is different from Bd samples
  mapped1 <- mapped %>% mutate(V2 = ifelse(V6 == '+', V2, V3-4))
  mapped1 <- mapped1[,c(1,2,4)]
  #count duplicated rows
  mapped1 <- setDT(mapped1)[,list(Count=.N),names(mapped1)]
  names(mapped1) <- c('Chrom', 'Site', 'Strand', 'Count')
  
  #create a dataframe to record rows of mapped1 can not be assigned to a UII
  unmapped <- data.frame(V1 = rep(0, nrow(mapped1)),
                         V2 = rep(0, nrow(mapped1)),
                         V3 = rep(0, nrow(mapped1)),
                         Count = rep(0, nrow(mapped1))
  )
  
  for (j in 1:nrow(mapped1)){
    #in this step will get rid of mapped reads which is not at theoretical sites
    ind_count <- which(counts.distri.Pf.all$Site==mapped1[j, ]$Site & counts.distri.Pf.all$Chrom == mapped1[j, ]$Chrom)
    if(length(ind_count)==0){
      unmapped[j, ] <- mapped1[j, ]
    }else{
      if(mapped1[j, ]$Strand == '+'){ind_count_all <- ind_count[1]}else{ind_count_all <- ind_count[2]}
      counts.distri.Pf.all[ind_count_all, (i+11)] <- mapped1[j, ]$Count
    }
  }
  output_for_unmapped <- './Output/unmapped_goodReads_remove_duplicates_markdup/'
  write.xlsx(unmapped, paste0(output_for_unmapped, samples2[i, 1], '_unmapped_goodReads.xlsx'), na.string='NA', keepNA=F)
  cat(paste('processing file', samples2[i, 1]))
  cat('\n')
}

#/mathspace/data01/mathbio_lab/Tnseq/202212_Novaseq/count_matrix
saveRDS(counts.distri.Pf.all, "./Input/Rdata_new/remove_duplicates_markdup/counts_distri_Pf_all.RData")
#optional: use function rowSums
##need to change and customized
counts.distri.Pf.all$Present_in_any_samples <- counts.distri.Pf.all %>% transmute(Present_in_any_samples  = Sample1+Sample2+Sample3+Sample4+Sample5+Sample6+Sample8)
counts.distri.Pf.all$Present_in_any_samples <- counts.distri.Pf.all$Present_in_any_samples[[1]]
counts.distri.Pf.all$Present_in_any_samples[which(counts.distri.Pf.all$Present_in_any_samples!=0)] <- 'yes'
counts.distri.Pf.all$Present_in_any_samples[which(counts.distri.Pf.all$Present_in_any_samples==0)] <- 'no'
counts.distri.Pf.all$Location <- counts.distri.Pf.all %>% transmute(Location = (Location1 + Location2))
counts.distri.Pf.all$Location <- counts.distri.Pf.all$Location[[1]]
table(counts.distri.Pf.all$Location)
counts.distri.Pf.all$Location[which(counts.distri.Pf.all$Location == 0)] <- 'intergenic'
counts.distri.Pf.all$Location[which(counts.distri.Pf.all$Location == 2)] <- 'intron'
counts.distri.Pf.all$Location[which(counts.distri.Pf.all$Location == 3)] <- 'exon'
#error in is.nan(tmp) : default method not implemented for type 'list', because transmute create list
table(counts.distri.Pf.all$Location)#this is the only the statistics for reads locations 

#/mathspace/data01/mathbio_lab/Tnseq/202212_Novaseq/count_matrix
write.xlsx(counts.distri.Pf.all, "./Output/reads_count_matrix_remove_duplicates_markdup/counts.distri.Pf.all.xlsx", na.string='NA', keepNA=F)

