library(tidyverse)
library(IRanges)
library(Biostrings) 
library(ShortRead)
library(openxlsx)
library(dplyr)
library(data.table)
library(rtracklayer)
library(parallel)

####This script is for generating countmatrix template and mapping reads to TTAA sites
#######################################################bedtools intersect -a -b -wa -wb > .bed !!!must use gtf file instead of gff file
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
num_cores <- 16
results <- mclapply(1:nrow(counts.distri.Pf.all), process_row, mc.cores = num_cores)

# Combine the results back into a single data frame
counts.distri.Pf.all <- bind_rows(results)


# loci in transcript but not in exon is in intron; introns = transcript - exons
#if location1+location2=0 ： intergenic 0+0
#if location1+location2=2 ： intron 0+2
#if location1+location2=3 ： exon 1+2

counts.distri.Pf.all$Location <- counts.distri.Pf.all %>% transmute(Location = (Location1 + Location2))
saveRDS(counts.distri.Pf.all, "/mathspace/data01/mathbio_lab/Tnseq/202410_API_Pf3D7_NextSeq/R_output/counts_distri_Pf_theo_only.RData")
table(counts.distri.Pf.all$Location)