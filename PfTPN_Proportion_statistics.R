library(tidyverse)
library(IRanges)
library(Biostrings) 
library(ShortRead)
library(openxlsx)
library(stringr)

getwd()

####################################################complete the statistics table
####################################################complete the statistics table
####################################################complete the statistics table
####################################################complete the statistics table
####################################################complete the statistics table
####################################################complete the statistics table


initial_reads <- read.table('/Users/sidaye/Documents/R/API_TnSeq/Input/statistics/Output/initial_reads.txt',fill = TRUE, header = F)
read_with_adaptors <- read.table('/Users/sidaye/Documents/R/API_TnSeq/Input/statistics/Output/read1_with_adaptor.txt',fill = TRUE, header = F)
Passing_filters <- read.table('/Users/sidaye/Documents/R/API_TnSeq/Input/statistics/Output/Passing_filters.txt',fill = TRUE, header = F)

statistics <- read.xlsx('/Users/sidaye/Documents/R/API_TnSeq/Input/statistics/Output/Proportion_statistics_all.xlsx', skipEmptyRows=F)

sample_order <- read.table('/Users/sidaye/Documents/R/API_TnSeq/Input/statistics/Output/sample_order.txt',fill = TRUE, header = F)
sample_order <- sample_order$V19
sample_order <- unlist(lapply(strsplit(sample_order,"/"),'[[',8))
sample_order <-gsub("_R1_trimmed.fastq.gz","",sample_order)

sample_index <- match(sample_order,colnames(statistics))
####!!!!need to change according to different error, since the starting fastq is merged, the the order need to be manually curated
#sample_order <- sample_order$V21

#########################customized part, it is up to whether need to exclude undetermined samples########
#########################customized part, it is up to whether need to exclude undetermined samples########
#sample_order <- lapply(strsplit(sample_order, '_L001_'), '[[', 1)

#will change Undetermined as S0 later on
#sample_order[121] <- "missing_S121"
#sample_order <-lapply(sample_order,str_sub,-6,-1)
#sample_order <- lapply(strsplit(unlist(sample_order), '_'), '[[', 2)

#sample_order <- as.numeric(gsub('S', '', sample_order))

#########################customized part, it is up to whether need to exclude undetermined samples########
#########################customized part, it is up to whether need to exclude undetermined samples########
###sample order is the order for initial_reads, read_with_adaptors and Passing_filters, it is determined by index order in Sean group's data
###Check the order in sample_order.txt and compare to the statistics table
####202410 batch
sample_order <-c(1, 3, 11, 19, 7, 15,  2, 4, 12, 20, 8, 16, 5, 13, 21, 9, 17, 6, 14, 22,10,18)
####202411 batch
sample_order <-c(1,7,6,3,2,5,4)
####202412 batch
sample_order <-c(1,7,6,3,2,5,4)
####202502 batch
sample_order <-c(1,9,8,2,3,6,7,4,5)
####202504 batch
sample_order <-c(2,16,15,6,4,8,7,14,12,1,18,17,5,3,13,11,10,9)
####202505 batch
sample_order <-c(2,7,5,9,8,17,15,1,6,4,16,14,12,11,18,3,10,13)
####202507 batch
sample_order <-sample_index
####202508 batch
sample_order <-sample_index


df <- data.frame(initial_reads = parse_number(initial_reads$V5),
                 read_with_adaptors = parse_number(read_with_adaptors$V5),
                 passing_filters = parse_number(Passing_filters$V5),
                 sample_order= sample_order)




df<- df[order(df$sample_order), ]


####!!!! please check the order manually, input initial_reads and read_with_adaptors as first and second rows

####Need to match the order of statistics tables
statistics[1,] <-df$initial_reads 
statistics[2,] <- df$read_with_adaptors

write.xlsx(statistics,"./Output/statistics202410.xlsx", rowNames=T)
write.xlsx(statistics,"./Output/statistics202411.xlsx", rowNames=T)
write.xlsx(statistics,"./Output/statistics202412.xlsx", rowNames=T)
write.xlsx(statistics,"./Output/statistics202502.xlsx", rowNames=T)
write.xlsx(statistics,"./Output/statistics202504.xlsx", rowNames=T)
write.xlsx(statistics,"./Output/statistics202505.xlsx", rowNames=T)
write.xlsx(statistics,"./Output/statistics202507.xlsx", rowNames=T)
write.xlsx(statistics,"./Output/statistics202508.xlsx", rowNames=T)
#Map to Pf
#mapping_to_Bd_sta <- read.table('/Users/sidaye/Documents/R/API_TnSeq/Input/statistics/Output/map_to_Bd_properlypaired.txt',fill = TRUE, header = F)
mapping_to_Pf_sta <- read.table('/Users/sidaye/Documents/R/API_TnSeq/Input/statistics/Output/map_to_Pf_properlypaired.txt',fill = TRUE, header = F)

mapping_to_Pf = (mapping_to_Pf_sta$V1)
#mapping_to_Pk = (mapping_to_Pk_sta$V1)/2

####actually the sample order is the same as "sample_order" above in countfiles for fastq filters or cutadapt
####we can check it in chimera:ls *.txt in mapping statistics folder
#sample_order， double check sample_order with mapping order
n <- 22
####202412 batch
n <- 7
####202502 batch
n <- 9
####202504 batch
n <- 18
####202505 batch
n <- 18
####202507 batch
n <- 35
####202508 batch
n <- 48

####need to check it in chimera:ls *.txt in mapping statistics folder
####this time, there are 17 Bd samples
####202412 batch
Pf_index <- c(1:22)
####202502 batch
Pf_index <- c(1:7)
####202504 batch
Pf_index <- c(1:9)
####202505 batch
Pf_index <- c(1:18)
####202507 batch
Pf_index <- c(1:35)
####202508 batch
Pf_index <- c(1:48)

#statistics[6,Bd_index] <- mapping_to_Bd
statistics[6,] <- mapping_to_Pf
#statistics[6,Pk_index] <- mapping_to_Pk

proportions_of_reads_Pf <- statistics[,Pf_index]
#proportions_of_reads_Pk <- statistics[,Pk_index]


write.xlsx(proportions_of_reads_Pf, "/Users/sidaye/Documents/R/API_TnSeq/Output/Proportion_statistics_Pf_202410.xlsx")
write.xlsx(proportions_of_reads_Pf, "/Users/sidaye/Documents/R/API_TnSeq/Output/Proportion_statistics_Pf_202411.xlsx")
write.xlsx(proportions_of_reads_Pf, "/Users/sidaye/Documents/R/API_TnSeq/Output/Proportion_statistics_Pf_202412.xlsx")
write.xlsx(proportions_of_reads_Pf, "/Users/sidaye/Documents/R/API_TnSeq/Output/Proportion_statistics_Pf_202502.xlsx")
write.xlsx(proportions_of_reads_Pf, "/Users/sidaye/Documents/R/API_TnSeq/Output/Proportion_statistics_Pf_202504.xlsx")
write.xlsx(proportions_of_reads_Pf, "/Users/sidaye/Documents/R/API_TnSeq/Output/Proportion_statistics_Pf_202505.xlsx")
write.xlsx(proportions_of_reads_Pf, "/Users/sidaye/Documents/R/API_TnSeq/Output/Proportion_statistics_Pf_202507.xlsx")
write.xlsx(proportions_of_reads_Pf, "/Users/sidaye/Documents/R/API_TnSeq/Output/Proportion_statistics_Pf_202508.xlsx")


#manipulate in excel, and reinput
#proportions_of_reads <- read.xlsx('./Input/Visualization/Proportions_of_reads.xlsx')
#proportions_of_reads <-proportions_of_reads[, -1]
#proportions_of_reads_Bd <- proportions_of_reads[,sort(Bd_index)]
#proportions_of_reads_Pk <- proportions_of_reads[,sort(Pk_index)]

#Info_sample_tranfor <- read.xlsx("./Output/Info_sample_tranformation_table.xlsx")
#names(proportions_of_reads_Bd) <- Info_sample_tranfor$Sample_Name[sort(Bd_index)]
#names(proportions_of_reads_Pk) <- Info_sample_tranfor$Sample_Name[sort(Pk_index)]

#write.xlsx(proportions_of_reads_Pk, "./Output/Proportion_statistics_Pk.xlsx")
#write.xlsx(proportions_of_reads_Bd, "./Output/Proportion_statistics_Bd.xlsx")

###########################################Proportions of PCR duplicates
###########################################Proportions of PCR duplicates
###########################################Proportions of PCR duplicates
###########################################Proportions of PCR duplicates
###########################################Proportions of PCR duplicates

####Be careful, order may be different########
#number of pairs of reads mapped to target genome
Pf.map.to.Pf <- read.table('./Input/statistics/Output/PCRdup_stat/Pf.all.positionsort.sta.txt',fill = TRUE, header = F)


#number of pairs of reads after removing the PCR duplicates
Pf.removed <- read.table('./Input/statistics/Output/PCRdup_stat/Pf.all.markdupremoved.sta.txt',fill = TRUE, header = F)
#Pf.removed  <- Pf.removed [c(1,3,2,5,4,7,6,9,8),]

#number of R1 extracted with MAPQ>=10 and properly paired after removing PCR duplicates
Pf.R1<- read.table('./Input/statistics/Output/PCRdup_stat/Pf.all.mappedR1.sta.txt',fill = TRUE, header = F)




#make use of the proportion statistics template, again; this time skipEmptyRows=T
tmp <- read.xlsx('./Input/statistics/Output/Proportion_statistics_all.xlsx', skipEmptyRows=T)

#Pf_index <- c(1,2,3,4,5,6,7,8,9)
Pf.tmp <- tmp[,Pf_index] 


Pf.tmp[1,] <- unlist(Pf.map.to.Pf)
Pf.tmp[2,] <- unlist(Pf.removed)
Pf.tmp[3,] <- unlist(Pf.R1)

rownames(Pf.tmp) <- c("number of reads in bam","number of reads after removing PCR duplicates in bam ","extracted properly paired R1 with MAPQ >= 30 without PCR duplicates")

Pf.tmp[nrow(Pf.tmp)+1,] <- Pf.tmp[1,] - Pf.tmp[2,]
rownames(Pf.tmp)[4] <- "absolute number of reads removed as PCR duplicates"

write.xlsx(Pf.tmp, "./Output/PCdup_sta_Pk_original202504.xlsx", rowNames=T)
write.xlsx(Pf.tmp, "./Output/PCdup_sta_Pk_original202505.xlsx", rowNames=T)
write.xlsx(Pf.tmp, "./Output/PCdup_sta_Pk_original202507.xlsx", rowNames=T)
write.xlsx(Pf.tmp, "./Output/PCdup_sta_Pk_original202508.xlsx", rowNames=T)


