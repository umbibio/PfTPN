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
getwd()
setwd('/Users/sidaye/Documents/R/API_TnSeq')


#####Step1:cm_Location_conversion_for_exon_extraction.R
#Convert the location columns of count matrix because of the palindrome structure of TTAA, and we need to use it the extract exons only for genes, this step is both required for edgeR model and CV_inverse model
#####Step2: Perturbation_analysis_Pk_countmatrix2DIGmatrix_preparation to turn cm into DIG.R(OPtional: remove TTAA sites locating at 99% of transcript)

######please note that the sites at 99%> transcript has not been filtered for count matrix2DIG_df transformation
######Please note that for binary model, one strand is exon, the other strand is intergenic, should include the other strand also

##################Mode2: DIG for CDS region(not include intron)##################
##################Mode2: DIG for CDS region(not include intron)##################
##################Mode2: DIG for CDS region(not include intron)##################
###Double check if includes total columns
#cm <- read.xlsx("./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon_UTR_annotation_202412_202502_202504_202505.xlsx")
cm <- read.xlsx("./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon_UTR_annotation_all.xlsx")

###Double-check if UTR annotation is included
table(cm$Location)
cm <- cm%>%dplyr::filter(Location=="exon")

################### 75 samples' name for essentialome only
#Sample_list <- colnames(cm[,12:(ncol(cm)-1)])#has no gene.description column and contains total column

Sample_list <- colnames(cm[,13:(ncol(cm))])#contaims gene.description column and has no total column
####Input total sample names
#Samples2 <- read.table("./Input/Sample_list_Pk_run123.txt", header = F)
####Input those genes' geneID of which have >=1 TTAA within CDS(exons)
Total.df2 <- data.frame(geneID=unique(cm$GeneID))###5343 pc genes including API/MITO genes
nrow(Total.df2) ###5556 genes have >1 TTAA
##To remove API/MITO genes, remained 5558 pc genes
#Total.df2 <- Total.df2 %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))


########Can also remove the TTAA locating at 99% of the transcript in order to increase confidence
getMatrix_DIG <- function(cm, Sample_list){
  cm2 <- cm%>%dplyr::select(all_of(Sample_list))
  cm2 <- cbind(cm$GeneID,cm2)
  colnames(cm2)[1] <- 'geneID'
  #turn the Count matrix with rows are sites into Count matrix with rows are genes
  #create empty dataframe 
  #No of samples for perturbation analysis
  No_samples <- length(Sample_list) 
  
  df <- as.data.frame(matrix(0, nrow=nrow(Total.df2), ncol=No_samples))
  df2 <- data.frame(geneID=Total.df2$geneID)
  colnames(df) <- Sample_list
  df <- cbind(df2, df)
  for (i in 2:ncol(df)){
    tmp <- cm2[,i]
    #convert data frame to data table 
    tmp2 <- data.frame(geneID=cm2$geneID,
                       sample=tmp)
    tmp3 <- setDT(tmp2)
    ob.gene.insertions <- tmp3[ ,list(sum=sum(sample)), by=geneID]
    Total.df <- left_join(df, ob.gene.insertions, by = "geneID")
    Total.df$sum[is.na(Total.df$sum)] <- 0
    df[,i] <- Total.df$sum
  }
  return(df)
}


DIG_df <- getMatrix_DIG(cm, Sample_list=Sample_list)

write.xlsx(DIG_df, "./Output/Perturbation_CDS/original_tables/DIG_Pf_all_CDS.xlsx")
