library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(doParallel)
library(edgeR)
library(limma)

#####Input golden list essential genes, current model is 22 essential genes in the golden list with removing the 8 outliers
#backgroundgenes <- read.xlsx("./Output/background_genelist51.xlsx")
backgroundgenes <- read.xlsx("./Output/background_genelist_54outof65curated.xlsx")

#####Input count matrix
cm <- read.xlsx("./Output/count_matrix/Pf_count_matrix_all.xlsx")
#cm <- read.xlsx("./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon_UTR_annotation_all.xlsx")
colnames(cm)
####Bg noise calculation based on each samples and removed bg noise for each sample separately
#names of all the samples or transfection pools
sample_list <- colnames(cm[,c(13:(ncol(cm)))]) ###Be careful if it includes the gene.description column
names <- sample_list
#number of samples or transfection pools
n<- length(sample_list)
n

########################Sites_level_Bg_noise calculation
Sites_level_Bg_noise <- function(matrix, n, names, backgroundgenes){
  df.final <- data.frame(TP=names,
                         Sites_level_noise=rep(NA, length(names)))
  for (i in 1:n){
    #convert data frame to data table 
    df <- as.data.frame(matrix)
    df2 <- df[,c(5,12+i)]
    df2 <- setDT(df2)
    colnames(df2)[2] <- 'Total'
    #To filter out backgroud genes' rows, this filtered df includes the genic region in bg genes, since the bg noise is used to removed every site‘s noise including intron and intergenic regions
    df2.filtered <- df2 %>% dplyr::filter(GeneID %in%  backgroundgenes$GeneID)
    #calculate the average sites level Bg noise for single transfection pools
    
    df.final$Sites_level_noise[i] <- sum(df2.filtered$Total)/nrow(df2.filtered)
    
  }
  return(df.final)
}

######Bg noise at sites level for each sample
Sites_level_Bg_noise <- Sites_level_Bg_noise(cm, n, names, backgroundgenes)
print(Sites_level_Bg_noise)

##############Remove Bg noise at sites level, make sure it is starting from 12
Remove_Sites_level_Bg <- function(countmatrix, Sites_level_Bg_noise, n){
  Sites_level_Bg_noise <- Sites_level_Bg_noise$Sites_level_noise
  for (i in 1:n){
    countmatrix[,12+i] <- ifelse(countmatrix[,12+i]< Sites_level_Bg_noise[i],0, countmatrix[,12+i]-Sites_level_Bg_noise[i])
    ###############round the counts, this will affect sites' raw read counts<=1 only, please note no need to round the counts after removing for DIG analysis
    #countmatrix[,12+i] <- round(countmatrix[,12+i]) 
  }
  return(countmatrix)
}
colnames(cm)
cm_Bg_removed_siteslevel <- Remove_Sites_level_Bg(cm, Sites_level_Bg_noise, n)

write.xlsx(cm_Bg_removed_siteslevel, "./Output/count_matrix/Pf_count_matrix_all_Bg_removed_siteslevel_final.xlsx")
#write.xlsx(cm_Bg_removed_siteslevel, "./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon_UTR_annotation_all_Bg_removed_siteslevel_final.xlsx")

#########################Perform CPM normalization for count matrix by cpm function in edgeR package
sites_ID <- cm_Bg_removed_siteslevel[,c(1:12)]
cm_matrix <- cm_Bg_removed_siteslevel[,c(13:ncol(cm_Bg_removed_siteslevel))]
cpm_cm_matrix <- cpm(cm_matrix)
cm_Bg_removed_siteslevel <- cbind(sites_ID, cpm_cm_matrix)
colnames(cm_Bg_removed_siteslevel)
write.xlsx(cm_Bg_removed_siteslevel, "./Output/count_matrix/Pf_count_matrix_all_Bg_removed_siteslevel_cpm_final.xlsx")
#write.xlsx(cm_Bg_removed_siteslevel, "./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon_UTR_annotation_all_Bg_removed_siteslevel_cpm_final.xlsx")


