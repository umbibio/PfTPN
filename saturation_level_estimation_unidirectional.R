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
library(cowplot)

getwd()
setwd('/Users/sidaye/Documents/R/API_TnSeq')

###Number of samples
#n <- 34
#WT
#n <- 19
#n <- 11
n <- 7
######For Pf, need to include 
cm_Pf <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix_bgremoved202412_202502_202504.xlsx") ###Should not be matrix attached overlaping rows 

#cm_Pf <- cm_Pf[,c(1:6,10,18,20,29,31,33,34)]
cm_Pf <- cm_Pf[,c(1:6,10,18,20,29,31,33,34)] ###7 samples after Day15###These indexes are for transposon matrix

cm_Pf$Total <- rowSums(cm_Pf[,7:(7+n-1)]) ###13 corresponds to matrix including gene.description column while 12 corresponds to matrix without

####To calculate total number of observed insertions after bg noise removal
cm_Pf$Total <- round(cm_Pf$Total)
nrow(cm_Pf %>% dplyr::filter(Total != 0)) ######Count>0 after bg noise removal
nrow(cm_Pf)

###Total number of genes
length(unique(append(cm_Pf$sense_geneID,cm_Pf$antisen_geneID)))-1
####Optional：we can set up a manually cutoff for Total column as well to remove bg further
cutoff <- 0
modified_cm_Pf  <- cm_Pf%>% dplyr::mutate(ad.total = ifelse(Total<=cutoff, 0,Total-cutoff))
################Coupon collection model##############
################Coupon collection model##############
################Coupon collection model##############
#########################################Unidirectional model,gene level saturation#########################################Upper bound, every site is non-essential
#########################################Unidirectional model,gene level saturation#########################################Upper bound, every site is non-essential
#########################################Unidirectional model,gene level saturation#########################################Upper bound, every site is non-essential
total.sites = 330136
total.genes = 5637

get_insertedgeneNo <- function(tmp){
  present_yes_matrix <- tmp
  sense_geneID_ls_yes <- present_yes_matrix %>% dplyr::filter(sense_geneID !='NA')
  sense_geneID_ls_yes <- sense_geneID_ls_yes$sense_geneID
  antisen_geneID_ls_yes <-present_yes_matrix %>% dplyr::filter(antisen_geneID !='NA') 
  antisen_geneID_ls_yes <- antisen_geneID_ls_yes$antisen_geneID
  #just to append all the geneID in count matrix and to see the number of unique geneIDs
  gene_ID_ls_yes <- append(sense_geneID_ls_yes,antisen_geneID_ls_yes) 
  ##Total acutual covered inserted genes
  TotalNo_actual_inserted_genes_yes <- length(unique(gene_ID_ls_yes))
  return(TotalNo_actual_inserted_genes_yes)
}


#######################x-axis is No of unique insertions
CouponC_model<- function(total.sites, tmp, max_n,ratio){
  ###Removing NA
  gene_list <- na.omit(unique(append(tmp$sense_geneID,tmp$antisen_geneID)))
  No.genes <- length(gene_list)
  if(ratio!=0){
    ###Randomly sample genes as essential genes based on the assuming ratio
    essential.genes.ind <- sample(No.genes, No.genes*ratio, replace = F)
    essential.genes <- gene_list[essential.genes.ind]
    tmp <- tmp[!((tmp$sense_geneID%in%essential.genes)|(tmp$antisen_geneID%in%essential.genes)),]
  }else{}
  ####shuffle all the sites randomly
  possible.sites <- sample(1:nrow(tmp))
  #all.sites <- rep(possible.sites, 10000000)
  ####The total number of trials, which is also the limit of x-axis
  nn <- seq(1, max_n, by=1000)
  #make an empty ss first with every element is 0, this is the y-axis of final plot to show how many genes are targeted by different number of trials
  ss <- rep(0, length(nn))
  #for loop for each trial
  count=1
  for (n in nn){
    r.s <- sample(possible.sites, n, replace = T)
    ####Only insertions in exons/CDS counts
    filtered_tmp <- tmp[r.s, ] %>% dplyr::filter(Assigned_location=='exon')
    ss[count] <- get_insertedgeneNo(filtered_tmp)
    nn[count] <-  length(unique(r.s))
    count=count+1
  }
  df <- data.frame(nn = nn, ss = ss)
  return(df)
}


set.seed(001)
Pf_model <- CouponC_model(total.sites=330136, modified_cm_Pf , max_n=500000, ratio=0)
#threshold_Pf <- Pf_model[grep(round(5637 * .95), Pf_model$ss)[1],]$nn 
threshold_Pf <- Pf_model[which(Pf_model$ss > 5637 * 0.95)[1], ]$nn
threshold_Pf

pp <- ggplot(Pf_model, aes(x= nn, y= ss)) + geom_point(color = '#595959') + 
  geom_smooth(method = "loess", span = 0.01, color = 'red')+ theme_cowplot() + 
  labs(x = "Unidirectional unique insertions", y="Targeted genes") +
  scale_x_continuous(limits = c(0, 350000), labels = scales::comma) + scale_y_continuous(limits = c(0, 6000), breaks=seq(0, 6000, 1000))+ 
  theme(plot.title = element_text(color="black", size=14), axis.text = element_text(size = 12),  axis.title=element_text(size=14)) +
  geom_vline(xintercept =threshold_Pf, color = '#cc0000', lty=2, size = 1.2)
pp


#########################################Unidirectional model,gene level saturation#########################################Upper bound, 40% gene is non-essential
#########################################Unidirectional model,gene level saturation#########################################Upper bound, 40% gene is non-essential
#########################################Unidirectional model,gene level saturation#########################################Upper bound, 40% gene is non-essential

set.seed(002)
Pf_model_lowerbound3_2 <- CouponC_model(total.sites=330136, modified_cm_Pf , max_n=500000, ratio=0.4)

threshold_Pf_lowerbound <- Pf_model_lowerbound3_2[which(Pf_model_lowerbound3_2$ss>=round(5637 * .6 *.95))[1],]$nn
threshold_Pf_lowerbound

final_plot <- pp +
  geom_point(data = Pf_model_lowerbound3_2, aes(x = nn, y = ss), color = '#595959') +
  geom_smooth(data = Pf_model_lowerbound3_2, color = 'blue',method = "loess", span = 0.01)+# You can choose the smoothing method you prefer
  geom_vline(xintercept =threshold_Pf_lowerbound, color = 'blue', lty=2, size = 1.2)

out.dir <- "./Output/Figures/"
ggsave(final_plot, filename = paste(out.dir,"Unidirectional_v2", '.pdf',sep = ""), width = 6,height = 3, dpi = 400)
