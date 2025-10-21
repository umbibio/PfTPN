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
n <- 34
#WT
n <- 19
n <- 11
n <- 7
#######################Count matrix should be both exon/geneID conversion version, conversion_for_exon matrix only repeat exon/exon_geneID, not include intron/intron_geneID
cm_Pf <- read.xlsx("./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon202412_202502_202504_Bg_removed_siteslevel.xlsx")
cm_Pf <- cm_Pf[,c(1:12,13,14,16,18,20,22,24,26,27,29,30,31,32,35,37,39,40,43,45)]###19 WT in total including the Day0 and Day8/9
cm_Pf <- cm_Pf[,c(1:12,14,16,22,24,26,31,32,35,37,39,40)]###11 WT in total excluding the Day0 and Day8/9
cm_Pf <- cm_Pf[,c(1:12,16,24,26,35,37,39,40)]###7 WT in total after Day15

cm_Pf$Total <- rowSums(cm_Pf[,13:(13+n-1)]) ###13 corresponds to matrix including gene.description column while 12 corresponds to matrix without

####To calculate total number of observed insertions after bg noise removal
cm_Pf$Total <- round(cm_Pf$Total)
nrow(cm_Pf %>% dplyr::filter(Total != 0)) ######Count>0 after bg noise removal
####249865 for 34 samples of WT+Mev 
nrow(cm_Pf)
###Total number of genes
length(unique(cm_Pf$GeneID))-1
####Optional：we can set up a manually cutoff for Total column as well to remove bg further
cutoff <- 0
modified_cm_Pf  <- cm_Pf%>% dplyr::mutate(ad.total = ifelse(Total<=cutoff, 0,Total-cutoff))
################Coupon collection model##############
################Coupon collection model##############
################Coupon collection model##############
#########################################Bidirectional model,gene level saturation#########################################Upper bound, every site is non-essential
#########################################Bidirectional model,gene level saturation#########################################Upper bound, every site is non-essential
#########################################Bidirectional model,gene level saturation#########################################Upper bound, every site is non-essential
get_insertedgeneNo_cm <- function(tmp){
  ###-1 is to remove the count of NA
  TotalNo_actual_inserted_genes_yes <- length(unique(tmp$GeneID))-1
  return(TotalNo_actual_inserted_genes_yes)
}

###Mode=1: count insertions within introns/exons(CDS) as targeting genes
###Mode=2: count insertions within exons(CDS) only as targeting genes.
CouponC_model_cm<- function(total.sites, tmp, max_n, ratio, mode){
  ###Removing NA
  gene_list <- na.omit(unique(tmp$GeneID))
  No.genes <- length(gene_list)
  if(ratio!=0){
    ###Randomly sample genes as essential genes based on the assuming ratio
    essential.genes.ind <- sample(No.genes, No.genes*ratio, replace = F)
    essential.genes <- gene_list[essential.genes.ind]
    tmp <- tmp[!tmp$GeneID%in%essential.genes,]
  }else{}
  #all.sites <- rep(possible.sites, 10000000)
  if(mode==1){
    nn <- seq(1, max_n, by=1000)
    ss <- rep(0, length(nn))
    possible.sites <- sample(1:nrow(tmp))
    count=1
    for (n in nn){
      r.s <- sample(possible.sites, n, replace = T)
      ss[count] <- get_insertedgeneNo_cm(tmp[r.s, ])
      nn[count] <-  length(unique(r.s))
      count=count+1
    }
    df <- data.frame(nn = nn, ss = ss)
    return(df)
  }else if(mode==2){
    nn <- seq(1, max_n, by=1000)
    ss <- rep(0, length(nn))
    possible.sites <- sample(1:nrow(tmp))
    count=1
    for (n in nn){
      r.s <- sample(possible.sites, n, replace = T)
      ####Only insertions in exons/CDS counts
      filtered_tmp <- tmp[r.s, ] %>% dplyr::filter(Location=='exon')
      ss[count] <- get_insertedgeneNo_cm(filtered_tmp)
      nn[count] <-  length(unique(r.s))
      count=count+1
    }
    df <- data.frame(nn = nn, ss = ss)
    return(df)
  }else{
    print("mode error: needs to be 1 or 2")
  }
}

######Optional: convert geneID for both exon and intron####
#modified_cm <- function(tmp){
#  index_NA <- which(is.na(tmp$GeneID))
#  for (i in 1:length(index_NA)){
#    tmp$GeneID[index_NA[i]]=ifelse(index_NA[i] %% 2 == 0,tmp$GeneID[index_NA[i]-1],tmp$GeneID[index_NA[i]+1])
#  }
#  return(tmp)
#}
######Optional: convert geneID for both exon and intron####


#set.seed(001)
#Pf_modelBi <- CouponC_model_cm(total.sites=660272, modified_cm_Pf, max_n=1000000, ratio=0, mode=1)


set.seed(0012)
Pf_modelBi2 <- CouponC_model_cm(total.sites=660272, modified_cm_Pf, max_n=1000000, ratio=0, mode=2)
threshold_Pf <- Pf_modelBi2[which(Pf_modelBi2$ss>= (5637 * .95))[1],]$nn
threshold_Pf ##63734

pp <- ggplot(Pf_modelBi2, aes(x= nn, y= ss)) + geom_point(color = '#595959') + 
  geom_smooth(method = "loess", span = 0.01, color = 'red')+ theme_cowplot() + 
  labs(x = "Bidirectional unique insertions", y="Targeted genes") +
  scale_x_continuous(limits = c(0, 600000), labels = scales::comma) + scale_y_continuous(limits = c(0, 6000), breaks=seq(0, 6000, 1000))+ 
  theme(plot.title = element_text(color="black", size=14), axis.text = element_text(size = 12),  axis.title=element_text(size=14)) +
  geom_vline(xintercept =threshold_Pf, color = '#cc0000', lty=2, size = 1.2)


#########################################Bidirection model#########################################Lower bound, 40% genes is essential
#########################################Bidirection model#########################################Lower bound, 40% genes is essential
#########################################Bidirection model#########################################Lower bound, 40% genes is essential
##############Randomly pick up the essential sites, and block or remove it###############
##############Randomly pick up the essential sites, and block or remove it###############
##############Randomly pick up the essential sites, and block or remove it###############

#set.seed(002)
#Pf_model_lowerbound3 <- CouponC_model_cm(total.sites=160126*2, modified_cm_Pf, max_n=1000000,ratio=0.4,mode = 1)

set.seed(0022)
Pf_model_lowerbound3_2 <- CouponC_model_cm(total.sites=660272, modified_cm_Pf, max_n=1000000,ratio=0.4,mode = 2)

threshold_Pf_lowerbound <- Pf_model_lowerbound3_2[which(Pf_model_lowerbound3_2$ss>=round(5637 * .6 *.95))[1],]$nn
threshold_Pf_lowerbound

final_plot <- pp +
  geom_point(data = Pf_model_lowerbound3_2, aes(x = nn, y = ss), color = '#595959') +
  geom_smooth(data = Pf_model_lowerbound3_2, color = 'blue',method = "loess", span = 0.01)+# You can choose the smoothing method you prefer
  geom_vline(xintercept =threshold_Pf_lowerbound, color = 'blue', lty=2, size = 1.2)

#print(final_plot)
out.dir <- "./Output/Figures/"
ggsave(final_plot, filename = paste(out.dir,"saturation_CDS2", '.pdf',sep = ""), width = 6,height = 3, dpi = 400)
###6x4:2
###8X4
#write.table(Pf_modelBi,"./Output/Saturation_level_estimation/upperbound_genelevel_df_max_n1000000_mode1_including_intron.txt", row.names = F, col.names = F, sep = "\t")
#write.table(Pf_model_lowerbound3,"./Output/Saturation_level_estimation/lowerbound_genelevel_df_max_n1000000_ratio04_mode1_including_intron.txt", row.names = F, col.names = F, sep = "\t")
write.table(Pf_modelBi2,"./Output/Saturation_level_estimation/upperbound_genelevel_df_max_n1000000_mode2_exononly.txt", row.names = F, col.names = F, sep = "\t")
write.table(Pf_model_lowerbound3_2,"./Output/Saturation_level_estimation/lowerbound_genelevel_df_max_n1000000_ratio04_mode2_exononly.txt", row.names = F, col.names = F, sep = "\t")

Pf_modelBi <- read.table("./Output/Saturation_level_estimation/upperbound_genelevel_df_max_n1000000_mode1_including_intron.txt")
Pf_model_lowerbound3 <- read.table("./Output/Saturation_level_estimation/lowerbound_genelevel_df_max_n1000000_ratio04_mode1_including_intron.txt")
Pf_modelBi2 <- read.table("./Output/Saturation_level_estimation/upperbound_genelevel_df_max_n1000000_mode2_exononly.txt")
colnames(Pf_modelBi2) <- c("nn","ss")
Pf_model_lowerbound3_2 <- read.table("./Output/Saturation_level_estimation/lowerbound_genelevel_df_max_n1000000_ratio04_mode2_exononly.txt")
colnames(Pf_model_lowerbound3_2) <- c("nn","ss")