library(grid)
library(EnhancedVolcano)
library(hrbrthemes)
library(viridis)
library(gghalves)
library(PupillometryR)
library(cowplot)
library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)

getwd()
#####To check the statistics of Bd
gtf <- read.table('./Input/Genome/PlasmoDB-63_Pfalciparum3D7.gtf', sep = '\t')
Total_gene_ls <-  gsub(' ', '', gsub(';', '', lapply(strsplit(gtf$V9, 'gene_id'), '[[', 2)))
Total_gene_ls <- unique(Total_gene_ls)
length(Total_gene_ls)###5720 genes in total
Total.df2 <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_WT_afterDay14_final_removeCDS.xlsx")
###############################################
nrow(Total.df2) ####5238
Total.df_pcgenes <- Total.df2 %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
No_Total.df_pcgenes <- nrow(Total.df_pcgenes)
print(No_Total.df_pcgenes) #####5208#####

TTAA0_pcgenes <- Total.df_pcgenes%>%dplyr::filter(Theo.num.unique.insertions==0) ####means No. of TTAA within exons only, excluding introns
nrow(TTAA0_pcgenes)#138
#5285-5208=77 #####77 protein coding genes(No API/MIT) has 0 TTAA within CDS

sum(TTAA0_pcgenes$geneID%in%test$`unique(append(tm_covered$sense_geneID, tm_covered$antisen_geneID))`)

################Input original transposon matrix count matrix###########
tcm <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix_all.xlsx")
tcm_WT <- tcm[,c(1:6,10,16,18,20,25,26,29,31,33,34,44,45,48,51,54,55,60,62,64,66,
                                                         68,72,74,76,78,80,84,86,88,90,94,96,98,101,107,110,114,115,123,
                                                         126,129,136,138,140)] ### WT samples after Day14
#WT_samplename <- colnames(tcm_WT)[7:ncol(tcm_WT)]
#write.table(WT_samplename, "./Output/WT_samplename.txt", row.names = F, col.names = F, quote = F)
#test <- read.table("./Output/WT_samplename.txt")

TTAAhits.greater.than.transcript99.TTAA.ID <- read.xlsx("./Output/TTAA_greater_than_99transcript/Pf.TTAAhits.greater.than.cds99.TTAA.ID.xlsx")
nrow(TTAAhits.greater.than.transcript99.TTAA.ID)
TTAAhits.greater.than.transcript99.TTAA.ID <- TTAAhits.greater.than.transcript99.TTAA.ID%>%dplyr::select(TTAA.ID,greater.than.transcript99)
table(TTAAhits.greater.than.transcript99.TTAA.ID$greater.than.transcript99)
colnames(TTAAhits.greater.than.transcript99.TTAA.ID)[2] <- "TTAA.ID.yes"
###essential_genelist removing 99% sites

essential_geneslist <- read.xlsx("./Output/background_genelist_54outof65curated.xlsx")
merged_df_all <- Total.df2
nrow(merged_df_all)

df <- data.frame(geneID=merged_df_all$geneID,
                 Product.Description=merged_df_all$Product.Description,
                 Theo.num.unique.insertions=merged_df_all$Theo.num.unique.insertions)
dim(df)

#####To calculate observed insertions within the CDS/Transcript of genes, and rule out 99% transcript TTAA sites
colnames(tcm)[grep('antisense_geneID',colnames(tcm))] <- 'antisen_geneID'
ob.counts <- function(transposon_count_matrix, TTAAhits.greater.than.transcript99.TTAA.ID, df, mode){
  transposon_count_matrix$Total <- rowSums(transposon_count_matrix%>%dplyr::select(contains('TPN')))
  transposon_count_matrix <- transposon_count_matrix%>% mutate (Present_in_any_samples=ifelse(Total>0, "yes","no"))
  present_yes_matrix <- transposon_count_matrix %>% dplyr::filter(Present_in_any_samples == 'yes')
  if (mode==1){
    exon_matrix <- present_yes_matrix %>% dplyr::filter(Assigned_location == 'exon'|Assigned_location == 'exon&intron')
  }else{
    exon_matrix <- present_yes_matrix %>% dplyr::filter(Assigned_location == 'exon'|Assigned_location == 'intron')
  }
  
  ####Rule out 99% transcript TTAA sites
  exon_matrix$TTAA.ID <- paste(exon_matrix$Chrom, exon_matrix$Site, sep = ":")
  exon_matrix.modified <- left_join(exon_matrix, TTAAhits.greater.than.transcript99.TTAA.ID, by = 'TTAA.ID')
  table(exon_matrix.modified$TTAA.ID.yes) #3404
  #nrow(exon_matrix) - #3404
  nrow(exon_matrix.modified) - 3404
  #filter out '1' labelled rows on greater.than.transcript99 column
  exon_matrix <- exon_matrix.modified[is.na(exon_matrix.modified$TTAA.ID.yes), ]
  nrow(exon_matrix) #3274
  
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
  df_all <- left_join(df,ob.gene.insertions, by = 'geneID')
  df_all$sum[is.na(df_all$sum)] <- 0
  colnames(df_all)[4] <- "sum.observed.insertions"
  return(df_all)
}

####CDS only####
####CDS only####
####CDS only####
df_all <- ob.counts(transposon_count_matrix=tcm_WT,TTAAhits.greater.than.transcript99.TTAA.ID, df,mode=1)
####CDS only####
####CDS only####
####CDS only####
nrow(df_all)


####Need to remove all api/mito genes and remain genes in 14 chromosomes
#df_all <- df_all[grepl("PKNH_", df_all$geneID),]
No_gene_with_TTAA <- df_all %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE)& !grepl("STR", geneID, fixed = TRUE))%>%nrow()
df_all <- df_all %>% dplyr::filter(!grepl("API", geneID, fixed = TRUE) & !grepl("MIT", geneID, fixed = TRUE))
nrow(df_all)###5205

mean_bg <- round(mean(df_all[df_all$geneID%in% essential_geneslist$GeneID,]$sum.observed.insertions))
print(mean_bg)
mean_bg <-2
##mean of bg within CDS is 71
##mean of bg within CDS is 75

#1model="protein_coding_genes"
#2mode="lncRNA"
#3mode="all"
No_0TTAA_genes <- 77
##For double check
No_Total.df_pcgenes <- 5285


#No of protein coding genes with 0 TTAA site. Count TTAA within exons and introns both
#No_0TTAA_genes <- 96
##For double check
#print(No_gene_with_TTAA2)
#No_Total.df_pcgenes-No_gene_with_TTAA2


data_plot <- function(df_all, mode){
  if(mode=="protein_coding_genes"){
    selected_df <- df_all%>% dplyr::filter((grepl("PF3D7", geneID)))
  }else if(mode=="lncRNA"){
    selected_df <- df_all%>% dplyr::filter((grepl("STRG", geneID)))
  }else{
    selected_df <-df_all
  }
  selected_df <- selected_df%>%dplyr::mutate(Theo.num.unique.insertions=ifelse(Theo.num.unique.insertions > 10,11,Theo.num.unique.insertions))
  df_plot <- selected_df %>% group_by(Theo.num.unique.insertions)%>%summarise(
    above_bg = sum(sum.observed.insertions > mean_bg),
    below_bg = sum(sum.observed.insertions <= mean_bg))
  
  if(mode=="protein_coding_genes"){
    df0 <- data.frame(Theo.num.unique.insertions=0,
                      above_bg=0,
                      below_bg=No_0TTAA_genes
    )
  }else if(mode=="lncRNA"){
    df0 <- data.frame(Theo.num.unique.insertions=0,
                      above_bg=0,
                      below_bg=70
    )
  }else{
    df0 <- data.frame(Theo.num.unique.insertions=0,
                      above_bg=0,
                      below_bg=No_0TTAA_genes+70)
  }
  df_plot_all <- rbind(df0, df_plot)
  #df_plot_all <- rbind(df_plot)
  return(df_plot_all)
}

data_summary <- data_plot(df_all,mode="protein_coding_genes")
data_summary$Total_genes <- data_summary$above_bg+data_summary$below_bg
data_summary$Proportion <- round(data_summary$Total_genes/sum(data_summary$Total_genes),2)
# Create a data frame with counts
colnames(data_summary)[1] <- "Category"
data_summary$Category[12] <- ">10"
data_summary$Total <- data_summary$above_bg+data_summary$below_bg
colnames(data_summary)
##delete Total_genes column
data_summary <- data_summary[,c(1,2,3,6)]

######Turn the dataframe into long format for plotting
data_summary <- data_summary%>% pivot_longer(-c(Category, Total), names_to = 'Category2', values_to = 'num')
data_summary$Prop <- data_summary$num/data_summary$Total
#data_summary$Total[c(FALSE, TRUE)] <- NA
data_summary$Category <- factor(data_summary$Category, levels=unique(data_summary$Category))

#data_summary <- data_summary%>%mutate(Category2=ifelse(Category2=="above_bg","above 0","below 0"))

ggplot(data_summary, aes(x = Category, y = num, fill = Category2)) +
  geom_bar(stat = "identity", position="stack", alpha = 1) +
  labs(title = "", x = "No. of TTAA within CDS", y = "No. of genes") +theme_cowplot()+
  ylim(c(0,2500))+
  scale_fill_manual(values = c("above_bg" = "#FFD586", "below_bg" = "#9AC9DB")) +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(0.2, 1), legend.justification = c(0, 1),legend.text = element_text(size = 18))+
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 18))+
  geom_point(data = data_summary, aes(x = Category, y = Prop*2000, color = Category2, group=Category2), size=3, pch=19) +
  geom_smooth(data = data_summary, aes(x = Category, y = Prop * 2000, color = Category2, group = Category2),
              method = "loess", se = FALSE, span = 1.2) +
  scale_y_continuous(breaks = c(0, 500, 1000, 1500, 2000), labels = c(0, 500, 1000, 1500, 2000),
                     sec.axis = sec_axis(~./2000, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), name = "Proportion"))+
  scale_color_manual(values = c("above_bg" = "#FFA500", "below_bg" = "#2C62D4"))

###Statistics for paper:
sum(data_summary$num)
###77 pcgenes has 0 TTAA
### pcgenes has >=1 TTAA
sum(unique(data_summary$Total)[6:12]) #4839 pc genes has >=5 TTAA sites


##8X4 inches
# Break the y-axis
# coord_cartesian(ylim = c(0, 400))

#how many with less than 5 sites that are non-essential that can be called with reasonable confidence? 
less5df <- merged_df_all%>% dplyr::filter(!grepl("API", GeneID.Pk_H, fixed = TRUE) & !grepl("MIT", GeneID.Pk_H, fixed = TRUE)&Theo.num.unique.insertions<5&Theo.num.unique.insertions>0)
less5df_confidence <- less5df%>%dplyr::filter(HMS>0.8)

########################
tcm <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13.xlsx")
TTAA_ID <- read.xlsx("./Output/transposon_matrix/all/158734TTAA_site_ID.xlsx")

filter_TTAA_ID <- function(countmatrix,TTAA_ID){
  countmatrix$ID <- paste(countmatrix$Chrom, countmatrix$Site, sep=":")
  countmatrix2 <- countmatrix%>%dplyr::filter(ID%in%TTAA_ID$ID)
  return(countmatrix2)
}
#cm<- filter_TTAA_ID(countmatrix=cm, TTAA_ID=TTAA_ID)

tcm2<- filter_TTAA_ID(countmatrix=tcm, TTAA_ID=TTAA_ID)
nrow(tcm2)
sum(tcm2$Total>0)
sum(tcm2$Total)
