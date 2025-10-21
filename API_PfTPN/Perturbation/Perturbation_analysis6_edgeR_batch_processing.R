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
library(EnhancedVolcano)
library(ggvenn)
library(VennDiagram)
library(ggVennDiagram)
getwd()
setwd('/Users/sidaye/Documents/R/API_TnSeq')
####Do not expect too much genes have differential insertions.

#Day20_M_0505 and Day20_M_0504 need to be excluded, since the data quality is super low

###################please note that the somple index in Perturb_Comparison_Info_all is the index in DIG_df
####Pipe in DIG CDS df or Genic df
DIG_df <- read.xlsx("./Output/Perturbation_CDS/original_tables/DIG_Pf_all_CDS_removeBg.xlsx")
dim(DIG_df)
colnames(DIG_df)
#Perturb_Comparison_Info_all <- read.xlsx("./Input/Perturbation_CDS/Perturbation_Comparison_Info.xlsx")


DIG <- function(Tnseq.c1, Tnseq.c2){
  ####Tnseq.c1=WT group
  ####Tnseq.c2=treated group
  Expr.c1 <- data.frame(Tnseq.c1[,2:ncol(Tnseq.c1)])
  Expr.c2 <- data.frame(Tnseq.c2[,2:ncol(Tnseq.c2)])
  Expr.c1.c2 <- cbind(Expr.c1, Expr.c2)
  rownames(Expr.c1.c2) <- Tnseq.c1[,1]
  
  ## Remove rows with low counts
  CPM  <- cpm(Expr.c1.c2)
  ##need to remove low counts, remove low counts can increase FDR
  keep <-  rowSums(CPM > 0) >= 0
  Expr.c1.c2 <- Expr.c1.c2[keep, ]
  print(paste('genes kept:', length(which(keep == T))))
  gene.id <- rownames(Expr.c1.c2)
  Group  <- factor(c(rep("0", ncol(Expr.c1)), rep("1", ncol(Expr.c2))))
  #determine which groups are compared to each other based on the design matrix used in the differential expression analysis
  design <- model.matrix(~Group)
  
  dge <- DGEList(counts=Expr.c1.c2, group = Group, genes = gene.id)
  #Perform TMM normalization
  dge <- calcNormFactors(dge)
  plotMDS(dge)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit <- glmFit(dge, design)
  #Likelihood Ratio Test (LRT) ,assumes a specific distribution for the count data (negative binomial)
  fit <- glmLRT(fit, coef = 2)
  #Compute adjusted p-values (e.g., using Benjamini-Hochberg correction)
  tab <- topTags(fit,n=Inf,adjust.method="BH")$table
  #remove missing values
  tab <- na.omit(tab)
  return(tab)
}

####Tnseq.c1=WT group
####Tnseq.c2=treated group
###Day8, 9
Tnseq.c1 <- DIG_df[,c(1,7,16,32,34,65,77,87,115,128)]
Tnseq.c2 <- DIG_df[,c(1,8,17,33,35,66,78,88,117,130)]

###Day14, 15, 16
Tnseq.c1 <- DIG_df[,c(1,3,11,20,21,39,40,67,69,71,79,81,89,96,102,118,121)]
Tnseq.c2 <- DIG_df[,c(1,4,10,22,23,41,42,68,70,72,80,82,90,98,104,120,123)]

###Day17, 18
Tnseq.c1 <- DIG_df[,c(1,13,24,26,43,46,59,73,83,91,105,124)]
Tnseq.c2 <- DIG_df[,c(1,12,25,27,44,47,54,60,74,84,92,107,126)]

###Day 19, 20, 26
Tnseq.c1 <- DIG_df[,c(1,5,15,28,29,49,50,55,61,75,85,133,63,135)]
Tnseq.c2 <- DIG_df[,c(1,6,14,30,51,56,62,76,86,134,64,136)] #52 need to be excluded,bad sample

####Merge all 
Tnseq.c1 <- DIG_df[,c(1,3,5,7,11,13,15,16)]
Tnseq.c2 <- DIG_df[,c(1,4,6,8,10,12,14,17)]


####or excluding Day8 samples
Tnseq.c1 <- DIG_df[,c(1,3,5,11,13,15)]
Tnseq.c2 <- DIG_df[,c(1,4,6,10,12,14)]



tab <- DIG(Tnseq.c1, Tnseq.c2)
tab$minus_log10Pvalue <- -log10(tab$PValue)
tab$minus_log10FDR <- -log10(tab$FDR)
tab <- tab[order(tab$logFC),]
tab$Bottom_rank <- seq(from=1, to=nrow(tab), by=1)

Product.description <- read.csv("./Input/5720_total_Pf_product_description.csv")
Product.description <- Product.description[,c(1,3,4)]
colnames(Product.description)[1] <- "genes"
tab2 <- left_join(tab, Product.description, by="genes")
####Add 2 columns show the total counts for treated and untreated counts after removing background noise
Tnseq.c1$Untreated_counts <- rowSums(Tnseq.c1[,c(2:ncol(Tnseq.c1))])
Tnseq.c1$Untreated_counts <- as.numeric(Tnseq.c1$Untreated_counts)
class(Tnseq.c1$Untreated_counts)
Tnseq.c1$Untreated_counts <- round(Tnseq.c1$Untreated_counts,digits = 6)
Tnseq.c2$Treated_counts <- rowSums(Tnseq.c2[,c(2:ncol(Tnseq.c2))])
Tnseq.c2$Treated_counts <- as.numeric(Tnseq.c2$Treated_counts)
class(Tnseq.c2$Treated_counts)
Tnseq.c2$Treated_counts <- round(Tnseq.c2$Treated_counts,digits = 6)
colnames(tab2)[1] <- "geneID"

tab3 <- left_join(tab2, Tnseq.c1[,c(1, ncol(Tnseq.c1))],by="geneID")
tab3 <- left_join(tab3, Tnseq.c2[,c(1, ncol(Tnseq.c2))],by="geneID")

#API_asssociated <- read.xlsx("./Input/ref_gene_list.xlsx", sheet=1)
#tab3 <- tab3%>%mutate(Category=ifelse(geneID%in%API_asssociated$Gene.ID,1,0))

#nrow(tab3%>%filter(FDR<=0.05))/nrow(tab3)
#nrow(tab3%>%filter(FDR<=0.05 & Category==1))/nrow(tab3%>%filter(Category==1))


#test <- matrix(c(nrow(tab3%>%filter(FDR<=0.05 & Category==1)), 
#                 nrow(tab3%>%filter(FDR>0.05 & Category==1)), 
#                 nrow(tab3%>%filter(FDR<=0.05 & Category!=1)),
#                 nrow(tab3%>%filter(FDR>0.05 & Category!=1))), nrow = 2, byrow = T,
#                      dimnames =
#                        list(c("API", "Non-API"),
#                             c("Diff", "Non-Diff")))
#print(test)
#fisher.test(test)
#library(vcd)

#mosaic(test, shade = TRUE, legend = TRUE, main = "Fisher's Exact Test: API-associated vs. All genes")


tab3$Gene.Name.or.Symbol <- ifelse(tab3$Gene.Name.or.Symbol=="N/A",tab3$geneID,tab3$Gene.Name.or.Symbol)
#tab2_API <- tab3 %>% dplyr::filter(geneID%in%API_asssociated$Gene.ID)
write.xlsx(tab3, "./Output/Perturbation_CDS/trending_analysis/Differential_insertions_merged_Day8_9_202510.xlsx")
write.xlsx(tab3, "./Output/Perturbation_CDS/trending_analysis/Differential_insertions_merged_Day14_15_16_202510.xlsx")
write.xlsx(tab3, "./Output/Perturbation_CDS/trending_analysis/Differential_insertions_merged_Day17_18_202510.xlsx")
write.xlsx(tab3, "./Output/Perturbation_CDS/trending_analysis/Differential_insertions_merged_Day19_20_26_202510.xlsx")
#write.xlsx(tab2_API, "./Output/Differential_insertions_merged_API_associated_Day14_15_202506.xlsx")

#C <- read.xlsx("./Output/Differential_insertions_merged_API_associated.xlsx")
#D <- read.xlsx("./Output/Differential_insertions_merged_API_associated_excludingDay8.xlsx")

#A <- read.xlsx("./Output/Differential_insertions_merged_API_associated.xlsx")
#A2 <- A%>%dplyr::filter(minus_log10FDR>1&logFC>0.585)
#A3 <- A%>%dplyr::filter(minus_log10FDR>1&logFC<(-0.585))
#B <- read.xlsx("./Output/Differential_insertions_merged_API_associated_excludingDay8.xlsx")
#B2 <- B%>%dplyr::filter(minus_log10FDR>1&logFC>0.585)
#B3 <- B%>%dplyr::filter(minus_log10FDR>1&logFC<(-0.585))

#venn11 <- list(WithDay8=A2$geneID, WithoutDay8=B2$geneID)
#venn22 <- list(WithDay8=A3$geneID, WithoutDay8=B3$geneID)
#color_vectors <- c("#9F7FBF","#87BFBF")
#p1 <- venn.diagram(
#  x = venn22,
#  category.names = c('',''),
#  filename = NULL,
  #fill = c("#FFFF00", "#0000FF", "#FF0000"),
#  fill = color_vectors,
#  alpha = 0.5,
#  cex = 1.5,
#  fontfamily = "sans", ###for numbers
  #cat.fontfamily = "sans",  # Set font family for category names to sans-serif
  ###Set.names size
#  cat.cex = 1.5,
#  cat.default.pos = "outer",
#  cat.fontfamily = "sans",
#  cat.dist = c(0.05, 0.05),
  #####0 means 12 o'clock
#  cat.pos = c(210, 150),
  #print.mode = c( "raw","percent"),
#  print.mode = c( "raw"),
#  direct.area = F,
#  disable.logging=T
#)
#grid.draw(p1)
#4X5.5
dev.off()

C <- read.xlsx("./Output/Differential_insertions_merged_Day20_202506.xlsx")
#D <- read.xlsx("./Output/Differential_insertions_merged_API_associated_Day20_202506.xlsx")

##############Volcano plot####################
EnhancedVolcano(tab3,
                title = "",
                subtitle = "",
                lab = tab3$Gene.Name.or.Symbol,
                pCutoff = 0.05,         # Significance threshold for FDR
                FCcutoff = 2,           # Threshold for fold change
                x = 'logFC',
                y = 'FDR',
                legendLabels = c('NS', 'Fold Change', 'FDR', 'FDR & Fold Change'),
                drawConnectors = TRUE,    # Draw lines connecting labels to points
                widthConnectors = 0.5,    # Width of the connector lines
                colConnectors = 'black',   # Color of the connector lines
                ylab = "-log10(FDR)",
                xlab = "log2FC",
                max.overlaps = 20,
                ylim=c(0,3))

EnhancedVolcano(tab2_API,
                title = "",
                subtitle = "",
                lab = tab2_API$Gene.Name.or.Symbol,
                pCutoff = 0.05,         # Significance threshold for FDR
                FCcutoff = 1.5,           # Threshold for fold change
                x = 'logFC',
                y = 'FDR',
                legendLabels = c('NS', 'Fold Change', 'FDR', 'FDR & Fold Change'),
                drawConnectors = TRUE,    # Draw lines connecting labels to points
                widthConnectors = 0.5,    # Width of the connector lines
                colConnectors = 'black',   # Color of the connector lines
                ylab = "-log10(FDR)",
                xlab = "log2FC",
                ylim=c(0,3))



