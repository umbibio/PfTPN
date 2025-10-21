library(tidyverse)
library(IRanges)
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
library(mixtools)
library(scales)
library(colorspace)
library(cowplot)
library(rtracklayer)
library(ggvenn)
library(ggExtra) 
library(patchwork)
##########How many CDS has no insertions across all the samples############
MIS <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_WT_afterDay15_final.xlsx")
MIS2 <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_MEV_afterDay15_final.xlsx")


MIS$sum.observed.insertions <- round(MIS$sum.observed.insertions)
sum(MIS$sum.observed.insertions==0)
tb <- as.data.frame(table(MIS$sum.observed.insertions))
tb2 <-tb%>%dplyr::mutate(Var1 = as.numeric(as.character(Var1)),Var1=ifelse(Var1 > 10,11,Var1))
tb3 <- tb2%>%group_by(Var1)%>%summarise(n=sum(Freq))
colnames(tb3)[1] <- "Category"
tb3$Category[12] <-">10"

tb3$Category <- factor(tb3$Category, levels=unique(tb3$Category))
ggplot(tb3, aes(x =Category, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue",alpha = 1) +
  geom_text(aes(label = n), vjust = -0.3, size = 4) +
  labs(title = " ",
       x = "Genes with CDS has insertions",
       y = "Number") + theme_cowplot()

MIS2 <- MIS2%>%dplyr::select(geneID, MIS)
colnames(MIS2)[2] <- "MIS2"
MISvsMIS2 <- left_join(MIS,MIS2,by="geneID")
nrow(MISvsMIS2)
colnames(MISvsMIS2)
#####Mode2: Discrepancy analysis
MISvsMIS2 <- MISvsMIS2%>%dplyr::filter(((MIS>0.93)|(MIS<0.14))&((MIS2>0.85)|(MIS2<0.13)))
nrow(MISvsMIS2)

Pf_all1 <- MISvsMIS2%>%dplyr::filter(MIS>0.93)
Pf_all0 <- MISvsMIS2%>%dplyr::filter(MIS<0.14)

Pk_all1 <-  MISvsMIS2%>%dplyr::filter(MIS2>0.85)
Pk_all0 <-  MISvsMIS2%>%dplyr::filter(MIS2<0.13)

venn11 <- list(Pf.MIS.MEV=Pk_all1$geneID, Pf.MIS.WT=Pf_all1$geneID)
venn00 <- list(Pf.MIS.MEV=Pk_all0$geneID, Pf.MIS.WT=Pf_all0$geneID)

p11 <- ggvenn(
  venn11, 
  fill_color = c("#FF5575",  "#6299FF"),
  stroke_size = 0.5, set_name_size = 6,text_size = 5)
p11

p00 <- ggvenn(
  venn00, 
  fill_color = c("#FF5575",  "#6299FF"),
  stroke_size = 0.5, set_name_size = 6,text_size = 5)
p00

p11+p00

discrepant1 <- MISvsMIS2%>%dplyr::filter(MIS>0.93&MIS2<0.13)
discrepant2 <- MISvsMIS2%>%dplyr::filter(MIS<0.14&MIS2>0.85)

gene.table <- read.xlsx("./Input/ref_gene_list.xlsx", sheet = 1)

venn22 <- list(discrepant1=discrepant1$geneID, discrepant2=discrepant2$geneID, Api_associated=gene.table$Gene.ID)


p22 <- ggvenn(
  venn22, 
  fill_color = c("#d39f8c",  "#b3cfd8","#ffe599"),
  stroke_size = 0.5, set_name_size = 6,text_size = 5)
p22

geneName <- intersect(discrepant2$geneID, gene.table$Gene.ID)
geneName <- discrepant1$geneID

discrepant2 <- MISvsMIS2%>%dplyr::filter(MIS<0.14&MIS2>0.85)
discrepant2_M <- discrepant2 %>% dplyr::filter(discrepant2$geneID%in%geneName)

