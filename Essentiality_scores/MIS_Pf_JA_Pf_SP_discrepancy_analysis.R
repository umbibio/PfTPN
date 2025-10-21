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

MIS <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_WTonly_afterDay15all.xlsx")
MIS_Sean <- MIS %>% dplyr::select('geneID','MIS')
colnames(MIS_Sean) <- c("GeneID.Pf_3D7","Pf.MIS2")

Pf_MIS_MFS <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/MIS_MFS_Pf.xlsx", startRow=2)
Pf_MIS_MFS <- Pf_MIS_MFS%>%dplyr::select('Gene_ID', 'Product.description', 'Gene.Identification', 'MIS', 'MFS','transcript.length')
colnames(Pf_MIS_MFS) <- c("GeneID.Pf_3D7","Pf_3D7.Product.description","Pf.phenotype","Pf.MIS","Pf.MFS","Pf.transcript.length")

#####Mode1: Direct comparison
MISvsMIS2 <- left_join(MIS_Sean,Pf_MIS_MFS, by='GeneID.Pf_3D7')

write.xlsx(MISvsMIS2,'./Output/MIS/MISvsMIS2.xlsx')
#####Mode2: Discrepancy analysis
MISvsMIS2 <- MISvsMIS2%>%dplyr::filter(((Pf.phenotype=='Mutable in CDS'&Pf.MIS>0.93)|(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.29))&
                                      ((Pf.MIS2>0.74)|(Pf.MIS2<0.07)))

Pf_all1 <- MISvsMIS2%>%dplyr::filter(Pf.phenotype=='Mutable in CDS'&Pf.MIS>0.93)
Pf_all0 <- MISvsMIS2%>%dplyr::filter(Pf.phenotype=='Non - Mutable in CDS'&Pf.MIS<0.29)

Pk_all1 <-  MISvsMIS2%>%dplyr::filter(Pf.MIS2>0.74)
Pk_all0 <-  MISvsMIS2%>%dplyr::filter(Pf.MIS2<0.07)

venn11 <- list(Pf.MIS.SP=Pk_all1$GeneID.Pf_3D7, Pf.MIS.JA=Pf_all1$GeneID.Pf_3D7)
venn00 <- list(Pf.MIS.SP=Pk_all0$GeneID.Pf_3D7, Pf.MIS.JA=Pf_all0$GeneID.Pf_3D7)

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

MISvsMIS2$diff <- MISvsMIS2$Pf.MIS2-MISvsMIS2$Pf.MIS

####To make the table####
Discrepancy_essential <- data.frame(geneID=unique(append(append(Pf_all0$geneID,Pb_all0$geneID),Pk_all0$geneID)),
                                    labels=NA)
# Label geneIDs
for (i in 1:length(Discrepancy_essential$geneID)) {
  if (Discrepancy_essential$geneID[i] %in% Reduce(intersect, venn00)) {
    Discrepancy_essential$labels[i] <- "PkPbPf shared"
  } else if (Discrepancy_essential$geneID[i] %in% Pb_all0$geneID & Discrepancy_essential$geneID[i] %in% Pk_all0$geneID & !( Discrepancy_essential$geneID[i] %in% Pf_all0$geneID)) {
    Discrepancy_essential$labels[i] <- "PkPb shared"
  } else if (Discrepancy_essential$geneID[i] %in% Pk_all0$geneID & Discrepancy_essential$geneID[i] %in% Pf_all0$geneID & !( Discrepancy_essential$geneID[i] %in% Pb_all0$geneID)) {
    Discrepancy_essential$labels[i] <- "PkPf shared"
  } else if  (Discrepancy_essential$geneID[i] %in% Pb_all0$geneID & Discrepancy_essential$geneID[i] %in% Pf_all0$geneID & !( Discrepancy_essential$geneID[i] %in% Pk_all0$geneID)) {
    Discrepancy_essential$labels[i] <- "PbPf shared"
  } else if (!(Discrepancy_essential$geneID[i] %in% Pb_all0$geneID) & Discrepancy_essential$geneID[i] %in% Pf_all0$geneID & !( Discrepancy_essential$geneID[i] %in% Pk_all0$geneID)) {
    Discrepancy_essential$labels[i] <- "Pf specific"
  } else if (Discrepancy_essential$geneID[i] %in% Pb_all0$geneID & !(Discrepancy_essential$geneID[i] %in% Pf_all0$geneID) & !( Discrepancy_essential$geneID[i] %in% Pk_all0$geneID)) {
    Discrepancy_essential$labels[i] <- "Pb specific"
  } else {
    Discrepancy_essential$labels[i] <- "Pk specific"
  }
}
table(Discrepancy_essential$labels)

Discrepancy_dispensable <- data.frame(geneID=unique(append(append(Pf_all1$geneID,Pb_all1$geneID),Pk_all1$geneID)),
                                      labels=NA)
# Label geneIDs
for (i in 1:length(Discrepancy_dispensable$geneID)) {
  if (Discrepancy_dispensable$geneID[i] %in% Reduce(intersect, venn11)) {
    Discrepancy_dispensable$labels[i] <- "PkPbPf shared"
  } else if (Discrepancy_dispensable$geneID[i] %in% Pb_all1$geneID & Discrepancy_dispensable$geneID[i] %in% Pk_all1$geneID & !( Discrepancy_dispensable$geneID[i] %in% Pf_all1$geneID)) {
    Discrepancy_dispensable$labels[i] <- "PkPb shared"
  } else if (Discrepancy_dispensable$geneID[i] %in% Pk_all1$geneID & Discrepancy_dispensable$geneID[i] %in% Pf_all1$geneID & !( Discrepancy_dispensable$geneID[i] %in% Pb_all1$geneID)) {
    Discrepancy_dispensable$labels[i] <- "PkPf shared"
  } else if  (Discrepancy_dispensable$geneID[i] %in% Pb_all1$geneID & Discrepancy_dispensable$geneID[i] %in% Pf_all1$geneID & !( Discrepancy_dispensable$geneID[i] %in% Pk_all1$geneID)) {
    Discrepancy_dispensable$labels[i] <- "PbPf shared"
  } else if (!(Discrepancy_dispensable$geneID[i] %in% Pb_all1$geneID) & Discrepancy_dispensable$geneID[i] %in% Pf_all1$geneID & !( Discrepancy_dispensable$geneID[i] %in% Pk_all1$geneID)) {
    Discrepancy_dispensable$labels[i] <- "Pf specific"
  } else if (Discrepancy_dispensable$geneID[i] %in% Pb_all1$geneID & !(Discrepancy_dispensable$geneID[i] %in% Pf_all1$geneID) & !( Discrepancy_dispensable$geneID[i] %in% Pk_all1$geneID)) {
    Discrepancy_dispensable$labels[i] <- "Pb specific"
  } else {
    Discrepancy_dispensable$labels[i] <- "Pk specific"
  }
}
table(Discrepancy_dispensable$labels)

Discrepancy_essential2 <- left_join(Discrepancy_essential,orthologs_1on1_filtered3,by="geneID")
Discrepancy_dispensable2<- left_join(Discrepancy_dispensable,orthologs_1on1_filtered3,by="geneID")


Pf_MIS_MFS <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/MIS_MFS_Pf.xlsx", startRow=2)
Pf_MIS_MFS2 <- Pf_MIS_MFS[,c(2,7,8)]
colnames(Pf_MIS_MFS2)[1] <- "geneID"
Total.df2 <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_WTonly_afterDay15all.xlsx")
colnames(Total.df2)
colnames(MISvsMIS2)
colnames(MISvsMIS2)[1] <- "geneID"
test <- left_join(MISvsMIS2,Total.df2,by="geneID")
test2 <- left_join(test,Pf_MIS_MFS2,by="geneID")
test3 <- test2[,c(1:13,20,21)]
write.xlsx(test3,'./Output/MIS/Discrepant_MISvsMIS2.xlsx')
