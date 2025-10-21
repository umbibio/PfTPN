library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(doParallel)
library(cowplot)
#####!!!please do not use exon converted version of tranposon matrix
Pk.tcm <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202311_Novaseq/Output/transposon_matrix/all/transposon_count_matrix75essentialomeonly_run13_bgremoved_cpm.xlsx")
Pk.tcm$Total <- round(Pk.tcm$Total)
totalgene <- length(unique(append(Pk.tcm$sense_geneID,Pk.tcm$antisen_geneID)))-1
####Number of genes with TTAA
#######Number of genes with 0, 1, 2, +3 unique insertions covered##############
NoGeneUniqueInsertionCov <- function(tcm, totalgene){
  tcm2 <- tcm%>%dplyr::filter(Present_in_any_samples=="yes")
  tcm2 <- tcm2%>%dplyr::filter(Assigned_location=="exon")
  sense_geneID <- as.data.frame(table(tcm2$sense_geneID))
  antisen_geneID<- as.data.frame(table(tcm2$antisen_geneID))
  all <- rbind(sense_geneID,antisen_geneID)
  all <- all[!is.na(all$Var1), ]
  zeroFreq <- totalgene-length(unique(all$Var1))
  freq_1 <- sum(all$Freq == 1)
  freq_2 <- sum(all$Freq == 2)
  freq_3 <- sum(all$Freq == 3)
  freq_gt3 <- sum(all$Freq> 3)
  dt <- data.frame(geneNo=c("0","1","2","3",">3"),
                   Freq=c(zeroFreq,freq_1 <- sum(all$Freq == 1),freq_2 <- sum(all$Freq == 2),freq_3 <- sum(all$Freq == 3),freq_gt3 <- sum(all$Freq> 3))
    
  )
  return(dt)
}
Pk_sta <- NoGeneUniqueInsertionCov(Pk.tcm, totalgene)
Pk_sta$geneNo <- factor(Pk_sta$geneNo, levels = c("0", "1", "2", "3", ">3"))
pk <- ggplot(Pk_sta, aes(x = geneNo, y = Freq)) +
  geom_bar(stat = "identity", fill = "#9F7FBF") +
  labs(
    title = "",
    x = "Unique insertions within CDS",
    y = "Number of genes"
  ) + theme_cowplot()+ylim(0, 5000)
ggsave(filename = "./Output/Figures/Pk_IIE.pdf",
       plot = pk, 
       width = 4, height = 4, 
       dpi = 300)


#####################Pf.SP######################Can not be exon converted transposon matrix###############
Pf.SP.tcm <- read.xlsx("./Output/transposon_matrix/all/transposon_count_matrix_202412_202502_202504_202505_bgremoved_cpm.xlsx")
Pf.SP.tcm_WT <- Pf.SP.tcm[,c(1:6,10,18,20,29,31,33,34,48,51,54,55)] ### WT samples after Day15
Pf.SP.tcm_WT <- Pf.SP.tcm[,c(1:6,10,11,17,18,19,20,29,30,31,32,33,34,35,36,48,49,50,51,52,53,54,55,56,57)] ### WT samples after Day15
Pf.SP.tcm_WT$Total <- rowSums(Pf.SP.tcm_WT[,7:ncol(Pf.SP.tcm_WT)])
Pf.SP.tcm_WT$Total <- round(Pf.SP.tcm_WT$Total)
Pf.SP.tcm_WT <- Pf.SP.tcm_WT%>%mutate(Present_in_any_samples=ifelse(Total==0,"no","yes"))
table(Pf.SP.tcm_WT$Present_in_any_samples)

totalgene_Pf.SP <- length(unique(append(Pf.SP.tcm_WT$sense_geneID,Pf.SP.tcm_WT$antisen_geneID)))-1
Pf.SP_sta <- NoGeneUniqueInsertionCov(Pf.SP.tcm_WT, totalgene_Pf.SP)
Pf.SP_sta$geneNo <- factor(Pf.SP_sta$geneNo, levels = c("0", "1", "2", "3", ">3"))

pf.SP <- ggplot(Pf.SP_sta, aes(x = geneNo, y = Freq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "",
    x = "Unique insertions within CDS",
    y = "Number of genes"
  ) + theme_cowplot()+ylim(0, 5000)
pf.SP
ggsave(filename = "./Output/Figures/Pf_SP_IIE.pdf",
       plot = pf.SP, 
       width = 4, height = 4, 
       dpi = 300)

#####################Pf.JA######################
Pf_MIS_MFS <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/202212_Novaseq/Input/MIS_MFS_Pf.xlsx", startRow=2)
Pf.SP.JA <- as.data.frame(table(Pf_MIS_MFS$`#.of.insertions.in.CDS`))
zeroFreq <- Pf.SP.JA$Freq[Pf.SP.JA$Var1==0]
freq_1 <- Pf.SP.JA$Freq[Pf.SP.JA$Var1==1]
freq_2 <- Pf.SP.JA$Freq[Pf.SP.JA$Var1==2]
freq_3 <- Pf.SP.JA$Freq[Pf.SP.JA$Var1==3]
freq_gt3 <- sum(Pf.SP.JA$Freq[Pf.SP.JA$Var1[5:nrow(Pf.SP.JA)]])
Pf_JA_sta <- data.frame(geneNo=c("0","1","2","3",">3"),
                 Freq=c(zeroFreq,freq_1 ,freq_2 ,freq_3 ,freq_gt3))

Pf_JA_sta$geneNo <- factor(Pf_JA_sta$geneNo, levels = c("0", "1", "2", "3", ">3"))

pf.JA <- ggplot(Pf_JA_sta, aes(x = geneNo, y = Freq)) +
  geom_bar(stat = "identity", fill = "#87BFBF") +
  labs(
    title = "",
    x = "Unique insertions within CDS",
    y = "Number of genes"
  ) + theme_cowplot()+ylim(0, 5000)
ggsave(filename = "./Output/Figures/Pf_JA_IIE.pdf",
       plot = pf.JA, 
       width = 4, height = 4, 
       dpi = 300)

Pk_sta$Group<- "Pk"
Pf.SP_sta$Group <- "Pf.SP"
Pf_JA_sta$Group <- "Pf.JA"
Pf_combined <- rbind(Pk_sta, Pf.SP_sta, Pf_JA_sta)
Pf_combined$geneNo <- factor(Pf_combined$geneNo, levels = c("0", "1", "2", "3", ">3"))                 

combined <- ggplot(Pf_combined, aes(x = geneNo, y = Freq, fill = Group,group = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width=0.8) +
  scale_fill_manual(values = c("Pk"="#9F7FBF","Pf.JA" = "#87BFBF", "Pf.SP" = "#4682B4")) +
  labs(
    title = "",
    x = "Unique insertions within CDS",
    y = "Number of genes",
    fill = "Group"
  ) +
  ylim(0, 5000) +
  theme_cowplot()+theme(legend.position = c(0.7, 0.9))

ggsave(filename = "./Output/Figures/Pf_JA_Pf_SP_Pk_IIE.pdf",
       plot = combined, 
       width = 6, height = 4, 
       dpi = 300)

