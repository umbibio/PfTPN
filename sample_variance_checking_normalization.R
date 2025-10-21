library(tidyverse)
library(IRanges)
library(openxlsx)
library(dplyr)
library(data.table)
library(rtracklayer)
library(parallel)
library(UpSetR)
library(cowplot)
library(limma)


cm <- read.xlsx("./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon202412_202502_202504.xlsx")
colnames(cm)
cm_ex <- cm[,c(12:45)]

cm_ex[cm_ex == 0] <- NA

####Notice: loess normalization will make counts as negative, add pseudo count to make sure every count is non-negative
####Since normalize function in limma is applied to log-expression values, so it needs be transformed into log values and retransformed back to count
#MFS_cm3 <- normalizeCyclicLoess(MFS_cm2_cpm)
#add pseudo count 1 to avoid infinity after log transformation
cm_log2 <- log2(as.matrix(cm_ex)+1)
cm_after <- normalizeBetweenArrays(cm_log2, method = "cyclicloess")
cm_after <- 2^cm_after

cm_after[cm_after  == 0] <- NA

cm_ex <- as.data.frame(cm_ex)
cm_ex2<- cm_ex %>% 
  pivot_longer(everything(), names_to = "Samples", values_to = "Reads")

cm_after <- as.data.frame(cm_after)
cm_after2<- cm_after %>% 
  pivot_longer(everything(), names_to = "Samples", values_to = "Reads")

boxplot_before <- cm_ex2%>%
  ggplot() +
  aes(y = Reads, 
      x = Samples,
      fill = Samples) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5)+ylim(0, 150)+
  xlab("Samples") +
  ylab("Counts per site") +
  ggtitle("") + theme_bw() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),legend.background = element_blank())+theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 18,color = "black",family='sans'),
        axis.title = element_text(size = 20,color = "black",family='sans'),
        legend.text = element_text(size = 12))+
  theme(
    panel.border = element_rect(color = "black", fill = NA))+ylim(c(0,15))


boxplot_after <- cm_after2%>%
  ggplot() +
  aes(y = Reads, 
      x = Samples,
      fill = Samples) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5)+ylim(0, 150)+
  xlab("Samples") +
  ylab("Counts per site") +
  ggtitle("") + theme_bw() + theme(
    plot.title = element_text(color="black", size=14, face="bold"), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),legend.background = element_blank())+theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+theme(panel.grid = element_blank())+
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 18,color = "black",family='sans'),
        axis.title = element_text(size = 20,color = "black",family='sans'),
        legend.text = element_text(size = 12))+
  theme(
    panel.border = element_rect(color = "black", fill = NA))+ylim(c(0,15))

out.dir <- "./Output/Figures/"
ggsave(filename = paste(out.dir,"Boxplot_before_normalization", '.pdf',sep = ""), plot=boxplot_before,width = 10,height = 4, dpi = 300)
ggsave(filename = paste(out.dir,"Boxplot_between_arrays_log2transformation_cyclic", '.pdf',sep = ""), plot=boxplot_after,width = 10,height = 4, dpi = 300)

cm_after[is.na(cm_after)] <- 0
cm_after_tab <- cbind(cm[,1:11], cm_after)
write.xlsx(cm_after_tab, "./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon202412_202502_202504_after_loessnormalization.xlsx")

