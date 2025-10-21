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
library(colorspace)
library(cowplot)
library(rtracklayer)
library(ggrepel)
library(edgeR)
library(cowplot)



Total.df2 <- read.xlsx("./Output/HMS/HMS_202510.xlsx")
Total.df2 <- Total.df2[order(Total.df2$HM), ]             # sort from low to high by HM
Total.df2$geneIndex <- seq_len(nrow(Total.df2))   

p.MIS <- ggplot(Total.df2, aes(x=geneIndex, y= HM)) +
  geom_point(aes(colour = HM)) +
  labs(x = "Rank-ordered genes", y="HMS")+
  scale_colour_gradient2(low = muted("blue"), mid = "white",
                         limits = c(0, 1),
                         breaks = c(0.00, 0.50, 1.00),
                         high = "red" , midpoint = 0.5,  name = "HMS")+
  ggtitle('Hybrid model score (HMS)') + scale_x_continuous(breaks=seq(0, 5000, 2500))+ylim(c(0,1))

p.MIS+theme(
  plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
  legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
  axis.text = element_text(size = 12),  axis.title=element_text(size=14), legend.background = element_blank())+theme_cowplot()

p.MIS.background <- ggplot(Total.df2, aes(x=geneIndex, y= HM)) +
  geom_point(color='grey') +
  labs(x = "Rank-ordered genes", y="HMS")+
  ggtitle('Hybrid model score (HMS)') + scale_x_continuous(breaks=seq(0, 5000, 2500)) + theme(
    plot.title = element_text(color="black", size=14, face="bold"), legend.position = c(0.15, 0.8), 
    legend.key = element_rect(fill = "transparent", colour = "transparent"), legend.text = element_text(size=12),
    axis.text = element_text(size = 16),  axis.title=element_text(size=16), legend.background = element_blank())+ylim(c(0,1))

wronggenes <- read.xlsx("./Input/WrongGenes_PfJAMIS.xlsx")
wronggenes.df <- Total.df2[Total.df2$geneID %in% wronggenes$Gene_ID, ]
geneName.df <- wronggenes%>%dplyr::select(Gene_ID,GeneName, wrongtype)
colnames(geneName.df)[1] <- "geneID"
wronggenes.df <- left_join(wronggenes.df, geneName.df, by="geneID")
wronggenes1 <- wronggenes.df%>% dplyr::filter(wrongtype=="KO.essential")
wronggenes2 <- wronggenes.df%>% dplyr::filter(wrongtype=="KO.dispensable")

Pf.HMS.SP <- p.MIS.background + 
  geom_point(data=Total.df2[Total.df2$geneID %in% wronggenes2$geneID, ], aes(x=geneIndex, y=HM),
             shape = 21, colour = "black", fill = "#2e832c", size = 5, stroke = 0.5)+
  geom_point(data=Total.df2[Total.df2$geneID %in% wronggenes1$geneID, ], aes(x=geneIndex, y=HM),
             shape = 21, colour = "black", fill = "red", size = 5, stroke = 0.5) +
  geom_text_repel(data=wronggenes.df, aes(x=geneIndex, y=HM, label=GeneName, color=wrongtype), size=4, box.padding = unit(0.6, "lines"),
                  segment.linetype=2,
                  max.overlaps = Inf,
                  show.legend=F,
                  nudge_x=100,
                  nudge_y=-0.03,
                  #0 indicates left alignment, 0.5 indicates center alignment, and 1 indicates right alignment. Adjusting hjust allows you to control how the labels are positioned horizontally relative to the data points.
                  hjust = 0,
                  min.segment.length = 0,
                  force=1,
                  fontface="italic",
                  family="sans")+
  scale_color_manual(values = c("KO.essential" = "red", "KO.dispensable" = "#2e832c")) 
Pf.HMS.SP 
ggsave(filename = "./Output/Figures/13Wronggenes_PfSP_final.pdf", plot=Pf.MIS.SP, width = 4,height = 4, dpi = 300)
