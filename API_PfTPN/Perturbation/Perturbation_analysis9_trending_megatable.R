library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggVennDiagram)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(scales)
library(ggpmisc)
library(edgeR)
library(venneuler)
library(VennDiagram)
library(grid)
library(EnhancedVolcano)
library(patchwork)

Total.df2 <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_WTonly_afterDay15all.xlsx")

###########merge megatable for drug perturbation###########
###########merge megatable for drug perturbation###########
###########merge megatable for drug perturbation###########
#megatable_drugs
input.dir.edgeR <- "./Output/Perturbation_CDS/trending_analysis/"
#input.dir.vc <- "./Output/Perturbation_CDS/Perturbation_analysis_all_model3_siteslevel_Bg_removed_Variation_coefficient/VC_batch_processing_CPM_normalized/SetA_and_SetB/"
#list.files can list the name of files 
count.files <- list.files(input.dir.edgeR)
#count.files2 <- list.files(input.dir.vc)

Mev.files <- count.files[grep("Differential_insertions", count.files)]

####order the names
Mev.files <- Mev.files[c(4,1,2,3)]

merge_megatable <- function(files, Total.df2, input.dir1){
  megatable <-Total.df2
  for(i in 1:length(files)){
    file <- read.xlsx(paste0(input.dir1,files[i]))
    file <- file[,c(1:9)]
    colnames(file)[2:9] <- paste0(colnames(file)[2:9], "_time",i)
    megatable <- left_join(megatable, file, by="geneID")
  }
  return(megatable)
}

Mev_megatable <- merge_megatable(Mev.files,Total.df2, input.dir.edgeR)
write.xlsx(Mev_megatable,"./Output/Perturbation_CDS/megatable/mev_megatable_202510.xlsx")

#######################Input megatables###############
Mev_megatable <- read.xlsx("./Output/Perturbation_CDS/megatable/mev_megatable_202510.xlsx")


###########To fit a linear regression model on multiple time points#########
time <- 1:4

get_slope_p <- function(row) {
  y <- as.numeric(row[c("logFC_time1","logFC_time2","logFC_time3","logFC_time4")])
  df <- data.frame(logFC = y, time = time)
  
  # remove rows with NA
  df <- na.omit(df)
  
  if (nrow(df) < 3) {
    return(c(slope = NA, pval = NA))  # not enough points
  }
  
  fit <- lm(logFC ~ time, data = df)
  coef_summary <- summary(fit)$coefficients
  slope <- coef_summary["time", "Estimate"]
  pval  <- coef_summary["time", "Pr(>|t|)"]
  return(c(slope = slope, pval = pval))
}

results <- t(apply(Mev_megatable, 1, get_slope_p))
results_df <- data.frame(geneID = Mev_megatable$geneID, results)

head(results_df)

results_df_unique <- results_df[!duplicated(results_df$geneID), ]
Mev_megatable2 <- left_join(Mev_megatable,results_df_unique,by="geneID")

Mev_megatable_up <- Mev_megatable2%>%dplyr::filter(slope>0&pval<0.05)
Mev_megatable_down <- Mev_megatable2%>%dplyr::filter(slope<0&pval<0.01)
write.xlsx(Mev_megatable2,"./Output/Perturbation_CDS/megatable/mev_megatable_202510_all.xlsx")
################Downstream is just for record#####################
################Downstream is just for record#####################
################Downstream is just for record#####################
###Go to Log2FC_trending_visualization.R########
######################Input any genes you what to check
geneName <- c("PF3D7_1322700","PF3D7_1350500","PF3D7_0621000","PF3D7_1005000","PF3D7_0616400","PF3D7_1239500","PF3D7_0608100","PF3D7_0806100","PF3D7_0711000")
geneName <- c("PF3D7_1214200","PF3D7_1340000","PF3D7_0219800")
extract_day <- function(time_str) {
  gsub("time", "", time_str)
}


#####Input of apicoplast-associated genes and known essential genes
gene.table <- read.xlsx("./Input/ref_gene_list.xlsx", sheet = 4) #sheet4 is known essential genes but no apicoplast
geneName <- gene.table$X1

gene.table <- read.xlsx("./Input/ref_gene_list.xlsx", sheet = 1) #sheet1 is WT Pf api-associated genes
gene.table <- gene.table%>%dplyr::filter(Category=="ESSENTIAL ")
gene.table <- gene.table%>%dplyr::filter(Category=="DISPENSABLE")
gene.table <- gene.table%>%dplyr::filter(Category=="UNKNOWN")
geneName <- gene.table$Gene.ID
######Only for drugs
trending_plot <- function(Mev_megatable, geneName){
  megatable <- Mev_megatable %>% filter(geneID %in% geneName)
  #####Be careful about the names of columns, since the names have been changed many times!!!!!
  df <- megatable %>% dplyr::select(contains("logFC"))
  #df2 <- megatable %>% dplyr::select(contains("mean_log2_FC_sites"))
  #df2 <- megatable %>% dplyr::select(contains("mean(log2FC_sites)"))
  #df3 <- megatable %>% dplyr::select(contains("mean_FC_sites"))
  df_plot <- data.frame(geneID=rep(megatable$geneID,each=ncol(df)),
                        Time=rep(unlist(lapply(strsplit(colnames(df), '_'), '[[', 2)), nrow(df)),
                        #Time=rep(unlist(lapply(strsplit(colnames(df), '_'), '[[', 3)), nrow(df)), #when Human vs rhesus needs to be 3
                        #cond=rep(unlist(lapply(strsplit(colnames(df), '_day'), '[[', 1)), nrow(df)),
                        logFC_edgeR=c(as.vector(t(df))))
  
  #df_plot$log2_mean_FC_sites <-  log2(df_plot$mean_FC_sites)
  df_plot$timeNumber <- as.numeric(sapply(df_plot$Time, extract_day))
  df_plot <- df_plot %>%
    arrange(geneID, timeNumber)
  df_plot$geneID <- factor(df_plot$geneID, levels = unique(df_plot$geneID))
  df_plot$Time <- factor(df_plot$Time, levels = unique(df_plot$Time))
  #df_plot$cond <- factor(df_plot$cond, levels = unique(df_plot$cond))
  return(df_plot)
}


df_plot <- trending_plot(Mev_megatable=Mev_megatable2, geneName)


#df_plot<- df_plot %>%arrange(geneID, match(Time, c("day3", "day4", "day6", "day9","day15")))
#df_plot$Time <- factor(df_plot$Time, levels = unique(df_plot$Time))
df_plot$Time <- factor(df_plot$Time, levels = c("time1", "time2", "time3","time4"))
#######For edgeR gene level LogFC model#######
#######For edgeR gene level LogFC model#######
#######For edgeR gene level LogFC model#######

df_plot$logFC_edgeR[is.na(df_plot$logFC_edgeR)] <- 0
p1 <- ggplot(df_plot,aes(x = Time, y=logFC_edgeR, group = geneID))+
  geom_line(size=1) +
  geom_point(size=1)+
  theme_bw() +
  facet_wrap(geneID~., scales = 'free')+theme(strip.background = element_rect(colour = "black", fill = "white",
                                                                              size=1))+
  labs(y = "log2FC_edgeR")
plot(p1)



# 4) plot WITHOUT selectLab
EnhancedVolcano(
  Mev_megatable2,
  title = "",
  subtitle = "",
  lab = tab3$plotLabel,      # <- use the custom labels
  pCutoff = 0.05,
  FCcutoff = 2,
  x = "slope",
  y = "pval",
  legendLabels = c("NS", "Fold Change", "FDR", "FDR & Fold Change"),
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "black",
  ylab = "-log10(FDR)",      # (label text only; your y is PValue)
  xlab = "log2FC",
  max.overlaps = 200,        # or Inf to avoid dropping labels
  ylim = c(0, 3)
)




