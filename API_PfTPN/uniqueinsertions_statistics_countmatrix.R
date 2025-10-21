library(tidyverse)
library(IRanges)
library(Biostrings) 
library(ShortRead)
library(openxlsx)
library(dplyr)
library(data.table)
library(rtracklayer)
library(parallel)
library(UpSetR)
library(cowplot)

#cm <- read.xlsx("./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon202412.xlsx") ###Double check if includes total columns
cm <- read.xlsx("./Output/count_matrix/Pf_count_matrix_all_Bg_removed_siteslevel.xlsx")

####WT samples after Day15
TM <- c(10,18,20,29,31,33,34,48,51,54,55,60,62,64,66,
        68,72,76,78,80,86,88,90,96,98,101,107,110,114,115,
        126,129,136,138,140)
####Mev samples after Day15
TM <- c(11,17,19,30,32,35,49,52,56,57,59,61,63,65,67,69,
  73,77,79,81,87,89,91,97,98,109,112,118,119,128,
  131,137,139,141)

TM2 <- TM+6

CM_WT<- cm[,c(1:12,TM2)]

UniqueIn <- function(cm){
  cm$Total <- rowSums(cm%>%dplyr::select(contains('TPN')))
  cm$Total <- round(cm$Total)
  num <- sum(cm$Total!=0)
  #sum(cm$Total != 0, na.rm = TRUE)
  return(num)
}
All_uniquein <- UniqueIn(cm)
WTafterD15_uniquein <- UniqueIn(CM_WT)
print(WTafterD15_uniquein)
sum(cm$Total != 0, na.rm = TRUE)
sum(is.na(CM_WT$Total))

cm <- cbind(cm1,cm2[,12:ncol(cm2)])
cm <- cbind(cm,cm3[,12:ncol(cm3)])

before
cm3<- cm[,c(13:18,20:ncol(cm))]

#cm$total <- rowSums(cm%>%dplyr::select(contains("pTN9")))
#table(cm$Present_in_any_samples)
#sum(cm$total)
###Final usable reads per samples
#colSums(cm%>%dplyr::select(contains("pTN9")))

mat_present <- cm[, 12:ncol(cm)]
colSums(mat_present)
#mat_present <- cm[, 13:87]
###Unique insertions per samples
function_createuniquelist <- function(x){
  for (i in 1:ncol(x)){
    #index <- which(x[,i]!=0)
    index <- is.na(x[,i])
    x[,i][index] <- 0
    index2 <- which(x[,i]!=0)
    x[,i][index2] <- 1
  }
  return(x)
}



mat_final <- function_createuniquelist(mat_present)
colSums(mat_final) #unique insertions coverages
upset(mat_final, empty.intersections = "on", nsets = 16, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")
####WT

mat_WT <- mat_final[,c(1,2,4,6,8,10,12,14,15)]
mat_WT$total_WT <- rowSums(mat_WT)
mat_Mev <- mat_final[,c(3,5,7,9,11,13,16)]
mat_Mev$total_Mev <- rowSums(mat_Mev)

mat_final_all <- mat_final 
mat_final_all$total <- rowSums(mat_final)

######exclude D0 to calculate the unique insertions coverage per sample
mat_final_all2 <- mat_final[,c(2:7,9:ncol(mat_final))]
mean(colSums(mat_final_all2))

df3 <- data.frame(
  Sample = c("All_WT","All_Mev","All", "Pf"),
  unique =c(sum(mat_WT$total_WT!=0),sum(mat_Mev$total_Mev!=0),sum(mat_final_all$total!=0), 38173*2)
)

ggplot(df3, aes(x = reorder(Sample, unique), y = unique)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label= unique), vjust=-0.5, size=5)+
  labs(title = " ",
       x = "Sample",
       y = "Unique insertions") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  scale_y_continuous(labels = scales::comma)  # Format y-axis numbers

####After bg noise removal, contains overlaped rows
cm <- read.xlsx("./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon202412_202502.xlsx")
cm <- read.xlsx("./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon202412_202502_Bg_removed_siteslevel.xlsx")
cm <- cm[c(1:660272),] ###overlaps=5061*2
mat_present <- cm[, 13:ncol(cm)]
mat_present <- cm[, 12:(ncol(cm)-1)] ###before bg correction
dim(mat_present)
#mat_present[,] <- round(mat_present)
mat_present$Total <- rowSums(mat_present)
mat_present$Total <- ifelse(mat_present$Total<0.8,0,mat_present$Total)

mat_WT <- mat_present[,c(1,2,4,6,8,10,12,14,15)]
mat_WT$total_WT <- rowSums(mat_WT)
mat_WT$total_WT <- ifelse(mat_WT$total_WT<0.8,0,mat_WT$total_WT)
mat_Mev <- mat_present[,c(3,5,7,9,11,13,16)]
mat_Mev$total_Mev  <- rowSums(mat_Mev)
mat_Mev$total_Mev <- ifelse(mat_Mev$total_Mev<0.8,0,mat_Mev$total_Mev)

df3 <- data.frame(
  Sample = c("All_WT","All_Mev","All", "Pf"),
  unique =c(sum(mat_WT$total_WT>0),sum(mat_Mev$total_Mev!=0),sum(mat_present$Total!=0), 38173*2)
)

ggplot(df3, aes(x = reorder(Sample, unique), y = unique)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label= unique), vjust=-0.5, size=5)+
  labs(title = " ",
       x = "Sample",
       y = "Unique insertions") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  scale_y_continuous(labels = scales::comma)  # Format y-axis numbers


df3 <- data.frame(
  Sample = c("PkTPN","PfTPN","PfTPN2.0"),
  unique =c(101515,1170,37576)
)

ggplot(df3, aes(x = reorder(Sample, unique), y = unique)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label= unique), vjust=-0.5, size=5)+
  labs(title = " ",
       x = "Sample",
       y = "Unique insertions per sample") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  scale_y_continuous(labels = scales::comma)  # Format y-axis numbers

#####Usable reads comparisons
df4 <- data.frame(
  Sample = c("PkTPN","PfTPN","PfTPN2.0"),
  usable_reads =c(1209198,31417,169791)
)

ggplot(df4, aes(x = reorder(Sample, usable_reads), y = usable_reads)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label= usable_reads), vjust=-0.5, size=5)+
  labs(title = " ",
       x = "Sample",
       y = "Usable reads per sample") +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  scale_y_continuous(labels = scales::comma)  # Format y-axis numbers


##########How many CDS has no insertions across all the samples############
MIS <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_WTonly_afterDay15all.xlsx")
MIS <- read.xlsx("./Output/MIS/MIS_readyforplot_bgremoved_WT_alltimepoints.xlsx")


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


