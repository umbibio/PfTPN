library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyverse)
library(readxl)
library(rtracklayer)
library(UpSetR)

#####Sheet 1 is gene list for essential genes accosiated with Apicoplast
gene.table <- read.xlsx("./Input/ref_gene_list.xlsx", sheet = 1)
table(gene.table$Category)
######To merge with the Pk, Pf and Pb essentialome
scores <- read.xlsx("/Users/sidaye/Documents/R/Tnseq/webapp/HMS_MFS_regression_trending_results_pcgenes_loess_normalization2.xlsx")
colnames(gene.table)[1] <- "GeneID.Pf_3D7"



merged_table <- left_join(gene.table, scores, by="GeneID.Pf_3D7")
table(merged_table$Category)
#DISPENSABLE  ESSENTIAL      UNKNOWN 
#17          38          47 

merged_table2 <- merged_table%>%dplyr::filter(HMS<0.26&Pf.MIS<0.5&Pb.Relative.Growth.Rate<0.2) #25 essential genes among Pk, Pf and Pb essentialomes
write.xlsx(merged_table2,"./Output/backgroup_gene_list.xlsx")

######Extract info from gff file
gff_file <- "./Input/New_Plasmodium_GFF.gff"
# Read the GFF file
gff_data <- import(gff_file, format = "gff")
######Extract CDS info
gff_data2 <- gff_data[gff_data$type=="CDS"]
######Extract the genes of interest
gff_data3 <- gff_data2[gff_data2$gene_id%in%merged_table$GeneID.Pf_3D7]


#Define the path to  Excel file
file_path <- "./Input/D24 Site Count Data-Trial 2,All Replicates.xlsx"
# Get the names of all sheets
sheet_names <- excel_sheets(file_path)
# Read all sheets into a list of data frames
all_sheets <- lapply(sheet_names, function(sheet) read_excel(file_path, sheet = sheet,col_names=F))
# Optionally, name the list elements by sheet names
names(all_sheets) <- sheet_names
# View the result (prints the first few rows of each sheet)
names(all_sheets)

n <- length(names(all_sheets))

for(i in 1:n){
  API_sheet <- all_sheets[[i]]
  colnames(API_sheet) <- c("V1","V2","V3","V4")
  API_sheet <- API_sheet%>%mutate(V3=ifelse(V3=="F","+","-"))
  API_sheet2 <- API_sheet%>%mutate(V1=V1,
                                   V2=V2,
                                   V3=V3,
                                   V4=V4)
  colnames(API_sheet2)[4] <- names(all_sheets)[i]
  all_sheets[[i]] <- API_sheet2
  
}

#API_sheet3 <- API_sheet2[,c(1,2,5,3,4)]
# Step 2: Merge all tables on V1, V2, V3
merged_table <- Reduce(function(x, y) merge(x, y, by = c("V1", "V2", "V3"), all = TRUE), all_sheets)


#To check distribution
tmp1 <- tmp1 %>% dplyr::mutate(UII = c(1: nrow(tmp1)))
histo.freq <- as.data.frame(table(tmp1[, 4]))
#rule out the row with 0 reads
histo.freq <- histo.freq[-1, ]
histo.freq$Freq <- as.numeric(histo.freq$Freq)
histo.freq <- histo.freq[1:60,]
histo.freq.figure <- ggplot(histo.freq, aes(x=Var1, y=Freq)) + 
  xlab('Number of reads mapped to UII') + 
  geom_bar(stat='identity') +
  ylab('Count')+
  ggtitle(samplesNamePk[i,])+
  ylim(0, 20000)
print(histo.freq.figure)






#To check complexity

#merged_table$Total <- rowSums(merged_table[,(4:ncol(merged_table))])
#####To check complexity among replicates and different treatments
#####Upset, transposon count matrix 
mat <- merged_table
#filter out rows that not present in any samples
#index_present <- which(mat$Present_in_any_samples == 'yes')
#mat_present <- mat[index_present,]

#only remain samples counts
#number of samples
mat_present <- mat[, 4:ncol(mat)]


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
#mat_final$Total <- rowSums(mat_final)
#mat_final2 <- mat_final[mat_final$Total!=0,]
#names(mat_final) <- c( "PkTPN23_Day14", "PkTPN23_Day14_REPEAT_index_cycles10", "PkTPN18_Day9_preLib_prep_SbfI_NotI_digested", "PkTPN18_Day9_No_SbfI_NotI_digestion", "PkTPN18_Day_Day9_SbfI_NotI_digested", "PkTPN18_Day9_NotI_digested", "merge_PkTPN23_Day14")
upset(mat_final, empty.intersections = "on", nsets = 11, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")


#####Have double check in IGV how to transform the loci and match them in IGV with theoretical TTAA#####
merged_table2 <- merged_table%>%mutate(V4=ifelse(V3=="+",V2+4-1, V2-4+1))
merged_table2 <- merged_table2[,c(1,2,ncol(merged_table2),4:ncol(merged_table2)-1)]
merged_table3 <- merged_table2%>%mutate(start=ifelse(V3=="+", V2, V4))
merged_table3 <- merged_table3%>%mutate(end=ifelse(V3=="+", V4, V2))
merged_table3 <- merged_table3%>%arrange(V1,start,end)
merged_table3 <- merged_table3[,c(1,2,3,4)]
merged_table4 <- merged_table3%>%dplyr::filter(grepl("Pf3D7", V1))
merged_table5 <- data.frame(V1=merged_table4$V1,
                            V2=merged_table4$V2,
                            V3=merged_table4$V4,
                            V4="Insert_site",
                            V5=0,
                            V6=merged_table4$V3)
options(scipen = 999)
merged_table6 <- merged_table5 %>%
  mutate(
    V2_temp = ifelse(V6 == "-", V3, V2),  # Temporarily store V3 in V2 when V6 == "-"
    V3 = ifelse(V6 == "-", V2, V3),       # Replace V3 with V2 when V6 == "-"
    V2 = V2_temp                          # Replace V2 with V3 stored in V2_temp
  ) %>%
  select(-V2_temp)  # Remove the temporary column
merged_table6$V2 <- as.numeric(merged_table6$V2)
merged_table6$V3 <- as.numeric(merged_table6$V3)
###gtf format to bed format
merged_table7 <- merged_table6%>%dplyr::mutate(V2=ifelse(V6=="+",V2-1,V2))
merged_table7 <- merged_table7%>%dplyr::mutate(V3=ifelse(V6=="+",V3-1,V3))

write.table(merged_table7, './Output/Original_countmatrixmerged_insertionsite_QC.bed', col.names = F, quote = F, row.names = F, sep = '\t')


####################library complexity, Upset checking######################
cm <- read.xlsx("./Output/count_matrix/Pf_count_matrix_202412_202502_202504_202505.xlsx") 
cm2 <- cm[,c(12:(ncol(cm)-1))]
cm2 <- cm2[,c(17,18,35,36)] ### D0
cm2 <- cm2[,c(19,20,21,22,38,39,40,41)] ### D14
cm2 <- cm2[,c(23,24,42,43)] ### D17
cm2 <- cm2[,c(25,26,45,46)] ### D18
cm2 <- cm2[,c(27,28,29,30,48,49,50,51)] ### D20

cm2 <- cm2[,c(44,47,52)] ### Acet samples
colnames(cm2)

cm2$Total <- rowSums(cm2)
cm2 <- cm2[cm2$Total!=0,]
cm2 <- cm2[,-ncol(cm2)]
mat_final <- function_createuniquelist(cm2)
#head(cm)
####Not include total column
cm <- cm%>%filter(Present_in_any_samples=="yes")
mat_final2 <- function_createuniquelist(cm[,c(12:ncol(cm))])
####Include total column
#mat_final2 <- function_createuniquelist(cm[,c(12:(ncol(cm)-1))])
colnames(mat_final2)
colnames(mat_final2) <- paste0(colnames(mat_final2), "_PE")
mat_final_PE <- mat_final2


upset(mat_final, empty.intersections = "on", nsets = 16, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")

upset(mat_final1112, empty.intersections = "on", nsets = 10, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")
upset(mat_final, empty.intersections = "on", nsets = 16, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")

mat_final3 <- mat_final2[,c(1,2,15:22)]
upset(mat_final3, empty.intersections = "on", nsets = 10, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")

mat_final4 <- mat_final2[,c(1:14)]
upset(mat_final4, empty.intersections = "on", nsets = 14, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")

mat_final_D0 <- mat_final2[,c(1:2)]
upset(mat_final_D0, empty.intersections = "on", nsets = 2, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")

mat_final_D12 <- mat_final2[,c(3:6)]
upset(mat_final_D12, empty.intersections = "on", nsets = 4, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")

mat_final_D18 <- mat_final2[,c(7:14)]
upset(mat_final_D18, empty.intersections = "on", nsets = 8, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")

mat_final_D20 <- mat_final2[,c(15:22)]
upset(mat_final_D20, empty.intersections = "on", nsets = 8, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")

###############202411 batch##########
upset(mat_final2, empty.intersections = "on", nsets = 16, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")

mat_final3 <- mat_final2[,c(1,2,23)]
upset(mat_final3, empty.intersections = "on", nsets = 3, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")

mat_final4 <- mat_final2[,c(23:29)]
upset(mat_final4, empty.intersections = "on", nsets = 7, text.scale = c(1.5, 1.5, 1.3, 1, 1.5, 1.5),show.numbers = "no")
#table(cm$D0_Rep1)
cm$total <- rowSums(cm[,c(12:ncol(cm))])
hist(cm$total[cm$total!=0], nclass=1000)
quantile(cm$total)
