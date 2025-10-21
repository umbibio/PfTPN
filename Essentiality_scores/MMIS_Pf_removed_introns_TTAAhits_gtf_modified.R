library(tidyverse)
library(openxlsx)
library(dplyr)
library(data.table)
library(ggplot2)
library(GenomicRanges)
library(GenomicFeatures)
library(syntenet)
##V4 is original TTAA_ID
##V2 is modified TTAA_ID after removing introns
##V10 is modified loci of start of transcript after removing introns
##V11 is modified loci of end of transcript after removing introns
##V13 is the strandness of genes

TTAAhits_R_gtf_include_contigs <- read.table('./Output/TTAAhitsPf/TTAAhits_R_gtf_include_contigs.bed', sep = '\t')


##use gff2GRangesList() to read GFF/GFF3/GTF files in a directory as a GRangesList object
# Path to directory containing FASTA files
gff_dir <- "./Input/Genome/Pf_3D758_gff"
dir(gff_dir) # see the contents of the directory
# Read all gff files in gff_dir and create a GRangesList object
gr <- gff2GRangesList(gff_dir)
# Convert GrangesList to Granges
gr <- gr$`PlasmoDB-63_Pfalciparum3D7`

#Extract all transcript rows or mRNA rows
#transcripts <- gr[gr$type== "transcript"]
transcripts <- gr[gr$type== "mRNA"]
transcripts  <- sort(transcripts )
transcripts_df <- as.data.frame(transcripts)
transcripts_df <- data.frame(chr=transcripts_df$seqnames,
                             transcript_start=transcripts_df$start,
                             transcript_end=transcripts_df$end,
                             trans_width=transcripts_df$width,
                             trans_strand=transcripts_df$strand,
                             gene_id=transcripts_df$gene_id)

#Extract all exons rows
exons <- gr[gr$type== "exon"]
exons <- sort(exons)

exons_df <- as.data.frame(exons)
exons_df <- exons_df %>% group_by(gene_id)%>%mutate(exon_num=n(),
                                                    exon.ord = ifelse(strand  == '+', 1:n(), seq(n(), 1, by = -1)),
                                                    multiple.exon = ifelse(n() > 1, T, F),
                                                    exon_ID=paste(gene_id,"exon",exon.ord,sep = "_"),
                                                    exon_loci=paste(gene_id,start,end, sep = ":"))

exons_intron_df <- exons_df %>% ungroup()%>%
  group_by(gene_id) %>%
  mutate(
    intron_start = lag(end) + 1,
    intron_end = lead(start) - 1
  )%>% mutate(intron_end=lag(intron_end)) %>% mutate(intron_width=intron_end - intron_start + 1) %>%
  mutate(intron.ord = ifelse(strand == '+', lag(exon.ord), exon.ord)) 

exons_intron_df$intron.ord[is.na(exons_intron_df$intron_width)] <- NA

exons_intron_df <- left_join(exons_intron_df,transcripts_df, by="gene_id")


length(unique(exons_df$gene_id))
genelist <- unique(exons_df$gene_id)


#####To stitched the exons by two categories:1 genes have only 1 exon,2 genes have multiple exons

removed_introns_TTAAhits_gtf_exons_extracted_modified <- function(TTAAhits_R_gtf_include_contigs, exons_intron_df, transcripts_df){
  exons_df_1exon <- exons_intron_df %>%dplyr::filter(multiple.exon==F)%>% ungroup()
  ###All genes
  genelist_1exon_all <- exons_df_1exon$gene_id
  exons_df_mutiexon <- exons_intron_df %>%dplyr::filter(multiple.exon==T)%>% ungroup()
  ###All genes
  genelist_mutiexon_all <- unique(exons_df_mutiexon$gene_id)
  TTAAhits_R_gtf_include_contigs$V15 <- lapply(strsplit(TTAAhits_R_gtf_include_contigs$V15,"gene_id "),"[[",2)
  TTAAhits_R_gtf_include_contigs$V15 <-gsub(";","",TTAAhits_R_gtf_include_contigs$V15)
  TTAAhits_R_gtf_exon <-TTAAhits_R_gtf_include_contigs%>%dplyr::filter(V9=="exon")
  TTAAhits_R_gtf_exon <- TTAAhits_R_gtf_exon%>%mutate(exon_loci=paste(V15,V10,V11, sep = ":"))
  TTAAhits_R_gtf_exon2 <- left_join(TTAAhits_R_gtf_exon,exons_intron_df, by="exon_loci")
  #####intersect TTAA hit exons
  gene_TTAA <- unique(TTAAhits_R_gtf_exon$V15)
  genelist_1exon <- intersect(gene_TTAA, genelist_1exon_all)
  genelist_mutiexon <- intersect(gene_TTAA, genelist_mutiexon_all)
  ####label V4 as original TTAA_ID
  TTAAhits_R_gtf_exon2$V4 <- paste(TTAAhits_R_gtf_exon2$V1 , paste(TTAAhits_R_gtf_exon2$V2, TTAAhits_R_gtf_exon2$V3, sep= "-"), sep = ":")
  combined_df <- data.frame()
  for (i in genelist_mutiexon){
    df_single_gene <- exons_intron_df%>% dplyr::filter(gene_id==i)
    TTAAhits_R_gtf_exon2_single_gene <- TTAAhits_R_gtf_exon2%>%dplyr::filter(gene_id==i)
    transcripts_df_single_gene <- transcripts_df%>%dplyr::filter(gene_id==i)
    intron_single_gene <- exons_intron_df[exons_intron_df$gene_id==i,]
    ####get exon ord numbers within single gene hit by TTAA
    exon_ord_all <- unique(TTAAhits_R_gtf_exon2_single_gene$exon.ord)
    strand_single <- as.character(TTAAhits_R_gtf_exon2_single_gene$V13)
    if(strand_single[1] == "+"){
      #####should not use unique to remove NA, since some introns may have same length
      intron_width <- na.omit(intron_single_gene$intron_width)
      #intron_width <- intron_width[-1]
      all_intron_width <- sum(intron_width)
      new_df <- data.frame()
      for (j in exon_ord_all){
        TTAA_single_exon <-  TTAAhits_R_gtf_exon2_single_gene[TTAAhits_R_gtf_exon2_single_gene$exon.ord==j,]
        TTAA_single_exon$accu_introns <- rep(sum(intron_width[1:(j-1)]), nrow(TTAA_single_exon))
        ###can not use next, since next is for next j not next line of code,try() function can tolerate error and move to next line of code
        if(j==1){TTAA_single_exon$V2 <- TTAA_single_exon$V2}else{TTAA_single_exon$V2 <- TTAA_single_exon$V2-TTAA_single_exon$accu_introns}
        TTAA_single_exon <- TTAA_single_exon%>%mutate(V3=V2+4)
        new_df <- rbind(new_df,TTAA_single_exon)
      }
      new_df$V9 <- rep("modi_trans",nrow(new_df))
      new_df$V10 <- new_df$transcript_start
      new_df$V11 <- new_df$transcript_end-all_intron_width
    }else{
      #####should not use unique to remove NA, since some introns may have same length
      intron_width <- na.omit(intron_single_gene$intron_width)
      #intron_width <- intron_width[-1]
      ####reverse the order
      intron_width <- rev(intron_width)
      all_intron_width <- sum(intron_width)
      new_df <- data.frame()
      for (j in exon_ord_all){
        TTAA_single_exon <-  TTAAhits_R_gtf_exon2_single_gene[TTAAhits_R_gtf_exon2_single_gene$exon.ord==j,]
        TTAA_single_exon$accu_introns <- rep(sum(intron_width[1:(j-1)]), nrow(TTAA_single_exon))
        if(j==1){TTAA_single_exon$V3 <-TTAA_single_exon$V3}else{TTAA_single_exon$V3 <- TTAA_single_exon$V3+TTAA_single_exon$accu_introns}
        TTAA_single_exon <- TTAA_single_exon%>%mutate(V2=V3-4)
        new_df <- rbind(new_df,TTAA_single_exon)
      }
      new_df$V9 <- rep("modi_trans",nrow(new_df))
      new_df$V10 <- new_df$transcript_start+all_intron_width
      new_df$V11 <- new_df$transcript_end
    }
    
    combined_df <- rbind(combined_df, new_df)
  }
  TTAAhits_R_gtf_only1exon <- TTAAhits_R_gtf_exon2[TTAAhits_R_gtf_exon2$gene_id%in% genelist_1exon,]
  TTAAhits_R_gtf_only1exon$accu_introns <- 0
  TTAAhits_R_gtf_only1exon$V9 <- rep("modi_trans",nrow(TTAAhits_R_gtf_only1exon))
  dim(TTAAhits_R_gtf_only1exon)
  cat(paste('TTAAhits_R_gtf_only1exon:', nrow(TTAAhits_R_gtf_only1exon)))
  cat('\n')
  cat(paste('No. of multiexon_gene:', length(unique(TTAAhits_R_gtf_only1exon$gene_id))))
  cat('\n')
  cat(paste('TTAAhits_R_gtf_multiexon:', nrow(combined_df)))
  cat('\n')
  cat(paste('No. of multiexon_gene:', length(unique(combined_df$gene_id))))
  cat('\n')
  combined_df <- rbind(TTAAhits_R_gtf_only1exon, combined_df)
  combined_df <- combined_df%>%arrange(V1,V2,V3, gene_id, exon.ord)
  cat(paste('all:', nrow(combined_df)))
  cat('\n')
  cat(paste('No. of all_gene:', length(unique(combined_df$gene_id))))
  return(combined_df)
}

TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info <- removed_introns_TTAAhits_gtf_exons_extracted_modified(TTAAhits_R_gtf_include_contigs, exons_intron_df, transcripts_df)

write.xlsx(TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info,"./Output/TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info.xlsx")

test <- read.xlsx("./Output/TTAAhits_R_gtf_include_contigs_exons_extracted_modified_all_info.xlsx")


