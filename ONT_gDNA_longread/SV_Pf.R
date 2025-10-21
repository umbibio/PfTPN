library(data.table)
library(stringr)
library(rtracklayer)
library(GenomicRanges)

vcf_df <- fread("./Input/long_read_ONT_gDNA/Pf3D7_longread.sniffles.vcf", skip = "#CHROM")
write.xlsx(vcf_df,"./Output/long_read_ONT_gDNA/original_106_SV_identified_without_filtering.xlsx")
nrow(vcf_df)
head(vcf_df)
vcf_df <- vcf_df[!grepl("API|MIT", vcf_df$`#CHROM`), ]
nrow(vcf_df)
vcf_filtered <- vcf_df[
  FILTER == "PASS" &
    grepl("PRECISE", INFO) &
    !grepl("IMPRECISE", INFO)
]
nrow(vcf_filtered)


vcf_df[, END := as.numeric(str_extract(INFO, "(?<=END=)[0-9]+"))]
vcf_df[, START := as.numeric(POS) - 1]
# Extract SVTYPE as well (optional)
vcf_df[, SVTYPE := str_extract(INFO, "(?<=SVTYPE=)[A-Z]+")]


vcf_sniffles2bed <- function(vcf_df) {
  # Load required package
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required. Please install it with install.packages('stringr')")
  }
  library(stringr)
  
  # Extract END from INFO (numeric)
  vcf_df[, END := as.numeric(str_extract(INFO, "(?<=END=)[0-9]+"))]
  
  # Compute 0-based START coordinate
  vcf_df[, START := as.numeric(POS) - 1]
  
  # Extract SVTYPE (e.g. DEL, INS, DUP, INV, BND)
  vcf_df[, SVTYPE := str_extract(INFO, "(?<=SVTYPE=)[A-Z]+")]
  
  bed_df <- vcf_df[, .(`#CHROM`, START, END, SVTYPE, FILTER)]
  bed_df$Len <- bed_df$END-bed_df$START
  # Return modified data.table
  return(bed_df)
}


bed_df_filtered <- vcf_sniffles2bed(vcf_filtered)
table(bed_df_filtered$SVTYPE)
bed_df_filtered_DEL <- bed_df_filtered%>%dplyr::filter(SVTYPE=="DEL")
bed_df_filtered_DEL2 <- bed_df_filtered_DEL%>%dplyr::filter(Len>100)
nrow(bed_df_filtered_DEL2)
fwrite(bed_df_filtered_DEL , file = "./Output/long_read_ONT_gDNA/sniffles_sv_deletion.bed", sep = "\t", col.names = FALSE)
fwrite(bed_df_filtered, file = "./Output/long_read_ONT_gDNA/sniffles_sv_PASS.bed", sep = "\t", col.names = FALSE)

bed_df_filtered_DEL_loci<- bed_df_filtered_DEL[, .(`#CHROM`, START, END)]
fwrite(bed_df_filtered_DEL_loci , file = "./Output/long_read_ONT_gDNA/sniffles_sv_deletion_loci_only.bed", sep = "\t", col.names = FALSE)

