library(readxl)  
library(openxlsx)
library(bedtoolsr)
library(rtracklayer)
library(dplyr)




#########################
#                       #
# Step 1: Read in files #
#                       #
#########################
# read in count matrix
cm_file <- "./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon_all.xlsx"
cm <- read_xlsx(cm_file)
class(cm)


# read in bed file with ttaa sites, as a dataframe
TTAA_file <- "./Output/TTAAhitsPf/PF3D7_theo_TTAA_R_modified_final.bed"
TTAA_sites <- read.table(TTAA_file)
class(TTAA_sites)

# read in gff file, using rtracklayer, which turns the gff file into a GRanges object
gff_file <- "./Input/Genome/PlasmoDB-63_Pfalciparum3D7.gff"
gff <- import(gff_file)
class(gff)



###########################
#                         #
# Step 2: Rename cols     #
#                         #
###########################


# TTAA_file col names
colnames(TTAA_sites)[1:6] <- c("chr", "start", "end","name","score","strand")

# create key
TTAA_sites$TTAA_key <- paste0(TTAA_sites$chr,":",TTAA_sites$start)
#collection locations
TTAA_loc <- unique(TTAA_sites[,c("chr","start","end","TTAA_key")])



# gff file

feature <- subset(gff, type %in% c("five_prime_UTR","three_prime_UTR"))

feature_df <- data.frame(
  chr = as.character(seqnames(feature)),
  start = start(feature)-1,
  end=end(feature),
  Feature = as.character(feature$type),
  geneID= sub("\\..*$","",unlist(mcols(feature)$Parent)),
  stringsAsFactors = FALSE
)

features_bed <- data.frame(
  V1 = feature_df$chr,
  V2 = feature_df$start,
  V3 = feature_df$end,
  V4 = feature_df$geneID, 
  V5 = feature_df$Feature
)



#################################
#                               #
# Step 3: intersect by gene     #
#                               #
#################################

# finds where genes intersect with ttaa sites by feature... 
hits <- bedtoolsr::bt.intersect(a = TTAA_loc, b=features_bed, wa = TRUE, wb=TRUE)

### label hits by col name... 
colnames(hits)[1:9] <- c("chr", "start","end", "TTAA_key", "bed_chr","bed_start", "bed_end","geneID", "Feature")



# annotate by gene not by site, 
annotate_by_gene <- hits %>%
  # group by key and gene id, so same key on different genes aren't connected
  group_by(TTAA_key, geneID) %>%
  summarise(
    Location_hit = case_when(
      # only need to update exon to UTR regions if applicable
      any(Feature =="five_prime_UTR") ~ "UTR5",
      any(Feature =="three_prime_UTR") ~"UTR3",
      TRUE ~ NA_character_
    ),
    .groups ="drop"
  ) %>%
  rename(GeneID=geneID)


new_cm <- cm %>%
  # join using TTAA_Key and geneID to keep overlapping genes separated
  dplyr::mutate(TTAA_key=paste0(Chrom,":",Site)) %>%
  dplyr::left_join(annotate_by_gene,by=c("TTAA_key","GeneID")) %>%
  dplyr::mutate(
    Location= ifelse(!is.na(Location_hit), Location_hit,Location)
  )%>%
  dplyr::select(names(cm))


write.xlsx(new_cm, "./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon_UTR_annotation_all_120225.xlsx")



