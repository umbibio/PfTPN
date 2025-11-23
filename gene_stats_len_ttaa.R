library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(writexl)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)

# Read in files
genome_dir   <- "./Input/Genome"
input_dir    <- "./Input"
bg_dir       <- file.path(input_dir, "background_genes_kritika")
output_dir   <- "./Output"

gtf_file     <- file.path(genome_dir, "PlasmoDB-63_Pfalciparum3D7.gtf")
fasta_file   <- file.path(genome_dir, "PlasmoDB-63_Pfalciparum3D7_Genome.fasta")

cm_file      <- file.path(output_dir, "PfTPN_count_matrix_renamed.xlsx")
meta_file    <- file.path(output_dir, "PfTPN_sample_metadata.xlsx")

kritika_file <- file.path(bg_dir, "Pf high confidence essential,dispensable gene list.xlsx")

# Kritika background lists (essential vs dispensable)
dispensable <- read_xlsx(kritika_file, sheet = 1) %>%
  mutate(class = "dispensable")

essential   <- read_xlsx(kritika_file, sheet = 2) %>%
  mutate(class = "essential")

background_genes <- bind_rows(dispensable, essential) %>%
  transmute(GeneID = Gene, class, Description)

# 0 for essential and 1 for dispensable
background_genes$class_lab <- ifelse(background_genes$class == 'dispensable', 1, 0)

# calculate CDS length + TTAA counts per gene
# Import GTF and keep CDS features
gtf <- import(gtf_file)
gtf_cds <- gtf[gtf$type == "CDS"] ## some genes only have exon annotation not CDS: PF3D7_0220400
gtf_cds_df <- as.data.frame(gtf_cds) 

# split into list of cds per gene and the reduce to non-overlapping genomic ranges
# using the reduce function from Granges and convert back to Granges object 
gene_cds_gr <- split(gtf_cds, gtf_cds_df$gene_id) %>%
  lapply(reduce) %>%
  GRangesList()

gene_lengths <- tibble(
  GeneID    = names(gene_cds_gr), ## gene names
  ## width calculates width of each entry (exon), sum all
  cds_length = sapply(gene_cds_gr, function(gr) sum(width(gr)))
)

# Read TTAA bed file and fix chromosome names
ttaa_bed_file <- file.path(output_dir, "PF3D7_theo_TTAA_R_modified_final.bed")

ttaa_gr <- import(ttaa_bed_file, format = "BED")


# Clean seqlevels: drop everything after first space
#old_levels <- seqlevels(ttaa_gr)
#new_levels <- sub(" .*", "", old_levels)  # "Pf3D7_01_v3 | organism=..." -> "Pf3D7_01_v3"
#names(new_levels) <- old_levels

#ttaa_gr <- renameSeqlevels(ttaa_gr, new_levels)

# Sanity check: now they should match
print(intersect(seqlevels(gtf_cds), seqlevels(ttaa_gr)))

# Count TTAA sites in CDS per gene
# for each gene overlap the cds with the bed of TTAA locations  
ttaa_cds_counts <- sapply(gene_cds_gr, function(gr) sum(countOverlaps(gr, ttaa_gr)))

gene_TTAA <- tibble(
  GeneID     = names(ttaa_cds_counts),
  n_TTAA_cds = as.integer(ttaa_cds_counts)
)

# Combined genome stats table
gene_stats <- left_join(gene_lengths, gene_TTAA, by = "GeneID")

head(gene_stats)
summary(gene_stats$n_TTAA_cds)

write_xlsx(
  gene_stats,
  file.path(genome_dir, "Pf3D7_gene_ttaa_len_summary_stats.xlsx")
)

write_xlsx(
  background_genes,
  file.path(genome_dir, "kritika_background_genes_kz.xlsx")
)

