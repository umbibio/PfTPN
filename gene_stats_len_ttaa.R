library(readxl)
library(dplyr)
library(writexl)
library(rtracklayer)
library(GenomicRanges)
# library(Biostrings)  # only needed if you later use the FASTA

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
genome_dir   <- "./Input/Genome"
input_dir    <- "./Input"
bg_dir       <- file.path(input_dir, "background_genes_kritika")
output_dir   <- "./Output"

gtf_file     <- file.path(genome_dir, "PlasmoDB-63_Pfalciparum3D7.gtf")
# fasta_file   <- file.path(genome_dir, "PlasmoDB-63_Pfalciparum3D7_Genome.fasta")

cm_file      <- file.path(output_dir, "PfTPN_count_matrix_renamed.xlsx")
meta_file    <- file.path(output_dir, "PfTPN_sample_metadata.xlsx")

kritika_file <- file.path(bg_dir, "Pf high confidence essential,dispensable gene list.xlsx")

# ============================================================
# 1) Kritika background genes
# ============================================================
dispensable <- read_xlsx(kritika_file, sheet = 1) %>%
  mutate(class = "dispensable")

essential   <- read_xlsx(kritika_file, sheet = 2) %>%
  mutate(class = "essential")

background_genes <- bind_rows(dispensable, essential) %>%
  transmute(
    GeneID      = Gene,
    class,
    Description
  ) %>%
  mutate(
    # 0 = essential, 1 = dispensable
    class_lab = if_else(class == "dispensable", 1L, 0L)
  )

write_xlsx(
  background_genes,
  file.path(genome_dir, "kritika_background_genes_kz.xlsx")
)

# ============================================================
# 2) CDS ranges and gene lengths
# ============================================================
gtf      <- import(gtf_file)
gtf_cds  <- gtf[gtf$type == "CDS"]   # note: some genes only have exon annotation (e.g. PF3D7_0220400)
gtf_cds_df <- as.data.frame(gtf_cds)

# CDS per gene, reduced to non-overlapping ranges
gene_cds_gr <- split(gtf_cds, gtf_cds_df$gene_id) %>%
  lapply(reduce) %>%
  GRangesList()

gene_lengths <- tibble(
  GeneID     = names(gene_cds_gr),
  cds_length = sapply(gene_cds_gr, function(gr) sum(width(gr)))
)

# ============================================================
# 3) TTAA BED and CDS-relative positions
# ============================================================
ttaa_bed_file <- file.path(output_dir, "PF3D7_theo_TTAA_R_modified_final.bed")
ttaa_gr       <- import(ttaa_bed_file, format = "BED")

# Sanity check
print(intersect(seqlevels(gtf_cds), seqlevels(ttaa_gr)))

## TTAA sites in CDS, with position relative to CDS start

# 3.1 Unlist CDS and tag GeneID
gene_cds_unlist <- unlist(gene_cds_gr)
gene_ids_unlist <- rep(names(gene_cds_gr), elementNROWS(gene_cds_gr))
mcols(gene_cds_unlist)$GeneID <- gene_ids_unlist

# 3.2 CDS offsets per segment (distance from CDS 5' to segment start)
cds_offset_list <- lapply(gene_cds_gr, function(gr) {
  strand_g <- as.character(unique(strand(gr)))
  if (length(strand_g) != 1L) stop("Gene on multiple strands?")
  
  if (strand_g == "+") {
    ord <- order(start(gr))
  } else if (strand_g == "-") {
    # for - strand, CDS "start" is at highest genomic coordinate
    ord <- order(end(gr), decreasing = TRUE)
  } else {
    stop("Unknown strand: ", strand_g)
  }
  
  w         <- width(gr)[ord]
  offs_ord  <- c(0L, head(cumsum(w), -1L))     # cumulative length before each segment
  offs_full <- integer(length(gr))
  offs_full[ord] <- offs_ord
  offs_full
})

cds_offset_unlist <- unlist(cds_offset_list, use.names = FALSE)
mcols(gene_cds_unlist)$cds_offset <- cds_offset_unlist

# 3.3 Overlaps TTAA ↔ CDS segments (strand-aware)
hits <- findOverlaps(ttaa_gr, gene_cds_unlist, ignore.strand = FALSE)

q <- queryHits(hits)      # TTAA index
s <- subjectHits(hits)    # CDS segment index

ttaa_chr   <- as.character(seqnames(ttaa_gr))[q]
ttaa_start <- start(ttaa_gr)[q]
ttaa_end   <- end(ttaa_gr)[q]

seg_start  <- start(gene_cds_unlist)[s]
seg_end    <- end(gene_cds_unlist)[s]
seg_strand <- as.character(strand(gene_cds_unlist))[s]
seg_gene   <- mcols(gene_cds_unlist)$GeneID[s]
seg_offset <- mcols(gene_cds_unlist)$cds_offset[s]

# 3.4 Position of TTAA relative to CDS start (1-based)
ttaa_cds_pos <- ifelse(
  seg_strand == "+",
  seg_offset + (ttaa_start - seg_start + 1L),
  seg_offset + (seg_end   - ttaa_end   + 1L)
)

# 3.5 Build TTAA table with fractional CDS position
gene_len_vec <- gene_lengths$cds_length
names(gene_len_vec) <- gene_lengths$GeneID
this_len <- gene_len_vec[seg_gene]

ttaa_cds_tbl <- tibble(
  GeneID        = seg_gene,
  chrom         = ttaa_chr,
  ttaa_start    = ttaa_start,
  ttaa_end      = ttaa_end,
  strand        = seg_strand,
  cds_length    = this_len,
  ttaa_cds_pos  = ttaa_cds_pos,
  ttaa_cds_frac = ttaa_cds_pos / this_len   # 0–1 position along CDS
)

# Per-gene TTAA counts
gene_TTAA <- ttaa_cds_tbl %>%
  count(GeneID, name = "n_TTAA_cds")

# Final gene_stats: one row per TTAA (per gene)
gene_stats <- ttaa_cds_tbl %>%
  left_join(gene_TTAA, by = "GeneID") %>%
  mutate(
    Site_ID = paste0(chrom, ":", ttaa_start, "-", ttaa_end)
  )

# quick checks
head(ttaa_cds_tbl)
head(gene_stats)
summary(gene_stats$n_TTAA_cds)

write_xlsx(
  gene_stats,
  file.path(genome_dir, "Pf3D7_gene_ttaa_rel_loc_len_summary_stats.xlsx")
)

# ============================================================
# 4) Merge with TN-seq count matrix
# ============================================================
cm          <- read_xlsx(cm_file)
sample_meta <- read_xlsx(meta_file)  # currently unused; handy if you later filter by treatment/day

rm_cols <- c(
  "R1_pointing_downsteam", "R1_pointing_upstream",
  "Location1", "Location2",
  "Run_repeats", "Theoretical"
)

# Drop unused meta columns but keep all sample count columns + core meta
cm_single <- cm %>%
  dplyr::select(-all_of(rm_cols))

# meta used to define a unique genomic insertion site (TTAA)
meta_keep       <- c("Chrom", "Site")
meta_for_counts <- c(meta_keep, "GeneID", "gene.description", "Location", "Present_in_any_samples")
count_cols      <- setdiff(names(cm_single), meta_for_counts)

# Collapse to one row per Chrom:Site, summing across strand / repeats
cm_single_sum <- cm_single %>%
  dplyr::group_by(dplyr::across(all_of(meta_keep))) %>%
  dplyr::summarise(
    dplyr::across(all_of(count_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # TTAA motif is 4 bp; original Site is start-1, so shift by +1
    ttaa_start = Site + 1L,
    ttaa_end   = Site + 4L,
    Site_ID    = paste0(Chrom, ":", ttaa_start, "-", ttaa_end),
    Present_in_any_samples = dplyr::if_else(
      rowSums(dplyr::across(all_of(count_cols))) > 0,
      "yes", "no"
    )
  ) %>%
  dplyr::relocate(Present_in_any_samples, .after = Site)

cat("cm_single_sum dim:", paste(dim(cm_single_sum), collapse = " x "), "\n")
cat("Unique Site_ID in cm_single_sum:", length(unique(cm_single_sum$Site_ID)), "\n")

cat("gene_stats dim:", paste(dim(gene_stats), collapse = " x "), "\n")
cat("Unique Site_ID in gene_stats:", length(unique(gene_stats$Site_ID)), "\n")

# Overlap sanity checks
overlap_cm_in_gs <- sum(cm_single_sum$Site_ID %in% gene_stats$Site_ID)
overlap_gs_in_cm <- sum(gene_stats$Site_ID %in% cm_single_sum$Site_ID)

cat("Sites in cm_single_sum that are in gene_stats:", overlap_cm_in_gs, "\n")
cat("Sites in gene_stats that are in cm_single_sum:", overlap_gs_in_cm, "\n")

# Keep only sites that fall in CDS (present in gene_stats), one row per Site_ID
cm_cds_counts <- cm_single_sum %>%
  dplyr::filter(Site_ID %in% gene_stats$Site_ID) %>%
  dplyr::group_by(Site_ID) %>%
  dplyr::summarise(
    dplyr::across(all_of(count_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

cat("cm_cds_counts dim:", paste(dim(cm_cds_counts), collapse = " x "), "\n")

# Merge counts into gene_stats (TTAA per gene; counts replicated if TTAA belongs to >1 gene)
gene_stats_counts <- gene_stats %>%
  dplyr::left_join(cm_cds_counts, by = "Site_ID")

cat("gene_stats_counts dim:", paste(dim(gene_stats_counts), collapse = " x "), "\n")
cat("Any NA GeneID after merge? ", sum(is.na(gene_stats_counts$GeneID)), "\n")

# Reorder columns: meta first, then per-sample counts
meta_cols_gs <- c(
  "GeneID", "chrom", "ttaa_start", "ttaa_end", "strand",
  "cds_length", "ttaa_cds_pos", "ttaa_cds_frac", "n_TTAA_cds", "Site_ID"
)
meta_cols_gs  <- intersect(meta_cols_gs, names(gene_stats_counts))
count_cols_gs <- setdiff(names(gene_stats_counts), meta_cols_gs)

gene_stats_counts <- gene_stats_counts %>%
  dplyr::select(dplyr::all_of(meta_cols_gs), dplyr::all_of(count_cols_gs))

write_xlsx(
  gene_stats_counts,
  file.path(genome_dir, "Pf3D7_gene_ttaa_rel_loc_len_counts_per_sample.xlsx")
)
