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
#bg_dir       <- file.path(input_dir, "background_genes_kritika")
bg_dir       <- file.path(input_dir)
output_dir   <- "./Output"

gtf_file     <- file.path(genome_dir, "PlasmoDB-63_Pfalciparum3D7.gtf")
# fasta_file   <- file.path(genome_dir, "PlasmoDB-63_Pfalciparum3D7_Genome.fasta")

cm_file      <- file.path(output_dir, "PfTPN_count_matrix_renamed_UTR_overlaps_corrected_120225.xlsx")
meta_file    <- file.path(output_dir, "PfTPN_sample_metadata_UTR_overlaps_corrected_120225.xlsx")

#kritika_file <- file.path(bg_dir, "Pf high confidence essential,dispensable gene list.xlsx")
bg_genes_files <- file.path(bg_dir, "pf_pk_all_verified_genes.xlsx")
# ============================================================
# 1) Kritika background genes
# ============================================================
#dispensable <- read_xlsx(kritika_file, sheet = 1) %>%
#  mutate(class = "dispensable")


dispensable <- read_xlsx(bg_genes_files, sheet = "pf_dispensable") %>%
    mutate(class = "dispensable")
  
essential   <- read_xlsx(bg_genes_files, sheet = "pf_essential") %>%
  mutate(class = "essential")

background_genes <- bind_rows(dispensable, essential) %>%
  mutate(
    # 0 = essential, 1 = dispensable
    class_lab = if_else(class == "dispensable", 1L, 0L)
  )

write_xlsx(
  background_genes,
  file.path(genome_dir, "bg_genes_single_list_pf.xlsx")
)

# ============================================================
# 2) CDS ranges and gene lengths
# ============================================================
gtf      <- import(gtf_file)
gtf_cds  <- gtf[gtf$type == "CDS"]   # note: some genes only have exon annotation (e.g. PF3D7_0220400). -->  I dont understand this comment
gtf_cds_df <- as.data.frame(gtf_cds)

# add in for exons... will filter later for genes that have cds regions in them
gtf_exon  <- gtf[gtf$type == "exon"]  
gtf_exon_df <- as.data.frame(gtf_exon)


## look at gene 
gene <- "PF3D7_0220400"
gtf_gene <- gtf[mcols(gtf)$gene_id ==gene]
gene_PF3D7_0220400 <- as.data.frame(gtf_gene)

cm_gene <- dplyr::filter(cm,GeneID==gene)

genes_check <- c(
  "PF3D7_0220400", # missing. --> none of the genes have UTRs 
  "PF3D7_0424300", 
  "PF3D7_1039100",
  "PF3D7_1149500",
  "PF3D7_1417400",
  "PF3D7_0206800",  # overlapping genes, large intersection of exons
  "PF3D7_0206700"
)

gtf_list <- setNames(lapply(genes_check, \(g) as.data.frame(gtf[mcols(gtf)$gene_id==g])),genes_check)
cm_list <- lapply(genes_check, \(g) dplyr::filter(cm,GeneID==g))

### before 
cds_prefilter <- split(gtf_cds, gtf_cds_df$gene_id)
old_len <- setdiff(genes_check, names(cds_prefilter))
cat("genes not contained in gtf_cds but in verified genes:\n", old_len)


# CDS per gene, reduced to non-overlapping ranges.  --> I don't understand this comment
gene_cds_gr <- split(gtf_cds, gtf_cds_df$gene_id) %>%
  lapply(reduce) %>%
  GRangesList()


### check diff
all_genes <- unique(gtf$gene_id)
genes_with_exon <-  unique(mcols(gtf[gtf$type=="exon"])$gene_id)
cds_genes <- unique(mcols(gtf[gtf$type=="CDS"])$gene_id)

cat("total genes in gtf: ", length(all_genes), "\n")
cat("genes with exon row: ", length(genes_with_exon), "\n")
cat("genes with cds row:", length(cds_genes),"\n")  # 5312 genes in file... 
cat("genes no cds row: ", length(setdiff(genes_with_exon, cds_genes)), "\n")


### check number of exon rows 
cat("CDS row count:", sum(gtf$type=="CDS"), "\n",    ## returns true for past two genes, but cannot see those in gtf?
    "exon row count:", sum(gtf$type=="exon"), "\n",
    "transcript row count:", sum(gtf$type=="transcript"), "\n"
)


## see if list of genes is in genes, has exon, has cds
gene_check_status <- data.frame(
  GeneID=genes_check, 
  in_genes = genes_check %in% all_genes,
  has_exon = genes_check %in% genes_with_exon,
  has_cds = genes_check %in% cds_genes
  
)

gene_check_status


### compared to IGV, cds-> exon only, exon/transcript -> include UTR region as well
gene_with_cds <- "PF3D7_0206800"
with_cds_table <- as.data.frame(gtf[mcols(gtf)$gene_id ==gene_with_cds])
View(with_cds_table)



### compared to IGV -> this gene has no UTRs so everything is labeled as an exon? 
gene_without_cds <- "PF3D7_1149500"
without_cds_table <-  as.data.frame(gtf[mcols(gtf)$gene_id == gene_without_cds])
View(without_cds_table)







##########################
#
# if gene has CDS -> implied UTR, want to use CDS region
# else no CDS -> implied no UTR, use exon ? 
#
# So, what I am doing here, is picking up the genes that have no cds regions and then concat the two lists 
#
###########################





# assume genes in this list do not have UTRs....
genes_exons_only <- setdiff(unique(gtf_exon_df$gene_id), unique(gtf_cds_df$gene_id))

gene_exon_only_gr <- split(
  gtf_exon[mcols(gtf_exon)$gene_id %in% genes_exons_only],
  gtf_exon_df$gene_id[gtf_exon_df$gene_id %in% genes_exons_only]) %>%
  lapply(reduce) %>%
  GRangesList()




################
#
# next, combine GRanges Lists, make sure that the number of genes in each makes sense 
#
##################


gene_final_exons_gr <- c(gene_cds_gr, gene_exon_only_gr)

cat("CDS genes (has UTR): ", length(gene_cds_gr), "\n")
cat("exon only genes (no UTR): ", length(gene_exon_only_gr), "\n")
cat("combined total: ",length(gene_final_exons_gr), "\n")
cat("confirm no genes duplicated: ", any(duplicated(names(gene_final_exons_gr))))





gene_lengths <- tibble(
  GeneID     = names(gene_final_exons_gr),
  cds_length = sapply(gene_final_exons_gr, function(gr) sum(width(gr)))
)










########################
#
# replaced all instances of gene_cds_gr with gene_final_exons_gr
#
# also add TTAAhitsPf/ to path directory of ttaa_bed_file... 
#
# otherwise, everything is the same 
#
########################



# ============================================================
# 3) TTAA BED and CDS-relative positions
# ============================================================

### changed path directory here to include TTAAhitsPf/
ttaa_bed_file <- file.path(output_dir, "TTAAhitsPf/PF3D7_theo_TTAA_R_modified_final.bed")
ttaa_gr       <- import(ttaa_bed_file, format = "BED")

# Sanity check
print(intersect(seqlevels(gene_final_exons_gr), seqlevels(ttaa_gr)))

## TTAA sites in CDS, with position relative to CDS start

# 3.1 Unlist CDS and tag GeneID
gene_cds_unlist <- unlist(gene_final_exons_gr)
gene_ids_unlist <- rep(names(gene_final_exons_gr), elementNROWS(gene_final_exons_gr))
mcols(gene_cds_unlist)$GeneID <- gene_ids_unlist

# 3.2 CDS offsets per segment (distance from CDS 5' to segment start)
cds_offset_list <- lapply(gene_final_exons_gr, function(gr) {
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
  file.path(genome_dir, "Pf3D7_gene_ttaa_rel_loc_len_summary_stats_1225_mcs.xlsx")
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
  file.path(genome_dir, "Pf3D7_gene_ttaa_rel_loc_len_counts_per_sample_1225_mcs.xlsx")
)
