## ---------------------------------------------------------------
## Generate theoretical TTAA insertion sites (fixed chr names)
## ---------------------------------------------------------------

library(tidyverse)
library(IRanges)
library(Biostrings)
library(rtracklayer)

## Paths
fasta_path <- "./Input/Genome/PlasmoDB-63_Pfalciparum3D7_Genome.fasta"
raw_bed    <- "./Output/PF3D7_theo_TTAA_R.bed"
final_bed  <- "./Output/PF3D7_theo_TTAA_R_modified_final.bed"

## 1. Load genome FASTA and clean chromosome names
Genome.fasta_all <- readDNAStringSet(fasta_path)

cat("FASTA names before cleaning:\n")
print(head(names(Genome.fasta_all)))

# Strip everything after first space or "|" so names become:
#   Pf3D7_01_v3, Pf3D7_02_v3, ...
names(Genome.fasta_all) <- sub("([[:space:]]|\\|).*", "", names(Genome.fasta_all))

cat("\nFASTA names after cleaning:\n")
print(head(names(Genome.fasta_all)))

## 2. Find all TTAA motifs
TTAA <- "TTAA"
Genome_theo_TTAA <- vmatchPattern(TTAA, Genome.fasta_all, fixed = TRUE)

## 3. Convert to a single GRanges of unique genomic TTAA sites
ttaa_gr <- as(Genome_theo_TTAA, "GRanges")

# Metadata for BED
strand(ttaa_gr)        <- "*"   # neutral strand for the raw BED
mcols(ttaa_gr)$name    <- TTAA
mcols(ttaa_gr)$score   <- 0L

## 4. Export raw single-strand BED (one row per genomic site)
export.bed(ttaa_gr, con = raw_bed)

## 5. Duplicate for + and − strand (Sida-style final BED)
ttaa_plus  <- ttaa_gr; strand(ttaa_plus)  <- "+"
ttaa_minus <- ttaa_gr; strand(ttaa_minus) <- "-"

ttaa_both <- c(ttaa_plus, ttaa_minus)

# This BED has: chrom, start, end, name, score, strand
export.bed(ttaa_both, con = final_bed)

cat("\nTTAA unique genomic sites:", length(ttaa_gr), "\n")
cat("Rows in final (+/-) BED:", length(ttaa_both), "\n")
cat("Wrote:\n  ", raw_bed, "\n  ", final_bed, "\n")
