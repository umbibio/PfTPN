library(dplyr)
library(readxl)
library(stringr)
library(ggplot2)
library(tidyr)
library(ggrepel)
## Read in meta data and count matrix
input_dir    <- "./Input"
bg_dir       <- file.path(input_dir, "background_genes_kritika")
output_dir   <- "./Output"
genome_dir   <- "./Input/Genome"

cm_file      <- file.path(output_dir, "PfTPN_count_matrix_renamed.xlsx")
meta_file    <- file.path(output_dir, "PfTPN_sample_metadata.xlsx")
background_genes_file <- file.path(genome_dir, "kritika_background_genes_kz.xlsx")

kritika_file <- file.path(bg_dir, "Pf high confidence essential,dispensable gene list.xlsx")

cm <- read_xlsx(cm_file)
sample_meta <- read_xlsx(meta_file)

# Columns that are NOT samples
meta_cols <- c(
  "Chrom","Site","R1_pointing_downsteam","R1_pointing_upstream",
  "GeneID","gene.description","Location","Location1","Location2",
  "Present_in_any_samples","Run_repeats","Theoretical"
)

sample_cols <- setdiff(colnames(cm), meta_cols)

# ctrl_samples <- sample_meta %>%
#   filter(treatment == "Control", !is.na(day), day >= 14) %>%
#   pull(new_name)


ctrl_samples <- sample_meta %>%
  filter( !is.na(day), day >= 14) %>%
  pull(new_name)

length(ctrl_samples)
print(ctrl_samples)

## look into exons only for control samples
cm_ctrl <- cm %>% dplyr::filter(Location == "exon") %>% 
  dplyr::select(GeneID, all_of(ctrl_samples))


# aggregate counts across all sites (total per gene per sample)
cm_gene_lev_ctrl <- cm_ctrl %>%
  group_by(GeneID) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

head(cm_gene_lev_ctrl)

# Join gene stats (CDS length + TTAA counts)
gene_stats <- read_xlsx("./Input/Genome/Pf3D7_gene_ttaa_len_summary_stats.xlsx")

# calculate features. Inner join becasue soem genes do not have CDS annotation
# and some genes have no insertion profile

gene_table <- inner_join(gene_stats, cm_gene_lev_ctrl, by = "GeneID") 


# background essential/dispensable labels. Some background genes have no insertion numbers
background_genes <- read_xlsx(background_genes_file)

background_genes_table <- inner_join(background_genes, gene_table , by = 'GeneID') 

dis.ind <- background_genes_table$class_lab == 1
ess.ind <- background_genes_table$class_lab == 0

quantile(background_genes_table$cds_length, probs = 0.1)

# mark outliers
background_genes_table <- background_genes_table %>%
  mutate(total_insertions_gene = rowSums(across(contains("TC")), na.rm = TRUE))
Qs <- background_genes_table %>%
  group_by(class) %>%
  summarise(
    q1 = quantile(total_insertions_gene, 0.15, na.rm = TRUE),
    q2 = quantile(total_insertions_gene, 0.85, na.rm = TRUE),
    .groups = "drop"
  )

background_genes_table <- background_genes_table %>% group_by(class) %>% 
  mutate(outlier = ifelse(class == 'dispensable' & total_insertions_gene < Qs$q1[Qs$class == 'dispensable'], 1,
                          ifelse(class == 'essential' & total_insertions_gene > Qs$q2[Qs$class == 'essential'], 1, 0)))

background_genes_table.h <- background_genes_table %>% dplyr::filter(outlier == 0)
ggplot(background_genes_table.h, aes(x = total_insertions_gene, fill = class)) + 
  geom_histogram() + geom_density() + facet_grid(class~.) + theme_bw()


# High-confidence background genes only
train_df <- background_genes_table.h %>%
  filter(class %in% c("essential", "dispensable")) %>%
  mutate(
    y = ifelse(class == "essential", 0L, 1L)
  )

# Build predictor matrix X
X <- train_df %>% ungroup() %>% 
  dplyr::select(contains('TC')) %>%
  as.matrix()

y <- train_df$y

# Remove rows with NAs
keep <- complete.cases(X, y)
X <- X[keep, , drop = FALSE]
y <- y[keep]

train_df <- train_df[keep, ]


## PCA on samples: transpose X
# PCA
pc <- prcomp(t(X), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pc$x[,1],
  PC2 = pc$x[,2],
  sample = colnames(X)
)

# Clean labels
pca_df$label_clean <- gsub('pTN9_|PTPN_', '', gsub('TC_', '', gsub("(_Control.*)$", "", pca_df$sample)))

# Plot
library(ggplot2)
library(ggrepel)

ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = label_clean), size = 3) +
  theme_bw(base_size = 14) +
  ggtitle("PCA on samples (Control ≥ day14)")


ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = sample),
    size = 3,
    max.overlaps = Inf   # ← allow unlimited overlaps
  ) +
  theme_bw(base_size = 14) +
  ggtitle("PCA on samples (Control ≥ day14)")




## Re do glmnet
library(dplyr)
library(readxl)
library(glmnet)
library(ggplot2)
library(ggrepel)


bad_ctrl_samples <- c(
  "TC_pTN9_day20_Control",
  "TC_PTPN_day14_Control_b0505_rep2",
  "TC_PTPN_day14_Control_b0505_rep2",
  "TC_PTPN_day17_Mev_b0505",
  "TC_PTPN_day20_Control_b0505_rep2"
)

ctrl_samples_clean <- setdiff(ctrl_samples, bad_ctrl_samples)
ctrl_samples_clean


# exon-only control counts, strand-specific insertions
cm_ctrl <- cm %>%
  filter(Location == "exon") %>%
  dplyr::select(GeneID, all_of(ctrl_samples_clean))

# per-gene insertion summary
gene_insertions <- cm_ctrl %>%
  group_by(GeneID) %>%
  summarise(
    n_sites = n(),  # number of (strand-specific) TTAA rows in matrix
    total_insertions = sum(across(all_of(ctrl_samples_clean)), na.rm = TRUE),
    n_sites_hit = sum(rowSums(across(all_of(ctrl_samples_clean)), na.rm = TRUE) > 0),
    prop_sites_hit = n_sites_hit / n_sites,
    .groups = "drop"
  )

# CDS length + TTAA per gene
gene_stats <- read_xlsx("./Input/Genome/Pf3D7_gene_ttaa_len_summary_stats.xlsx")

gene_table <- inner_join(gene_stats, gene_insertions, by = "GeneID") %>%
  mutate(
    mean_insertions_per_TTAA = total_insertions / n_TTAA_cds,
    Msg_score = log10((total_insertions + 1) / n_TTAA_cds)
  )

background_genes <- read_xlsx(background_genes_file)

background_genes_table <- inner_join(background_genes, gene_table, by = "GeneID")

# mark outlier genes by total insertions within each class
Qs <- background_genes_table %>%
  group_by(class) %>%
  summarise(
    q1 = quantile(total_insertions, probs = 0.15),
    q2 = quantile(total_insertions, probs = 0.85),
    .groups = "drop"
  )

background_genes_table <- background_genes_table %>%
  group_by(class) %>%
  mutate(
    outlier = case_when(
      class == "dispensable" &
        total_insertions < Qs$q1[Qs$class == "dispensable"] ~ 1L,
      class == "essential" &
        total_insertions > Qs$q2[Qs$class == "essential"] ~ 1L,
      TRUE ~ 0L
    ),
    class_lab = ifelse(class == "essential", 0L, 1L)
  ) %>%
  ungroup()

background_genes_table.h <- background_genes_table %>% filter(outlier == 0)


train_df <- background_genes_table.h %>%
  filter(class %in% c("essential", "dispensable")) %>%
  mutate(y = ifelse(class == "essential", 0L, 1L))

X <- train_df %>%
  dplyr::select(Msg_score, cds_length, n_TTAA_cds,
         prop_sites_hit, total_insertions, mean_insertions_per_TTAA) %>%
  as.matrix()

y <- train_df$y

keep <- complete.cases(X, y)
X <- X[keep, , drop = FALSE]
y <- y[keep]
train_df <- train_df[keep, ]

cvfit <- cv.glmnet(
  X, y,
  family = "binomial",
  alpha  = 1,
  nfolds = 10,
  type.measure = "deviance"
)

lambda_min  <- cvfit$lambda.min
lambda_min
coef(cvfit, s = "lambda.min")

# posterior P(dispensable) on training bg genes
train_df$prob_dis <- as.numeric(
  predict(cvfit, newx = X, s = "lambda.min", type = "response")
)


X_all <- gene_table %>%
  dplyr::select(Msg_score, cds_length, n_TTAA_cds,
         prop_sites_hit, total_insertions, mean_insertions_per_TTAA) %>%
  as.matrix()

keep_all <- complete.cases(X_all)
prob_all <- rep(NA_real_, nrow(X_all))
prob_all[keep_all] <- predict(
  cvfit,
  newx = X_all[keep_all, , drop = FALSE],
  s    = "lambda.min",
  type = "response"
)

gene_table$prob_dis_glmnet <- prob_all


df_plot <- gene_table %>%
  arrange(prob_dis_glmnet) %>%
  mutate(rank = row_number()) %>%
  left_join(background_genes_table.h %>% dplyr::select(GeneID, class),
            by = "GeneID") %>%
  mutate(
    color = case_when(
      class == "essential"   ~ "red",
      class == "dispensable" ~ "green",
      TRUE                   ~ "grey80"
    ),
    is_bg = !is.na(class)
  )

bg_eval <- df_plot %>%
  filter(is_bg) %>%
  mutate(
    true = ifelse(class == "essential", 0L, 1L),
    pred = ifelse(prob_dis_glmnet >= 0.5, 1L, 0L)
  )

conf_mat <- table(Predicted = bg_eval$pred, True = bg_eval$true)
conf_mat

## quick diagnostic plot
plot(
  df_plot$rank,
  df_plot$prob_dis_glmnet,
  pch = 20,
  cex = 0.6,
  col = df_plot$color,
  xlab = "Ranked Genes",
  ylab = "Posterior P(dispensable)",
  main = "glmnet scores (outlier samples removed)"
)

points(
  df_plot$rank[df_plot$is_bg],
  df_plot$prob_dis_glmnet[df_plot$is_bg],
  pch = 21,
  bg  = df_plot$color[df_plot$is_bg],
  col = "black",
  lwd = 1.2,
  cex = 1.6
)

legend(
  "topleft",
  legend = c("Essential (BG)", "Dispensable (BG)", "Other"),
  col    = c("red", "green", "grey80"),
  pt.bg  = c("red", "green", "grey80"),
  pch    = 21,
  pt.cex = c(1.6, 1.6, 0.8),
  bty    = "n"
)

# JA false calls overlay (optional, same glmnet fit)

ko_table <- read_xlsx("./Input/JA_falsecalls/Gene KO table.xlsx", sheet = 4) %>%
  transmute(
    GeneID   = Gene_ID,
    class_JA = KO,
    geneName = geneName
  )

JA_df_plot <- gene_table %>%
  arrange(prob_dis_glmnet) %>%
  mutate(rank = row_number()) %>%
  left_join(ko_table, by = "GeneID") %>%
  mutate(
    class_simple = case_when(
      class_JA == "Essential" ~ "Essential",
      class_JA == "Disp"      ~ "Dispensable",
      TRUE                    ~ "Other"
    ),
    is_JA = !is.na(class_JA)
  )

class_colors <- c(
  "Essential"   = "red",
  "Dispensable" = "green4",
  "Other"       = "grey70"
)

ggplot(JA_df_plot, aes(rank, prob_dis_glmnet)) +
  geom_point(aes(color = class_simple), size = 1.3, alpha = 0.8) +
  scale_color_manual(values = class_colors) +
  geom_point(
    data = JA_df_plot %>% filter(is_JA),
    aes(rank, prob_dis_glmnet, fill = class_simple),
    shape = 21, size = 3, color = "black", stroke = 1
  ) +
  scale_fill_manual(values = class_colors, guide = "none") +
  geom_text_repel(
    data = JA_df_plot %>% filter(is_JA),
    aes(label = geneName),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.25,
    max.overlaps = Inf
  ) +
  labs(
    x = "Ranked Genes",
    y = "Posterior P(dispensable)",
    title = "glmnet Posterior Scores (clean controls)",
    color = "JA class"
  ) +
  theme_bw(base_size = 13)
