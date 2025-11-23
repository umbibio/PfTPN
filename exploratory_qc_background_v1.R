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

ctrl_samples <- sample_meta %>%
  filter(treatment == "Control", !is.na(day), day >= 14) %>%
  pull(new_name)

length(ctrl_samples)
print(ctrl_samples)

## look into exons only for control samples
cm_ctrl <- cm %>% dplyr::filter(Location == "exon") %>% 
  dplyr::select(GeneID, all_of(ctrl_samples))


# Total insertions per gene
gene_insertions <- cm_ctrl %>%
  group_by(GeneID) %>%
  summarise(
    n_sites = n(),  # number of TTAA sites for this gene in the matrix
    # total insertions across all control samples and all sites
    total_insertions = sum(across(all_of(ctrl_samples)), na.rm = TRUE),
    # at the site level: is this TTAA hit in any control sample?
    n_sites_hit = sum(
      rowSums(across(all_of(ctrl_samples)), na.rm = TRUE) > 0
    ),
    prop_sites_hit = n_sites_hit / n_sites,
    .groups = "drop"
  )

head(gene_insertions)

# Join gene stats (CDS length + TTAA counts)

gene_stats <- read_xlsx("./Input/Genome/Pf3D7_gene_ttaa_len_summary_stats.xlsx")

# calculate features. Inner join becasue soem genes do not have CDS annotation
# and some genes have no insertion profile

gene_table <- inner_join(gene_stats, gene_insertions, by = "GeneID") %>%
  mutate(mean_insertions_per_TTAA = total_insertions / n_TTAA_cds,
         Msg_score = log10((total_insertions+1) / (n_TTAA_cds)))


# background essential/dispensable labels. Some background genes have no insertion numbers
background_genes <- read_xlsx(background_genes_file)

background_genes_table <- inner_join(background_genes, gene_table , by = 'GeneID') 

## Do a quick outlier detection
ggplot(background_genes_table, aes(x = cds_length, y = total_insertions, colour = class)) + 
  geom_point() + geom_smooth(method = 'lm') + theme_bw()

ggplot(background_genes_table, aes(x = cds_length, fill = class)) + 
  geom_histogram() + facet_grid(class~.) + theme_bw()

dis.ind <- background_genes_table$class_lab == 1
ess.ind <- background_genes_table$class_lab == 0

quantile(background_genes_table$cds_length, probs = 0.1)

# mark outliers
Qs <- background_genes_table %>% group_by(class) %>% 
  summarise(q1 = quantile(total_insertions, probs = 0.15),
            q2 = quantile(total_insertions, probs = 0.85))

background_genes_table <- background_genes_table %>% group_by(class) %>% 
  mutate(outlier = ifelse(class == 'dispensable' & total_insertions < Qs$q1[Qs$class == 'dispensable'], 1,
                          ifelse(class == 'essential' & total_insertions > Qs$q2[Qs$class == 'essential'], 1, 0)))

background_genes_table.h <- background_genes_table %>% dplyr::filter(outlier == 0)
ggplot(background_genes_table.h, aes(x = total_insertions, fill = class)) + 
  geom_histogram() + geom_density() + facet_grid(class~.) + theme_bw()


# Plot distributions for essential vs dispensable

plot_density <- function(df, variable, xlab) {
  ggplot(df, aes(x = .data[[variable]], color = class, fill = class)) +
    geom_density(alpha = 0.3) +
    theme_bw() +
    labs(x = xlab, y = "Density", title = paste("Distribution of", xlab))
}

p1 <- plot_density(background_genes_table.h, "n_TTAA_cds", "TTAA sites in CDS")
p2 <- plot_density(background_genes_table.h, "total_insertions", "Total Insertions (Control ≥ Day14)")
p3 <- plot_density(background_genes_table.h, "Msg_score", "Msg Score")

library(patchwork)
p1|p2|p3

## Base R density plots
dis.ind <- background_genes_table.h$class_lab == 1
ess.ind <- background_genes_table.h$class_lab == 0

vals_ess <- background_genes_table.h$Msg_score[ess.ind]
vals_dis <- background_genes_table.h$Msg_score[dis.ind]

d1 <- density(vals_ess, na.rm = TRUE)
d2 <- density(vals_dis, na.rm = TRUE)

## determine common x/y limits
xlim <- range(c(d1$x, d2$x))
ylim <- range(c(d1$y, d2$y))

## main density plot
plot(d1,
     col = "red",
     lwd = 2,
     xlim = xlim,
     ylim = ylim,
     xlab = "MSg (Control ≥ Day14)",
     ylab = "Density",
     main = "Density of MSg: Essential (red) vs Dispensable (green)")

lines(d2, col = "darkgreen", lwd = 2)

## add rugs (raw points along x-axis)
rug(vals_ess, col = rgb(1, 0, 0, 0.5), lwd = 1.5)
rug(vals_dis, col = rgb(0, 0.6, 0, 0.5), lwd = 1.5)

## legend
legend("topright",
       legend = c("Essential", "Dispensable"),
       col = c("red", "darkgreen"),
       lwd = 2,
       bty = "n")



## Boxplots
ggplot(background_genes_table.h, aes(x = class, y = total_insertions, fill = class)) +
  geom_boxplot(
    width = 0.5,
    outlier.size = 1.5,
    outlier.alpha = 0.6
  ) +
  ## Jittered points (raw data)
  geom_jitter(
    width = 0.20,
    alpha = 0.55,
    size = 1.2,
    show.legend = FALSE
  ) +
 
  scale_y_continuous(trans = "log10") +
  
  labs(
    x = "",
    y = "Total Insertions (log10 scale)",
    title = "Distribution of Insertions per Gene\nEssential vs Dispensable"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )



ggplot(background_genes_table.h, aes(x = class, y = Msg_score, fill = class)) +
  geom_boxplot(
    width = 0.5,
    outlier.size = 1.5,
    outlier.alpha = 0.6
  ) +
  ## Jittered points (raw data)
  geom_jitter(
    width = 0.20,
    alpha = 0.55,
    size = 1.2,
    show.legend = FALSE
  ) +
  
  #scale_y_continuous(trans = "log10") +
  
  labs(
    x = "",
    y = "MSg",
    title = "Distribution of MSg\nEssential vs Dispensable"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )


## glmnet
library(glmnet)
# High-confidence background genes only
train_df <- background_genes_table.h %>%
  filter(class %in% c("essential", "dispensable")) %>%
  mutate(
    y = ifelse(class == "essential", 0L, 1L)
  )

# Build predictor matrix X
X <- train_df %>% ungroup() %>% 
  dplyr::select(Msg_score, cds_length, n_TTAA_cds, prop_sites_hit,
         total_insertions, mean_insertions_per_TTAA) %>%
  as.matrix()

y <- train_df$y

# Remove rows with NAs
keep <- complete.cases(X, y)
X <- X[keep, , drop = FALSE]
y <- y[keep]
train_df <- train_df[keep, ]


## fit the path of lambda values
cvfit <- cv.glmnet(
  X, y,
  family = "binomial",
  alpha  = 1,        # LASSO
  nfolds = 10,
  type.measure = "deviance"
)

# Inspect CV curve
plot(cvfit)

lambda_min <- cvfit$lambda.min
lambda_1se <- cvfit$lambda.1se
lambda_min; lambda_1se

coef(cvfit, s = "lambda.min")
coef(cvfit, s = "lambda.1se")

train_df$prob_dis <- as.numeric(
  predict(cvfit, newx = X, s = "lambda.min", type = "response")
)

head(train_df %>% dplyr::select(GeneID, class, y, prob_dis))


## predict on all genes
X_all <- gene_table %>%
  dplyr::select(Msg_score, cds_length, n_TTAA_cds,prop_sites_hit,
         total_insertions, mean_insertions_per_TTAA) %>%
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
  mutate(color = case_when(
    class == "essential" ~ "red",
    class == "dispensable" ~ "green",
    TRUE ~ "grey80"
  ),
  is_bg = ifelse(is.na(class), FALSE, TRUE))
## ---- Plot ----

plot(
  df_plot$rank,
  df_plot$prob_dis_glmnet,
  pch = 20,
  cex = 0.6,
  col = df_plot$color,
  xlab = "Ranked Genes",
  ylab = "Posterior P(Essential)",
  main = "glmnet Posterior Essentiality Scores"
)

## Add background genes on TOP with bold styling
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
  col = c("red", "green", "grey80"),
  pt.bg = c("red", "green", "grey80"),
  pch = 21,
  pt.cex = c(1.6, 1.6, 0.8),
  bty = "n"
)


# confusion matrix
bg_eval <- df_plot %>%
  filter(!is.na(class)) %>%           # keep only background genes
  mutate(
    true = ifelse(class == "essential", 0, 1),
    pred = ifelse(prob_dis_glmnet >= 0.5, 1, 0)
  )

conf_mat <- table(
  Predicted = bg_eval$pred,
  True      = bg_eval$true
)

conf_mat


## Compare to JA false calls
## Load John Adams KO table
ko_table <- readxl::read_xlsx("./Input/JA_falsecalls/Gene KO table.xlsx", sheet = 4) %>%
  transmute(
    GeneID = Gene_ID,
    class = KO,
    geneName = geneName
  )

JA_plot <- left_join(ko_table, gene_table, by = 'GeneID')



JA_df_plot <- gene_table %>%
  arrange(prob_dis_glmnet) %>%
  mutate(rank = row_number()) %>%
  left_join(ko_table %>% dplyr::select(GeneID, class),
            by = "GeneID") %>%
  mutate(color = case_when(
    class == "Essential" ~ "red",
    class == "Disp" ~ "green",
    TRUE ~ "grey80"
  ),
  is_JA = ifelse(is.na(class), FALSE, TRUE))

plot(
  JA_df_plot$rank,
  JA_df_plot$prob_dis_glmnet,
  pch = 20,
  cex = 0.6,
  col = JA_df_plot$color,
  xlab = "Ranked Genes",
  ylab = "Posterior P(Essential)",
  main = "glmnet Posterior Essentiality Scores"
)

## Add background genes on TOP with bold styling
points(
  JA_df_plot$rank[JA_df_plot$is_JA],
  JA_df_plot$prob_dis_glmnet[JA_df_plot$is_JA],
  pch = 21,
  bg  = JA_df_plot$color[JA_df_plot$is_JA],
  col = "black",
  lwd = 1.2,
  cex = 1.6
)

legend(
  "topleft",
  legend = c("Essential (BG)", "Dispensable (BG)", "Other"),
  col = c("red", "green", "grey80"),
  pt.bg = c("red", "green", "grey80"),
  pch = 21,
  pt.cex = c(1.6, 1.6, 0.8),
  bty = "n"
)


## JA
library(dplyr)
library(ggplot2)
library(ggrepel)

## Build plotting table
JA_df_plot <- gene_table %>%
  arrange(prob_dis_glmnet) %>%
  mutate(rank = row_number()) %>%
  left_join(
    ko_table %>% 
      dplyr::select(GeneID, class, geneName),
    by = "GeneID"
  ) %>%
  mutate(
    class_simple = case_when(
      class == "Essential" ~ "Essential",
      class == "Disp"      ~ "Dispensable",
      TRUE                 ~ "Other"
    ),
    is_JA = !is.na(class)
  )

## Color palette for classes
class_colors <- c(
  "Essential"   = "red",
  "Dispensable" = "green4",
  "Other"       = "grey70"
)

## Base plot
p <- ggplot(JA_df_plot, aes(rank, prob_dis_glmnet)) +
  geom_point(aes(color = class_simple), size = 1.3, alpha = 0.8) +
  scale_color_manual(values = class_colors) +
  labs(
    x = "Ranked Genes",
    y = "Posterior P(Essential)",
    title = "glmnet Posterior Essentiality Scores",
    color = "Class"
  ) +
  theme_bw(base_size = 13)

## Add highlighted JA false-call genes (class-colored, black outline)
p <- p +
  geom_point(
    data = JA_df_plot %>% filter(is_JA),
    aes(rank, prob_dis_glmnet, fill = class_simple),
    shape = 21,
    size = 3.0,
    color = "black",
    stroke = 1.0
  ) +
  scale_fill_manual(values = class_colors, guide = "none") +
  geom_text_repel(
    data = JA_df_plot %>% filter(is_JA),
    aes(label = geneName),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.25,
    max.overlaps = Inf
  )

p

