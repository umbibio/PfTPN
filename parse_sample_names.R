library(readxl)   ## faster than openxlsx for reading
library(openxlsx)
library(dplyr)
library(stringr)
library(purrr)
library(writexl)
library(glue)

# 1) Read count matrix ---------------------------------------------------

cm_file <- "./Output/count_matrix/Pf_count_matrix_Location_conversion_for_exon_UTR_annotation_all.xlsx"
cm <- read_xlsx(cm_file)

# columns that are NOT samples
meta_cols <- c(
  "Chrom","Site","R1_pointing_downsteam","R1_pointing_upstream",
  "GeneID","gene.description","Location","Location1","Location2",
  "Present_in_any_samples","Run_repeats","Theoretical"
)

sample_cols <- setdiff(names(cm), meta_cols)


# 2) Build sample metadata from column names -----------------------------

sample_meta <- tibble(raw_name = sample_cols) %>%
  mutate(
    # High-level group
    group = case_when(
      str_detect(raw_name, "^AM\\d+")   ~ "AM",
      str_detect(raw_name, "^SM\\d+")   ~ "SM",
      str_detect(raw_name, "^Exp\\d+")  ~ "Exp",
      str_detect(raw_name, "^pTN9")     ~ "pTN9",
      str_detect(raw_name, "^Day\\d")   ~ "PTPN",
      str_detect(raw_name, "^D0_exp")   ~ "Plasmid_D0",
      TRUE                              ~ "Other"
    ),
    
    # Experiment number (for AM/SM/Exp batches)
    exp_number = case_when(
      group %in% c("AM", "SM") ~ as.integer(str_match(raw_name, "^[A-Z]{2}(\\d+)_")[, 2]),
      group == "Exp"           ~ as.integer(str_match(raw_name, "^Exp(\\d+)_")[, 2]),
      TRUE                     ~ NA_integer_
    ),
    
    # Day of collection
    day = case_when(
      str_detect(raw_name, "Day(\\d+)") ~ as.integer(str_match(raw_name, "Day(\\d+)")[, 2]),
      str_detect(raw_name, "_D(\\d+)")  ~ as.integer(str_match(raw_name, "_D(\\d+)")[, 2]),
      TRUE                              ~ NA_integer_
    ),
    
    # Treatment condition (single pass, robust)
    treatment = case_when(
      str_detect(raw_name, "Acet")                 ~ "Acet",
      str_detect(raw_name, "Mev")                  ~ "Mev",
      str_detect(raw_name, "Day\\d+_M(_|_0|_T)")   ~ "Mev",   # Day14_M_TPN etc.
      str_detect(raw_name, "D\\d+M_")              ~ "Mev",   # AM2_D17M_0507_TPN etc.
      TRUE                                         ~ "Control"
    ),
    
    # Batch; allow optional _2_ etc. between batch and TPN
    batch = str_match(raw_name, "_(05\\d{2})(?:_\\d+)?_TPN$")[, 2],
    
    # Person (optional)
    person = case_when(
      group == "AM" ~ "Alexis",
      group == "SM" ~ "Sebastian",
      TRUE          ~ NA_character_
    ),
    
    # Time-course: everything except plasmid D0
    is_timecourse = group != "Plasmid_D0"
  ) %>%
  group_by(group, exp_number, day, treatment, batch) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  mutate(
    # optional label pieces
    exp_label   = if_else(is.na(exp_number), "", glue("_exp{exp_number}")),
    batch_label = if_else(is.na(batch),       "", glue("_b{batch}")),
    rep_label   = if_else(replicate <= 1,     "", glue("_rep{replicate}")),
    
    # final standardized names
    new_name = case_when(
      # plasmid controls
      group == "Plasmid_D0" ~
        glue("Plasmid_D0_{treatment}{batch_label}{rep_label}"),
      
      # time-course samples (PTPN, pTN9, AM, SM, Exp, …)
      is_timecourse ~
        glue("TC_{group}{exp_label}_day{day}_{treatment}{batch_label}{rep_label}"),
      
      # fallback (shouldn't really be used with current data)
      TRUE ~ glue("{group}_{treatment}{batch_label}{rep_label}")
    ),
    
    new_name = as.character(new_name)
  )

# Quick sanity checks
print(sample_meta, n = nrow(sample_meta))
any(duplicated(sample_meta$new_name))


# 3) Rename columns in the count matrix ---------------------------------

# map raw_name -> new_name
name_map <- setNames(sample_meta$new_name, sample_meta$raw_name)

cm_renamed <- cm %>%
  rename_with(~ name_map[.x], .cols = all_of(sample_cols))


# 4) Write outputs -------------------------------------------------------

write_xlsx(sample_meta,   "./Output/PfTPN_sample_metadata.xlsx")
write_xlsx(cm_renamed,    "./Output/PfTPN_count_matrix_renamed.xlsx")
