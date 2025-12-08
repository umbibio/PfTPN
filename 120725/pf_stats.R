library(readxl)
library(dplyr)
library(ggplot2)


#
# 1. Read in files
# 2. Break control samples up by day
# 3. Get data for each day
# 4. Read in verified genes
# 5. Create data frames for data/verified genes
#




###
#
# Read in Files 
#
###

#Chrom	Site	R1_pointing_downsteam	R1_pointing_upstream	GeneID	gene.description	Location	Location1	Location2	Present_in_any_samples	Run_repeats	Theoretical	TC_pTN9_day0_Control	TC_pTN9_day15_Control	TC_pTN9_day15_Mev	TC_pTN9_day20_Control	TC_pTN9_day20_Mev	TC_pTN9_day8_Control	TC_pTN9_day8_Mev	Plasmid_D0_Control	TC_PTPN_day14_Mev	TC_PTPN_day14_Control	TC_PTPN_day17_Mev	TC_PTPN_day17_Control	TC_PTPN_day20_Mev	TC_PTPN_day20_Control	TC_PTPN_day8_Control	TC_PTPN_day8_Mev	Plasmid_D0_Control_b0504	Plasmid_D0_Control_b0504_rep2	TC_PTPN_day14_Control_b0504	TC_PTPN_day14_Control_b0504_rep2	TC_PTPN_day14_Mev_b0504	TC_PTPN_day14_Mev_b0504_rep2	TC_PTPN_day17_Control_b0504	TC_PTPN_day17_Mev_b0504	TC_PTPN_day18_Control_b0504	TC_PTPN_day18_Mev_b0504	TC_PTPN_day20_Control_b0504	TC_PTPN_day20_Control_b0504_rep2	TC_PTPN_day20_Mev_b0504	TC_PTPN_day20_Mev_b0504_rep2	TC_PTPN_day8_Control_b0504	TC_PTPN_day8_Mev_b0504	TC_PTPN_day9_Control_b0504	TC_PTPN_day9_Mev_b0504	Plasmid_D0_Control_b0505	Plasmid_D0_Control_b0505_rep2	TC_PTPN_day14_Acet	TC_PTPN_day14_Control_b0505	TC_PTPN_day14_Control_b0505_rep2	TC_PTPN_day14_Mev_b0505	TC_PTPN_day14_Mev_b0505_rep2	TC_PTPN_day17_Control_b0505	TC_PTPN_day17_Mev_b0505	TC_PTPN_day18_Acet	TC_PTPN_day18_Control_b0505	TC_PTPN_day18_Mev_b0505	TC_PTPN_day20_Acet	TC_PTPN_day20_Control_b0505	TC_PTPN_day20_Control_b0505_rep2	TC_PTPN_day20_Mev_b0505	TC_PTPN_day20_Mev_b0505_rep2	TC_PTPN_day9_Acet	TC_AM_exp2_day17_Mev_b0507	TC_AM_exp2_day19_Control_b0507	TC_AM_exp2_day19_Mev_b0507	TC_AM_exp3_day28_Control_b0507	TC_AM_exp3_day28_Mev_b0507	TC_AM_exp8_day17_Control_b0507	TC_AM_exp8_day17_Mev_b0507	TC_AM_exp8_day20_Control_b0507	TC_AM_exp8_day20_Mev_b0507	TC_AM_exp8_day26_Control_b0507	TC_AM_exp8_day26_Mev_b0507	TC_AM_exp8_day8_Control_b0507	TC_AM_exp8_day8_Mev_b0507	TC_SM_exp5_day16_Control_b0507	TC_SM_exp5_day16_Mev_b0507	TC_SM_exp6_day14_Control_b0507	TC_SM_exp6_day14_Mev_b0507	TC_SM_exp6_day16_Control_b0507	TC_SM_exp6_day16_Mev_b0507	TC_SM_exp6_day18_Control_b0507	TC_SM_exp6_day18_Mev_b0507	TC_SM_exp6_day20_Control_b0507	TC_SM_exp6_day20_Mev_b0507	TC_SM_exp6_day8_Control_b0507	TC_SM_exp6_day8_Mev_b0507	TC_SM_exp7_day14_Control_b0507	TC_SM_exp7_day14_Mev_b0507	TC_SM_exp7_day16_Control_b0507	TC_SM_exp7_day16_Mev_b0507	TC_SM_exp7_day18_Control_b0507	TC_SM_exp7_day18_Mev_b0507	TC_SM_exp7_day20_Control_b0507	TC_SM_exp7_day20_Mev_b0507	TC_SM_exp7_day8_Control_b0507	TC_SM_exp7_day8_Mev_b0507	TC_Exp_exp3_day15_Control_b0508	TC_Exp_exp3_day15_Mev_b0508	TC_Exp_exp3_day17_Control_b0508	TC_Exp_exp3_day17_Mev_b0508	TC_Exp_exp3_day28_Control_b0508	TC_Exp_exp3_day28_Mev_b0508	TC_Exp_exp5_day14_Acet_b0508	TC_Exp_exp5_day16_Control_b0508	TC_Exp_exp5_day16_Acet_b0508	TC_Exp_exp5_day16_Mev_b0508	TC_Exp_exp5_day18_Acet_b0508	TC_Exp_exp5_day20_Acet_b0508	TC_Exp_exp6_day14_Acet_b0508	TC_Exp_exp6_day16_Control_b0508	TC_Exp_exp6_day16_Acet_b0508	TC_Exp_exp6_day16_Mev_b0508	TC_Exp_exp6_day18_Control_b0508	TC_Exp_exp6_day18_Acet_b0508	TC_Exp_exp6_day18_Mev_b0508	TC_Exp_exp6_day20_Acet_b0508	TC_Exp_exp6_day40_Control_b0508	TC_Exp_exp6_day40_Control_b0508_rep2	TC_Exp_exp6_day40_Acet_b0508	TC_Exp_exp6_day40_Acet_b0508_rep2	TC_Exp_exp6_day40_Mev_b0508	TC_Exp_exp6_day40_Mev_b0508_rep2	TC_Exp_exp6_day8_Control_b0508	TC_Exp_exp6_day8_Acet_b0508	TC_Exp_exp6_day8_Mev_b0508	TC_Exp_exp7_day14_Control_b0508	TC_Exp_exp7_day14_Acet_b0508	TC_Exp_exp7_day14_Mev_b0508	TC_Exp_exp7_day16_Control_b0508	TC_Exp_exp7_day16_Acet_b0508	TC_Exp_exp7_day16_Mev_b0508	TC_Exp_exp7_day18_Control_b0508	TC_Exp_exp7_day18_Acet_b0508	TC_Exp_exp7_day18_Mev_b0508	TC_Exp_exp7_day20_Acet_b0508	TC_Exp_exp7_day8_Control_b0508	TC_Exp_exp7_day8_Acet_b0508	TC_Exp_exp7_day8_Mev_b0508	TC_Exp_exp8_day17_Control_b0508	TC_Exp_exp8_day17_Mev_b0508	TC_Exp_exp8_day20_Control_b0508	TC_Exp_exp8_day20_Mev_b0508	TC_Exp_exp8_day26_Control_b0508	TC_Exp_exp8_day26_Mev_b0508
pf_count_matrix_file <- read_xlsx("Output/PfTPN_count_matrix_renamed_UTR_overlaps_corrected_120225.xlsx")

# raw_name	group	exp_number	day	treatment	batch	person	is_timecourse	replicate	exp_label	batch_label	rep_label	new_name
pf_metadata <- read_xlsx("Output/PfTPN_sample_metadata_UTR_overlaps_corrected_120225.xlsx")

# GeneID	chrom	ttaa_start	ttaa_end	strand	cds_length	ttaa_cds_pos	ttaa_cds_frac	n_TTAA_cds	Site_ID	TC_pTN9_day0_Control	TC_pTN9_day15_Control	TC_pTN9_day15_Mev	TC_pTN9_day20_Control	TC_pTN9_day20_Mev	TC_pTN9_day8_Control	TC_pTN9_day8_Mev	Plasmid_D0_Control	TC_PTPN_day14_Mev	TC_PTPN_day14_Control	TC_PTPN_day17_Mev	TC_PTPN_day17_Control	TC_PTPN_day20_Mev	TC_PTPN_day20_Control	TC_PTPN_day8_Control	TC_PTPN_day8_Mev	Plasmid_D0_Control_b0504	Plasmid_D0_Control_b0504_rep2	TC_PTPN_day14_Control_b0504	TC_PTPN_day14_Control_b0504_rep2	TC_PTPN_day14_Mev_b0504	TC_PTPN_day14_Mev_b0504_rep2	TC_PTPN_day17_Control_b0504	TC_PTPN_day17_Mev_b0504	TC_PTPN_day18_Control_b0504	TC_PTPN_day18_Mev_b0504	TC_PTPN_day20_Control_b0504	TC_PTPN_day20_Control_b0504_rep2	TC_PTPN_day20_Mev_b0504	TC_PTPN_day20_Mev_b0504_rep2	TC_PTPN_day8_Control_b0504	TC_PTPN_day8_Mev_b0504	TC_PTPN_day9_Control_b0504	TC_PTPN_day9_Mev_b0504	Plasmid_D0_Control_b0505	Plasmid_D0_Control_b0505_rep2	TC_PTPN_day14_Acet	TC_PTPN_day14_Control_b0505	TC_PTPN_day14_Control_b0505_rep2	TC_PTPN_day14_Mev_b0505	TC_PTPN_day14_Mev_b0505_rep2	TC_PTPN_day17_Control_b0505	TC_PTPN_day17_Mev_b0505	TC_PTPN_day18_Acet	TC_PTPN_day18_Control_b0505	TC_PTPN_day18_Mev_b0505	TC_PTPN_day20_Acet	TC_PTPN_day20_Control_b0505	TC_PTPN_day20_Control_b0505_rep2	TC_PTPN_day20_Mev_b0505	TC_PTPN_day20_Mev_b0505_rep2	TC_PTPN_day9_Acet	TC_AM_exp2_day17_Mev_b0507	TC_AM_exp2_day19_Control_b0507	TC_AM_exp2_day19_Mev_b0507	TC_AM_exp3_day28_Control_b0507	TC_AM_exp3_day28_Mev_b0507	TC_AM_exp8_day17_Control_b0507	TC_AM_exp8_day17_Mev_b0507	TC_AM_exp8_day20_Control_b0507	TC_AM_exp8_day20_Mev_b0507	TC_AM_exp8_day26_Control_b0507	TC_AM_exp8_day26_Mev_b0507	TC_AM_exp8_day8_Control_b0507	TC_AM_exp8_day8_Mev_b0507	TC_SM_exp5_day16_Control_b0507	TC_SM_exp5_day16_Mev_b0507	TC_SM_exp6_day14_Control_b0507	TC_SM_exp6_day14_Mev_b0507	TC_SM_exp6_day16_Control_b0507	TC_SM_exp6_day16_Mev_b0507	TC_SM_exp6_day18_Control_b0507	TC_SM_exp6_day18_Mev_b0507	TC_SM_exp6_day20_Control_b0507	TC_SM_exp6_day20_Mev_b0507	TC_SM_exp6_day8_Control_b0507	TC_SM_exp6_day8_Mev_b0507	TC_SM_exp7_day14_Control_b0507	TC_SM_exp7_day14_Mev_b0507	TC_SM_exp7_day16_Control_b0507	TC_SM_exp7_day16_Mev_b0507	TC_SM_exp7_day18_Control_b0507	TC_SM_exp7_day18_Mev_b0507	TC_SM_exp7_day20_Control_b0507	TC_SM_exp7_day20_Mev_b0507	TC_SM_exp7_day8_Control_b0507	TC_SM_exp7_day8_Mev_b0507	TC_Exp_exp3_day15_Control_b0508	TC_Exp_exp3_day15_Mev_b0508	TC_Exp_exp3_day17_Control_b0508	TC_Exp_exp3_day17_Mev_b0508	TC_Exp_exp3_day28_Control_b0508	TC_Exp_exp3_day28_Mev_b0508	TC_Exp_exp5_day14_Acet_b0508	TC_Exp_exp5_day16_Control_b0508	TC_Exp_exp5_day16_Acet_b0508	TC_Exp_exp5_day16_Mev_b0508	TC_Exp_exp5_day18_Acet_b0508	TC_Exp_exp5_day20_Acet_b0508	TC_Exp_exp6_day14_Acet_b0508	TC_Exp_exp6_day16_Control_b0508	TC_Exp_exp6_day16_Acet_b0508	TC_Exp_exp6_day16_Mev_b0508	TC_Exp_exp6_day18_Control_b0508	TC_Exp_exp6_day18_Acet_b0508	TC_Exp_exp6_day18_Mev_b0508	TC_Exp_exp6_day20_Acet_b0508	TC_Exp_exp6_day40_Control_b0508	TC_Exp_exp6_day40_Control_b0508_rep2	TC_Exp_exp6_day40_Acet_b0508	TC_Exp_exp6_day40_Acet_b0508_rep2	TC_Exp_exp6_day40_Mev_b0508	TC_Exp_exp6_day40_Mev_b0508_rep2	TC_Exp_exp6_day8_Control_b0508	TC_Exp_exp6_day8_Acet_b0508	TC_Exp_exp6_day8_Mev_b0508	TC_Exp_exp7_day14_Control_b0508	TC_Exp_exp7_day14_Acet_b0508	TC_Exp_exp7_day14_Mev_b0508	TC_Exp_exp7_day16_Control_b0508	TC_Exp_exp7_day16_Acet_b0508	TC_Exp_exp7_day16_Mev_b0508	TC_Exp_exp7_day18_Control_b0508	TC_Exp_exp7_day18_Acet_b0508	TC_Exp_exp7_day18_Mev_b0508	TC_Exp_exp7_day20_Acet_b0508	TC_Exp_exp7_day8_Control_b0508	TC_Exp_exp7_day8_Acet_b0508	TC_Exp_exp7_day8_Mev_b0508	TC_Exp_exp8_day17_Control_b0508	TC_Exp_exp8_day17_Mev_b0508	TC_Exp_exp8_day20_Control_b0508	TC_Exp_exp8_day20_Mev_b0508	TC_Exp_exp8_day26_Control_b0508	TC_Exp_exp8_day26_Mev_b0508
pf_stats_per_site_counts <- read_xlsx("Input/Genome/Pf3D7_gene_ttaa_rel_loc_len_counts_per_sample_1225_mcs.xlsx")
# GeneID	chrom	ttaa_start	ttaa_end	strand	cds_length	ttaa_cds_pos	ttaa_cds_frac	n_TTAA_cds	Site_ID
pf_stats <- read_xlsx("Input/Genome/Pf3D7_gene_ttaa_rel_loc_len_counts_per_sample_1225_mcs.xlsx")

# GeneID	cds_length	n_TTAA_cds
pf_stats_per_gene <- read_xlsx("Input/Genome/Pf3D7_gene_ttaa_len_summary_stats_1225.xlsx")

# Columns that are NOT samples
meta_cols <- c(
  "Chrom","Site","R1_pointing_downsteam","R1_pointing_upstream",
  "GeneID","gene.description","Location","Location1","Location2",
  "Present_in_any_samples","Run_repeats","Theoretical"
)

sample_cols <- setdiff(colnames(pf_total_per_site), meta_cols)

###########
#
# Break up samples by day
#
###########


ctrl_samples_d8 <- sample_meta %>%
  filter(treatment == "Control", !is.na(day), day %in% c(8,9)) %>%
  pull(new_name)


ctrl_samples_d14 <- sample_meta %>%
  filter(treatment == "Control", !is.na(day), day %in% c(14,15)) %>%
  pull(new_name)


ctrl_samples_d16 <- sample_meta %>%
  filter(treatment == "Control", !is.na(day), day >= 16, day <= 19) %>%
  pull(new_name)


ctrl_samples_d20 <- sample_meta %>%
  filter(treatment == "Control", !is.na(day), day >= 20, day <= 25) %>%
  pull(new_name)

ctrl_samples_d26 <- sample_meta %>%
  filter(treatment == "Control", !is.na(day), day >=26) %>%
  pull(new_name)




## probably can turn into a function that does this for each, but going to repeat everything for now
## look into exons only for control samples

###########
#
# Get stats for each day... 
#
###########

######p 8 
pf_8 <- pf_total_per_site %>% dplyr::filter(Location == "exon") %>% 
  dplyr::select(GeneID, all_of(ctrl_samples_d8))  %>%
  dplyr::mutate(site_total= rowSums(across(all_of(ctrl_samples_d8)), na.rm=TRUE)) %>%
  group_by(GeneID) %>%
  summarise(
    TTAA_sites = n(),
    n_sites_hit = sum(site_total >0, na.rm=TRUE),
    max_total= max(site_total,na.rm =TRUE),
    Yes_Ratio = n_sites_hit/TTAA_sites, 
    
    Total_Per_Gene = sum(site_total, na.rm =TRUE),
    Total = Total_Per_Gene / pmax(TTAA_sites,1),
    MSG = log10((Total_Per_Gene+1)/pmax(Total,1))
  ) 


pf8_stats <- pf_stats_per_site_counts %>% 
  dplyr::select(GeneID, ttaa_cds_frac, all_of(ctrl_samples_d8))  %>%
  dplyr::mutate(site_total= rowSums(across(all_of(ctrl_samples_d8)), na.rm=TRUE)) %>%
  dplyr::group_by(GeneID) %>%
  dplyr::summarise(
    first_hit= ifelse(
      any(site_total>0),
      min(ttaa_cds_frac[site_total >0]), 
                          1.0), 
    first_substantial_hit= ifelse(any(site_total>3),
                                  min(ttaa_cds_frac[site_total >3]), 
                                      1.0),
    .groups="drop"
  )%>% mutate(day_num=8)

  





##### p14 
pf_14 <- pf_total_per_site %>% dplyr::filter(Location == "exon") %>% 
  dplyr::select(GeneID, all_of(ctrl_samples_d14))  %>%
  dplyr::mutate(site_total= rowSums(across(all_of(ctrl_samples_d14)), na.rm=TRUE)) %>%
  group_by(GeneID) %>%
  summarise(
    TTAA_sites = n(),
    n_sites_hit = sum(site_total >0, na.rm=TRUE),
    max_total= max(site_total,na.rm =TRUE),
    Yes_Ratio = n_sites_hit/TTAA_sites, 
    
    Total_Per_Gene = sum(site_total, na.rm =TRUE),
    Total = Total_Per_Gene / pmax(TTAA_sites,1),
    MSG = log10((Total_Per_Gene+1)/pmax(Total,1))
  )



pf14_stats <- pf_stats_per_site_counts %>% 
  dplyr::select(GeneID, ttaa_cds_frac, all_of(ctrl_samples_d14))  %>%
  dplyr::mutate(site_total= rowSums(across(all_of(ctrl_samples_d14)), na.rm=TRUE)) %>%
  dplyr::group_by(GeneID) %>%
  dplyr::summarise(
    first_hit= ifelse(
      any(site_total>0),
      min(ttaa_cds_frac[site_total >0]), 
      1.0), 
    first_substantial_hit= ifelse(any(site_total>3),
                                  min(ttaa_cds_frac[site_total >3]), 
                                  1.0),
    .groups="drop"
  )%>% mutate(day_num=14)








pf_16 <- pf_total_per_site %>% dplyr::filter(Location == "exon") %>% 
  dplyr::select(GeneID, all_of(ctrl_samples_d16))  %>%
  dplyr::mutate(site_total= rowSums(across(all_of(ctrl_samples_d16)), na.rm=TRUE)) %>%
  group_by(GeneID) %>%
  summarise(
    TTAA_sites = n(),
    n_sites_hit = sum(site_total >0, na.rm=TRUE),
    max_total= max(site_total,na.rm =TRUE),
    Yes_Ratio = n_sites_hit/TTAA_sites, 
    
    Total_Per_Gene = sum(site_total, na.rm =TRUE),
    Total = Total_Per_Gene / pmax(TTAA_sites,1),
    MSG = log10((Total_Per_Gene+1)/pmax(Total,1))
  )





pf16_stats <- pf_stats_per_site_counts %>% 
  dplyr::select(GeneID, ttaa_cds_frac, all_of(ctrl_samples_d16))  %>%
  dplyr::mutate(site_total= rowSums(across(all_of(ctrl_samples_d16)), na.rm=TRUE)) %>%
  dplyr::group_by(GeneID) %>%
  dplyr::summarise(
    first_hit= ifelse(
      any(site_total>0),
      min(ttaa_cds_frac[site_total >0]), 
      1.0), 
    first_substantial_hit= ifelse(any(site_total>3),
                                  min(ttaa_cds_frac[site_total >3]), 
                                  1.0),
    .groups="drop"
  )%>% mutate(day_num=16)




pf_20 <- pf_total_per_site %>% dplyr::filter(Location == "exon") %>% 
  dplyr::select(GeneID, all_of(ctrl_samples_d20))  %>%
  dplyr::mutate(site_total= rowSums(across(all_of(ctrl_samples_d20)), na.rm=TRUE)) %>%
  group_by(GeneID) %>%
  summarise(
    TTAA_sites = n(),
    n_sites_hit = sum(site_total >0, na.rm=TRUE),
    max_total= max(site_total,na.rm =TRUE),
    Yes_Ratio = n_sites_hit/TTAA_sites, 
    
    Total_Per_Gene = sum(site_total, na.rm =TRUE),
    Total = Total_Per_Gene / pmax(TTAA_sites,1),
    MSG = log10((Total_Per_Gene+1)/pmax(Total,1))
  )


pf20_stats <- pf_stats_per_site_counts %>% 
  dplyr::select(GeneID, ttaa_cds_frac, all_of(ctrl_samples_d20))  %>%
  dplyr::mutate(site_total= rowSums(across(all_of(ctrl_samples_d20)), na.rm=TRUE)) %>%
  dplyr::group_by(GeneID) %>%
  dplyr::summarise(
    first_hit= ifelse(
      any(site_total>0),
      min(ttaa_cds_frac[site_total >0]), 
      1.0), 
    first_substantial_hit= ifelse(any(site_total>3),
                                  min(ttaa_cds_frac[site_total >3]), 
                                  1.0),
    .groups="drop"
  )%>% mutate(day_num=20)



pf_26 <- pf_total_per_site %>% dplyr::filter(Location == "exon") %>% 
  dplyr::select(GeneID, all_of(ctrl_samples_d26))  %>%
  dplyr::mutate(site_total= rowSums(across(all_of(ctrl_samples_d26)), na.rm=TRUE)) %>%
  group_by(GeneID) %>%
  summarise(
    TTAA_sites = n(),
    n_sites_hit = sum(site_total >0, na.rm=TRUE),
    max_total= max(site_total,na.rm =TRUE),
    Yes_Ratio = n_sites_hit/TTAA_sites, 
    
    Total_Per_Gene = sum(site_total, na.rm =TRUE),
    Total = Total_Per_Gene / pmax(TTAA_sites,1),
    MSG = log10((Total_Per_Gene+1)/pmax(Total,1))
  )




pf26_stats <- pf_stats_per_site_counts %>% 
  dplyr::select(GeneID, ttaa_cds_frac, all_of(ctrl_samples_d26))  %>%
  dplyr::mutate(site_total= rowSums(across(all_of(ctrl_samples_d26)), na.rm=TRUE)) %>%
  dplyr::group_by(GeneID) %>%
  dplyr::summarise(
    first_hit= ifelse(
      any(site_total>0),
      min(ttaa_cds_frac[site_total >0]), 
      1.0), 
    first_substantial_hit= ifelse(any(site_total>3),
                                  min(ttaa_cds_frac[site_total >3]), 
                                  1.0),
    .groups="drop"
  )%>% mutate(day_num=26)



##########
#
# read in verified genes
#
#########

essentiality_path <- "112025_files/pk_pf_comparison/essential_gene_lists/pf_pk_all_verified_genes.xlsx"

# read pf essential genes
pf_essential <- read_xlsx(essentiality_path, sheet="pf_essential")


# filter pf_total_per_site
pf_essential_ids <- pf_essential$GeneID

pf_essential_list <- pf_total_per_site %>%
  dplyr::filter(GeneID %in% pf_essential_ids)


nrow(pf_essential_list)
head(pf_essential_list)


pf_dispensable <-read_xlsx(essentiality_path, sheet="pf_dispensable")
pf_dispensable_ids <- pf_dispensable$GeneID

pf_dispensable_list <- pf_total_per_site %>%
  dplyr::filter(GeneID %in% pf_dispensable_ids)


nrow(pf_dispensable_list)
head(pf_dispensable_list)




#
#
# create data frames
#
#


pf_all <- bind_rows(
  pf_8 %>% mutate(day_bin="Day8", day_num = 8),
  pf_14 %>% mutate(day_bin="Day14", day_num = 14),
  pf_16 %>% mutate(day_bin="Day16", day_num = 16),
  pf_20 %>% mutate(day_bin="Day20", day_num = 20),
  pf_26 %>% mutate(day_bin="Day26", day_num = 26)
  ) %>%
  mutate(day_bin = factor(day_bin, levels = c("Day8","Day14","Day16","Day20","Day26")))

pf_stats_all <- dplyr::bind_rows(
  pf8_stats, pf14_stats, pf16_stats, pf20_stats, pf26_stats
)

pf_verified_all <- pf_all %>%
  mutate(class = case_when(
    GeneID %in% pf_essential_ids ~ "Essential",
    GeneID %in% pf_dispensable_ids ~ "Dispensable"
  )) %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::mutate(class=factor(class, levels=c("Essential", "Dispensable")))





pf_verified_all_by_site <- pf_stats_all %>%
  mutate(class = case_when(
    GeneID %in% pf_essential_ids ~ "Essential",
    GeneID %in% pf_dispensable_ids ~ "Dispensable"
  )) %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::mutate(class=factor(class, levels=c("Essential", "Dispensable")))



#
#
# plots
#
#



#
#
#. Log of Number of Sites hit over time
#
#



p <- ggplot(pf_verified_all, aes(x=day_num, y = log1p(n_sites_hit), group = GeneID)) +
  geom_line(alpha=0.3, linewidth=0.3, color="grey70")+
  geom_point() +
  geom_smooth(aes(group=1),method="lm",color="cyan", linewidth=1, level=0.95)+
  facet_wrap(~ class) +
  scale_x_continuous(breaks = c(8,14,16,20,26))+
  #scale_color_manual(values = c(Mean="red", Median="cyan",name =""))+
  theme(strip.text=element_text(size=14), 
        plot.title = element_text(size=20)
        )+
  labs(title="Log of Number of Sites Hit",x="Day", y ="Number of Sites Hit")

p


w <- wilcox.test(log1p(n_sites_hit) ~ class, data= pf_verified_all)
w

t <- t.test(log1p(n_sites_hit) ~ class, data=pf_verified_all)
t



days <- c(8,14,16,20,26)

#############################
#
# Log of Number of Sites Hit 
#
#############################


p_box <- ggplot(pf_verified_all, aes(x=day_num, y=log1p(n_sites_hit), fill=class)) +
  geom_boxplot(aes(group=day_num), outlier.shape=NA) +
  geom_jitter(width=0.85,size=0.6, alpha=0.15, color="grey20")+
  geom_smooth(aes(group=1),method="lm",color="black", linewidth=1, level=0.95)+
  facet_wrap(~ class) +
  guides(fill="none")+ # hides redundant legend. 
  scale_x_continuous(breaks=days)+
  theme(strip.text=element_text(size=14),         plot.title = element_text(size=20))+
  labs(title="Log of Number of Sites Hit",x="Day", y ="Log of Number of Sites Hit")

p_box


### need to center at day8... 
pf_verified_all <- pf_verified_all %>%
  dplyr::mutate(day_centered_at_8 = day_num-8)

# gets y-intercept at d8
m <- lm(log1p(n_sites_hit) ~ class * day_centered_at_8, data=pf_verified_all)
anova(m)
summary(m)










#############################
#
# Ratio of Number of Sites Hit 
#
#############################




p_box_ratio_sites_hit <- ggplot(pf_verified_all, aes(x=day_num, y=Yes_Ratio, fill=class)) +
  geom_boxplot(aes(group=day_num), outlier.shape=NA) +
  geom_jitter(width=0.85,size=0.6, alpha=0.15, color="grey20")+
  geom_smooth(aes(group=1),method="lm",color="black", linewidth=1, level=0.95)+
  facet_wrap(~ class) +
  guides(fill="none")+ # hides redundant legend. 
  scale_x_continuous(breaks=days)+
  theme(strip.text=element_text(size=14)    ,    plot.title = element_text(size=20))+
  labs(title="Percent of TTAA Sites Hit",x="Day", y ="Percent of Sites Hit of Sites Hit")

p_box_ratio_sites_hit


### need to center at day8... 
pf_verified_all <- pf_verified_all %>%
  dplyr::mutate(day_centered_at_8 = day_num-8)

# gets y-intercept at d8
m <- lm(Yes_Ratio ~ class * day_centered_at_8, data=pf_verified_all)
anova(m)
summary(m)



#############################
#
# Ratio of Number of Sites Hit without day 8
#
#############################


p_box_ratio_sites_hit_no_d8 <- ggplot(dplyr::filter(pf_verified_all, day_num!=8), aes(x=day_num, y=Yes_Ratio, fill=class)) +
  geom_boxplot(aes(group=day_num), outlier.shape=NA) +
  geom_jitter(width=0.85,size=0.6, alpha=0.15, color="grey20")+
  geom_smooth(aes(group=1),method="lm",color="black", linewidth=1, level=0.95)+
  facet_wrap(~ class) +
  guides(fill="none")+ # hides redundant legend. 
  scale_x_continuous(breaks=days)+
  theme(strip.text=element_text(size=14)   ,     plot.title = element_text(size=20))+
  labs(title="Percent of Sites Hit Over Time (Excluding D8)",x="Day", y ="Percent of Sites Hit")

p_box_ratio_sites_hit_no_d8



### need to center at day1... 
pf_verified_all <- pf_verified_all %>%
  dplyr::mutate(day_centered_at_14 = day_num-14)



# gets y-intercept at d8
m <- lm(Yes_Ratio ~ class * day_centered_at_14, data=dplyr::filter(pf_verified_all, day_num!=8))
anova(m)
summary(m)






#############################
#
# site with max total... 
#
#############################




p_box_most_TTAA_hits <- ggplot(pf_verified_all, aes(x=day_num, y=log10(pmax(max_total,1)), fill=class)) +
  geom_boxplot(aes(group=day_num), outlier.shape=NA) +
  geom_jitter(width=0.85,size=0.6, alpha=0.15, color="grey20")+
  geom_smooth(aes(group=1),method="lm",color="black", linewidth=1, level=0.95)+
  facet_wrap(~ class) +
  guides(fill="none")+ # hides redundant legend. 
  scale_x_continuous(breaks=days)+
  theme(strip.text=element_text(size=14),        plot.title = element_text(size=20))+
  labs(title="Site with Most TTAA Hits",x="Day", y ="Site with Most TTAA Hits")

p_box_most_TTAA_hits


### need to center at day8... 
pf_verified_all <- pf_verified_all %>%
  dplyr::mutate(day_centered_at_8 = day_num-8)

# gets y-intercept at d8
m <- lm(log10(pmax(max_total,1)) ~ class * day_centered_at_8, data=pf_verified_all)
anova(m)
summary(m)




###########
#
# MSG
#
###########


p_box_MSG <- ggplot(pf_verified_all, aes(x=day_num, y=MSG, fill=class)) +
  geom_boxplot(aes(group=day_num), outlier.shape=NA) +
  geom_jitter(width=0.85,size=0.6, alpha=0.15, color="grey20")+
  geom_smooth(aes(group=1),method="lm",color="black", linewidth=1, level=0.95)+
  facet_wrap(~ class) +
  guides(fill="none")+ # hides redundant legend. 
  scale_x_continuous(breaks=days)+
  theme(strip.text=element_text(size=14),         plot.title = element_text(size=20))+
  labs(title="MSG score",x="Day", y ="MSG score")

p_box_MSG


### need to center at day8... 
pf_verified_all <- pf_verified_all %>%
  dplyr::mutate(day_centered_at_8 = day_num-8)

# gets y-intercept at d8
m <- lm(MSG ~ class * day_centered_at_8, data=pf_verified_all)
anova(m)
summary(m)




trend_by_gene_MSG <- pf_verified_all %>%
  group_by(GeneID, class) %>%
  group_modify(~{
    fit <- lm(MSG~ day_centered_at_8, data=.x)
    s <- summary(fit)$coefficients
    tibble(
      beta = s["day_centered_at_8", "Estimate"],
      p=s["day_centered_at_8","Pr(>|t|)"],
      se = s["day_centered_at_8", "Std. Error"],
      n=nrow(.x)
    )
  })%>%
  ungroup() %>%
  mutate(p_adj=p.adjust(p,"BH")) %>%
  arrange(beta)

trend_by_gene_MSG <- trend_by_gene_MSG %>%
  mutate(rank = row_number())

x <- ggplot(trend_by_gene_MSG, aes(x=rank, y=beta, color=class))+
  geom_point()+
  geom_hline(yintercept = 0, linetype=2)+
  labs(title= "Slope of MSG Score vs Day")
x
wilcox.test(beta ~ class, data=trend_by_gene_MSG)
   

fit_w <- lm(beta ~ class, data=trend_by_gene_MSG)
summary(fit_w)


anova(lm(beta~ class, data= trend_by_gene_MSG))



###########
#
# Total per Gene
#
###########

p_box_tot <- ggplot(pf_verified_all, aes(x=day_num, y=log(pmax(Total_Per_Gene,1)), fill=class)) +
  geom_boxplot(aes(group=day_num), outlier.shape=NA) +
  geom_jitter(width=0.85,size=0.6, alpha=0.15, color="grey20")+
  geom_smooth(aes(group=1),method="lm",color="black", linewidth=1, level=0.95)+
  facet_wrap(~ class) +
  guides(fill="none")+ # hides redundant legend. 
  scale_x_continuous(breaks=days)+
  theme(strip.text=element_text(size=14),        plot.title = element_text(size=20))+
  labs(title="Total per Gene",x="Day", y ="Total Per Gene")

p_box_tot


### need to center at day8... 
pf_verified_all <- pf_verified_all %>%
  dplyr::mutate(day_centered_at_8 = day_num-8)

# gets y-intercept at d8
m <- lm(log(pmax(Total_Per_Gene,1)) ~ class * day_centered_at_8, data=pf_verified_all)
anova(m)
summary(m)






###########
#
# First location of CDS disruption 
#
###########
p_first_loc <- ggplot(pf_verified_all_by_site, aes(x=day_num, y=first_hit, fill=class)) +
  geom_boxplot(aes(group=day_num), outlier.shape=NA) +
  geom_jitter(width=0.85,size=0.6, alpha=0.15, color="grey20")+
  geom_smooth(aes(group=1),method="lm",color="black", linewidth=1, level=0.95)+
  facet_wrap(~ class) +
  guides(fill="none")+ # hides redundant legend. 
  scale_x_continuous(breaks=days)+
  theme(strip.text=element_text(size=14)  ,    plot.title = element_text(size=20))+
  labs(title="First Disruption in CDS",x="Day", y ="First Disruption in CDS")

p_first_loc


### need to center at day8... 
pf_verified_all_by_site <- pf_verified_all_by_site %>%
  dplyr::mutate(day_centered_at_8 = day_num-8)

# gets y-intercept at d8
m <- lm(first_hit ~ class * day_centered_at_8, data=pf_verified_all_by_site)
anova(m)
summary(m)


#####
#
# First substantial disruption 
#
#####


###########
#
# First location of CDS disruption 
#
###########
p_first_substantial_loc <- ggplot(pf_verified_all_by_site, aes(x=day_num, y=first_substantial_hit, fill=class)) +
  geom_boxplot(aes(group=day_num), outlier.shape=NA) +
  geom_jitter(width=0.85,size=0.6, alpha=0.15, color="grey20")+
  geom_smooth(aes(group=1),method="lm",color="black", linewidth=1, level=0.95)+
  facet_wrap(~ class) +
  guides(fill="none")+ # hides redundant legend. 
  scale_x_continuous(breaks=days)+
  theme(strip.text=element_text(size=14)  ,    plot.title = element_text(size=20))+
  labs(title="First Disruption in CDS >3 Hits ",x="Day", y ="First Substantial Disruption in CDS")

p_first_substantial_loc


### need to center at day8... 
pf_verified_all_by_site <- pf_verified_all_by_site %>%
  dplyr::mutate(day_centered_at_8 = day_num-8)

# gets y-intercept at d8
m <- lm(first_substantial_hit ~ class * day_centered_at_8, data=pf_verified_all_by_site)
anova(m)
summary(m)





