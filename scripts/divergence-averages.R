### Sonia Cebrián Camisón
## Proyecto Culex SNPs
# 28 junio 2025
# Script to calculate average of divergence measures along the genome

# Load libraries
library(dplyr)
library(readr)

#list files
pi_files  <- list.files(pattern = "pipiens1234_all_N90_chr[1,2,3]_win10000bp_pi.txt")
dxy_files <- list.files(pattern = "pipiens1234_all_N90_chr[1,2,3]_win10000bp_dxy.txt")
fst_files <- list.files(pattern = "pipiens1234_all_N90_chr[1,2,3]_win10000bp_fst.txt")
da_files <- list.files(pattern = "pipiens1234_all_N90_chr[1,2,3]_win10000bp_Da.txt")

# Load and bind all files
dxy_all <- bind_rows(lapply(dxy_files, read.table, header = TRUE))
fst_all <- bind_rows(lapply(fst_files, read.table, header = TRUE))
pi_all  <- bind_rows(lapply(pi_files, read.table, header = TRUE))
da_all  <- bind_rows(lapply(da_files, read.table, header = TRUE))

# DXY comparisons
dxy_all <- dxy_all %>%
  filter(!is.na(avg_dxy)) %>%
  mutate(comparison = case_when(
    (pop1 == "mol" & pop2 == "pip") | (pop1 == "pip" & pop2 == "mol") ~ "molestus vs pipiens",
    (pop1 == "mol" & pop2 == "mix") | (pop1 == "mix" & pop2 == "mol") ~ "molestus vs admixed",
    (pop1 == "pip" & pop2 == "mix") | (pop1 == "mix" & pop2 == "pip") ~ "pipiens vs admixed"
  ))

fst_all <- fst_all %>%
  filter(!is.na(avg_wc_fst)) %>%
  mutate(
  avg_wc_fst = ifelse(avg_wc_fst < 0, 0, avg_wc_fst), 
  comparison = case_when(
    (pop1 == "mol" & pop2 == "pip") | (pop1 == "pip" & pop2 == "mol") ~ "molestus vs pipiens",
    (pop1 == "mol" & pop2 == "mix") | (pop1 == "mix" & pop2 == "mol") ~ "molestus vs admixed",
    (pop1 == "pip" & pop2 == "mix") | (pop1 == "mix" & pop2 == "pip") ~ "pipiens vs admixed"
  ))

da_all <- da_all %>%
  filter(!is.na(avg_da)) %>%
  mutate(
    comparison = case_when(
      (pop1 == "mol" & pop2 == "pip") | (pop1 == "pip" & pop2 == "mol") ~ "molestus vs pipiens",
      (pop1 == "mol" & pop2 == "mix") | (pop1 == "mix" & pop2 == "mol") ~ "molestus vs admixed",
      (pop1 == "pip" & pop2 == "mix") | (pop1 == "mix" & pop2 == "pip") ~ "pipiens vs admixed"
    ))

pi_all <- pi_all %>%
  filter(!is.na(avg_pi)) %>%
  mutate(comparison = case_when(
    pop == "mol" ~ "molestus",
    pop == "pip" ~ "pipiens",
    pop == "mix" ~ "admixed"
  ))

# DXY
dxy_stats <- dxy_all %>%
  group_by(comparison) %>%
  summarise(
    stat = "DXY",
    mean = mean(avg_dxy, na.rm = TRUE),
    sd   = sd(avg_dxy, na.rm = TRUE),
    se   = sd / sqrt(n()),
    median = median(avg_dxy, na.rm = TRUE),
    n = n()
  )

# DXY
da_stats <- da_all %>%
  group_by(comparison) %>%
  summarise(
    stat = "DA",
    mean = mean(avg_da, na.rm = TRUE),
    sd   = sd(avg_da, na.rm = TRUE),
    se   = sd / sqrt(n()),
    median = median(avg_da, na.rm = TRUE),
    n = n()
  )

# FST
fst_stats <- fst_all %>%
  group_by(comparison) %>%
  summarise(
    stat = "FST",
    mean = mean(avg_wc_fst, na.rm = TRUE),
    sd   = sd(avg_wc_fst, na.rm = TRUE),
    se   = sd / sqrt(n()),
    median = median(avg_wc_fst, na.rm = TRUE),
    n = n()
  )

# π
pi_stats <- pi_all %>%
  group_by(comparison) %>%
  summarise(
    stat = "π",
    mean = mean(avg_pi, na.rm = TRUE),
    sd   = sd(avg_pi, na.rm = TRUE),
    se   = sd / sqrt(n()),
    median = median(avg_pi, na.rm = TRUE),
    n = n()
  )
dxy_stats

combined_stats <- bind_rows(dxy_stats, fst_stats, pi_stats, da_stats)
# Save to file
write.table(combined_stats, "genomewide_summary_stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)



###AVERAGES IN PEAKS VERSUS OUTSIDE PEAKS
# Peak regions in Mb
fst_all <- fst_all %>%
  mutate(window_mid = (window_pos_1 + window_pos_2) / 2 / 1e6)  # in Mb
peak_chr1_fst <- subset(fst_all, chromosome == "1" & window_mid >= 75 & window_mid <= 100)
peak_chr3_fst <- subset(fst_all, chromosome == "3" & window_mid >= 150 & window_mid <= 160)

# Combine all peak regions
fst_peaks <- bind_rows(peak_chr1_fst, peak_chr3_fst)

# Background = everything not in the peaks
fst_background <- anti_join(fst_all, fst_peaks)

# Summary for peaks
fst_peak1_stats <- peak_chr1_fst %>%
  group_by(comparison) %>%
  summarise(
    stat   = "FST",
    region = "peak1",
    mean   = mean(avg_wc_fst, na.rm = TRUE),
    sd     = sd(avg_wc_fst, na.rm = TRUE),
    se     = sd / sqrt(n()),
    median = median(avg_wc_fst, na.rm = TRUE),
    n      = n()
  )

fst_peak3_stats <- peak_chr3_fst %>%
  group_by(comparison) %>%
  summarise(
    stat   = "FST",
    region = "peak3",
    mean   = mean(avg_wc_fst, na.rm = TRUE),
    sd     = sd(avg_wc_fst, na.rm = TRUE),
    se     = sd / sqrt(n()),
    median = median(avg_wc_fst, na.rm = TRUE),
    n      = n()
  )
# Summary for background
fst_background_stats <- fst_background %>%
  group_by(comparison) %>%
  summarise(
    stat   = "FST",
    region = "background",
    mean   = mean(avg_wc_fst, na.rm = TRUE),
    sd     = sd(avg_wc_fst, na.rm = TRUE),
    se     = sd / sqrt(n()),
    median = median(avg_wc_fst, na.rm = TRUE),
    n      = n()
  )

# Combine both
fst_region_stats <- bind_rows(fst_peak1_stats, fst_peak3_stats, fst_background_stats)
write.table(fst_region_stats, "regions_fst_summary_stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Peak regions in Mb
da_all <- da_all %>%
  mutate(window_mid = (window_pos_1 + window_pos_2) / 2 / 1e6)  # in Mb
peak_chr1_da <- subset(da_all, chromosome == "1" & window_mid >= 75 & window_mid <= 100)
peak_chr3_da <- subset(da_all, chromosome == "3" & window_mid >= 150 & window_mid <= 160)

# Combine all peak regions
da_peaks <- bind_rows(peak_chr1_da, peak_chr3_da)

# Background = everything not in the peaks
da_background <- anti_join(da_all, da_peaks)

# Summary for peaks
da_peak1_stats <- peak_chr1_da %>%
  group_by(comparison) %>%
  summarise(
    stat   = "da",
    region = "peak1",
    mean   = mean(avg_da, na.rm = TRUE),
    sd     = sd(avg_da, na.rm = TRUE),
    se     = sd / sqrt(n()),
    median = median(avg_da, na.rm = TRUE),
    n      = n()
  )

da_peak3_stats <- peak_chr3_da %>%
  group_by(comparison) %>%
  summarise(
    stat   = "da",
    region = "peak3",
    mean   = mean(avg_da, na.rm = TRUE),
    sd     = sd(avg_da, na.rm = TRUE),
    se     = sd / sqrt(n()),
    median = median(avg_da, na.rm = TRUE),
    n      = n()
  )
# Summary for background
da_background_stats <- da_background %>%
  group_by(comparison) %>%
  summarise(
    stat   = "da",
    region = "background",
    mean   = mean(avg_da, na.rm = TRUE),
    sd     = sd(avg_da, na.rm = TRUE),
    se     = sd / sqrt(n()),
    median = median(avg_da, na.rm = TRUE),
    n      = n()
  )


# Combine both
da_region_stats <- bind_rows(da_peak1_stats, da_peak3_stats, da_background_stats)
write.table(da_region_stats, "regions_da_summary_stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

