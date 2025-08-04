# sonia cebrian camison
# 1/agosto/2025
# estructura culex pipiens
# this script is to detect fst windows with the more divergence between pipiens and molestus 
#and to plot a representative manhattan plot about that data


# ---- Load libraries ----
library(dplyr)
library(ggplot2)
library(readr)
library(scales)
library(tidyr)

# ---- Parameters ----


# ---- 1. Load FST files ----
fst1 <- read.table("pipiens1234_all_N90_chr1_win10000bp_fst.txt", header = TRUE)
fst2 <- read.table("pipiens1234_all_N90_chr2_win10000bp_fst.txt", header = TRUE)
fst3 <- read.table("pipiens1234_all_N90_chr3_win10000bp_fst.txt", header = TRUE)

# ---- 2. Filter for pipiens vs molestus and calculate window midpoints (Mb) ----

fst1 <- fst1 %>%
  filter((pop1 == "pip" & pop2 == "mol") | (pop1 == "mol" & pop2 == "pip")) %>%
  mutate(window_mid = (window_pos_1 + window_pos_2) / 2 / 1e6)
fst2 <- fst2 %>%
  filter((pop1 == "pip" & pop2 == "mol") | (pop1 == "mol" & pop2 == "pip")) %>%
  mutate(window_mid = (window_pos_1 + window_pos_2) / 2 / 1e6)
fst3 <- fst3 %>%
  filter((pop1 == "pip" & pop2 == "mol") | (pop1 == "mol" & pop2 == "pip")) %>%
  mutate(window_mid = (window_pos_1 + window_pos_2) / 2 / 1e6)


#--- 3. subset areas outside the inversion ----
fst1 <- subset(fst1, window_mid < 85.5 | window_mid > 102.5)
fst3 <- subset(fst3, window_mid < 3 | window_mid > 19)
fst3 <- subset(fst3, window_mid < 148.5 | window_mid > 162)

# ---4. merge files

fst_all <- rbind(fst1,fst2,fst3)

# calculate 1% highest fst values
threshold <- quantile(fst1$avg_wc_fst, 0.99, na.rm = TRUE)
top1 <- subset(fst1, avg_wc_fst >= threshold)
top2 <- subset(fst2, avg_wc_fst >= threshold)
top3 <- subset(fst3, avg_wc_fst >= threshold)
