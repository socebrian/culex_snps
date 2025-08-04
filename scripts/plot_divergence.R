### Sonia Cebrián Camisón
## Proyecto Culex SNPs
# 1 junio 2025
# Script to plot DXY, FST, and π along chromosome or genome

# ---- Load libraries ----
library(ggplot2)
library(dplyr)
library(tidyverse)
library(patchwork)

plot_title <- "Divergence along chromosome 1"
fig <- "DA_divergence_pipmol_chr1.png"

fst_title <- "FST divergence along chromosome 2"
fig_fst <- "fst_pipmol_chrom2.png"

title_median <- "FST divergence along chromosome 3 (pipiens vs molestus)"
fig_median <- "fst_only_pipmol_chrom3.png"

plot_date_title <- "Divergence in chromosome 2 (μ =4.85e-9, ,  g =20 days)"
fig2 <- "line_div-diff_chr2_pipmol.png"
dxy_summ_file <- "div_times_pipmol_pipiens_ch2.csv"


inv_start <- 4.8
inv_end   <- 17.5
inv_start2 <- 150
inv_end2 <- 160
avg_div <- 0.0096

# ---- 1. Load Pixy output ----
pi  <- read.table("pipiens1234_all_N90_chr3_win10000bp_pi.txt", header = TRUE)
dxy <- read.table("pipiens1234_all_N90_chr3_win10000bp_dxy.txt", header = TRUE)
fst <- read.table("pipiens1234_all_N90_chr3_win10000bp_fst.txt", header = TRUE)
da <- read.table("pipiens1234_all_N90_chr3_win10000bp_Da.txt", header = TRUE)

try <- subset(fst, fst$avg_wc_fst < 0)
min(try$avg_wc_fst)
#6000 de 29.9mil chrom 1, 17mil de 55mil en ch2 y  9mil de 46mil en chrom3
#valores hasta -0.09

# ---- 2. Add window midpoint (in Mb) ----
pi  <- pi  %>% mutate(window_mid = (window_pos_1 + window_pos_2) / 2 / 1e6)
dxy <- dxy %>% mutate(window_mid = (window_pos_1 + window_pos_2) / 2 / 1e6)
fst <- fst %>% mutate(window_mid = (window_pos_1 + window_pos_2) / 2 / 1e6)
da <- da %>% mutate(window_mid = (window_pos_1 + window_pos_2) / 2 / 1e6)

# ---- 3. Extract 3 comparisons for DXY and FST ----
dxy_all <- bind_rows(
  dxy %>% filter((pop1 == "mol" & pop2 == "pip") | (pop1 == "pip" & pop2 == "mol")) %>% mutate(comparison = "molestus vs pipiens"),
  dxy %>% filter((pop1 == "mol" & pop2 == "mix") | (pop1 == "mix" & pop2 == "mol")) %>% mutate(comparison = "molestus vs admixed"),
  dxy %>% filter((pop1 == "pip" & pop2 == "mix") | (pop1 == "mix" & pop2 == "pip")) %>% mutate(comparison = "pipiens vs admixed")
)

fst_all <- bind_rows(
  fst %>% filter((pop1 == "mol" & pop2 == "pip") | (pop1 == "pip" & pop2 == "mol")) %>% mutate(comparison = "molestus vs pipiens"),
  fst %>% filter((pop1 == "mol" & pop2 == "mix") | (pop1 == "mix" & pop2 == "mol")) %>% mutate(comparison = "molestus vs admixed"),
  fst %>% filter((pop1 == "pip" & pop2 == "mix") | (pop1 == "mix" & pop2 == "pip")) %>% mutate(comparison = "pipiens vs admixed")
)

da_all <- bind_rows(
  da %>% filter((pop1 == "mol" & pop2 == "pip") | (pop1 == "pip" & pop2 == "mol")) %>% mutate(comparison = "molestus vs pipiens"),
  da %>% filter((pop1 == "mol" & pop2 == "mix") | (pop1 == "mix" & pop2 == "mol")) %>% mutate(comparison = "molestus vs admixed"),
  da %>% filter((pop1 == "pip" & pop2 == "mix") | (pop1 == "mix" & pop2 == "pip")) %>% mutate(comparison = "pipiens vs admixed")
)

# ---- 4. Format for plotting ----
# Filter valid rows and label by population
pi_long <- pi %>%
  filter(!is.na(avg_pi)) %>%
  mutate(comparison = case_when(
    pop == "mol" ~ "molestus",
    pop == "pip" ~ "pipiens",
    pop == "mix" ~ "admixed"
  )) %>%
  select(window_mid, value = avg_pi, comparison) %>%
  mutate(stat = "π")


dxy_long <- dxy_all %>%
  select(window_mid, value = avg_dxy, comparison) %>%
  mutate(stat = "DXY")

fst_long <- fst_all %>%
  select(window_mid, value = avg_wc_fst, comparison) %>%
  mutate(stat = "FST")

da_long <- da_all %>%
  select(window_mid, value = avg_da, comparison) %>%
  mutate(stat = "DA")

plot_data <- bind_rows(pi_long, dxy_long, fst_long, da_long)

# ---- 5. Define inversion region (in Mb) ----
#movido arriba

# ---- 6. Plot ----# Combine all
# Plot

# Filter by stat
plot_pi  <- plot_data %>% filter(stat == "π")
plot_dxy_fst_da <- plot_data %>% filter(stat %in% c("DXY", "FST", "DA"))

# Plot DXY + FST
p1 <- ggplot(plot_dxy_fst_da, aes(x = window_mid, y = value, color = comparison)) +
  geom_smooth(method = "loess", span = 0.005, se = FALSE, linewidth = 0.7) +
  annotate("rect", xmin = inv_start, xmax = inv_end, ymin = -Inf, ymax = Inf, fill = "#44B7C2", alpha = 0.1, color = NA) +
  annotate("rect", xmin = inv_start2, xmax = inv_end2, ymin = -Inf, ymax = Inf, fill = "#44B7C2", alpha = 0.1, color = NA) +
  facet_grid(stat ~ ., scales = "free_y", switch = "y") +
  labs(x = NULL, y = NULL, color = " ") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text.y = element_text(angle = 0),
    panel.spacing.y = unit(1, "lines"),
    legend.position = "none"
  ) +
  scale_color_manual(
    values = c(
      "molestus vs pipiens"  = "#cf1259",
      "molestus vs admixed" = "#D38217",
      "pipiens vs admixed"  = "#034B7B"
    )
  )

p1
# Plot π separately
p2 <- ggplot(plot_pi, aes(x = window_mid, y = value, color = comparison)) +
  geom_smooth(method = "loess", span = 0.005, se = FALSE, linewidth = 0.7) +
 annotate("rect", xmin = inv_start, xmax = inv_end, ymin = -Inf, ymax = Inf,  fill = "#44B7C2", alpha = 0.1, color = NA) +
  annotate("rect", xmin = inv_start2, xmax = inv_end2, ymin = -Inf, ymax = Inf, fill = "#44B7C2", alpha = 0.1, color = NA) +
  facet_grid(stat ~ ., scales = "free_y", switch = "y") +
  labs(x = "Genomic position (Mb)", y = NULL, color = " ") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text.y = element_text(angle = 0),
    panel.spacing.y = unit(1, "lines"),
    legend.position = "right"
  ) +
  scale_color_manual(
    values = c(
      "molestus" = "#FFAE4A",
      "pipiens"  = "#44B7C2",
      "admixed"  = "grey50"
    )
  )

# Combine with patchwork
p_ch3 <- p1 / p2 + plot_layout(heights = c(3, 1))

# First, remove individual legends (we’ll add a single global one)
p_ch1_clean <- p_ch1 + theme(legend.position = "none") + plot_annotation(title = "Chromosome 1")
p_ch2_clean <- p_ch2 + theme(legend.position = "none") + plot_annotation(title = "Chromosome 2")
p_ch3_clean <- p_ch3 + theme(legend.position = "none") + plot_annotation(title = "Chromosome 3")
p_ch3_clean
# Combine horizontally, adjust layout
combined_plot <- (p_ch1_clean | p_ch2_clean | p_ch3_clean) +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

# Print combined plot
print(combined_plot)

# Save to file
ggsave("divergence_figure_provisional.png", plot = combined_plot, width = 14.5, height = 10, dpi = 300)

# Create the plot using points

p_fst_point <- ggplot(fst_all, aes(x = window_mid, y = avg_wc_fst)) +
  geom_point(aes(color = comparison), size = 1, alpha = 0.5) +
  annotate("rect", xmin = inv_start, xmax = inv_end, ymin = -Inf, ymax = Inf,
           fill = "#44B7C2", alpha = 0.1, color = NA) +
  annotate("rect", xmin = inv_start2, xmax = inv_end2, ymin = -Inf, ymax = Inf,
           fill = "#44B7C2", alpha = 0.1, color = NA) +
  labs(
    x = "Genomic position (Mb)",
    y = "FST",
    color = "Comparison",
    title = fst_title
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  scale_color_manual(
    values = c(
      "molestus vs pipiens"  = "#cf1259",
      "molestus vs admixed" = "#D38217",
      "pipiens vs admixed"  = "#034B7B"
    )
  )

print(p_fst_point)

ggsave(fig_fst, plot = p_fst_point, width = 10, height = 6, dpi = 300)

## solocomarsion mol-pip con su media
#calculate median
median_pipmol_fst <- fst_all %>%
  filter(comparison == "molestus vs pipiens") %>%
  summarise(median_fst = median(avg_wc_fst, na.rm = TRUE)) %>%
  pull(median_fst)
#create column sayin if its over or under median
fst_all <- fst_all %>%
  mutate(fst_class_pipmol = ifelse(avg_wc_fst > median_pipmol_fst, "above", "below"))

# Step 1: Filter for only "pipiens vs molestus" (regardless of order)
fst_filtered <- fst_all %>%
  filter((pop1 == "pip" & pop2 == "mol") | (pop1 == "mol" & pop2 == "pip"))

# Step 4: Plot
fst_median <- ggplot(fst_filtered, aes(x = window_mid, y = avg_wc_fst)) +
  geom_point(aes(color = fst_class_pipmol), size = 1, alpha = 0.8, shape = 18) +  # diamond shape like your example
 annotate("rect", xmin = inv_start, xmax = inv_end, ymin = -Inf, ymax = Inf,
          fill = "#44B7C2", alpha = 0.1, color = NA) +
 annotate("rect", xmin = inv_start2, xmax = inv_end2, ymin = -Inf, ymax = Inf,
           fill = "#44B7C2", alpha = 0.1, color = NA) +
  scale_color_manual(
    values = c("below" = "gray80", "above" = "#cf1259"),
    guide = "none"  # hides legend
  ) +
  labs(
    title = title_median,
    x = "Genomic position (Mb)",
    y = "FST"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
  )

print(fst_median)
ggsave(fig_median, plot = fst_median, width = 10, height = 6, dpi = 300)

#----------------- Divergence difference plots--------------------#
# Define comparison labels based on pop1 and pop2 values
dxy_inv <- dxy %>%
  filter(window_mid >= inv_start, window_mid <= inv_end)

#dxy_summary <- dxy_inv %>%
dxy_summary <- dxy %>%
  mutate(comparison = case_when(
    (pop1 == "mol" & pop2 == "pip") | (pop1 == "pip" & pop2 == "mol") ~ "molestus vs pipiens",
    (pop1 == "mol" & pop2 == "mix") | (pop1 == "mix" & pop2 == "mol") ~ "molestus vs admixed",
    (pop1 == "pip" & pop2 == "mix") | (pop1 == "mix" & pop2 == "pip") ~ "pipiens vs admixed"
  )) %>%
  group_by(comparison) %>%
  summarise(mean_dxy = mean(avg_dxy, na.rm = TRUE),
            sd_dxy = sd(avg_dxy, na.rm = TRUE),
            n = n(),
            se_dxy = sd_dxy / sqrt(n))


# Define parameter sets: mutation rate and generational time (based on Haba et al 2025)
params <- data.frame(
  label = c("Best-guess", "Minimum", "Maximum"),
  mu = c(4.85e-9, 8e-9, 1e-9),
  g  = c(0.0548, 0.0548, 0.0822)
)

# Expand grid to compute all combinations
expanded <- merge(dxy_summary, params)

# Calculate divergence time
expanded <- expanded %>%
  mutate(
    generations = mean_dxy / (2 * mu),
    years = generations * g
  )

# View output
print(expanded)

# Add second axis (divergence time in generations)
# Set μ  and gfor plotting secondary axis (Best-guess)

g <- 0.0548   # generation time in years
mu <- 4.85e-9  # mutation rate


p2 <- ggplot(dxy_summary, aes(x = comparison, y = mean_dxy)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_dxy - se_dxy, ymax = mean_dxy + se_dxy),
                width = 0.2) +
  # geom_hline(yintercept = 2000, color = "red")+
   labs(
    x = NULL,
    title = plot_date_title
  ) +
  scale_y_continuous(
    name = "Mean DXY",
    sec.axis = sec_axis(
      transform = ~ . / (2 * mu) * g,
      name = "Divergence time (years)"
    )
  ) +
  theme_minimal(base_size = 14)+
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(p2)
ggsave(fig2, plot = p2, width = 6, height = 6, dpi = 300)
write.csv(expanded, dxy_summ_file, row.names = FALSE)


## MEANS OF ESTIMATES FOR RESULTS
# Mean DXY genome-wide (not restricted to inversion)
# DXY
dxy_means <- dxy_all %>%
  group_by(comparison) %>%
  summarise(
    stat = "DXY",
    mean_value = mean(avg_dxy, na.rm = TRUE),
    sd_value   = sd(avg_dxy, na.rm = TRUE),
    se_value   = sd_value / sqrt(n()),
    min_value  = min(avg_dxy, na.rm = TRUE),
    max_value  = max(avg_dxy, na.rm = TRUE),
    n          = n()
  )

# FST
fst_means <- fst_all %>%
  filter(!is.na(avg_wc_fst), avg_wc_fst >= 0) %>%  # remove negative FSTs
  group_by(comparison) %>%
  summarise(
    stat = "FST",
    mean_value = mean(avg_wc_fst, na.rm = TRUE),
    sd_value   = sd(avg_wc_fst, na.rm = TRUE),
    se_value   = sd_value / sqrt(n()),
    min_value  = min(avg_wc_fst, na.rm = TRUE),
    max_value  = max(avg_wc_fst, na.rm = TRUE),
    n          = n()
  )

# π
pi_means <- pi_long %>%
  group_by(comparison) %>%
  summarise(
    stat = "π",
    mean_value = mean(value, na.rm = TRUE),
    sd_value   = sd(value, na.rm = TRUE),
    se_value   = sd_value / sqrt(n()),
    min_value  = min(value, na.rm = TRUE),
    max_value  = max(value, na.rm = TRUE),
    n          = n()
  )

print(pi_means)

combined_means <- bind_rows(dxy_means, fst_means, pi_means)

write.table(combined_means, "summary_divergence_stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

