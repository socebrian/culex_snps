## Sonia cebrian camison
## proyecto culex snps
# 7 febrero 2024
# Script to plot PCA analysis

#!/usr/bin/env Rscript

# Load libraries
library(tidyverse)

###############################################################################
# 1) Specify your sample names
###############################################################################
molestus <- c("pmo_SRR12125693", "pmo_SRR12125694", "pmo_SRR12125695", "pmo_SRR12125696", "pmo_SRR12125697")  # <-- EDIT to match your sample prefixes
pipiens <- c("ppi_SRR12125688", "ppi_SRR12125689", "ppi_SRR12125691", "ppi_SRR12125692", "ppi_SRR12125699", "ppi_SRR12125700", "ppi_SRR12125701")

species <-"pipiens"
samples <- pipiens
chr <- 3
###############################################################################
# 2) Read the GLOBAL coverage distribution files for each sample
###############################################################################
# Mosdepth global distribution files typically end in ".mosdepth.global.dist.txt"
# and have two columns: coverage, fraction
global_data_list <- lapply(samples, function(s) {
  file_gdist <- paste0(s, ".md.mosdepth.global.dist.txt")  # e.g., "sampleA.mosdepth.global.dist.txt"
  
  
  # Adjust `col_names` if your file's columns differ
  df <- read_tsv(file_gdist, 
                 col_names = c("chrom", "coverage", "fraction"),
                 comment = "#") %>%
    mutate(sample = s)
  return(df)
})

# Combine into one data frame
global_data <- bind_rows(global_data_list)

###############################################################################
# 3) Read the REGIONS coverage files (e.g., ".regions.bed.gz") for each sample
###############################################################################
# These files typically have columns: chrom, start, end, mean_coverage
regions_data_list <- lapply(samples, function(s) {
  file_regions <- paste0(s, ".md.regions.bed.gz")  # or .bed if uncompressed
  
  # If still gzipped, read it with read_tsv() directly (if supported) or use pipe("zcat ...")
  df <- read_tsv(file_regions, 
                 col_names = c("chrom", "start", "end", "mean_coverage"),
                 comment = "#") %>%
    mutate(sample = s,
           mid = (start + end) / 2)  # midpoint of each region for plotting
  return(df)
})

# Combine all samples into one data frame
regions_data <- bind_rows(regions_data_list)

# We assume you have only three chromosomes of interest: chr1, chr2, chr3
#chroms_of_interest <- c("1", "2", "3")

###############################################################################
# 4) PLOT 1: Global Coverage Distribution (all samples in one plot)
###############################################################################
p <- paste0("coverage_distribution_",species, ".png")
dist <- ggplot(global_data, aes(x=coverage, y=fraction, color=sample)) +
  geom_line(size=1, alpha=0.5) +
  theme_bw() +
  labs(
    title=paste0("Global Coverage Distribution Cx. pipiens ",species),
    x="Coverage",
    y="Proportion of Genome >= Coverage"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )+
  scale_x_continuous(limits = c(0, 140), breaks = seq(0, 140, 10))+
  scale_y_continuous(limits = c(0, 1.00), breaks = seq(0, 1.00, 0.10))
dist

ggsave(p, device = "png", units = "in", plot=dist, width=5, height=4)


###############################################################################
# 5) PLOT 2,3,4: Coverage Along Each Chromosome (three separate plots)
###############################################################################
# We'll loop over the three chromosomes, generate one PNG per chromosome
#for (chr in chroms_of_interest) {
  # Subset data to just this chromosome
  chr_data <- regions_data %>% filter(chrom == chr)
  # Construct output file name
  f <- paste0("coverage_",species, "_", chr, ".png")
  
  # Make plot


fig <- ggplot(chr_data, aes(x=mid/1e6, y=mean_coverage, color=sample)) +
    geom_line(size=0.8, alpha=0.5) +
    theme_bw() +
    labs(
      title = paste("Coverage Along chromosome", chr, "(All Samples)"),
      x = "Position Mbp (midpoint of window)",
      y = "Mean Depth"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )+
  coord_cartesian(ylim=c(0, 50))

fig
 
ggsave(f, device = "png", units = "in", plot=fig, width=7, height=4)
   
     