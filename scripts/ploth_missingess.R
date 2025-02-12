## Sonia cebrian camison
## proyecto culex snps
# 11/02/2025
# Script to plot Missing data per individuo of a vcf file

# Load ggplot2
library(ggplot2)

# If df is not yet in your environment, read it in:
df <- read.table("2202pipiens123_pmolestus_3_filteredGT_pruned.imiss ", header = TRUE, stringsAsFactors = FALSE)
sp <- read.table("sample-species.txt ", header = TRUE, stringsAsFactors = FALSE)


merged <- merge(df, sp, by.x = "INDV", by.y = "id")
merged$species <- as.factor(merged$species)
str(merged)

chrom <- 3
f <- paste0("pipiens123_pmol_MISSind_", chrom, ".png")

# Create a scatter plot of F_MISS*100 vs INDV
fig <- ggplot(merged, aes(x = INDV, y = F_MISS * 100, colour = species)) +
  geom_point(size = 2) + 
  labs(
    x = "Sample (INDV)",
    y = "Missing Data (%)",
    title = paste0("Percentage of Missing Data per Sample. Chromosome ", chrom)
  ) +
  theme_bw() +
  # Rotate x-axis labels if needed
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
        legend.position = "bottom",
        legend.direction = "horizontal"  # optional, makes legend items line up horizontally
  )
fig
ggsave(f, device = "png", units = "in", plot=fig, width=7, height=4)

