# sonia cebrian camison
# 1/agosto/2025
# estructura culex pipiens
# this script is to plot enrichment analyses


# Load libraries
library(ggplot2)
library(readr)
library(dplyr)

# Load enrichment results (significant terms only)
BP1q <- load_enrich("1q_enrich_BP_significant.tsv", "1q", "BP")
MF1q <- load_enrich("1q_enrich_MF_significant.tsv","1q", "MF")
CC1q <- load_enrich("1q_enrich_CC_significant.tsv", "1q", "CC")

BP3p <- load_enrich("3p_enrich_BP_significant.tsv", "3p", "BP")
MF3p <- load_enrich("3p_enrich_MF_significant.tsv","3p", "MF")
CC3p <- load_enrich("3p_enrich_CC_significant.tsv","3p", "CC")

BP3q <- load_enrich("3q_enrich_BP_significant.tsv", "3q", "BP")
MF3q <- load_enrich("3q_enrich_MF_significant.tsv", "3q", "MF")
CC3q <- load_enrich("3q_enrich_CC_significant.tsv", "3q", "CC")
# Optional: clean long GO term names
enrich_df$Term <- gsub(" \\(.*\\)", "", enrich_df$Term)

# Create GeneRatio column (Significant / Annotated)
enrich_df <- enrich_df %>%
  mutate(GeneRatio = Significant / Annotated)

# Plot
ggplot(enrich_df, aes(x = GeneRatio, y = reorder(Term, GeneRatio))) +
  geom_point(aes(size = Significant, color = Fisher)) +
  scale_color_gradient(low = "red", high = "blue", name = "p-value") +
  labs(
    title = "GO Enrichment Dotplot (Biological Process)",
    x = "Gene Ratio (Significant / Annotated)",
    y = "GO Term",
    size = "Gene Count"
  ) +
  theme_minimal(base_size = 12)

# Function to load and tag enrichment file
load_enrich <- function(file, inversion, category) {
  df <- read_tsv(file, show_col_types = FALSE) %>%
    mutate(
      inversion = inversion,
      category = category,
      Term = gsub(" \\(.*\\)", "", Term),  # Clean long names
      GeneRatio = Significant / Annotated
    )
  return(df)
}

# Load each enrichment file (you must adapt paths and labels)
BP1q <- read_tsv("1q_enrich_BP_significant.tsv")
MF1q <- read_tsv("1q_enrich_MF_significant.tsv")
CC1q <- read_tsv("1q_enrich_CC_significant.tsv")

BP3p <- read_tsv("3p_enrich_BP_significant.tsv")
MF3p <- read_tsv("3p_enrich_MF_significant.tsv")
CC3p <- read_tsv("3p_enrich_CC_significant.tsv")

BP3q <- read_tsv("3q_enrich_BP_significant.tsv")
MF3q <- read_tsv("3q_enrich_MF_significant.tsv")
CC3q <- read_tsv("3q_enrich_CC_significant.tsv")

# Combine all
all_enrich <- bind_rows(BP1q, BP3p, BP3q, CC1q, CC3p, CC3q, MF1q, MF3p, MF3q)

# Plot
combined_plot <- ggplot(all_enrich, aes(x = GeneRatio, y = reorder(Term, GeneRatio))) +
  geom_point(aes(size = Significant, color = Fisher)) +
  scale_color_gradient(low = "red", high = "blue", name = "p-value") +
  facet_grid(category ~ inversion, scales = "free_y", space = "free_y") +
  labs(
    title = "GO Enrichment Dotplot",
    x = "Gene Ratio (Significant / Annotated)",
    y = "GO Term",
    size = "Gene Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))+
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA))


ggsave("inversion_enrichment2.png", plot = combined_plot, width = 10, height = 15, dpi = 300)
