## Sonia cebrian camison
## proyecto culex snps
# 27 enero 2024
# Script to plot PCA analysis

library(ggplot2)

# read eigenvalues
eigenvals <- read.table("2202pipiens123_3_f1n_BDQMC_pruned_PCA.eigenval",
                        col.names = c("eigenval"))
# convert eigenvalues to percentage
percents <- eigenvals$eigenval/sum(eigenvals$eigenval)*100

# read eigenvector
eigenvec <- read.table("2202pipiens123_1_f1n_BDQMC_pruned_PCA.eigenvec", sep="\t", header=T)
eigenvec$batch <- as.factor(eigenvec$batch)
eigenvec$id <- as.factor(eigenvec$id)
eigenvec <- eigenvec[eigenvec$id != "PP1122",]
eigenvec <- eigenvec[eigenvec$id != "PP1001",]
#variables
scree_plot_file <- "mod_eigenval_2202pipiens123_3_f1n_BDQMC.pdf"
pca_plot_file <- "mod_pca12_2202pipiens123_1_f1n_BDQMC.pdf"
pc.x <- eigenvec$X
pc.y <- eigenvec$X.1
x_label <- paste0("PC1-", round(percents[1], 2), "%")
y_label <- paste0("PC2-", round(percents[2], 2), "%")

# scree-plot of eigenvalues
(scree <- ggplot() +
  geom_line(aes(x = as.numeric(row.names(eigenvals)), y = percents)) +
  geom_point(aes(x = as.numeric(row.names(eigenvals)), y = percents),
    fill = "grey", shape = 21, size = 2.5) +
  xlab(paste0("Number of Principal Component")) +
  ylab(paste0("Percentage of variance Explained")) +
  scale_x_continuous(n.breaks = 10) +
  theme_bw())

ggsave(filename = scree_plot_file, plot = scree)

#### this first part is from https://github.com/Enricobazzi/EBD-PopGenIntro/blob/main/practica2/practica2.md
###it was adapted to plot my data

(pca <- ggplot() +
  geom_point(
    data = eigenvec,
    aes(x = pc.x, y = pc.y, fill = prov, shape = batch),
    size  = 2.5,
    color = "black"   # <--- Give an explicit outline color
  ) +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_fill_brewer(palette = "Set1") +
  xlab(x_label) +
  ylab(y_label) +
  theme_bw() +
  # Override the legend so the shapes in the *fill* legend also appear filled
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, color = "black")),
    shape = guide_legend(override.aes = list(color = "black"))
  ))

ggsave(filename = pca_plot_file, plot = pca)
