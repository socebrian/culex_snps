## Sonia cebrian camison
## proyecto culex snps
# 27 enero 2024
# Script to plot PCA analysis

library(ggplot2)

# read eigenvalues
eigenvals <- read.table("2202pipiens123_pmolestus_3_filteredGT_pruned_PCA.eigenval",
                        col.names = c("eigenval"))
# convert eigenvalues to percentage
percents <- eigenvals$eigenval/sum(eigenvals$eigenval)*100

# read eigenvector
eigenvec <- read.table("2202pipiens123_pmolestus_3_filteredGT_pruned_PCA.eigenvec", sep="\t", header=T)
#eigenvec$batch <- as.factor(eigenvec$batch)
eigenvec$species <- as.factor(eigenvec$species)
eigenvec$id <- as.factor(eigenvec$id)

#eigenvec <- eigenvec[eigenvec$id != "PP1122",]
#eigenvec <- eigenvec[eigenvec$id != "PP1001",]
#variables
f1 <- "2202pipiens123_pmol_3.pdf"
f2 <- "pca12_2202pipiens123_pmol_3.png"
pc.x <- eigenvec$X
pc.y <- eigenvec$X.1
x_label <- paste0("PC1-", round(percents[1], 2), "%")
y_label <- paste0("PC2-", round(percents[2], 2), "%")
maint <- "Chromosome 3"

# scree-plot of eigenvalues
(fig1 <- ggplot() +
  geom_line(aes(x = as.numeric(row.names(eigenvals)), y = percents)) +
  geom_point(aes(x = as.numeric(row.names(eigenvals)), y = percents),
    fill = "grey", shape = 21, size = 2.5) +
  xlab(paste0("Number of Principal Component")) +
  ylab(paste0("Percentage of variance Explained")) +
  scale_x_continuous(n.breaks = 10) +
  theme_bw())

ggsave(f1, device = "png", units = "in", plot=fig1, width=5, height=4)

#### this first part is from https://github.com/Enricobazzi/EBD-PopGenIntro/blob/main/practica2/practica2.md
###it was adapted to plot my data

(fig2 <- ggplot() +
  geom_point(
    data = eigenvec, alpha = 0.5,
    aes(x = pc.x, y = pc.y, colour = species) #, shape = prov),
    #size  = 2.5 #,
    #color = "black"   # <--- Give an explicit outline color
  ) +
  #scale_shape_manual(values = c(22, 24, 25, 21)) +
  #scale_fill_brewer(palette = "Set1") +
  labs(title= maint)+
  xlab(x_label) +
  ylab(y_label) +
  theme_bw() +
  # Override the legend so the shapes in the *fill* legend also appear filled
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, color = "black")),
    shape = guide_legend(override.aes = list(color = "black"))
  ))

ggsave(f2, device = "png", units = "in", plot=fig2, width=5, height=4)
