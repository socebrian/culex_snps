## sonia cebrian camison
## culex snps project
## 10/5/2025
## scebrian27@gmail.com /sonia.cebrian@ebd.csic.es

#this iscript is to plot phylogenetic trees from mitochondrials equence allignments

library(ape)
library(ggtree)
library(ggplot2)

tree <- read.tree("mitoc_pipiens1234_tree.contree")
rooted_tree <- root(tree, outgroup = "Anopheles_cruzii", resolve.root = TRUE)
ultra_tree <- chronos(rooted_tree)


p <- ggtree(ultra_tree)

# Now extract the data frame from the plot
tree_data <- p$data

# Plot again, using data filtering inline
# Define the tree plot object
p <- ggtree(ultra_tree) +
  geom_tiplab(size = 2.5, align = TRUE, linesize = 0.3) +
  geom_text2(data = subset(p$data, !isTip),
             aes(label = label),
             size = 2.5,  # same as tip labels now
             hjust = -0.2) +
  theme_tree()

# Export as a tall, narrow PNG at high resolution
ggsave("ultrametric_tree_labeled.png",
       plot = p,
       width = 8,        # inches
       height = 40,      # inches — very tall
       dpi = 300,
       units = "in")
