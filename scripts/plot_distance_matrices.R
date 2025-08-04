##sonia cebrian camison
##proyecto culex snps
#01-30-2025
#this script is to plot genetic vs geographi distances in culex pipiens sequencing batches 1 2 and 3

# Load necessary libraries
install.packages("vegan")
library(vegan)  # For Mantel test
library(ggplot2)
library(reshape2)  # To reshape matrices into long format

## geogrtaphic matrix ----
# Load the geographic distance matrix
geo_dist <- read.csv("geo-distances-pipiens1234.csv", header=T, sep=",")
#transform in matrix
geo_matrix <- dcast(geo_dist, InputID ~ TargetID, value.var="Distance")
#fix format and colnames
rownames(geo_matrix) <- geo_matrix$InputID
geo_matrix$InputID <- NULL 
geo_matrix <- as.matrix(geo_matrix)

str(geo_matrix)
dim(geo_matrix)
#get lowe triangle indices
lower_tri <- lower.tri(geo_matrix)
str(lower_tri)

# Extract pairwise distances as vectors
geo_vec <- geo_matrix[lower_tri]   # Geographic distance

# Load the genetic distance matrix----
gen_dist <- read.table("pipiens1234_1-IBS.mdist")
f <- "pipiens1234_1-IBS.png"
tit <- "Genetic vs Geographic Distance (whole genome, MAF>3%)"

#setting first column as colnames (just when needed)----
#colnames(gen_dist) <- gen_dist[1,]
#gen_dist <- gen_dist[-1, ]
#head(gen_dist)
#setting first row as rownames
rownames(gen_dist) <- gen_dist[,1]
gen_dist <- gen_dist[,-1]

#turn into matrix----
gen_dist <- as.matrix(apply(gen_dist, 2, as.numeric))
rownames(gen_dist) <- colnames(gen_dist)  # Restore row names in gen_dist

#check
str(gen_dist)
dim(gen_dist)

identical(rownames(geo_matrix), colnames(geo_matrix))  # Should return TRUE
identical(rownames(gen_dist), colnames(gen_dist))  # Should return TRUE

#using lower tri generated in geo matrix part of the script
gen_vec <- gen_dist[lower_tri]     # Genetic distance

df <- data.frame(GeographicDistance = geo_vec, GeneticDistance = gen_vec)

fig <- ggplot(df, aes(x = GeographicDistance, y = GeneticDistance)) +
  geom_point(alpha=0.6, color="blue") +  # Scatterplot points
  geom_smooth(method="lm", color="red", se=FALSE) +  # Regression line
  labs(title=tit,
       x="Geographic Distance (meters)",
       y="Genetic Distance (1 - IBS)") +
  
  theme_bw()

fig

ggsave(f, device = "png", units = "in", plot=fig, width=5, height=4)

#