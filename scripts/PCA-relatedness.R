## Sonia Cebrián Camisón
## Proyecto Culex SNPs
# 11 mayo 2025
# Script to run PCA robust to related samples

# Load required Bioconductor libraries
library(SNPRelate)
library(GENESIS)
library(GWASTools)

#-------------------------------------------------------------
# Step 1: Convert PLINK files to GDS format (only run once)
#-------------------------------------------------------------
 snpgdsBED2GDS(
     bed.fn = "pipiens1234_all_PCA.bed",
     bim.fn = "pipiens1234_all_PCA.bim",
     fam.fn = "pipiens1234_all_PCA.fam",
     out.gdsfn = "pipiens1234_all_PCA.gds"
 )

#-------------------------------------------------------------
# Step 2: Calculate KING-robust kinship matrix using SNPRelate
#-------------------------------------------------------------
genofile <- snpgdsOpen("pipiens1234_all_PCA.gds")  # open GDS file for SNPRelate
king <- snpgdsIBDKING(genofile)
kinmat <- king$kinship

#------------------------------------------------------#
#Step 2.1: do average kinship cpefficient before closing the genofile
# Convert upper triangle of kinship matrix to a vector (excluding self-comparisons)
kin_values <- kinmat[upper.tri(kinmat)]

# Calculate mean and standard deviation
mean_kinship <- mean(kin_values, na.rm = TRUE)
sd_kinship <- sd(kin_values, na.rm = TRUE)

# Print results
cat("Mean kinship coefficient:", mean_kinship, "\n")
cat("Standard deviation:", sd_kinship, "\n")
# ------------------------------------------------------------

snpgdsClose(genofile)  # close the file to avoid conflicts
#showfile.gds(closeall = TRUE) #in case you delete the object before closing the file

#-------------------------------------------------------------
# Step 3: Prepare data for GENESIS and run PCAir
#-------------------------------------------------------------
gds <- GdsGenotypeReader("pipiens1234_all_PCA.gds")  # open GDS for GWASTools
genoData <- GenotypeData(gds)

#sample ids as matrix row and column names
sample_ids <- king$sample.id
dimnames(kinmat) <- list(sample_ids, sample_ids)

# Run PCAir using the kinship matrix
pca <- pcair(genoData, kinobj = kinmat, kin.thresh = 0.125)

#-------------------------------------------------------------
# Step 4: Plot PCA results
#-------------------------------------------------------------
plot(pca$vectors[,1], pca$vectors[,2], col = 'blue',
     xlab = "PC1", ylab = "PC2", main = "PCAir: PC1 vs PC2")

#-------------------------------------------------------------
# Step 5: export eigenvec and eigenval to plot with my better script
#-------------------------------------------------------------
eigenvec <- as.data.frame(pca$vectors)


eigenval <- pca$values
write.csv(eigenvec, "pcair_pipiens1234_eigenvec.csv", row.names = T)
write.csv(data.frame(PC = 1:length(eigenval), Eigenvalue = eigenval), 
          "pcair_pipiens1234_eigenval.csv", row.names = FALSE)


#-------------------------------------------------------------
# Step 6: to know excluded samples
#-------------------------------------------------------------

# Unrelated individuals (used for PCA)
unrelated <- pca$unrels

# Related individuals (projected onto PCs, but not used to compute them)
related <- pca$rels
related
