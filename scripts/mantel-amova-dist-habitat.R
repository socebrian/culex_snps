##sonia cebrian camison
##proyecto culex snps
#21-07-2025
#this script is to do numeric test of the relation between genetic distances, geographic distances and habitat type
#working directory: SNPCulex/gen-geo-distances
# Load necessary libraries

library(vegan)  # For Mantel test
library(adegenet) #for AMOVA?
library(ggplot2)
library(reshape2)
library(pegas) #for AMOVA
# Load the  distances files
geo_dist <- read.csv("geo-distances-pipiens1234.csv", header=T, sep=",")
gen_dist <- read.table("pipiens1234_1-IBS.mdist")
metadata <- read.table("metadata-pipiens-217samples.txt", sep="\t", header=T)
habitat <- metadata[,c(1,4)]

## geogrtaphic matrix ----
#transform in matrix
geo_matrix <- dcast(geo_dist, InputID ~ TargetID, value.var="Distance")
#fix format and colnames
rownames(geo_matrix) <- geo_matrix$InputID
geo_matrix$InputID <- NULL 
geo_matrix <- as.matrix(geo_matrix)

str(geo_matrix)
dim(geo_matrix)

#genetic distance matrix----
#turn into matrix
gen_dist <- as.matrix(apply(gen_dist, 2, as.numeric))
rownames(gen_dist) <- colnames(gen_dist)  # Restore row names in gen_dist
#check
str(gen_dist)
dim(gen_dist)

identical(rownames(geo_matrix), colnames(geo_matrix))  # Should return TRUE
identical(rownames(gen_dist), colnames(gen_dist))  # Should return TRUE

#mantel test: gen vs geo
mantel_result <- mantel(gen_dist, geo_matrix, method = "pearson", permutations = 999)
dim(gen_dist)
dim(geo_matrix)
mantel_result
################ mantel partial y AMOVA para tipo de habitat ################
#pone mismo orden que en dist geneticas
orden <- rownames(as.matrix(gen_dist))
habitat <- habitat[match(orden, habitat$ID), ]

# Codifica el tipo de hábitat como numérico
# urbano = 0, natural = 1 (por ejemplo)
hab_num <- as.numeric(as.factor(habitat$habitat))  # urbano = 1, natural = 2
str(hab_num)

# Crea matriz de diferencias binarias (1 si distinto tipo de hábitat, 0 si igual)
habitat_dist <- as.dist(outer(hab_num, hab_num, FUN = function(x, y) as.numeric(x != y)))
attr(habitat_dist, "Size")        # número de muestras

gen_dist <- as.dist(gen_dist)
geo_matrix <- as.dist(geo_matrix)
attr(gen_dist, "Size")
attr(geo_matrix, "Size")
str(geo_matrix)
# Partial Mantel: ¿diferencias genéticas explicadas por hábitat, controlando por distancia geográfica?
partial_mantel <- mantel.partial(
  gen_dist,
  habitat_dist,
  geo_matrix,
  method = "pearson",
  permutations = 999
)

print(partial_mantel)



#############################################################
# VIOLINPLOT habitat vs dist genetica
##############################################################
str(gen_dist)
str(gen_mat)
gen_mat <- as.matrix(gen_dist) # convertir a matriz
dim(gen_mat)
dist_values <- gen_mat[lower.tri(gen_mat)]

tipo_habitat <- habitat$habitat
# Crear una matriz lógica: TRUE si las muestras comparten hábitat, FALSE si no
same_habitat <- outer(tipo_habitat, tipo_habitat, FUN = "==")

# Extraer los valores del triángulo inferior de la matriz (evita repeticiones y diagonal)
same_values <- same_habitat[lower.tri(same_habitat)]

# Crear un vector de grupos: "Mismo hábitat" vs "Distinto hábitat"
habitat_group <- ifelse(same_values, "Same habitat", "Different habitat")

habitat_group


# Crear dataframe para graficar
df <- data.frame(
  distancia = dist_values,
  grupo = habitat_group
)

# Graficar
fig <- ggplot(df, aes(x = grupo, y = distancia, fill = grupo)) +
  geom_violin(trim = FALSE, alpha = 0.5) +                                # violín
  geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.6) +           # caja
  stat_summary(fun = "mean", geom = "point", shape = 20, size = 3, 
               color = "black") +                                        # media anotada
  labs(x = "Habitat comparison", 
       y = "Genetic distance",
       title = " ") +
  scale_fill_manual(values = c("tomato", "lightblue")) +
  theme_minimal() +
  theme(legend.position = "none")
fig

ggsave
ggsave("gen-dist-habitat.png", device = "png", units = "in", plot=fig, width=5, height=5)


#########################################
# AMOVA
############################################
str(gen_dist) #need a distance matrix
str(metadata)
metadata <- metadata[match(orden, metadata$ID), ]
metadata$location <- as.factor(metadata$location)
metadata$habitat <- as.factor(metadata$habitat)
metadata$admixture <- as.factor(metadata$admixture)

str(metadata)

amova_result <- amova(gen_dist ~ habitat, data = metadata, nperm = 999)
amova_result

amova_admix <- amova(gen_dist ~ admixture, data = metadata, nperm = 999)
amova_admix
######################################
#AMOVA de todos los factores a 
###############################

#primero comprobamos correlaciones entre variables

amova_hab_loc <- amova(gen_dist ~ habitat/location, data= metadata, nperm=999)
amova_hab_loc
amova_3fac <- amova(gen_dist ~ habitat/admixture/location, data= metadata, nperm=999)
amova_3fac
amova_eco_loc <- amova(gen_dist ~ location/admixture, data= metadata, nperm=999)
amova_eco_loc


#ade4 

labels(gen_dist)[1:10]
metadata$ID[1:10]
str(gen_dist)
