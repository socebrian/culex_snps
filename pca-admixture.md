# Proyec: SNPS Culex - ECOGEN
# pipieline for first analyses: PCA, local PCA, Admixture
## sonia cebrian camison scebrian27@gmail.com/sonia.cebrian@ebd.csic.es
### 13/1/2025


### 8. Linkage- desequilibrium and conventional PCA (culex env on talapas)

I am working with PLINK and sometimes it was gwtting confused by any FORMAT field that was not GT (it returned Error saying GT can not contain Floating values). So first thing I get rid of all FORMAT fields other than GT
```bash
# first let only GT in format so plink dont freak out
bcftools annotate -x ^FORMAT/GT pipiens123_1_1nBDQF_M.bcf  -Oz -o pipiens123_1_filteredGT.vcf.gz
tabix -p vcf stripped.vcf.gz
```
Then with this vcfs I run  <ld-prunning2.sh>. This script prune out SNPs under linkage desequilibrium, recode new vcf without the prunned out sites and then use this new vcf to perform a PCA.

#### ld-prunning2.sh
```bash

#set variables
input=${1} # file.bcf
prefix=${2} #pipiens123_ch1

# perform linkage pruning - i.e. identify prune sites 
echo "performing prunning on ${input}"

plink --bcf  $input --double-id \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.5 --out $prefix

echo "prunning complete"
echo "recoding vcf"

plink --bcf $input --set-missing-var-ids @:# \
  --extract ${prefix}.prune.in \
  --recode vcf-iid \
  --out ${prefix}_pruned

echo "performing pca from recoded vcf"

plink --vcf ${prefix}_pruned.vcf \
      --double-id \
      --set-missing-var-ids @:# \
      --make-bed \
      --pca \
      --out ${prefix}_PCA

# sbatch -n 1 --cpus-per-task=15 --mem-per-cpu=20G -t 12:00:00 --job-name=prune_1 prune-linkage.sh pipiens123_1_filter1n_qd_M.bcf pipienschr1
```
#### 

#### How many SNPS after prunning?
* chrom 1: 2167 removed, leaving 10949.
* chrom 2: 2479 removed, leaving 12893.
* chrom 3: 2919 removed, leaving 14790

This generates two files per chromosome:
pipiens123_ch1_PCA.eigenval           
pipiens123_ch1_PCA.eigenvec

I downloaded them, added info on localities, provinces and batch to the eigenvectors data and plotted the results on R. I got a problem with my Rstudio and it closed the script <pipiens123_pca_filtrosmal.R> without saving, so now it is corrupt. More or less it looked like that:

```r
library(ggplot2)

# read eigenvalues
eigenvals <- read.table("pca/data1.plink.pca.eigenval",
                        col.names = c("eigenval"))

# read eigenvectors
eigenvecs <- read.table("pca/data1.plink.pca.eigenvec",
                        header = TRUE)[,-1]

# convert eigenvalues to percentage
percents <- eigenvals$eigenval/sum(eigenvals$eigenval)*100

# scree-plot of eigenvalues
scree <- ggplot() +
  geom_line(aes(x = as.numeric(row.names(eigenvals)), y = percents)) +
  geom_point(aes(x = as.numeric(row.names(eigenvals)), y = percents),
    fill = "grey", shape = 21, size = 2.5) +
  xlab(paste0("Number of Principal Component")) +
  ylab(paste0("Percentage of variance Explained")) +
  scale_x_continuous(n.breaks = 10) +
  theme_bw()

ggsave(filename = "pca/data1.scree_plot.pdf", plot = scree)

#### this first part is from https://github.com/Enricobazzi/EBD-PopGenIntro/blob/main/practica2/practica2.md
###it was adapted to plot my data

pc12_prov <- ggplot() +
  geom_point(
    data = eigenvecs_1,
    aes(x = v1, y = v2, fill = prov, shape = batch),
    size  = 2.5,
    color = "black"   # <--- Give an explicit outline color
  ) +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_fill_brewer(palette = "Set1") +
  xlab(paste0("PC1 - ", round(percents_1[1], 2), "%")) +
  ylab(paste0("PC2 - ", round(percents_1[2], 2), "%")) +
  theme_bw() +
  # Override the legend so the shapes in the *fill* legend also appear filled
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, color = "black")),
    shape = guide_legend(override.aes = list(color = "black"))
  )

ggsave(filename = "R/pipiens123_pc12_prov_plot.pdf", plot = pc12_prov)
```

### 9. local PCA (vcf env on poppy)
NOTE: I used files before LD prunning instead of the prunned ones.

I followed <https://github.com/petrelharp/local_pca/tree/master/templated> to run <run_lostruct.R> and <summarize_run-Rmd>.
(Some packages may be neccesary to install following the error, I did so inside an R session)

As indicated, i placed my vcf files and indexes on /data
pipiens123_1_1nBDQF_M.bcf      
pipiens123_2_1nBDQF_M.bcf      
pipiens123_3_1nBDQF_M.bcf
pipiens123_1_1nBDQF_M.bcf.csi  
pipiens123_2_1nBDQF_M.bcf.csi  
pipiens123_3_1nBDQF_M.bcf.csi

Then, from my shell
```bash
# I used a .txt file instead of .tsv for population, is that okay?
Rscript run_lostruct.R -i ./data \
                       -I populations-mpileup.txt \ 
                       -t snp \
                       -s 1000 

```

And on the r session again
```r
templater::render_template(
  "summarize_run.Rmd",
  output = "lostruct_results/type_snp_size_1000_weights_none_jobid_488519/run_summary.html",
  change.rootdir = TRUE
)
```
This generates a summary that I renamed as <pipiens123_localPCA_results.html> 




### questions to that point
* use molestus from other studies in PCA to see if any clustering is due to biotypes. Where to introduce them? In the snp calling already? Later? Can i do a different VCF and then merge it with the one y already got?
* PCA (plink) vs Local PCA ¿cual hacer y por qué? 
* local PCA se hace tambien con biallelic positions?
* local PCA interpretation?
* any other plots to do with the results?

### 10. ADMIXTURE /STRUCTURE analyses
en principio structure es mas robusto cuando no hay una clara estrcutra poblacional pero es muy sensible a unaven sampling por localidades asi que tener en cuenta. Structure (ESPAÑA)

### 11. Calculate Fst and Ne
 cuantificar diferentiation Fst usando VCFtools y tamaños poblacionales efectivos (ESPAÑA)

## 12. Landscape genetics analyses
I have raster images with monthly values for 2022 of:
- Tmin
- Tmax
- rainfall
- soil moisture
- vapor pressure
- evapotranspiration
- over surface water
- NDVI

* landscape analysis. alternativa a circuit theory?

