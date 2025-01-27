Glossary to understand the name of the files:
pca12_0302pipiens123_1_f1n_BDQMC

* pca12 - the numbers correspond to the plotted principal components
* 0302 - refers to the variable filters applied to the vcf in order DP, GQ, MISSING - in this case, genotype depth was nnot filtered (0), GQ is at least 30 and missing data per loci is under 20% (0.2, thus the 2). 
* pipiens - species Cx pipiens
* _1_ - refers to the chromosome, in this examples, chromosome 1
* f1n - indicates that the vcf before DP, GQ and MISSING were already filter to keep only callable sites (without highly repeated areas and in which alternate allele had at least 1 read)
* BDQMC - indicates all the filters applied after f1n. B (biallelic), D (depth), Q (genotype quality), M (miising genotypes), C (minor allele count). All vcfs were aso filtered with a MAC (minor allele count) of 2 or more,  to keep only biallelic snps and also filtered loci depth above (mean DP * 2). Rest of filters are variable

## mod_*.pdf
In the original plots there were some obvious outliers, that were removed for better visualization.
PP1122 (seville, batch 2) and PP1001 (seville, batch 1).
PP1001 did not represent such an obvious outliers in the case of 2202 files corresponding to chromosomes 2 and 3, but it was in the case of chromosome 1, so I deleted it anyways to be constant through the visualization process.
