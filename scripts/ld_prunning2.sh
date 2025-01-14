#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --mem-per-cpu=5G
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scebrian27@gmail.com
#SBATCH --output=prune_ld_%j.out
#SBATCH --error=prune_ld_%j.err

#loasd modules
module load miniconda3
conda activate culex
module load plink

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

