#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --mem-per-cpu=5G
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scebrian27@gmail.com
#SBATCH --output=filter_M_%j.out
#SBATCH --error=filter_M_%j.err

#this script will take vcf files, and apply some basic filters to the putput bcf files. vcf files must be indexed.
#the goal of the filters is to clean data to explore missing data and distributions
input=${1} #
output=${2}

module load miniconda3
conda activate vcf

echo 'F_MISSING>0.2'

bcftools filter -e 'F_MISSING>=0.2'  -Ob -o "${output}_M.bcf" "$input"   

# sbatch -n 1 --cpus-per-task=15 --mem-per-cpu=20G -t 12:00:00 --job-name=M_1 variant-filter-M.sh pipiens123_1_filter1n_qd.vcf.gz pipiens123_whatever
