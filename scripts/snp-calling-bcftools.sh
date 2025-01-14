#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --mem-per-cpu=5G
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scebrian27@gmail.com
#SBATCH --output=snpcall_%j.out
#SBATCH --error=snpcall_%j.err

#HAY QUE MODIFICARLO PARA QUE EXPORTE EN BCF
#this script runs bcftools mpileup and call over all the bam files in a directory 
# when they are divided by regions and then merge the vcf files

module load miniconda3
conda activate snpbcf
module load samtools

thr=${SLURM_CPUS_PER_TASK}
chromosome=${1} #1, 2 or 3

# Paths and variables
dir="/home/scamison/kernlab/scamison/SNPcalling"
reference="/home/scamison/kernlab/scamison/data/idCulPipi.chrom.fa"
population_file="/home/scamison/kernlab/scamison/SNPcalling/populations-mpileup.txt"  # File with population information
 
# Chromosome BAM lists. list must contain the COMPLETE path to the bam files
list="${dir}/bam-list.txt" 	

# Run bcftools mpileup and call
bcftools mpileup -Ou  --threads $thr -f $reference -q 20 -Q 15 -b $list -r "$chromosome" --annotate "FORMAT/AD,FORMAT/DP,INFO/AD" | \
bcftools call -m  --threads $thr -G $population_file -Ob --format-fields "GQ,GP" -o ${dir}/pipiens123_${chromosome}_pop.bcf


# sbatch -n 1 --cpus-per-task=10 --mem-per-cpu=20G -t 12-00:00:00 --job-name=SNPcall2 snp-calling-bcftools.sh 2