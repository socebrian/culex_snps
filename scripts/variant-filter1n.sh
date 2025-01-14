#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --mem-per-cpu=5G
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scebrian27@gmail.com
#SBATCH --output=filter1n_%j.out
#SBATCH --error=filter1n_%j.err

#this script will take vcf files, filter repetitions and then filter also min allele count for the minor allele

module load miniconda3
conda activate vcf

# Get the chromosome number from the command line argument
chromosome=${1}
bcf=${2} #pipiens123

# Define file paths based on the chromosome number
repeat_bed="variant_filtering/repeats/s${chromosome}.cxpipiens.repeats.bed"
input_bcf="SNPcalling/${bcf}_${chromosome}_pop.bcf"
repeat_masked_vcf="variant_filtering/${bcf}_${chromosome}_repeatmasked.vcf.gz"
filtered_vcf="variant_filtering/${bcf}_${chromosome}_filter1n.vcf.gz"

# Print the file paths
echo "Repeat BED file: ${repeat_bed}"
echo "Input BCF file: ${input_bcf}"
echo "Repeat masked VCF file: ${repeat_masked_vcf}"
echo "Filtered VCF file: ${filtered_vcf}"

# Run bcftools commands
bcftools filter -M ^$repeat_bed --soft-filter repeat_masked_super${chromosome} --threads 30 -Oz -o $repeat_masked_vcf $input_bcf
bcftools view -f PASS,. --min-ac 1 --threads 30 -Oz -o $filtered_vcf $repeat_masked_vcf
