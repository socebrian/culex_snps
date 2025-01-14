#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --output=merge_bams_%A_%a.out
#SBATCH --error=merge_bams_%A_%a.err
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-5
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scebrian27@gmail.com


#this scripts merges all the sorted bam correspondign to the same sample but different lanes
module load miniconda3
conda activate culex
module load samtools

# Directory containing BAM files
dir=${1} #3.repetidas
input=${2} #3.list-to-merge-bams.txt

sample_list="/home/scamison/kernlab/scamison/${dir}/mark-duplicates/${input}"
bam_dir="/home/scamison/kernlab/scamison/${dir}/mark-duplicates"
out_dir="/home/scamison/kernlab/scamison/${dir}/merged"

# Read the sample ID for the current array task
sample=$(awk -v taskid="$SLURM_ARRAY_TASK_ID" '$1 == taskid {print $2}' "$sample_list")

# Debug: Print the sample variable
echo "Sample: $sample"

# Find all BAM files for the current sample
bam_files=$(ls "$bam_dir"/"${sample}"*.md.bam)

# Check if BAM files exist
if [ -z "$bam_files" ]; then
    echo "Error: Input BAM files not found for sample ${sample} at ${bam_dir}"
    exit 1
fi

# Merge the BAM files
output_bam="$out_dir/${sample}.merged.bam"
samtools merge "$output_bam" $bam_files

echo "Merged BAM file created: $output_bam"

##sbatch --array=1-5 -n 1 --cpus-per-task=1 --mem-per-cpu=2G -t 02:00:00 --job-name=3rep.merge  merge-bams-array.sh 3.repetidas  list-merge-bams.txt
