#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --mem-per-cpu=500M
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=scebrian27@gmail.com
#SBATCH --output=mark_dup_%A_%a.out
#SBATCH --error=mark_dup_%A_%a.err
#SBATCH --array=1-50

#this script takes a list of bam files and remove duplicates

# Load any necessary modules
module load miniconda3
conda activate vcf
module load racs-eb/1
module load samtools
module load picard

#variables
dir=${1} #3.repetidas
input=${2} #list of bam files

# Specify the path to the config file
sample_list="/home/scamison/kernlab/scamison/${dir}/bams/${input}"
bam_dir="/home/scamison/kernlab/scamison/${dir}/bams"
out_dir="/home/scamison/kernlab/scamison/${dir}/mark-duplicates"

# Read the sample ID for the current array task
sample=$(awk -v taskid="$SLURM_ARRAY_TASK_ID" '$1 == taskid {print $2}' "$sample_list")
in_bam="$bam_dir/${sample}.PR.sorted.rg.bam"
# Define output filenames
output_bam="$out_dir/${sample}.md.bam"
output_metrics="$out_dir/${sample}_marked_dup_metrics.txt"


# Debug: Print the sample variable
echo "Sample: $sample"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I=$in_bam \
    O=$output_bam \
    M=$output_metrics \
    REMOVE_DUPLICATES=true

##sbatch --array=1-5 -n 1 --cpus-per-task=10 --mem-per-cpu=20 -t 12:00:00 --job-name=3rep.markdup  mark-duplicates-array.sh 3.repetidas list-mark-duplicates.txt

