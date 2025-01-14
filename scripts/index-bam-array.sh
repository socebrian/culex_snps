#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --mem-per-cpu=500M
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=scebrian27@gmail.com
#SBATCH --output=bam_index_%A_%a.out
#SBATCH --error=bam_index_%A_%a.err
#SBATCH --array=1-50

#this script takes a list of bam files and index them giving info about the reads, lanes and batch 

# Load any necessary modules
module load miniconda3
conda activate culex
module load java/17
module load picard

#variables
dir=${1} #3.repetidas
input=${2} #list-index-output.txt

# Specify the path to the config file
config="/home/scamison/kernlab/scamison/${dir}/bams/${input}"


IFS=$'\t' read -r -a samples <<< $(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$config")
# Extract parameters
id=${samples[1]}
lane=${samples[2]}
r=${samples[3]}
ext=${samples[4]}

# Debug information
echo "Sample List: $config"
echo "Batch: $dir"
echo "ID: $id"
echo "Lane: $lane"
echo "R: $r"
echo "Ext: $ext"
echo "Executing: sbatch index-bam.sh $id $lane $r $dir"

# File paths
input_bam="/home/scamison/kernlab/scamison/${dir}/bams/${id}.${lane}.${r}.bam"
sorted_bam="/home/scamison/kernlab/scamison/${dir}/bams/${id}.${lane}.${r}.sorted.bam"
sorted_rg_bam="/home/scamison/kernlab/scamison/${dir}/bams/${id}.${lane}.${r}.sorted.rg.bam"

# Check if input BAM file exists
if [[ ! -f "$input_bam" ]]; then
    echo "Error: Input BAM file not found at $input_bam"
    exit 1
fi

# Sort BAM file
echo "Sorting BAM file: $input_bam"
samtools sort "$input_bam" -o "$sorted_bam"

# Add or replace read group
#id flag in picard should be only the mosquito ID so different lanes can be merged afterwards wihout problems messing with the metadata
echo "Adding or replacing read groups in BAM file: $sorted_bam"
picard AddOrReplaceReadGroups \
    I="$sorted_bam" \
    O="$sorted_rg_bam" \
    RGID="${id}" \ 
    RGLB="$dir" \
    RGPL=Illumina \
    RGPU="${id}.${lane}.${r}" \
    RGSM="$id" \
    VALIDATION_STRINGENCY=SILENT

# Index the sorted BAM file with read groups
echo "Indexing BAM file: $sorted_rg_bam"
samtools index "$sorted_rg_bam"

##sbatch --array=1-5 -n 1 --cpus-per-task=1 --mem-per-cpu=2G -t 02:00:00 --job-name=3rep.index  index-bam.sh 3.repetidas list-index-output.txt
#2gb para un ba son muchos, con 1g sobra seguramente
#esto esta hecho para cuando no era array cuidao