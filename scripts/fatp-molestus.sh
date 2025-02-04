#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --mem-per-cpu=20G
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=scebrian27@gmail.com
#SBATCH --output=fastpc_%A_%a.out
#SBATCH --error=fastp_%A_%a.err
#SBATCH --array=1-50

#this script will list all the files in a directory and write the output to a file
#specificly designed to list the files to perform fastp
#loasd modules
module load miniconda3
conda activate fastp

#variables
dir=${1} #3.repetidas
input=${2} #3.list-fastp-output.txt

# Put bash in strict mode
set -e
set -o pipefail

# Specify the path to the config file
config="/home/scamison/kernlab/scamison/${dir}/data/${input}"

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Extract the name for the current $SLURM_ARRAY_TASK_ID

# Run fastp
fastp \
    -i /home/scamison/kernlab/scamison/${dir}/data/${sample}_1.fq.gz -I /home/scamison/kernlab/scamison/${dir}/data/${sample}_2.fq.gz \
    -o /home/scamison/kernlab/scamison/${dir}/data/${sample}.R1.fastp.fq.gz -O /home/scamison/kernlab/scamison/${dir}/data/${sample}.R2.fastp.fq.gz \
    -h /home/scamison/kernlab/scamison/${dir}/fastp/${sample}.fastp.html -j /home/scamison/kernlab/scamison/${dir}/fastp/${sample}.fastp.json \
    --unpaired1 /home/scamison/kernlab/scamison/${dir}/data/${sample}.R1.unpaired.fq.gz --unpaired2 /home/scamison/kernlab/scamison/${dir}/data/${sample}.R2.unpaired.fq.gz \
    --failed_out /home/scamison/kernlab/scamison/${dir}/data/${sample}.failed.fq.gz \
    --trim_poly_g \
    --trim_poly_x \
    --length_required 30 \
    --correction \
    --detect_adapter_for_pe #\
   # --adapter_fasta /mnt/lustre/scratch/nlsas/home/csic/dbl/scc/data/polyG.fasta


    #sbatch --array=1-24 -n 1 --cpus-per-task=10 --mem-per-cpu=20G -t 12:00:00 fastp-molestus.sh molestus list-fastp-output.txt