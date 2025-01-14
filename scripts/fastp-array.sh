#!/bin/bash
#SBATCH --job-name=fastp_array
#SBATCH --output=fastp_array_%A_%a.out
#SBATCH --error=fastp_array_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scebrian27@gmail.com
#SBATCH --array=1-5

#variables
dir=${1} #3.repetidas
input=${2} #3.list-fastp-output.txt

# Specify the path to the config file
config="/mnt/lustre/scratch/nlsas/home/csic/dbl/scc/${dir}/data/${input}"

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Extract the name for the current $SLURM_ARRAY_TASK_ID
name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

##sbatch run-fastp-array.sh 3.repetidas/data list-fastp-output.txt
module load cesga/2020
module load miniconda3
source activate culex
module load fastp

fastp \
    -i /mnt/lustre/scratch/nlsas/home/csic/dbl/scc/${dir}/data/${sample}_1.fq.gz -I /mnt/lustre/scratch/nlsas/home/csic/dbl/scc/${dir}/data/${sample}_2.fq.gz \
    -o /mnt/lustre/scratch/nlsas/home/csic/dbl/scc/${dir}/data/${name}.R1.fastp.fq.gz -O /mnt/lustre/scratch/nlsas/home/csic/dbl/scc/${dir}/data/${name}.R2.fastp.fq.gz \
    -h /mnt/lustre/scratch/nlsas/home/csic/dbl/scc/${dir}/fastp/${name}.fastp.html -j /mnt/lustre/scratch/nlsas/home/csic/dbl/scc/${dir}/fastp/${name}.fastp.json \
    --unpaired1 /mnt/lustre/scratch/nlsas/home/csic/dbl/scc/${dir}/data/${name}.R1.unpaired.fq.gz --unpaired2 /mnt/lustre/scratch/nlsas/home/csic/dbl/scc/${dir}/data/${name}.R2.unpaired.fq.gz \
    --failed_out /mnt/lustre/scratch/nlsas/home/csic/dbl/scc/${dir}/data/${name}.failed.fq.gz \
    --trim_poly_g \
    --trim_poly_x \
    --length_required 30 \
    --correction \
    --detect_adapter_for_pe \
    --adapter_fasta /mnt/lustre/scratch/nlsas/home/csic/dbl/scc/data/polyG.fasta


    #sbatch --array=1-5 fastp-array.sh 3.repetidas list-fastp-output.txt