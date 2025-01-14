#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --mem-per-cpu=5G
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sonia.cebrian@ebd.csic.es
#SBATCH --output=qualimap_%j.out
#SBATCH --error=qualimap_%j.err

#this script runs qualimap over all *merged.bam files in a directory and the runs multibam qc over all the qualimap outputs
module load miniconda3
conda activate culex

dir=${1} #3.repetidas

# Paths and variables
input_dir="/home/scamison/kernlab/scamison/${dir}/merged"	
output_dir="/home/scamison/kernlab/scamison/${dir}/qualimaps"

## Loop through each merged BAM file in the input directory
for bam_file in "$input_dir"/*merged.bam; do
    # Check if the BAM file exists to avoid errors
    if [[ -f "$bam_file" ]]; then
        # Extract sample ID from filename
        sample_id=$(basename "$bam_file" .merged.bam)

        echo "Running qualimap for sample: ${sample_id}"
        qualimap bamqc -bam "$bam_file" -outdir "${output_dir}/${sample_id}" -outfile "${sample_id}.merged.bamqc.pdf"  -outformat PDF:HTML
    else
        echo "No BAM files found in ${input_dir}"
    fi
done

# sbatch --array=1-5 -n 1 --cpus-per-task=1 --mem-per-cpu=2G -t 02:00:00 --job-name=3rep.qualimap  qualimap.sh 3.repetidas
