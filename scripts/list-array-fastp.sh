#!/bin/bash
#SBATCH --output=output-list-array.txt
#SBATCH --mem-per-cpu=500M
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=sonia.cebrian@ebd.csic.es

#this script will list all the files in a directory and write the output to a file
#specificly designed to list the files to perform fastp

# Directory containing the files
dir=${1}
output=${2}

# Construct the full path for the output file
output_path="/mnt/lustre/scratch/nlsas/home/csic/dbl/scc/${dir}/${output}"

# Create or clear the output file and add column headers
echo "id sample name" > $output_path

# Initialize a counter for the ID column
id=1

# Iterate over all _1.fq.gz files in the directory
for file in /mnt/lustre/scratch/nlsas/home/csic/dbl/scc/${dir}/*_1.fq.gz; do
    # Extract the sample name
    sample=$(basename "$file" "_1.fq.gz")
    # Extract the name (first 6 characters + "." + L code)
    name="${sample:0:6}.${sample: -2}"
    # Write the information to the output file
    echo "${id} ${sample} ${name}" >> $output_path

    # Increment the ID counter
    id=$((id + 1))
done


#sbatch list-array-fastp.sh 3.repetidas/data list-fastp-output.txt