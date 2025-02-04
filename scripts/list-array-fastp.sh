#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --mem-per-cpu=500M
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=scebrian27@gmail.com
#SBATCH --output=list-fastp_%j.out
#SBATCH --error=list-fastp_%j.err

#this script will list all the files in a directory and write the output to a file
#specificly designed to list the files to perform fastp

# Directory containing the files
dir=${1} #4.enero/data
output=${2} #list-fastp.txt

# Construct the full path for the output file
output_path="/home/scamison/kernlab/scamison/${dir}/${output}"

# Create or clear the output file and add column headers
echo "id sample name" > $output_path

# Initialize a counter for the ID column
id=1

# Iterate over all _1.fq.gz files in the directory
for file in /home/scamison/kernlab/scamison/${dir}/*_1.fq.gz; do
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