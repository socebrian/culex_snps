#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --mem-per-cpu=500M
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=scebrian27@gmail.com
#SBATCH --output=list_index_%j.out
#SBATCH --error=list_index_%j.err


dir=${1} # 3.repetidas/bams
output=${2} # name of output file

# Directory containing BAM files
bam_dir="/home/scamison/kernlab/scamison/${dir}"
output_path="/home/scamison/kernlab/scamison/${dir}/${output}"

# Write the header line to the output file
echo -e "taskid\tsample" > $output_path

# Initialize a counter for the ID column
taskid=1

# Extract unique sample IDs (PP****) from the filenames and write to the sample list
ls "$bam_dir"/*.bam | awk -F'/' '{print $NF}' | awk -F'.' '{print $1}' | awk -F'-' '{print $1}' | sort | uniq | while read sample; do
    echo -e "${taskid}\t${sample}" >> $output_path
    taskid=$((taskid + 1))
done

echo "Sample list created: $output_path"

#sbatch list-merge-bams.sh 3.repetidas/bams list-merge-bams.txt