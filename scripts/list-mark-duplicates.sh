#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --mem-per-cpu=500M
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=scebrian27@gmail.com
#SBATCH --output=list_markd_%j.out
#SBATCH --error=list_markd_%j.err

# Directory and output file passed as arguments
dir=${1} # 1.noviembre/bams
output=${2} # list-mark-duplicates.txt

# Directory containing BAM files
bam_dir="/home/scamison/kernlab/scamison/${dir}"
output_path="/home/scamison/kernlab/scamison/${dir}/${output}"

# Write the header line to the output file
echo -e "taskid\tsample" > $output_path

# Initialize a counter for the ID column
taskid=1

# Extract unique BAM file paths and write to the sample list
ls "$bam_dir"/*.PR.sorted.rg.bam 2>/dev/null | awk -F'/' '{print $NF}' | sed 's/\.PR\.sorted\.rg\.bam//' | sort | uniq | while read sample; do
    echo -e "${taskid}\t${sample}" >> $output_path
    taskid=$((taskid + 1))
done

echo "Sample list created: $output_path"

# Example usage:
# sbatch list-mark-duplicates.sh 3.repetidas/bams list-mark-duplicates.txt
