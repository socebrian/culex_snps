#!/bin/bash
#SBATCH --account=kernlab
#SBATCH --partition=kern
#SBATCH --mem-per-cpu=500M
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=scebrian27@gmail.com
#SBATCH --output=check_files_%j.out
#SBATCH --error=check_files_%j.err
#SBATCH --array=1-50

# Directory containing .err and .out files
log_dir=${1}
report_file=${2}

# Initialize report file
echo "Job Error Report" > "$report_file"
echo "================" >> "$report_file"

# Check .err files for common errors
echo "Checking .err files for errors..." >> "$report_file"
for err_file in "$log_dir"/*.err; do
    if [[ -s "$err_file" ]] && grep -i -E "error|fail|abort|segmentation fault" "$err_file"; then
        job_id=$(basename "$err_file" | sed 's/.*-\([0-9]*\)\.err/\1/')
        echo "Error found in Job ID $job_id ($err_file):" >> "$report_file"
        grep -i -E "error|fail|abort|segmentation fault" "$err_file" >> "$report_file"
        echo "----------------------" >> "$report_file"
    fi
done

# Check .out files for exit code 1 (indicating error)
echo "Checking .out files for errors..." >> "$report_file"
for out_file in "$log_dir"/*.out; do
    if [[ -s "$out_file" ]] && grep -i "exit code 1" "$out_file"; then
        job_id=$(basename "$out_file" | sed 's/.*-\([0-9]*\)\.out/\1/')
        echo "Exit code 1 (Error) in Job ID $job_id ($out_file):" >> "$report_file"
        grep -i "exit code 1" "$out_file" >> "$report_file"
        echo "----------------------" >> "$report_file"
    fi
done

# Final message
if [ -s "$report_file" ]; then
    echo "Issues found. Check $report_file for details."
else
    echo "No issues found in .err or .out files."
fi

# Example of use
## sbatch -n 1 --mem-per-cpu=1G -t 00:30:00 --job-name=check-files check-files.sh path_to_logs report_file.txt
