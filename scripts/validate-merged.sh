#!/bin/bash
#this script is to check that all my merged bam files are sorted

# Define directories containing merged BAM files
DIRS=("1.noviembre/merged" "2.mayo/merged" "3.repetidas/merged" "4.enero/merged")  # Replace with your actual paths
OUTPUT_SUMMARY="all_merged_validation_summary.txt"
ERRORS_FILE="merged_validation_errors.txt"

# Clear any existing output files
> "$OUTPUT_SUMMARY"
> "$ERRORS_FILE"

# Loop through each directory
for dir in "${DIRS[@]}"; do
    echo "Scanning directory: $dir"

    # Loop through each merged BAM in the directory
    for bam in "$dir"/*.merged.bam; do
        [ -e "$bam" ] || continue  # Skip if no files matched

        base=$(basename "$bam" .merged.bam)
        summary_file="${dir}/${base}.summary.txt"

        echo "Validating: $bam"

        picard ValidateSamFile I="$bam" MODE=SUMMARY > "$summary_file" 2>&1

        # Append the full summary to the combined summary file
        {
            echo "===== Validation Summary: $bam ====="
            cat "$summary_file"
            echo -e "\n"
        } >> "$OUTPUT_SUMMARY"

        # Extract errors and append to error file with filenames
if grep -q '^ERROR' "$summary_file"; then
    echo "Errors in $bam:" >> "$ERRORS_FILE"
    grep '^ERROR' "$summary_file" >> "$ERRORS_FILE"
    echo "" >> "$ERRORS_FILE"
fi

    done
done

# Append errors to the end of the main summary
echo -e "\n\n===== ALL ERRORS EXTRACTED =====\n" >> "$OUTPUT_SUMMARY"
cat "$ERRORS_FILE" >> "$OUTPUT_SUMMARY"
