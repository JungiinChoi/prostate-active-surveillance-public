#!/usr/bin/env bash

rscript_dir=" /users/zwang3/PAS_clean/cleaned_code"
rscript_log_file_name="${HOME}/generated-files-cv/${job_name}_rscript_log.txt"
mkdir -p "${HOME}/generated-files-cv" # Make a directory if non-existent
touch "$rscript_log_file_name" # Create the file if non-existent
module load conda_R
echo "Starting ${job_name}..."
Rscript "${rscript_dir}/cv-to-run-single.R" 