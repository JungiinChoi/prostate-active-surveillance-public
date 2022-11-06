node_name='shared.q@compute-*'
mem_gb_free=10
n_core=1

for to_mask in 11 2 3 4 5 6 7 8 9 10; do
    job_name="cv${to_mask}"
    mkdir -p "${HOME}/generated-files-cv" # Make a directory if non-existent
    log_file_name="${HOME}/generated-files-cv/${job_name}_qsub_log.txt"
    qsub \
      -q "${node_name}" `# If you need a specific node` \
      -l mem_free="${mem_gb_free}G",h_vmem="${mem_gb_free}G" `# Specify memory requirement` \
      -pe local $n_core `# Parallel environment for multi-threading` \
      -N $job_name `# Give a human-readable name to the submitted job so that you can find it later` \
      -o $log_file_name `# Direct output messages` \
      -e $log_file_name `# Direct errors` \
      -m e -M zwang238@jh.edu `# Send an email when the job completes or aborts` \
      -v to_mask=${to_mask},job_name=${job_name} `# Assign variables to be passed to the bash script` \
      submit_single_cv.sh
done