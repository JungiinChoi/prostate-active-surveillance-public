node_name='shared.q@compute-*'
mem_gb_free=10
n_core=1
workdir="/users/jchoi/PAS"
mri_role="MRI_ST"
n_cv_fold=${1:-1} 

for ((to_mask=1;to_mask<=$n_cv_fold;to_mask++)) do
    job_name="cv${to_mask}_${n_cv_fold}_${mri_role}"
    #mkdir -p "${workdir}/generated-files-sh" # Make a directory if non-existent
    log_file_name="${workdir}/logs/${job_name}_qsub_log.txt"
    qsub \
      -q "${node_name}" `# If you need a specific node` \
      -l mem_free="${mem_gb_free}G",h_vmem="${mem_gb_free}G" `# Specify memory requirement` \
      -pe local $n_core `# Parallel environment for multi-threading` \
      -N $job_name `# Give a human-readable name to the submitted job so that you can find it later` \
      -o $log_file_name `# Direct output messages` \
      -e $log_file_name `# Direct errors` \
      # -m e -M jchoi177@jh.edu `# Send an email when the job completes or aborts` \
      -v K=${n_cv_fold},to_mask=${to_mask},mri_role=${mri_role},J=${J},workdir=${workdir},job_name=${job_name} `# Assign variables to be passed to the bash script` \
      submit_single_cv.sh
done
#sh submit_simulation.sh: submit for non-cv (K = 1)
#sh submit_simulation.sh K: submit for K fold cv 
