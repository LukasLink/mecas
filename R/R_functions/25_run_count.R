# run_count.R

run_count <- function(config_path, project_root_dir, cli_args){
  #-----------------------------------------------------------------------------
  # Run setup
  #-----------------------------------------------------------------------------
  
  config <- yaml::read_yaml(config_path)
  # Priority explicit overrides > Rmd params > config.yaml
  project_setup(
    project_root_dir = project_root_dir,
    config_path = config_path,
    use_old_suffix_construction = FALSE
  )

  #-----------------------------------------------------------------------------
  # Optional: handle use_only_these_controls
  #-----------------------------------------------------------------------------
  
  # If config$controls$use_only_these_controls contains entries, use them.
  # Otherwise, remove the object to preserve old downstream behavior.
  
  if (exists("use_only_these_controls_list") &&
      length(use_only_these_controls_list) == 0) {
    rm(use_only_these_controls_list)
  }
  #-----------------------------------------------------------------------------
  # Create Symlinks
  #-----------------------------------------------------------------------------
  manifest_output_path <- get_file_path(
    rds_output_folder,
    paste0("standardized_fastq_manifest_", file_info_suffix, ".tsv")
  )
  manifest <- prepare_fastq_inputs(
    fastq_dir = input_folder,
    fastq_name_table_file_path = fastq_name_table_xlsx,
    output_symlink_dir = fastq_symlinks_folder,
    manifest_output_path = manifest_output_path,
    strict_file_match = strict_file_match
  )
  log_info(paste("Using standardized FASTQ folder: ", fastq_symlinks_folder))
  log_info(paste("FASTQ manifest written to: ", manifest_output_path))
  
  #-----------------------------------------------------------------------------
  # Make slurm settings 
  #-----------------------------------------------------------------------------
  slurm_settings <- list(
    account    = slurm_account,
    qos        = slurm_qos,
    cpus       = slurm_cpus,
    mem        = slurm_mem,
    wall_time  = slurm_wall_time,
    partition  = slurm_partition,
    array      = slurm_array,
    email      = slurm_email
  )
  # slurm_settings <- list(
  #   account    = slurm_account,
  #   qos        = slurm_qos,
  #   cpus       = 1,
  #   mem        = "1g",
  #   wall_time  = "01:01:00",
  #   partition  = slurm_partition,
  #   array      = 1,
  #   email      = slurm_email
  # )
  #-----------------------------------------------------------------------------
  # Run QC filtering
  #-----------------------------------------------------------------------------
  if (opt$qc_filtering_run) {
    
    if (is.na(opt$qc_min_length) || is.null(opt$qc_min_length)) {
      opt$qc_min_length <- infer_qc_min_length(
        manifest = manifest,
        fastq_col = "symlink_file",
        n_lines = 10000
      )
    }
    
    manifest <- manifest %>%
      mutate(
        qc_filtered_paths = file.path(
          qc_filtered_folder,
          ensure_gz_suffix(symlink_file_basename)
        )
      )
    
    run_shell_step(
      step_name = "QC_filtering",
      script_path = file.path(project_root_dir, "shell", "QC_filtering.sh"),
      args = c(
        "--manifest", manifest_output_path,
        "--output-dir", qc_filtered_folder,
        "--min-qual", as.character(opt$qc_min_qual),
        "--qual-offset", as.character(opt$qc_qual_offset),
        "--min-length", as.character(opt$qc_min_length),
        "--use-modules", if (opt$use_modules) "true" else "false",
        "--seqtk-module", opt$seqtk_module
      ),
      slurm_settings = slurm_settings,
      machine = opt$machine,
      log_dir = log_folder
    )
    
  } else {
    
    manifest <- manifest %>%
      mutate(
        qc_filtered_paths = symlink_file
      )
  }
  #-----------------------------------------------------------------------------
  # Optional bcwithqc symlink creation goes here
  #-----------------------------------------------------------------------------
  if (read_counting == "bcwithqc") {
    
    manifest <- prepare_bcwithqc_inputs(
      manifest = manifest,
      output_symlink_dir = bcwithqc_symlinks_folder,
      overwrite_symlinks = TRUE,
      manifest_output_path = file.path(bcwithqc_symlinks_folder,
                                       "bcwithqc_symlink_manifest.tsv")
    )
    
    log_info(paste("Using bcwithqc FASTQ folder: ", bcwithqc_symlinks_folder))
  }
  #-----------------------------------------------------------------------------
  # Read counting
  #-----------------------------------------------------------------------------
  
  if (identical(opt$read_counting, "align_UMI_tools")) {
    
    # Make sure the manifest on disk contains qc_filtered_paths
    write.table(
      manifest,
      file = manifest_output_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    run_shell_step(
      step_name = "align_UMI_tools",
      script_path = file.path(project_root_dir, "shell", "align_UMI_tools.sh"),
      args = c(
        "--manifest", manifest_output_path,
        "--output-dir", opt$output_folder,
        "--star-index-folder", opt$star_index_folder,
        "--data-type", opt$data_type,
        "--umi-regex", opt$UMI_regex,
        "--threads", as.character(opt$slurm_cpus),
        "--use-modules", if (isTRUE(opt$use_modules)) "true" else "false",
        "--star-module", opt$star_module,
        "--samtools-module", opt$samtools_module,
        "--umi-tools-module", opt$umi_tools_module
      ),
      slurm_settings = slurm_settings,
      machine = opt$machine,
      log_dir = log_folder
    )
  }
  
  #-----------------------------------------------------------------------------
  # Get count_df_long 
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Optional Violin Plots
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # MAUDE
  #-----------------------------------------------------------------------------
  

  
  run_shell_step(
    step_name = "test_dummy",
    script_path = file.path(project_root_dir, "shell", "test_dummy.sh"),
    slurm_settings = slurm_settings,
    machine = "slurm", # LOREM -> change this to machine = machine in final setup
    log_dir = log_folder,
    extra_slurm_sbatch_lines = c("#SBATCH --constraint=avx512")
  )
}


