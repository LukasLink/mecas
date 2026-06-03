# run_count.R

run_count <- function(config_path, project_root_dir, cli_args){
  #-----------------------------------------------------------------------------
  # Run setup
  #-----------------------------------------------------------------------------
  message("Begining Project setup...")
  
  config <- yaml::read_yaml(config_path)
  # Priority explicit overrides > Rmd params > config.yaml
  project_setup(
    project_root_dir = project_root_dir,
    config_path = config_path,
    setup_mode = "count",
    use_old_suffix_construction = FALSE
  )
  run_prepare_inputs_stage <- should_run_stage(opt$start_with, "beginning")
  run_read_counting_stage  <- should_run_stage(opt$start_with, "read_counting")
  run_maude_stage          <- should_run_stage(opt$start_with, "MAUDE_analysis")
  run_plots_stage          <- should_run_stage(opt$start_with, "generate_plots")
  
  logger::log_info("Finished Project setup.")
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
  
  if (run_prepare_inputs_stage) {
    
    logger::log_info("Creating symlinks of input files...")
    
    manifest <- prepare_fastq_inputs(
      fastq_dir = input_folder,
      fastq_name_table_file_path = fastq_name_table_xlsx,
      output_symlink_dir = fastq_symlinks_folder,
      manifest_output_path = manifest_output_path,
      strict_file_match = strict_file_match
    )
    
    log_info("Finished Creating symlinks of input files.")
    
  } else {
    
    log_info("Skipping input preparation because start_with is: {opt$start_with}")
    
    manifest <- load_existing_tsv_or_stop(
      manifest_output_path,
      "standardized FASTQ manifest"
    )
  }
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
    email      = slurm_email)
  #-----------------------------------------------------------------------------
  # Run Read Counting
  #-----------------------------------------------------------------------------
  count_df_long <- run_read_counting(
    manifest = manifest,
    manifest_output_path = manifest_output_path,
    run_read_counting_stage = run_read_counting_stage,
    slurm_settings = slurm_settings,
    opt = opt,
    project_root_dir = project_root_dir,
    log_folder = log_folder
  )
  #-----------------------------------------------------------------------------
  # Optional Violin Plots
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Prepare Data for MAUDE
  #-----------------------------------------------------------------------------

  maude_counts_df <- run_prepare_data_for_MAUDE(
    count_df_long = count_df_long,
    opt = opt)
  
  #-----------------------------------------------------------------------------
  # pre MAUDE plots
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Run MAUDE
  #-----------------------------------------------------------------------------
  if (isTRUE(run_maude_stage)){log_info("Starting initial MAUDE run ...")}
  
  maude_results <- run_MAUDE(maude_counts_df = maude_counts_df,
                             opt = opt,
                             file_suffix = file_suffix,
                             run_maude_stage = run_maude_stage,
                             run_plots_stage = run_plots_stage)
  #-----------------------------------------------------------------------------
  # post MAUDE plots (Waterfall)
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Exporting MAUDE results + Combine Replicates
  #-----------------------------------------------------------------------------
  if (isTRUE(run_maude_stage)){
  export_df <- handle_auto_combine_replicates_and_export(
    file_info_suffix = file_info_suffix,
    opt = opt
    )
  }
  #-----------------------------------------------------------------------------
  # Finding Consensus Hits: MAUDE
  #-----------------------------------------------------------------------------
  if (isTRUE(opt$consensus_run)) {
    consensus_results <- run_consensus_call(
      maude_counts_df = maude_counts_df,
      opt = opt,
      run_maude_stage = run_maude_stage,
      run_plots_stage = run_plots_stage,
      n_reps = opt$consensus_n_reps,
      high_confidence_fdr = opt$consensus_high_confidence_FDR_threshold,
      high_confidence_hits_in_reps = opt$consensus_high_confidence_hits_in_X_reps,
      explorative_fdr = opt$consensus_explorative_FDR_threshold,
      explorative_hits_in_reps = opt$consensus_explorative_hits_in_X_reps
    )
  }
  
    
    
    
}


