# run_count.R

run_count <- function(config_path, project_root_dir, cli_args) {
  #-----------------------------------------------------------------------------
  # Run setup
  #-----------------------------------------------------------------------------
  message("Beginning Project setup...")
  
  config <- yaml::read_yaml(config_path)
  
  # Priority explicit overrides > Rmd params > config.yaml
  cfg <- project_setup(
    project_root_dir = project_root_dir,
    config_path = config_path,
    setup_mode = "count",
    use_old_suffix_construction = FALSE
  )
  
  run_prepare_inputs_stage <- should_run_stage(cfg$run$start_with, "beginning")
  run_read_counting_stage  <- should_run_stage(cfg$run$start_with, "read_counting")
  run_maude_stage          <- should_run_stage(cfg$run$start_with, "MAUDE_analysis")
  run_plots_stage          <- should_run_stage(cfg$run$start_with, "generate_plots")
  
  logger::log_info("Finished Project setup.")
  
  #-----------------------------------------------------------------------------
  # Optional: handle use_only_these_controls
  #-----------------------------------------------------------------------------
  # If cfg$controls$use_only_these_controls contains entries, use them.
  # If empty, downstream functions should treat this as NULL/character(0).
  #
  # NOTE:
  # Do not rm() objects anymore. This should now be handled from cfg.
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Create Symlinks
  #-----------------------------------------------------------------------------
  
  manifest_output_path <- get_file_path(
    cfg$paths$rds_output_folder,
    paste0("standardized_fastq_manifest_", cfg$suffix$file_info_suffix, ".tsv")
  )
  
  if (run_prepare_inputs_stage) {
    
    logger::log_info("Creating symlinks of input files...")
    
    manifest <- prepare_fastq_inputs(
      fastq_dir = cfg$paths$input_folder,
      fastq_name_table_file_path = cfg$files$fastq_name_table_xlsx,
      output_symlink_dir = cfg$paths$fastq_symlinks_folder,
      manifest_output_path = manifest_output_path,
      strict_file_match = cfg$files$strict_file_match
    )
    
    logger::log_info("Finished creating symlinks of input files.")
    
  } else {
    
    logger::log_info(
      "Skipping input preparation because start_with is: {cfg$run$start_with}"
    )
    
    manifest <- load_existing_tsv_or_stop(
      manifest_output_path,
      "standardized FASTQ manifest"
    )
  }
  #-----------------------------------------------------------------------------
  # Create reference genome
  #-----------------------------------------------------------------------------
  run_make_reference(
    config_path = config_path,
    project_root_dir = project_root_dir,
    cli_args = cli_args,
    cfg = cfg, 
    run_setup = FALSE
  )
  #-----------------------------------------------------------------------------
  # Run Read Counting
  #-----------------------------------------------------------------------------
  
  if (isTRUE(run_read_counting_stage) || isTRUE(run_maude_stage)) {
    
    count_df_long <- run_read_counting(
      manifest = manifest,
      manifest_output_path = manifest_output_path,
      run_read_counting_stage = run_read_counting_stage,
      cfg = cfg,
      project_root_dir = project_root_dir
    )
  }
  
  #-----------------------------------------------------------------------------
  # Run MAUDE
  #-----------------------------------------------------------------------------
  
  if (isTRUE(run_maude_stage)) {
    
    # Prepare data for MAUDE
    maude_counts_df <- run_prepare_data_for_MAUDE(
      count_df_long = count_df_long,
      cfg = cfg
    )
    
    logger::log_info("Starting initial MAUDE run ...")
    
    maude_results <- run_MAUDE(
      maude_counts_df = maude_counts_df,
      cfg = cfg,
      file_suffix = cfg$suffix$file_suffix,
      run_maude_stage = run_maude_stage,
      run_plots_stage = run_plots_stage
    )
    
    # Export MAUDE results + combine replicates
    export_df <- handle_auto_combine_replicates_and_export(
      file_info_suffix = cfg$suffix$file_info_suffix,
      maude_results = maude_results,
      cfg = cfg
    )
    
    #-----------------------------------------------------------------------------
    # Finding Consensus Hits: MAUDE
    #-----------------------------------------------------------------------------
    
    if (isTRUE(cfg$consensus$run)) {
      consensus_results <- run_consensus_call(
        maude_counts_df = maude_counts_df,
        cfg = cfg,
        run_maude_stage = run_maude_stage,
        run_plots_stage = run_plots_stage
      )
    }
  }
  
  #-----------------------------------------------------------------------------
  # Run plots
  #-----------------------------------------------------------------------------
  
  if (isTRUE(run_plots_stage)) {
    
    logger::log_info("Entering plotting stage...")
    
    run_create_plots(
      cfg = cfg,
      file_info_suffix = cfg$suffix$file_info_suffix
    )
    
    logger::log_info("Finished plotting stage...")
  }
  
  logger::log_info("DONE!")
  
  invisible(cfg)
}