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
  # Optional: handle use_only_these_controls
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
}