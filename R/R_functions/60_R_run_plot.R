# run_plot.R

run_plot <- function(config_path, project_root_dir, cli_args){
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
  # Run plot
  #-----------------------------------------------------------------------------
  run_create_plots(opt = opt, file_info_suffix = file_info_suffix)

  log_info("DONE!")  
}