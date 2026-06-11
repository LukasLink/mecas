# run_plot.R

run_plot <- function(config_path, project_root_dir, cli_args) {
  message("Beginning Project setup...")
  
  cfg <- project_setup(
    project_root_dir = project_root_dir,
    config_path = config_path,
    setup_mode = "plot",
    use_old_suffix_construction = FALSE
  )
  
  logger::log_info("Finished Project setup.")
  
  run_create_plots(
    cfg = cfg,
    file_info_suffix = cfg$suffix$file_info_suffix
  )
  
  logger::log_info("DONE!")
  
  invisible(cfg)
}