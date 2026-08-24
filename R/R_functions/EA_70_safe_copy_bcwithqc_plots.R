# R/R_functions/EA_70_safe_copy_bcwithqc_plots.R

# R/R_functions/EA_50_safe_copy_bcwithqc_plots.R

safe_copy_bcwithqc_plots <- function(cfg, output_subdir = "04_bcwithqc_QC_plots") {
  bcwithqc_dir <- cfg$paths$bcwithqc_output_folder
  plot_dir <- file.path(cfg$paths$plots_output_folder, output_subdir)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!dir.exists(bcwithqc_dir)) {
    logger::log_warn("bcwithqc output directory does not exist: {bcwithqc_dir}")
    return(invisible(character()))
  } else {
    log_info("Copying QC_metrics from bcwithqc outputs...")
  }
  
  sample_dirs <- list.dirs(bcwithqc_dir, recursive = FALSE, full.names = TRUE)
  copied_files <- character()
  
  for (sample_dir in sample_dirs) {
    sample_name <- basename(sample_dir)
    
    source_files <- list.files(
      sample_dir,
      pattern = "(\\.png$|summary\\.tsv$)",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    for (source_file in source_files) {
      destination_name <- paste0(sample_name, "_", basename(source_file))
      destination_file <- file.path(plot_dir, destination_name)
      
      # Avoid silently overwriting files if bcwithqc contains identically named
      # matching files in multiple subdirectories for the same sample.
      if (file.exists(destination_file)) {
        relative_path <- substring(source_file, nchar(sample_dir) + 2)
        relative_dir <- dirname(relative_path)
        
        if (!identical(relative_dir, ".")) {
          safe_relative_dir <- gsub("[^A-Za-z0-9_.-]+", "_", relative_dir)
          destination_name <- paste0(sample_name, "_", safe_relative_dir, "_", basename(source_file))
          destination_file <- file.path(plot_dir, destination_name)
        }
      }
      
      copied <- file.copy(source_file, destination_file, overwrite = TRUE)
      
      if (isTRUE(copied)) {
        copied_files <- c(copied_files, destination_file)
      } else {
        logger::log_warn("Failed to copy bcwithqc QC file: {source_file}")
      }
    }
  }
  
  if (length(copied_files) == 0) {
    logger::log_warn("No bcwithqc PNG or summary.tsv files were found to copy from: {bcwithqc_dir}")
  } else {
    logger::log_info("Copied {length(copied_files)} bcwithqc QC files to: {plot_dir}")
  }
  
  invisible(copied_files)
}