# R/cli_make_reference.R

#-------------------------------------------------------------------------------
# Identify paths
#-------------------------------------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  
  script_arg <- grep("^--file=", args, value = TRUE)
  
  if (length(script_arg) == 0) {
    stop("Could not find script path. Are you running the file with Rscript?")
  }
  
  script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
  dirname(script_path)
}

script_dir <- get_script_dir()
project_root_dir <- normalizePath(file.path(script_dir, ".."))

cli_args <- commandArgs(trailingOnly = TRUE)

is_snakemake <- "--snakemake" %in% cli_args
cli_args <- base::setdiff(cli_args, "--snakemake")

if (length(cli_args) < 1) {
  stop(
    "make_reference requires an input file.\n",
    "Old mode:       Rscript R/cli_make_reference_files.R <config.yaml>\n",
    "Snakemake mode: Rscript R/cli_make_reference_files.R <resolved_config.rds> --snakemake",
    call. = FALSE
  )
}

input_path <- normalizePath(cli_args[1], mustWork = TRUE)
extra_args <- cli_args[-1]
if (isTRUE(is_snakemake)) {
  #---------------------------------------------------------------------------
  # New Snakemake behavior
  #---------------------------------------------------------------------------
  
  bootstrap_path <- file.path(
    project_root_dir,
    "R",
    "R_snake_functions",
    "s00_bootstrap.R"
  )
  
  if (!file.exists(bootstrap_path)) {
    stop("Could not find bootstrap.R at:\n", bootstrap_path, call. = FALSE)
  }
  
  source(bootstrap_path)
  bootstrap_pipeline(project_root_dir)
  
  cfg <- tryCatch(
    readRDS(input_path),
    error = function(e) {
      stop("Could not load config.rds at ", input_path, "\nOriginal error: ", e$message)
    }
  )
  
  if (is.null(cfg$paths$log_files$make_reference)) {
    stop("No make_reference log file configured in cfg$paths$log_files.", call. = FALSE)
  }
  
  initialize_pipeline_logger(cfg$paths$log_files$make_reference)
  
  logger::log_info("Running make_reference_files in Snakemake mode.")
  logger::log_info("Resolved cfg: {input_path}")
  
  run_make_reference(
    config_path = cfg$paths$snake$resolved_config_yaml,
    project_root_dir = project_root_dir,
    cli_args = extra_args,
    run_setup = FALSE,
    cfg = cfg, 
    snakemake = TRUE) 
  
} else {
  #---------------------------------------------------------------------------
  # Old behavior
  #---------------------------------------------------------------------------
  
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop(
      "The R package 'yaml' is required to read config.yaml.\n",
      "Install it with: install.packages('yaml')",
      call. = FALSE
    )
  }
  
  source(file.path(project_root_dir, "R", "R_functions", "zzz_source_all.R"))
  
  logger::log_info("Running make_reference_files in legacy mode.")
  
  run_make_reference(
    config_path = input_path,
    project_root_dir = project_root_dir,
    cli_args = extra_args,
    snakemake = FALSE
  )
}