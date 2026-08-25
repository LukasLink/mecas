#!/usr/bin/env Rscript

# R/cli_B_infer_QC_filter_params.R

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

if (length(cli_args) < 1) {
  stop("QC filtering requires an input file <resolved_config.rds>.\n", call. = FALSE)
}

input_path <- normalizePath(cli_args[1], mustWork = TRUE)
extra_args <- cli_args[-1]

#-------------------------------------------------------------------------------
# Bootstrap
#-------------------------------------------------------------------------------
bootstrap_path <- file.path(
  project_root_dir,
  "R",
  "R_functions",
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
#-------------------------------------------------------------------------------
# Logging
#-------------------------------------------------------------------------------
if (is.null(cfg$paths$log_files$infer_QC_filter_params)) {
  stop("No QC_filter log file configured in cfg$paths$infer_QC_filter_params", call. = FALSE)
}

initialize_pipeline_logger(cfg$paths$log_files$infer_QC_filter_params)

logger::log_info("Starting infer_QC_filter_params...")

run_infer_QC_filter_params(cfg = cfg, project_root_dir = project_root_dir)

writeLines(
  paste("infer_QC_filter_params finished at", Sys.time()),
  cfg$paths$snake$done$infer_QC_filter_params
)
log_info("infer_QC_filter_params complete.")