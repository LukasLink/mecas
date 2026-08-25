#!/usr/bin/env Rscript

# R/cli_F_reads_per_umi_count.R

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

if (length(cli_args) != 2L) {
  stop(
    "Usage: cli_F_reads_per_umi_count.R <resolved_config.rds> <sample>.",
    call. = FALSE
  )
}

input_path <- normalizePath(cli_args[1], mustWork = TRUE)
sample <- cli_args[2]

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
if (is.null(cfg$paths$log_files$reads_per_umi_count)) {
  stop("No umis_as_sublibs log file configured in cfg$paths$log_files.", call. = FALSE)
}

initialize_pipeline_logger(cfg$paths$log_files$reads_per_umi_count)

logger::log_info("Starting umis_as_sublibs...")

run_get_one_sample_umis_per_sublib(cfg = cfg, sample = sample)

log_info("umis_as_sublibs complete.")