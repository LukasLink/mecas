#!/usr/bin/env Rscript

# R/cli_setup.R

# ==============================================================================
# Setup CLI for Snakemake pipeline
# ==============================================================================

args <- commandArgs(trailingOnly = FALSE)

# ------------------------------------------------------------------------------
# Identify this script path and project root
# ------------------------------------------------------------------------------

file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 0) {
  stop("Could not determine script path from commandArgs().", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)

# Assumes:
# project_root/
# └── R/
#     └── cli_setup.R
project_root_dir <- normalizePath(
  file.path(dirname(script_path), ".."),
  mustWork = TRUE
)

# ------------------------------------------------------------------------------
# Read optional config path argument
# ------------------------------------------------------------------------------

trailing_args <- commandArgs(trailingOnly = TRUE)

if (length(trailing_args) >= 1 && nzchar(trailing_args[1])) {
  config_path <- normalizePath(trailing_args[1], mustWork = TRUE)
} else {
  config_path <- file.path(project_root_dir, "config.yaml")
  
  if (!file.exists(config_path)) {
    stop(
      "No config path was provided and default config.yaml was not found:\n",
      config_path,
      call. = FALSE
    )
  }
  
  warning(
    "No config path provided. Using default config.yaml from project_root_dir:\n",
    config_path,
    call. = FALSE
  )
  
  config_path <- normalizePath(config_path, mustWork = TRUE)
}

# ------------------------------------------------------------------------------
# Source and run bootstrap
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# Run project setup
# ------------------------------------------------------------------------------

cfg <- project_setup(
  project_root_dir = project_root_dir,
  config_path = config_path,
  setup_mode = "setup",
  only_one_logger = FALSE
)

#-----------------------------------------------------------------------------
# Create Symlinks
#-----------------------------------------------------------------------------
logger::log_info("Creating symlinks of input files...")

strict_file_match_value <- cfg$paths$strict_file_match

log_info(
  paste0(
    "`cfg$paths$strict_file_match`: ",
    "typeof = ", typeof(strict_file_match_value),
    "; class = ", paste(class(strict_file_match_value), collapse = ", "),
    "; length = ", length(strict_file_match_value),
    "; value = ",
    paste(capture.output(dput(strict_file_match_value)), collapse = "")
  )
)


manifest <- prepare_fastq_inputs(
  fastq_dir = cfg$paths$input_folder,
  fastq_name_table_file_path = cfg$paths$fastq_name_table_xlsx,
  output_symlink_dir = cfg$paths$fastq_symlinks_folder,
  manifest_output_path = cfg$paths$manifest,
  strict_file_match = cfg$paths$strict_file_match
)
if (!file.exists(cfg$paths$manifest)) {
  stop("FASTQ manifest was not created: ", cfg$paths$manifest, call. = FALSE)
}
log_info("Finished creating symlinks of input files.")
# ------------------------------------------------------------------------------
# Save resolved cfg
# ------------------------------------------------------------------------------
saveRDS(cfg, cfg$paths$snake$resolved_config_rds)

yaml::write_yaml(
  cfg,
  cfg$paths$snake$resolved_config_yaml
)

writeLines(
  paste("Setup finished at", Sys.time()),
  cfg$paths$snake$done$setup
)

log_info("Setup complete.")
log_info("Project root: ", project_root_dir)
log_info("Config path: ", config_path)
log_info("Resolved cfg: ", cfg$paths$snake$resolved_config_rds)