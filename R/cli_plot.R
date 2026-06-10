# R/cli_plot.R
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
  stop("plot requires a config file to be provided: Rscript R/cli_count.R <config.yaml>", call. = FALSE)
}
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop(
    "The R package 'yaml' is required to read config.yaml.\n",
    "Install it with: install.packages('yaml')"
  )
}
config_path <- normalizePath(cli_args[1], mustWork = TRUE)

source(file.path(project_root_dir, "R", "R_functions", "zzz_source_all.R"))
#-------------------------------------------------------------------------------
# Run Plot function
#-------------------------------------------------------------------------------
run_plot(
  config_path = config_path,
  project_root_dir = project_root_dir,
  cli_args = cli_args[-1]
)