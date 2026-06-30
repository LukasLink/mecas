# R/R_snake_functions/zzz_snake_source_all.R

source_snake_functions <- function(project_root_dir) {
  stopifnot(is.character(project_root_dir), length(project_root_dir) == 1)
  
  project_root_dir <- normalizePath(project_root_dir, mustWork = TRUE)
  
  source_dirs <- c(
    file.path(project_root_dir, "R", "R_functions"),
    file.path(project_root_dir, "R", "R_snake_functions")
  )
  
  missing_dirs <- source_dirs[!dir.exists(source_dirs)]
  
  if (length(missing_dirs) > 0) {
    stop(
      "Could not find required R function directory/directories:\n",
      paste(missing_dirs, collapse = "\n"),
      call. = FALSE
    )
  }
  
  r_files <- unlist(lapply(source_dirs, function(this_dir) {
    list.files(this_dir, pattern = "\\.R$", full.names = TRUE)
  }))
  
  r_files <- sort(r_files)
  
  exclude_files <- c(
    "zzz_source_all.R",          # legacy source-all file
    "zzz_snake_source_all.R",    # this file
    "s00_bootstrap.R"            # bootstrap file
  )
  
  r_files <- r_files[!basename(r_files) %in% exclude_files]
  
  for (f in r_files) {
    source(f, local = .GlobalEnv)
  }
  
  invisible(r_files)
}