# R/R_functions/AA_10_setup_functions.R

library(fs)

#===============================================================================
# Setup block functions
#===============================================================================

# Create a directory
make_clean_dir <- function(base_path, sub_path = NULL) {
  full_path <- if (is.null(sub_path) || sub_path == "") base_path else fs::path(base_path, sub_path)
  fs::dir_create(full_path, recurse = TRUE)
  fs::path_norm(full_path)
}

# construct a file path
get_file_path <- function(folder_path, file_name) {
  fs::path_norm(fs::path(folder_path, file_name))
}
# Old homebrew versions of the two

# make_clean_dir <- function(base_path, sub_path) {
#   
#   full_path <- paste0(base_path, "/", sub_path)
#   
#   full_path <- gsub("///", "/", full_path, fixed = TRUE)
#   full_path <- gsub("//", "/", full_path, fixed = TRUE)
#   
#   if (!dir.exists(full_path)) {
#     dir.create(full_path, recursive = TRUE)
#   }
#   return(full_path)
# }
# get_file_path <- function(folder_path, file_name){
#   full_path <- paste0(folder_path, "/",file_name)
#   full_path <- gsub("///", "/", full_path, fixed = TRUE)
#   full_path <- gsub("//", "/", full_path, fixed = TRUE)
#   return(full_path)
# }
#===============================================================================
# logging functions
#===============================================================================
initialize_pipeline_logger <- function(log_file) {
  log_dir <- dirname(log_file)
  
  if (!dir.exists(log_dir)) {
    stop("Log directory does not exist: ", log_dir, call. = FALSE)
  }
  
  logger::log_appender(logger::appender_tee(log_file))
  logger::log_layout(
    logger::layout_glue_generator(
      format = "{format(time, '%Y-%m-%d %H:%M:%S')} [{level}] {msg}"
    )
  )
  
  globalCallingHandlers(
    message = function(m) {
      logger::log_info("MESSAGE: {conditionMessage(m)}")
    },
    warning = function(w) {
      logger::log_warn("WARNING: {conditionMessage(w)}")
    },
    error = function(e) {
      logger::log_error("ERROR: {conditionMessage(e)}")
    }
  )
  
  logger::log_info("----------------------------------------------------------")
  logger::log_info("----------------------------------------------------------")
  logger::log_info("Logger initialized. Find log file at: {log_file}")
  invisible(log_file)
}

stop_log <- function(..., call. = FALSE) {
  msg <- paste0(...)

  stop(msg, call. = call.)
}
#===============================================================================
# Helpers for deciding if count command should run a specific part of the pipeline
#===============================================================================

should_run_stage <- function(start_with, stage) {
  stage_order <- c("beginning", "read_counting", "MAUDE_analysis", "generate_plots")
  
  if (!start_with %in% stage_order) {
    stop_log(
      "Invalid `start_with`: ", start_with,
      "\nAllowed values are: ",
      paste(stage_order, collapse = ", ")
    )
  }
  
  if (!stage %in% stage_order) {
    stop_log(
      "Invalid pipeline stage: ", stage,
      "\nAllowed values are: ",
      paste(stage_order, collapse = ", ")
    )
  }
  
  match(stage, stage_order) >= match(start_with, stage_order)
}

load_existing_tsv_or_stop <- function(path, label) {
  if (!file.exists(path)) {
    stop_log(
      "Cannot resume pipeline: missing ", label, ".\n",
      "Expected file:\n  ", path,
      "\n\nEither run with an earlier `start_with` value or check your output folder/suffix settings."
    )
  }
  
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}
get_file_info_suffix <- function(file_suffix) {
  x <- sub("\\.rds$", "", file_suffix)
  x <- sub("^_", "", x)
  x
}

memory_used_gb <- function() {
  sum(gc()[, 2]) / 1024
}


memory_rss_gb <- function() {
  status <- readLines("/proc/self/status")
  rss_line <- grep("^VmRSS:", status, value = TRUE)
  
  if (length(rss_line) != 1) {
    return(NA_real_)
  }
  
  rss_kb <- as.numeric(
    sub("^VmRSS:\\s+([0-9]+).*", "\\1", rss_line)
  )
  
  rss_kb / 1024^2
}