# parse_coverage_file

#===============================================================================
# make coverage file
#===============================================================================

parse_coverage_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
  block_start_pattern <- "Checking wrong alignments for:"
  
  # find each sample block start, allowing logger/date prefix
  start_idx <- grep(block_start_pattern, lines, fixed = TRUE)
  
  if (length(start_idx) == 0) {
    stop(
      "No sample blocks found. Expected lines containing: 'Checking wrong alignments for:'",
      call. = FALSE
    )
  }
  
  # define block ends
  end_idx <- c(start_idx[-1] - 1, length(lines))
  
  parse_reads <- function(line) {
    if (is.na(line) || length(line) == 0) return(NA_real_)
    as.numeric(gsub(",", "", sub(".*?:\\s*([0-9,]+)\\s*\\(.*", "\\1", line)))
  }
  
  parse_perc <- function(line) {
    if (is.na(line) || length(line) == 0) return(NA_real_)
    as.numeric(sub(".*\\(([-0-9.]+)%\\).*", "\\1", line))
  }
  
  parse_coverage <- function(line) {
    if (is.na(line) || length(line) == 0) return(NA_real_)
    as.numeric(sub(".*sgRNA coverage:\\s*([-0-9.]+)%.*", "\\1", line))
  }
  
  out <- lapply(seq_along(start_idx), function(k) {
    block <- lines[start_idx[k]:end_idx[k]]
    
    sample_name <- sub(".*Checking wrong alignments for:\\s*", "", block[1])
    
    correct_line <- block[grep("Correct", block)][1]
    wrong_line   <- block[grep("Wrong", block)][1]
    cov_line     <- block[grep("sgRNA coverage:", block, fixed = TRUE)][1]
    
    correct_reads <- parse_reads(correct_line)
    correct_perc  <- parse_perc(correct_line)
    
    wrong_reads <- parse_reads(wrong_line)
    wrong_perc  <- parse_perc(wrong_line)
    
    coverage <- parse_coverage(cov_line)
    
    data.frame(
      sample_name   = sample_name,
      correct_reads = correct_reads,
      correct_perc  = paste0(sprintf("%.2f", correct_perc), "%"),
      wrong_reads   = wrong_reads,
      wrong_perc    = paste0(sprintf("%.2f", wrong_perc), "%"),
      coverage      = paste0(sprintf("%.2f", coverage), "%"),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::mutate(
        sample_name = sub(" ", "", sample_name)
      )
  })
  
  do.call(rbind, out)
}