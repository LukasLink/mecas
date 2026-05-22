infer_qc_min_length <- function(manifest,
                                fastq_col = "symlink_file",
                                n_lines = 10000) {
  
  if (!fastq_col %in% colnames(manifest)) {
    stop_log("Manifest does not contain column: ", fastq_col)
  }
  
  read_lengths <- lapply(seq_len(nrow(manifest)), function(i) {
    fastq_path <- manifest[[fastq_col]][i]
    sample_name <- manifest$pipeline_name[i]
    
    if (!file.exists(fastq_path)) {
      stop_log("FASTQ file does not exist for sample ", sample_name, ": ", fastq_path)
    }
    
    con <- if (grepl("\\.gz$", fastq_path)) {
      gzfile(fastq_path, open = "rt")
    } else {
      file(fastq_path, open = "rt")
    }
    
    on.exit(close(con), add = TRUE)
    
    lines <- readLines(con, n = n_lines, warn = FALSE)
    
    if (length(lines) < 4) {
      stop_log("FASTQ file appears too short or invalid for sample ", sample_name, ": ", fastq_path)
    }
    
    seq_lines <- lines[seq(2, length(lines), by = 4)]
    lengths <- nchar(seq_lines)
    
    most_common_length <- as.integer(names(sort(table(lengths), decreasing = TRUE)[1]))
    
    data.frame(
      sample_name = sample_name,
      fastq_path = fastq_path,
      inferred_read_length = most_common_length,
      stringsAsFactors = FALSE
    )
  })
  
  read_lengths <- do.call(rbind, read_lengths)
  
  unique_lengths <- unique(read_lengths$inferred_read_length)
  
  if (length(unique_lengths) != 1) {
    msg <- paste0(
      "Could not automatically determine a single QC minimum read length.\n",
      "Different input FASTQ files appear to have different dominant read lengths:\n\n",
      paste(
        paste0(
          read_lengths$sample_name,
          ": ",
          read_lengths$inferred_read_length
        ),
        collapse = "\n"
      ),
      "\n\nSet `qc_filtering.min_length` manually in config.yaml, ",
      "or set `qc_filtering.run: false` and perform QC filtering yourself."
    )
    
    stop_log(msg)
  }
  
  inferred <- unique_lengths[1]
  
  logger::log_info("Automatically inferred qc_min_length: {inferred} nt")
  
  inferred
}

ensure_gz_suffix <- function(x) {
  ifelse(grepl("\\.gz$", x), x, paste0(x, ".gz"))
}