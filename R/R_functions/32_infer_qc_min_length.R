infer_qc_min_length <- function(cfg,
                                manifest,
                                fastq_col = "symlink_file",
                                read_col = "read",
                                n_lines = 10000,
                                inferred_min_length_fraction = 0.80) {
  
  if (!fastq_col %in% colnames(manifest)) {
    stop_log("Manifest does not contain column: ", fastq_col)
  }
  
  read_types <- normalize_manifest_read_type(
    manifest = manifest,
    read_col = read_col
  )
  
  # A paired-end manifest should contain both R1 and R2
  if (
    any(read_types %in% c("R1", "R2")) &&
    !all(c("R1", "R2") %in% read_types)
  ) {
    stop_log(
      "Paired-end reads were detected, but the manifest does not contain ",
      "both R1 and R2 entries."
    )
  }
  
  read_lengths <- lapply(seq_len(nrow(manifest)), function(i) {
    fastq_path <- manifest[[fastq_col]][i]
    sample_name <- manifest$pipeline_name[i]
    read_type <- read_types[i]
    
    if (!file.exists(fastq_path)) {
      stop_log(
        "FASTQ file does not exist for sample ",
        sample_name,
        " (",
        read_type,
        "): ",
        fastq_path
      )
    }
    
    con <- if (grepl("\\.gz$", fastq_path)) {
      gzfile(fastq_path, open = "rt")
    } else {
      file(fastq_path, open = "rt")
    }
    
    on.exit(close(con), add = TRUE)
    
    lines <- readLines(con, n = n_lines, warn = FALSE)
    
    if (length(lines) < 4) {
      stop_log(
        "FASTQ file appears too short or invalid for sample ",
        sample_name,
        " (",
        read_type,
        "): ",
        fastq_path
      )
    }
    
    seq_lines <- lines[seq(2, length(lines), by = 4)]
    lengths <- nchar(seq_lines)
    
    most_common_length <- as.integer(
      names(sort(table(lengths), decreasing = TRUE)[1])
    )
    
    data.frame(
      sample_name = sample_name,
      read_type = read_type,
      fastq_path = fastq_path,
      dominant_read_length = most_common_length,
      stringsAsFactors = FALSE
    )
  })
  
  read_lengths <- do.call(rbind, read_lengths)
  
  present_read_types <- c("SE", "R1", "R2")
  present_read_types <- present_read_types[
    present_read_types %in% read_lengths$read_type
  ]
  
  inferred_lengths <- stats::setNames(
    integer(length(present_read_types)),
    present_read_types
  )
  
  for (current_read_type in present_read_types) {
    current_rows <- read_lengths[
      read_lengths$read_type == current_read_type,
      ,
      drop = FALSE
    ]
    
    unique_lengths <- unique(current_rows$dominant_read_length)
    
    if (length(unique_lengths) != 1L) {
      msg <- paste0(
        "Could not automatically determine a single QC minimum read length ",
        "for read type ",
        current_read_type,
        ".\n",
        "Different input FASTQ files have different dominant read lengths:\n\n",
        paste(
          paste0(
            current_rows$sample_name,
            " [",
            current_rows$read_type,
            "]: ",
            current_rows$dominant_read_length
          ),
          collapse = "\n"
        ),
        "\n\nSet `qc_filtering.min_length` manually in config.yaml, ",
        "or set `qc_filtering.run: false` and perform QC filtering yourself."
      )
      
      if (isTRUE(cfg$qc_filtering$run)) {
        stop_log(msg)
      } else {
        log_warn(
          "Error in QC filtering inference, which was ignored because ",
          "QC filtering is turned off."
        )
        log_warn(msg)
      }
    }
    
    dominant_read_length <- unique_lengths[[1]]
    
    inferred_min_length <- max(1L,floor(dominant_read_length * inferred_min_length_fraction))
    
    inferred_lengths[[current_read_type]] <- inferred_min_length
    
    logger::log_info(
      paste0(
        "Automatically inferred QC minimum length for ",
        "{current_read_type}: {inferred_min_length} nt ",
        "(dominant input length: {dominant_read_length} nt; ",
        "fraction: {inferred_min_length_fraction})."
      )
    )
  }
  
  inferred_lengths
}


normalize_manifest_read_type <- function(manifest, read_col = "read") {
  # Backwards compatibility with older single-end manifests
  if (!read_col %in% colnames(manifest)) {
    return(rep("SE", nrow(manifest)))
  }
  
  read_type <- toupper(trimws(as.character(manifest[[read_col]])))
  
  # Empty or NA means single-end
  read_type[is.na(read_type) | read_type == ""] <- "SE"
  
  invalid_read_types <- setdiff(unique(read_type), c("SE", "R1", "R2"))
  
  if (length(invalid_read_types) > 0) {
    stop_log(
      "Manifest contains invalid values in column `",
      read_col,
      "`: ",
      paste(invalid_read_types, collapse = ", "),
      ". Expected empty values, R1, or R2."
    )
  }
  
  read_type
}


