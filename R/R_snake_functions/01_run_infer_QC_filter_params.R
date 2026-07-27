run_infer_QC_filter_params <- function(cfg, project_root_dir) {
  
  manifest <- load_existing_tsv_or_stop(
    cfg$paths$manifest,
    "standardized FASTQ manifest"
  )
  
  manifest_read_types <- normalize_manifest_read_type(
    manifest = manifest,
    read_col = "read"
  )
  
  required_read_types <- c("SE", "R1", "R2")
  required_read_types <- required_read_types[
    required_read_types %in% manifest_read_types
  ]
  
  if (length(required_read_types) == 0L) {
    stop_log(
      "Could not determine a supported read type from the manifest. ",
      "Expected SE, R1, and/or R2."
    )
  }
  
  resolved_min_lengths <- switch(
    cfg$qc_filtering$min_length_state,
    
    infer = {
      inferred <- infer_qc_min_length(
        cfg = cfg,
        manifest = manifest,
        fastq_col = "symlink_file",
        read_col = "read",
        n_lines = 10000
      )
      
      inferred[required_read_types]
    },
    
    provided_single = {
      stats::setNames(
        rep(
          cfg$qc_filtering$min_length_single,
          length(required_read_types)
        ),
        required_read_types
      )
    },
    
    provided_R1R2 = {
      if ("SE" %in% required_read_types) {
        stop_log(
          "The manifest contains single-end data, but `qc_filtering.min_length` ",
          "was configured with separate R1/R2 values."
        )
      }
      
      c(
        R1 = cfg$qc_filtering$min_length_R1,
        R2 = cfg$qc_filtering$min_length_R2
      )[required_read_types]
    },
    
    stop_log(
      "Unknown qc_filtering.min_length_state: ",
      cfg$qc_filtering$min_length_state
    )
  )
  
  missing_lengths <- required_read_types[
    !required_read_types %in% names(resolved_min_lengths) |
      is.na(resolved_min_lengths[required_read_types])
  ]
  
  if (length(missing_lengths) > 0L) {
    stop_log(
      "No QC minimum length was resolved for required read type(s): ",
      paste(missing_lengths, collapse = ", ")
    )
  }
  
  get_min_length <- function(read_type) {
    if (!read_type %in% names(resolved_min_lengths)) {
      return("")
    }
    
    as.character(resolved_min_lengths[[read_type]])
  }
  
  
  writeLines(
    c(
      paste0("QC_FILTERING_RUN=", if (isTRUE(cfg$qc_filtering$run)) "true" else "false"),
      paste0("MANIFEST=", shQuote(cfg$paths$manifest)),
      paste0("QC_FILTERED_FOLDER=", shQuote(cfg$paths$qc_filtered_folder)),
      paste0("QC_MIN_QUAL=", shQuote(as.character(cfg$qc_filtering$min_qual))),
      paste0("QC_MIN_LENGTH_SE=", shQuote(get_min_length("SE"))),
      paste0("QC_MIN_LENGTH_R1=", shQuote(get_min_length("R1"))),
      paste0("QC_MIN_LENGTH_R2=", shQuote(get_min_length("R2")))
      ),
    cfg$paths$snake$QC_filtering_params_sh
  )
  
  logger::log_info("QC minimum-length mode: {cfg$qc_filtering$min_length_state}")
  
  logger::log_info(
    "Resolved QC minimum lengths: {values}",
    values = paste(
      paste0(names(resolved_min_lengths), "=", resolved_min_lengths),
      collapse = ", "
    )
  )
  
  logger::log_info(
    "QC filtering params written to: ",
    "{cfg$paths$snake$QC_filtering_params_sh}"
  )
  
  invisible(cfg)
}