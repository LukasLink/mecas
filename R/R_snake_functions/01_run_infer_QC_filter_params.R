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
  
  if (qc_min_length_is_missing(cfg$qc_filtering$min_length)) {
    cfg$qc_filtering$min_length <- infer_qc_min_length(
      cfg = cfg,
      manifest = manifest,
      fastq_col = "symlink_file",
      read_col = "read",
      n_lines = 10000
    )
  }
  
  resolved_min_lengths <- normalize_qc_min_lengths(
    min_length = cfg$qc_filtering$min_length,
    required_read_types = required_read_types
  )
  
  # Preserve a scalar for purely single-end data.
  # Store named values for paired-end or mixed data.
  if (identical(required_read_types, "SE")) {
    cfg$qc_filtering$min_length <- unname(
      resolved_min_lengths[["SE"]]
    )
  } else {
    cfg$qc_filtering$min_length <- as.list(
      resolved_min_lengths
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
      paste0(
        "QC_FILTERING_RUN=",
        if (isTRUE(cfg$qc_filtering$run)) "true" else "false"
      ),
      paste0(
        "MANIFEST=",
        shQuote(cfg$paths$manifest)
      ),
      paste0(
        "QC_FILTERED_FOLDER=",
        shQuote(cfg$paths$qc_filtered_folder)
      ),
      paste0(
        "QC_MIN_QUAL=",
        shQuote(as.character(cfg$qc_filtering$min_qual))
      ),
      paste0(
        "QC_QUAL_OFFSET=",
        shQuote(as.character(cfg$qc_filtering$qual_offset))
      ),
      paste0(
        "QC_MIN_LENGTH_SE=",
        shQuote(get_min_length("SE"))
      ),
      paste0(
        "QC_MIN_LENGTH_R1=",
        shQuote(get_min_length("R1"))
      ),
      paste0(
        "QC_MIN_LENGTH_R2=",
        shQuote(get_min_length("R2"))
      )
    ),
    cfg$paths$snake$QC_filtering_params_sh
  )
  
  saveRDS(
    cfg,
    cfg$paths$snake$resolved_config_rds
  )
  
  yaml::write_yaml(
    cfg,
    cfg$paths$snake$resolved_config_yaml
  )
  
  logger::log_info(
    "QC filtering params written to: ",
    "{cfg$paths$snake$QC_filtering_params_sh}"
  )
  
  invisible(cfg)
}