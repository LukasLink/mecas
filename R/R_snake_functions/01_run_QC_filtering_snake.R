run_QC_filtering_snake <- function(cfg, project_root_dir) {
  
  manifest <- load_existing_tsv_or_stop(
    cfg$paths$manifest,
    "standardized FASTQ manifest"
  )
  if (is.na(cfg$qc_filtering$min_length) || is.null(cfg$qc_filtering$min_length)) {
    cfg$qc_filtering$min_length <- infer_qc_min_length(
      cfg = cfg,
      manifest = manifest,
      fastq_col = "symlink_file",
      n_lines = 10000
    )
  }

  writeLines(
    c(
      paste0("QC_FILTERING_RUN=", if (isTRUE(cfg$qc_filtering$run)) "true" else "false"),
      paste0("MANIFEST=", shQuote(cfg$paths$manifest)),
      paste0("QC_FILTERED_FOLDER=", shQuote(cfg$paths$qc_filtered_folder)),
      paste0("QC_MIN_QUAL=", shQuote(as.character(cfg$qc_filtering$min_qual))),
      paste0("QC_QUAL_OFFSET=", shQuote(as.character(cfg$qc_filtering$qual_offset))),
      paste0("QC_MIN_LENGTH=", shQuote(as.character(cfg$qc_filtering$min_length)))
    ),
    cfg$paths$snake$QC_filtering_params_sh
  )
  
  saveRDS(cfg, cfg$paths$snake$resolved_config_rds)
  yaml::write_yaml(cfg, cfg$paths$snake$resolved_config_yaml)
  
  logger::log_info("QC filtering params written to: {cfg$paths$snake$QC_filtering_params_sh}")
  
  invisible(cfg)
}