# R/R_snake/CA_11_process_bcwithqc_data.R

#===============================================================================
# get_count_df_long
#===============================================================================


process_bcwithqc_data <- function(
    cfg,
    parent_dir = cfg$paths$bcwithqc_output_folder,
    manifest = cfg$manifest,
    merged_sgRNA_df = cfg$merged_sgRNA_df,
    data_type = cfg$counting$data_type,
    check_alignments = TRUE
) {
  
  if (!(data_type %in% c("reads", "umis"))) {
    stop("Invalid `data_type`. Expected `reads` or `umis`.", call. = FALSE)
  }
  
  if (!dir.exists(parent_dir)) {
    stop("bcwithqc output folder does not exist:\n  ", parent_dir, call. = FALSE)
  }
  
  if (!is.data.frame(manifest)) {
    stop("`manifest` must be a data frame.", call. = FALSE)
  }
  
  required_manifest_columns <- c(
    "pipeline_name",
    "bin_name",
    "sublibrary",
    "sample"
  )
  
  missing_manifest_columns <- setdiff(
    required_manifest_columns,
    colnames(manifest)
  )
  
  if (length(missing_manifest_columns) > 0) {
    stop(
      "The FASTQ manifest is missing columns required for bcwithqc count import:\n",
      paste0("  - ", missing_manifest_columns, collapse = "\n"),
      call. = FALSE
    )
  }
  
  # Paired-end samples have two manifest rows. At this stage, each pipeline_name
  # should be represented exactly once.
  pipeline_manifest <- manifest %>%
    dplyr::transmute(
      pipeline_name = trimws(as.character(pipeline_name)),
      bin_name = trimws(as.character(bin_name)),
      sublibrary = trimws(as.character(sublibrary)),
      sample_id = trimws(as.character(sample))
    ) %>%
    dplyr::distinct()
  
  conflicting_pipeline_names <- pipeline_manifest %>%
    dplyr::count(pipeline_name, name = "n_metadata_rows") %>%
    dplyr::filter(n_metadata_rows != 1)
  
  if (nrow(conflicting_pipeline_names) > 0) {
    stop(
      "Some `pipeline_name` values map to conflicting manifest metadata:\n",
      paste0(
        "  - ",
        conflicting_pipeline_names$pipeline_name,
        collapse = "\n"
      ),
      call. = FALSE
    )
  }
  
  
  actual_output_dirs <- list.dirs(
    parent_dir,
    recursive = FALSE,
    full.names = FALSE
  )
    
  unexpected_output_dirs <- setdiff(
    actual_output_dirs,
    pipeline_manifest$pipeline_name
  )
  
  if (length(unexpected_output_dirs) > 0) {
    logger::log_warn(
      "Ignoring bcwithqc output directories that are not present in the current manifest: {paste(unexpected_output_dirs, collapse = ', ')}"
    )
  }
  
  results <- vector(
    mode = "list",
    length = nrow(pipeline_manifest)
  )
  
  for (i in seq_len(nrow(pipeline_manifest))) {
    
    pipeline_name <- pipeline_manifest$pipeline_name[i]
    bin_name <- pipeline_manifest$bin_name[i]
    sublibrary <- pipeline_manifest$sublibrary[i]
    sample_id <- pipeline_manifest$sample_id[i]
    
    output_dir <- file.path(
      parent_dir,
      pipeline_name
    )
    
    if (!dir.exists(output_dir)) {
      stop(
        "Expected bcwithqc output directory does not exist for `",
        pipeline_name,
        "`:\n  ",
        output_dir,
        call. = FALSE
      )
    }
    
    matrix_folder_name <- if (data_type == "umis") {
      "raw_umis_bc_matrix"
    } else {
      "raw_reads_bc_matrix"
    }
    
    matrix_dir <- file.path(
      output_dir,
      matrix_folder_name
    )
    
    if (!dir.exists(matrix_dir)) {
      stop(
        "Expected bcwithqc matrix directory does not exist for `",
        pipeline_name,
        "`:\n  ",
        matrix_dir,
        call. = FALSE
      )
    }
    
    results[[i]] <- tryCatch(
      {
        read_bcwithqc_data(
          dir_path = matrix_dir,
          cfg = cfg,
          sgRNA_df = merged_sgRNA_df,
          sub_lib = sublibrary,
          name = pipeline_name,
          check_alignments = check_alignments
        ) %>%
          dplyr::rename(sgRNA = sgrna_id) %>%
          dplyr::mutate(
            bin_name = .env$bin_name,
            sublib = .env$sublibrary,
            sample = .env$sample_id,
            exp = paste(.env$sublibrary, .env$sample_id, sep = "_")
          )
      },
      error = function(e) {
        stop(
          "Failed to import bcwithqc counts for `",
          pipeline_name,
          "`:\n",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }
  
  count_df_long <- dplyr::bind_rows(results)
  
  if (nrow(count_df_long) == 0) {
    stop(
      "No bcwithqc count data were imported from:\n  ",
      parent_dir,
      call. = FALSE
    )
  }
  
  count_df_long
}


