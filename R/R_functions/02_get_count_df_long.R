#===============================================================================
# get_count_df_long
#===============================================================================

# -----------------------------------------------------------------------------
# Helpers for bcwithqc-derived count matrices
# -----------------------------------------------------------------------------
.prepare_bcwithqc_sgrna_annotation <- function(sgRNA_df) {
  required_cols <- c("seq", "sgrna_id", "sublib")
  missing_cols <- setdiff(required_cols, colnames(sgRNA_df))
  
  if (length(missing_cols) > 0) {
    stop(
      "`sgRNA_df` is missing required columns for bcwithqc import:\n",
      paste0("  - ", missing_cols, collapse = "\n"),
      call. = FALSE
    )
  }
  
  if ("type" %in% colnames(sgRNA_df)) {
    check_df <- sgRNA_df %>%
      dplyr::select(seq, sgrna_id, sublib, type)
  } else {
    logger::log_warn(
      paste0(
        "`merged_sgRNA_df` does not contain a `type` column. ",
        "Falling back to grep-based control annotation using sgRNA names. ",
        "Please update the library file to include a `type` column with values: ",
        "non_targeting_control, targeting_control, targeting."
      )
    )
    
    check_df <- sgRNA_df %>%
      dplyr::select(seq, sgrna_id, sublib) %>%
      dplyr::mutate(
        type = dplyr::case_when(
          grepl("^CONTROL_C_NONTARG_", sgrna_id) ~ "non_targeting_control",
          grepl("^CONTROL_C_", sgrna_id) ~ "targeting_control",
          TRUE ~ "targeting"
        )
      )
  }
  
  allowed_types <- c(
    "non_targeting_control",
    "targeting_control",
    "targeting"
  )
  
  invalid_types <- unique(check_df$type[!check_df$type %in% allowed_types])
  
  if (length(invalid_types) > 0) {
    stop(
      "`merged_sgRNA_df$type` contains invalid values:\n",
      paste0("  - ", invalid_types, collapse = "\n"),
      "\n\nAllowed values are:\n",
      paste0("  - ", allowed_types, collapse = "\n"),
      call. = FALSE
    )
  }
  
  return(check_df)
}

.apply_bcwithqc_alignment_checks_and_filter <- function(df,
                                                        sub_lib,
                                                        name,
                                                        cfg,
                                                        check_alignments = TRUE) {
  if (check_alignments == TRUE) {
    logger::log_info("Checking wrong alignments for: {name}")
    
    non_control_df <- df %>%
      dplyr::filter(type == "targeting", count > 0)
    
    correct_sum <- non_control_df %>%
      dplyr::filter(sublib == sub_lib) %>%
      dplyr::summarise(total = sum(count), .groups = "drop") %>%
      dplyr::pull(total)
    
    wrong_sum <- non_control_df %>%
      dplyr::filter(sublib != sub_lib) %>%
      dplyr::summarise(total = sum(count), .groups = "drop") %>%
      dplyr::pull(total)
    
    total_sum <- correct_sum + wrong_sum
    
    if (total_sum == 0) {
      logger::log_warn(
        "In file: {name} no reads/UMIs were detected. This is bad >.<"
      )
    } else if (total_sum < 1000) {
      logger::log_warn(
        paste0(
          "In file: {name} less than 1000 reads/UMIs were detected. ",
          "If you were not sequencing with low depth this is bad!"
        )
      )
    }
    
    logger::log_info("UMI/read count alignment stats:")
    
    if (total_sum > 0) {
      logger::log_info(
        "  Correct aligned to sgRNAs in sublibrary:     {correct_sum} ({round(100 * correct_sum / total_sum, 2)}%)"
      )
      logger::log_info(
        "  Wrong aligned to sgRNAs not in sublibrary:   {wrong_sum} ({round(100 * wrong_sum / total_sum, 2)}%)"
      )
      
      if (!is.na(wrong_sum / total_sum) && (wrong_sum / total_sum) > 0.5) {
        logger::log_warn(
          paste0(
            "More than 50% of reads aligned to sgRNAs from the wrong sublibrary. ",
            "Check whether the sublibraries in the library file and ",
            "fastq_name_table_xlsx have the same numbering."
          )
        )
      }
    } else {
      logger::log_info("  Correct aligned to sgRNAs in sublibrary:     0")
      logger::log_info("  Wrong aligned to sgRNAs not in sublibrary:   0")
    }
  }
  
  same_controls_in_all_sublibraries <-
    cfg$controls$same_controls_in_all_sublibraries
  
  if (same_controls_in_all_sublibraries == TRUE) {
    df <- df %>%
      dplyr::filter(
        sublib == sub_lib |
          type %in% c("non_targeting_control", "targeting_control")
      )
  } else {
    df <- df %>%
      dplyr::filter(sublib == sub_lib)
  }
  
  n_non_targeting_controls <- df %>%
    dplyr::filter(type == "non_targeting_control") %>%
    dplyr::pull(sgrna_id) %>%
    unique() %>%
    length()
  
  if (n_non_targeting_controls < 3) {
    logger::log_warn(
      paste0(
        "Only ", n_non_targeting_controls,
        " non-targeting controls detected after filtering. ",
        "MAUDE likely requires at least 3 non-targeting controls and may break downstream."
      )
    )
  } else if (n_non_targeting_controls < 50) {
    logger::log_warn(
      paste0(
        "Only ", n_non_targeting_controls,
        " non-targeting controls detected after filtering. ",
        "This may be statistically risky for normalization/modeling."
      )
    )
  }
  
  if (check_alignments == TRUE) {
    logger::log_info(
      "sgRNA coverage: {round(sum(df$count > 0) / nrow(df) * 100, 2)}%"
    )
    logger::log_info("--------------------------------------------")
  }
  
  df %>%
    dplyr::mutate(group_category = type) %>%
    dplyr::select(-sublib, -type)
}

read_bcwithqc_data <- function(
    dir_path,
    cfg,
    sgRNA_df = cfg$merged_sgRNA_df,
    sub_lib,
    name,
    check_alignments = TRUE
) {
  
  if (missing(sub_lib) || is.null(sub_lib) || length(sub_lib) != 1 || is.na(sub_lib) || !nzchar(sub_lib)) {
    stop("`sub_lib` must be provided to `read_bcwithqc_data()`.", call. = FALSE)
  }
  
  if (missing(name) || is.null(name) || length(name) != 1 || is.na(name) || !nzchar(name)) {
    stop("`name` must be provided to `read_bcwithqc_data()`.", call. = FALSE)
  }
  
  barcodes_path <- file.path(dir_path, "barcodes.tsv.gz")
  matrix_path <- file.path(dir_path, "matrix.mtx.gz")
  
  if (!file.exists(barcodes_path)) {
    stop("Missing bcwithqc barcode file for `", name, "`:\n  ", barcodes_path, call. = FALSE)
  }
  
  if (!file.exists(matrix_path)) {
    stop("Missing bcwithqc matrix file for `", name, "`:\n  ", matrix_path, call. = FALSE)
  }
  
  barcode_seqs <- readr::read_tsv(
    barcodes_path,
    col_names = FALSE,
    show_col_types = FALSE
  ) %>%
    dplyr::pull(1) %>%
    basename()
  
  mat <- Matrix::readMM(matrix_path)
  
  if (!inherits(mat, "CsparseMatrix")) {
    mat <- as(mat, "CsparseMatrix")
  }
  
  barcode_counts <- Matrix::colSums(mat)
  
  if (length(barcode_counts) != length(barcode_seqs)) {
    stop(
      "bcwithqc barcode/matrix size mismatch in `", name, "`: ",
      length(barcode_seqs), " barcodes but ",
      length(barcode_counts), " matrix columns.",
      call. = FALSE
    )
  }
  
  count_df <- tibble::tibble(
    seq = barcode_seqs,
    count = as.numeric(barcode_counts)
  )
  
  check_df <- .prepare_bcwithqc_sgrna_annotation(sgRNA_df)
  
  unmatched_count_sum <- count_df %>%
    dplyr::anti_join(check_df %>% dplyr::select(seq), by = "seq") %>%
    dplyr::summarise(total = sum(count), .groups = "drop") %>%
    dplyr::pull(total)
  
  if (!is.na(unmatched_count_sum) && unmatched_count_sum > 0) {
    logger::log_warn(
      "In file: {name} {unmatched_count_sum} reads/UMIs map to barcode sequences not present in `merged_sgRNA_df`."
    )
  }
  
  df <- check_df %>%
    dplyr::left_join(count_df, by = "seq") %>%
    dplyr::mutate(count = dplyr::coalesce(count, 0)) %>%
    dplyr::select(sgrna_id, count, sublib, type)
  
  .apply_bcwithqc_alignment_checks_and_filter(
    df = df,
    sub_lib = sub_lib,
    name = name,
    cfg = cfg,
    check_alignments = check_alignments
  )
}

process_bcwithqc_data <- function(
    cfg,
    parent_dir = cfg$paths$bcwithqc_output_folder,
    manifest = cfg$manifest,
    merged_sgRNA_df = cfg$merged_sgRNA_df,
    data_type = cfg$counting$data_type,
    skip_list = cfg$skip$files,
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
  
  is_skipped <- vapply(
    pipeline_manifest$pipeline_name,
    function(name) {
      any(vapply(
        skip_list,
        function(x) {
          if (startsWith(x, "re:")) {
            grepl(sub("^re:", "", x), name, perl = TRUE)
          } else {
            identical(x, name)
          }
        },
        logical(1)
      ))
    },
    logical(1)
  )
  
  if (any(is_skipped)) {
    logger::log_info(
      "Skipping bcwithqc outputs due to skip list: {paste(pipeline_manifest$pipeline_name[is_skipped], collapse = ', ')}"
    )
  }
  
  skipped_pipeline_names <- pipeline_manifest$pipeline_name[is_skipped]
  
  pipeline_manifest <- pipeline_manifest[!is_skipped, , drop = FALSE]
  
  actual_output_dirs <- list.dirs(
    parent_dir,
    recursive = FALSE,
    full.names = FALSE
  )
  
  expected_or_skipped_dirs <- c(
    pipeline_manifest$pipeline_name,
    skipped_pipeline_names
  )
  
  unexpected_output_dirs <- setdiff(
    actual_output_dirs,
    expected_or_skipped_dirs
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


