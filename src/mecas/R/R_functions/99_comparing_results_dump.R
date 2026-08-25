#===============================================================================
# get_count_df_long
#===============================================================================
read_file_to_df <- function(file_name,
                            cfg,
                            folder_path = cfg$paths$dedup_output_folder,
                            suffix_to_rm = "_dedup_idxstats.txt",
                            check_alignments = TRUE) {
  
  file_path <- file.path(folder_path, file_name)
  df <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  name <- sub(suffix_to_rm, "", file_name)
  sub_lib <- stringr::str_match(name, "^[A-Za-z]+_(L\\d+)_[^_]+")[, 2]
  exp <- stringr::str_replace(name, "^([^_]+_)", "")
  
  df <- df %>% 
    dplyr::rename(sgRNA = V1, count = V3, length = V2) %>%
    dplyr::select(-V4) %>%
    dplyr::mutate(exp = exp) %>%
    dplyr::filter(sgRNA != "*")
  
  merged_sgRNA_df <- cfg$merged_sgRNA_df
  
  if ("type" %in% colnames(merged_sgRNA_df)) {
    check_df <- merged_sgRNA_df %>%
      dplyr::select(sgrna_id, sublib, type)
  } else {
    logger::log_warn(
      paste0(
        "`merged_sgRNA_df` does not contain a `type` column. ",
        "Falling back to grep-based control annotation using sgRNA names. ",
        "Please update the library file to include a `type` column with values: ",
        "non_targeting_control, targeting_control, targeting."
      )
    )
    
    check_df <- merged_sgRNA_df %>%
      dplyr::select(sgrna_id, sublib) %>%
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
  
  joined_df <- df %>%
    dplyr::left_join(check_df, by = c("sgRNA" = "sgrna_id"))
  
  if (any(duplicated(joined_df$sgRNA))) {
    logger::log_warn(
      "Non-unique `sgRNA` after join; redoing join with 1:1 pairing by order."
    )
    
    df_idx <- df %>%
      dplyr::group_by(sgRNA) %>%
      dplyr::mutate(.pair_id = dplyr::row_number()) %>%
      dplyr::ungroup()
    
    check_idx <- check_df %>%
      dplyr::group_by(sgrna_id) %>%
      dplyr::mutate(.pair_id = dplyr::row_number()) %>%
      dplyr::ungroup()
    
    df <- df_idx %>%
      dplyr::left_join(
        check_idx,
        by = c("sgRNA" = "sgrna_id", ".pair_id")
      ) %>%
      dplyr::select(-.pair_id)
  } else {
    df <- joined_df
  }
  
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
      ) %>%
      dplyr::select(-sublib)
  } else {
    df <- df %>%
      dplyr::filter(sublib == sub_lib) %>%
      dplyr::select(-sublib)
  }
  
  n_non_targeting_controls <- df %>%
    dplyr::filter(type == "non_targeting_control") %>%
    dplyr::pull(sgRNA) %>%
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
  
  df <- df %>%
    dplyr::mutate(group_category = type) %>%
    dplyr::select(-type)
  
  return(df)
}

process_folder_files <- function(cfg,
                                 folder_path = cfg$paths$dedup_output_folder,
                                 skip_list = cfg$skip$files,
                                 check_alignments = TRUE) {
  
  data_type <- cfg$counting$data_type
  
  if (data_type == "umis") {
    files <- list.files(
      folder_path,
      pattern = "^[ILU]_L\\d+_[^_]+_dedup_idxstats\\.txt$",
      full.names = TRUE
    )
  } else if (data_type == "reads") {
    files <- list.files(
      folder_path,
      pattern = "^[ILU]_L\\d+_[^_]+_Aligned\\.sortedByCoord\\.out_idxstats\\.txt$",
      full.names = TRUE
    )
  } else {
    stop("Invalid data_type. data_type must be 'reads' or 'umis'", call. = FALSE)
  }
  
  combined_df <- dplyr::bind_rows(lapply(files, function(file) {
    file_name <- basename(file)
    
    if (data_type == "umis") {
      name <- sub("_dedup_idxstats.txt", "", file_name)
      matches <- stringr::str_match(
        file_name,
        "^([ILU])_L(\\d+)_([^_]+)_dedup_idxstats\\.txt$"
      )
      suffix_to_remove <- "_dedup_idxstats.txt"
    }
    
    if (data_type == "reads") {
      name <- sub("_Aligned.sortedByCoord.out_idxstats.txt", "", file_name)
      matches <- stringr::str_match(
        file_name,
        "^([ILU])_L(\\d+)_([^_]+)_Aligned\\.sortedByCoord\\.out_idxstats\\.txt$"
      )
      suffix_to_remove <- "_Aligned.sortedByCoord.out_idxstats.txt"
    }
    
    skip_this_file <- any(vapply(skip_list, function(x) {
      if (startsWith(x, "re:")) {
        grepl(sub("^re:", "", x), name, perl = TRUE)
      } else {
        identical(x, name)
      }
    }, logical(1)))
    
    if (skip_this_file) {
      logger::log_info("{name} skipped due to skip list")
      return(NULL)
    }
    
    if (is.na(matches[1])) {
      stop("Filename does not match expected pattern: ", file_name, call. = FALSE)
    }
    
    condition <- dplyr::case_when(
      matches[2] == "I" ~ "input",
      matches[2] == "L" ~ "lower",
      matches[2] == "U" ~ "upper"
    )
    
    sublib <- paste0("sublib_", matches[3])
    sample <- paste0("sample_", matches[4])
    
    df <- read_file_to_df(
      file_name = file_name,
      cfg = cfg,
      folder_path = folder_path,
      check_alignments = check_alignments,
      suffix_to_rm = suffix_to_remove
    )
    
    df <- df %>%
      dplyr::mutate(
        condition = condition,
        sublib = sublib,
        sample = sample
      ) %>%
      dplyr::select(-length)
    
    return(df)
  }))
  
  return(combined_df)
}
# -----------------------------------------------------------------------------
# Helpers for bcwithqc-derived count matrices
# -----------------------------------------------------------------------------
.normalise_sublib_to_L <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    grepl("^L\\d+$", x) ~ x,
    grepl("^sublib_?\\d+$", x, ignore.case = TRUE) ~ paste0("L", stringr::str_extract(x, "\\d+")),
    grepl("^ICS\\d+$", x, ignore.case = TRUE) ~ paste0("L", stringr::str_extract(x, "\\d+")),
    grepl("^\\d+$", x) ~ paste0("L", x),
    TRUE ~ x
  )
}

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
  
  check_df %>%
    dplyr::mutate(sublib = .normalise_sublib_to_L(sublib))
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

read_bcwithqc_data <- function(dir_path,
                               cfg,
                               sgRNA_df = cfg$merged_sgRNA_df,
                               sub_lib = NULL,
                               name = NULL,
                               check_alignments = TRUE) {
  if (is.null(name)) {
    name <- basename(dirname(dir_path))
  }
  
  if (is.null(sub_lib)) {
    sub_lib <- stringr::str_match(name, "^[A-Za-z]+_(L\\d+)_[^_]+$")[, 2]
  }
  
  if (is.na(sub_lib) || is.null(sub_lib)) {
    stop(
      "Could not infer sublibrary from bcwithqc folder name: ",
      name,
      call. = FALSE
    )
  }
  
  sub_lib <- .normalise_sublib_to_L(sub_lib)
  
  barcodes_path <- file.path(dir_path, "barcodes.tsv.gz")
  matrix_path <- file.path(dir_path, "matrix.mtx.gz")
  
  if (!file.exists(barcodes_path)) {
    stop("Missing bcwithqc barcode file: ", barcodes_path, call. = FALSE)
  }
  
  if (!file.exists(matrix_path)) {
    stop("Missing bcwithqc matrix file: ", matrix_path, call. = FALSE)
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
      "bcwithqc barcode/matrix size mismatch in ", name, ": ",
      length(barcode_seqs), " barcodes but ", length(barcode_counts),
      " matrix columns.",
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

process_bcwithqc_data <- function(cfg,
                                  parent_dir = cfg$paths$bcwithqc_output_folder,
                                  merged_sgRNA_df = cfg$merged_sgRNA_df,
                                  data_type = cfg$counting$data_type,
                                  skip_list = cfg$skip$files,
                                  check_alignments = TRUE) {
  if (!(data_type %in% c("reads", "umis"))) {
    stop("Invalid data_type. data_type must be 'reads' or 'umis'", call. = FALSE)
  }
  
  if (!dir.exists(parent_dir)) {
    stop("bcwithqc_output_folder does not exist: ", parent_dir, call. = FALSE)
  }
  
  subdirs <- list.dirs(parent_dir, recursive = FALSE, full.names = TRUE)
  
  if (length(subdirs) == 0) {
    logger::log_warn("No bcwithqc sample subdirectories found in: {parent_dir}")
    return(dplyr::bind_rows(list()))
  }
  
  results <- list()
  
  for (dir_path in subdirs) {
    folder_name <- basename(dir_path)
    
    matches <- stringr::str_match(folder_name, "^([A-Za-z]+)_(L\\d+)_([^_]+)$")
    
    if (any(is.na(matches))) {
      logger::log_warn(
        "Skipping bcwithqc folder with unexpected name: {folder_name}. Expected pattern examples: I_L1_1, L_L2_R1, or U_L10_5."
      )
      next
    }
    
    name <- folder_name
    
    is_skipped <- any(vapply(skip_list, function(x) {
      if (startsWith(x, "re:")) {
        grepl(sub("^re:", "", x), name, perl = TRUE)
      } else {
        identical(x, name)
      }
    }, logical(1)))
    
    if (is_skipped) {
      logger::log_info("{name} skipped due to skip list")
      next
    }
    
    condition_code <- matches[, 2]
    sub_lib <- matches[, 3]
    sample_id <- matches[, 4]
    
    condition <- dplyr::case_when(
      condition_code == "I" ~ "input",
      condition_code == "L" ~ "lower",
      condition_code == "U" ~ "upper",
      TRUE ~ condition_code
    )
    
    sublib <- paste0("sublib_", stringr::str_remove(sub_lib, "^L"))
    sample <- paste0("sample_", sample_id)
    exp <- stringr::str_replace(name, "^([^_]+_)", "")
    
    if (data_type == "umis") {
      matrix_dir <- file.path(dir_path, "raw_umis_bc_matrix")
    } else {
      matrix_dir <- file.path(dir_path, "raw_reads_bc_matrix")
    }
    
    if (!dir.exists(matrix_dir)) {
      logger::log_warn(
        "Skipping {name}: expected bcwithqc matrix directory does not exist: {matrix_dir}"
      )
      next
    }
    
    df <- tryCatch({
      read_bcwithqc_data(
        dir_path = matrix_dir,
        cfg = cfg,
        sgRNA_df = merged_sgRNA_df,
        sub_lib = sub_lib,
        name = name,
        check_alignments = check_alignments
      ) %>%
        dplyr::rename(sgRNA = sgrna_id) %>%
        dplyr::mutate(
          condition = condition,
          sublib = sublib,
          sample = sample,
          exp = exp
        )
    }, error = function(e) {
      logger::log_warn(
        "Skipping bcwithqc folder {folder_name} due to error: {conditionMessage(e)}"
      )
      return(NULL)
    })
    
    if (!is.null(df)) {
      results[[length(results) + 1]] <- df
    } else {
      logger::log_warn("{matrix_dir} yielded df -> NULL")
    }
  }
  
  dplyr::bind_rows(results)
}


normalize_count_df_long <- function(count_df_long,
                                    norm_method = "control_median",
                                    return_info = FALSE) {
  allowed_norm_methods <- c("control_median")
  
  if (!(norm_method %in% allowed_norm_methods)) {
    cat("ERROR: ", norm_method, "is not an allowed normalization method.\n")
    cat("Implemented normalization methods: ", allowed_norm_methods, "\n")
    stop()
  }
  
  pseudocount_added <- FALSE
  
  if (norm_method == "control_median") {
    
    norm_fac <- count_df_long %>%
      dplyr::filter(group_category %in% c(
        "targeting_control",
        "non_targeting_control",
        "kept_control"
      )) %>%
      dplyr::group_by(condition, sublib, sample) %>%
      dplyr::summarise(norm_factor = median(count), .groups = "drop")
    
    med_count <- median(count_df_long$count)
    
    count_df_long_continue <- count_df_long
    
    if (med_count == 0 | any(norm_fac$norm_factor == 0)) {
      cat("WARNING: at least one normalization factor is 0\n")
      cat("Median of all counts:", med_count, "\n")
      
      zero_nf <- norm_fac %>%
        dplyr::filter(norm_factor == 0)
      
      if (!(nrow(zero_nf) == 0)) {
        cat("Groups with norm_factor == 0:\n")
        print(norm_fac)
        cat("Adding global pseudocount (+1).\n")
        
        pseudocount_added <- TRUE
        
        count_df_long_plus_one <- count_df_long %>%
          dplyr::mutate(count = count + 1)
        
        norm_fac <- count_df_long_plus_one %>%
          dplyr::filter(group_category %in% c(
            "targeting_control",
            "non_targeting_control",
            "kept_control"
          )) %>%
          dplyr::group_by(condition, sublib, sample) %>%
          dplyr::summarise(norm_factor = median(count), .groups = "drop")
        
        med_count <- median(count_df_long_plus_one$count)
        
        count_df_long_continue <- count_df_long_plus_one
      }
      
      cat("--------------------------------------------\n")
    }
    
    return_df <- count_df_long_continue %>%
      dplyr::inner_join(norm_fac, by = c("condition", "sublib", "sample")) %>%
      dplyr::mutate(norm_count = (count * med_count) / norm_factor) %>%
      dplyr::mutate(count = round(norm_count, 2)) %>%
      dplyr::select(-c(norm_factor, norm_count))
  }
  
  return_df <- return_df %>%
    dplyr::mutate(count = round(count))
  
  if (isTRUE(return_info)) {
    return(list(
      count_df_long = return_df,
      pseudocount_added = pseudocount_added
    ))
  }
  
  return(return_df)
}



delete_bam_and_reads <- function(
    root,
    dry_run = TRUE
) {
  if (!dir.exists(root)) {
    stop("Directory does not exist: ", root)
  }
  
  all_files <- list.files(path = root,recursive = TRUE,full.names = TRUE,all.files = TRUE,no.. = TRUE)
  
  file_info <- file.info(all_files)
  all_files <- all_files[!file_info$isdir]
  
  targets <- all_files[
    grepl("\\.bam$", basename(all_files)) |
      basename(all_files) == "reads.tsv"
  ]
  
  if (length(targets) == 0) {
    message("No .bam or reads.tsv files found.")
    return(invisible(character()))
  }
  
  target_info <- file.info(targets)
  total_bytes <- sum(target_info$size, na.rm = TRUE)
  total_gb <- total_bytes / 1024^3
  
  if (dry_run) {
    message("Dry run: these files would be deleted:")
    print(targets)
    message("Total size that would be deleted: ",round(total_gb, 2)," GB")
    return(invisible(targets))
  }
  
  unlink(targets, force = FALSE)
  
  failed <- targets[file.exists(targets)]
  deleted <- setdiff(targets, failed)
  
  deleted_bytes <- sum(file.info(deleted)$size, na.rm = TRUE)
  
  if (length(failed) > 0) {
    warning("Some files could not be deleted:")
    print(failed)
    
    message("Deleted ",length(deleted)," of ",length(targets)," files.")
    message("Estimated total deleted: ",round((total_bytes - sum(file.info(failed)$size, na.rm = TRUE)) / 1024^3, 2)," GB")
  } else {
    message("Deleted ",length(targets)," files, freeing approximately ",round(total_gb, 2)," GB.")
  }
  
  invisible(targets)
}


compare_log_runtimes <- function(old_time,
                                 new_time,
                                 recursive = FALSE,
                                 tz = "UTC",
                                 strict = TRUE) {
  old_dirs <- as.character(unlist(old_time, use.names = FALSE))
  new_dirs <- as.character(unlist(new_time, use.names = FALSE))
  
  if (length(old_dirs) != length(new_dirs)) {
    stop("old_time and new_time must contain the same number of directories.")
  }
  
  if (any(!dir.exists(old_dirs))) {
    stop("These old_time directories do not exist: ",
         paste(old_dirs[!dir.exists(old_dirs)], collapse = ", "))
  }
  
  if (any(!dir.exists(new_dirs))) {
    stop("These new_time directories do not exist: ",
         paste(new_dirs[!dir.exists(new_dirs)], collapse = ", "))
  }
  
  pair_names <- names(old_time)
  
  if (is.null(pair_names) ||
      length(pair_names) != length(old_dirs) ||
      all(pair_names == "")) {
    pair_names <- names(new_time)
  }
  
  if (is.null(pair_names) ||
      length(pair_names) != length(old_dirs) ||
      all(pair_names == "")) {
    pair_names <- paste0("pair_", seq_along(old_dirs))
  }
  
  pair_names[pair_names == ""] <- paste0("pair_", which(pair_names == ""))
  
  timestamp_re <- "^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}"
  
  parse_log_runtime <- function(file) {
    lines <- readLines(file, warn = FALSE)
    
    matches <- regexpr(timestamp_re, lines, perl = TRUE)
    timestamp_strings <- regmatches(lines, matches)
    timestamp_strings <- timestamp_strings[nzchar(timestamp_strings)]
    
    timestamps <- as.POSIXct(
      timestamp_strings,
      format = "%Y-%m-%d %H:%M:%S",
      tz = tz
    )
    
    timestamps <- timestamps[!is.na(timestamps)]
    
    if (length(timestamps) == 0) {
      return(data.frame(
        first_timestamp = as.POSIXct(NA_character_, format = "%Y-%m-%d %H:%M:%S", tz = tz),
        last_timestamp = as.POSIXct(NA_character_, format = "%Y-%m-%d %H:%M:%S", tz = tz),
        elapsed_seconds = NA_real_,
        n_timestamps = 0L
      ))
    }
    
    first_ts <- min(timestamps)
    last_ts <- max(timestamps)
    
    data.frame(
      first_timestamp = first_ts,
      last_timestamp = last_ts,
      elapsed_seconds = as.numeric(difftime(last_ts, first_ts, units = "secs")),
      n_timestamps = length(timestamps)
    )
  }
  
  scan_log_dir <- function(dir) {
    files <- list.files(
      dir,
      pattern = "_(01_preprocess|03_count)\\.log$",
      full.names = TRUE,
      recursive = recursive
    )
    
    empty_ts <- as.POSIXct(character(), format = "%Y-%m-%d %H:%M:%S", tz = tz)
    
    if (length(files) == 0) {
      return(data.frame(
        sample_name = character(),
        log_type = character(),
        key = character(),
        path = character(),
        first_timestamp = empty_ts,
        last_timestamp = empty_ts,
        elapsed_seconds = numeric(),
        n_timestamps = integer(),
        stringsAsFactors = FALSE
      ))
    }
    
    base_names <- basename(files)
    
    sample_name <- sub("_(01_preprocess|03_count)\\.log$", "", base_names)
    
    log_type <- ifelse(
      grepl("_01_preprocess\\.log$", base_names),
      "preprocess",
      "count"
    )
    
    key <- paste(sample_name, log_type, sep = "__")
    
    if (any(duplicated(key))) {
      stop(
        "Duplicate log files found in directory: ", dir, "\nDuplicates: ",
        paste(unique(key[duplicated(key)]), collapse = ", ")
      )
    }
    
    runtimes <- do.call(rbind, lapply(files, parse_log_runtime))
    
    cbind(
      data.frame(
        sample_name = sample_name,
        log_type = log_type,
        key = key,
        path = files,
        stringsAsFactors = FALSE
      ),
      runtimes
    )
  }
  
  summarize_runtime_comparison <- function(old_runtime,
                                           new_runtime,
                                           diff_runtime,
                                           comparison_level) {
    keep <- !is.na(old_runtime) & !is.na(new_runtime) & !is.na(diff_runtime)
    
    old_runtime <- old_runtime[keep]
    new_runtime <- new_runtime[keep]
    diff_runtime <- diff_runtime[keep]
    
    improvement_percent <- ((old_runtime - new_runtime) / old_runtime) * 100
    
    if (length(diff_runtime) == 0) {
      return(data.frame(
        comparison_level = comparison_level,
        runtime_improvement_percent = NA_real_,
        n = 0L,
        
        mean_old_runtime_seconds = NA_real_,
        mean_new_runtime_seconds = NA_real_,
        mean_diff_seconds = NA_real_,
        sd_diff_seconds = NA_real_,
        max_diff_seconds = NA_real_,
        min_diff_seconds = NA_real_,
        
        mean_old_runtime_minutes = NA_real_,
        mean_new_runtime_minutes = NA_real_,
        mean_diff_minutes = NA_real_,
        sd_diff_minutes = NA_real_,
        max_diff_minutes = NA_real_,
        min_diff_minutes = NA_real_
      ))
    }
    
    data.frame(
      comparison_level = comparison_level,
      runtime_improvement_percent = mean(improvement_percent, na.rm = TRUE),
      n = length(diff_runtime),
      
      mean_old_runtime_seconds = mean(old_runtime),
      mean_new_runtime_seconds = mean(new_runtime),
      mean_diff_seconds = mean(diff_runtime),
      sd_diff_seconds = if (length(diff_runtime) > 1) sd(diff_runtime) else NA_real_,
      max_diff_seconds = max(diff_runtime),
      min_diff_seconds = min(diff_runtime),
      
      mean_old_runtime_minutes = mean(old_runtime) / 60,
      mean_new_runtime_minutes = mean(new_runtime) / 60,
      mean_diff_minutes = mean(diff_runtime) / 60,
      sd_diff_minutes = if (length(diff_runtime) > 1) sd(diff_runtime) / 60 else NA_real_,
      max_diff_minutes = max(diff_runtime) / 60,
      min_diff_minutes = min(diff_runtime) / 60
    )
  }
  all_file_pairs <- vector("list", length(old_dirs))
  
  for (i in seq_along(old_dirs)) {
    old_tbl <- scan_log_dir(old_dirs[i])
    new_tbl <- scan_log_dir(new_dirs[i])
    
    missing_in_new <- setdiff(old_tbl$key, new_tbl$key)
    missing_in_old <- setdiff(new_tbl$key, old_tbl$key)
    
    if (length(missing_in_new) > 0 || length(missing_in_old) > 0) {
      msg <- paste0(
        "Mismatched log files for directory pair '", pair_names[i], "'.",
        if (length(missing_in_new) > 0) {
          paste0("\nPresent only in old_time: ", paste(missing_in_new, collapse = ", "))
        } else "",
        if (length(missing_in_old) > 0) {
          paste0("\nPresent only in new_time: ", paste(missing_in_old, collapse = ", "))
        } else ""
      )
      
      if (strict) {
        stop(msg, call. = FALSE)
      } else {
        warning(msg, call. = FALSE)
      }
    }
    
    cmp <- merge(
      old_tbl,
      new_tbl,
      by = c("key", "sample_name", "log_type"),
      all = TRUE,
      suffixes = c("_old", "_new"),
      sort = FALSE
    )
    
    if (nrow(cmp) > 0) {
      cmp$dir_pair <- pair_names[i]
      cmp$old_dir <- old_dirs[i]
      cmp$new_dir <- new_dirs[i]
      
      cmp$time_diff_seconds <- cmp$elapsed_seconds_new - cmp$elapsed_seconds_old
      cmp$time_diff_minutes <- cmp$time_diff_seconds / 60
      cmp$time_diff_hours <- cmp$time_diff_seconds / 3600
      
      cmp <- cmp[, c(
        "dir_pair",
        "old_dir",
        "new_dir",
        "sample_name",
        "log_type",
        "path_old",
        "path_new",
        "first_timestamp_old",
        "last_timestamp_old",
        "elapsed_seconds_old",
        "n_timestamps_old",
        "first_timestamp_new",
        "last_timestamp_new",
        "elapsed_seconds_new",
        "n_timestamps_new",
        "time_diff_seconds",
        "time_diff_minutes",
        "time_diff_hours"
      )]
    }
    
    all_file_pairs[[i]] <- cmp
  }
  
  file_pairs <- do.call(rbind, all_file_pairs)
  rownames(file_pairs) <- NULL
  
  if (nrow(file_pairs) > 0) {
    file_pairs$log_type <- factor(file_pairs$log_type, levels = c("preprocess", "count"))
    file_pairs <- file_pairs[order(file_pairs$dir_pair, file_pairs$sample_name, file_pairs$log_type), ]
    file_pairs$log_type <- as.character(file_pairs$log_type)
    rownames(file_pairs) <- NULL
  }
  
  make_sample_pairs <- function(file_pairs) {
    if (nrow(file_pairs) == 0) {
      return(data.frame(
        dir_pair = character(),
        old_dir = character(),
        new_dir = character(),
        sample_name = character(),
        old_preprocess_seconds = numeric(),
        old_count_seconds = numeric(),
        old_total_seconds = numeric(),
        new_preprocess_seconds = numeric(),
        new_count_seconds = numeric(),
        new_total_seconds = numeric(),
        time_diff_seconds = numeric(),
        time_diff_minutes = numeric(),
        time_diff_hours = numeric(),
        stringsAsFactors = FALSE
      ))
    }
    
    groups <- split(
      file_pairs,
      paste(file_pairs$dir_pair, file_pairs$sample_name, sep = "___"),
      drop = TRUE
    )
    
    sample_tbl <- do.call(rbind, lapply(groups, function(x) {
      get_elapsed <- function(col, type) {
        value <- x[x$log_type == type, col]
        if (length(value) == 0) return(NA_real_)
        as.numeric(value[1])
      }
      
      old_preprocess <- get_elapsed("elapsed_seconds_old", "preprocess")
      old_count <- get_elapsed("elapsed_seconds_old", "count")
      
      new_preprocess <- get_elapsed("elapsed_seconds_new", "preprocess")
      new_count <- get_elapsed("elapsed_seconds_new", "count")
      
      old_total <- if (any(is.na(c(old_preprocess, old_count)))) {
        NA_real_
      } else {
        old_preprocess + old_count
      }
      
      new_total <- if (any(is.na(c(new_preprocess, new_count)))) {
        NA_real_
      } else {
        new_preprocess + new_count
      }
      
      time_diff <- new_total - old_total
      
      data.frame(
        dir_pair = x$dir_pair[1],
        old_dir = x$old_dir[1],
        new_dir = x$new_dir[1],
        sample_name = x$sample_name[1],
        old_preprocess_seconds = old_preprocess,
        old_count_seconds = old_count,
        old_total_seconds = old_total,
        new_preprocess_seconds = new_preprocess,
        new_count_seconds = new_count,
        new_total_seconds = new_total,
        time_diff_seconds = time_diff,
        time_diff_minutes = time_diff / 60,
        time_diff_hours = time_diff / 3600,
        stringsAsFactors = FALSE
      )
    }))
    
    rownames(sample_tbl) <- NULL
    sample_tbl[order(sample_tbl$dir_pair, sample_tbl$sample_name), ]
  }
  
  sample_pairs <- make_sample_pairs(file_pairs)
  
  preprocess_pairs <- file_pairs[file_pairs$log_type == "preprocess", , drop = FALSE]
  count_pairs <- file_pairs[file_pairs$log_type == "count", , drop = FALSE]
  
  summary_df <- rbind(
    summarize_runtime_comparison(
      old_runtime = preprocess_pairs$elapsed_seconds_old,
      new_runtime = preprocess_pairs$elapsed_seconds_new,
      diff_runtime = preprocess_pairs$time_diff_seconds,
      comparison_level = "preprocess_pairs"
    ),
    summarize_runtime_comparison(
      old_runtime = count_pairs$elapsed_seconds_old,
      new_runtime = count_pairs$elapsed_seconds_new,
      diff_runtime = count_pairs$time_diff_seconds,
      comparison_level = "count_pairs"
    ),
    summarize_runtime_comparison(
      old_runtime = sample_pairs$old_total_seconds,
      new_runtime = sample_pairs$new_total_seconds,
      diff_runtime = sample_pairs$time_diff_seconds,
      comparison_level = "sample_pairs"
    )
  )
  
  rownames(summary_df) <- NULL
  
  list(
    file_pairs = file_pairs,
    sample_pairs = sample_pairs,
    summary = summary_df
  )
}

build_conjugate_maude_summary <- function(
    conjugates_overlap,
    merged_sgRNA_df,
    DCA_dir,
    PA_dir,
    GAL_dir,
    mode = c("use_all", "use_only_confirmed_hits")
) {
  mode <- match.arg(mode)
  
  merged_sgRNA_df <- merged_sgRNA_df %>%
    mutate(
      entrez = if_else(
        is.na(entrez),
        sub("^CONTROL_C_(.*)_\\d+$", "\\1", sgrna_id),
        as.character(entrez)
      )
    )
  
  mean_na <- function(x) {
    x <- suppressWarnings(as.numeric(as.character(x)))
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }
  
  row_mean_na <- function(mat) {
    mat <- as.matrix(mat)
    storage.mode(mat) <- "numeric"
    out <- rowMeans(mat, na.rm = TRUE)
    out[is.nan(out)] <- NA_real_
    out
  }
  
  to_logical_flag <- function(x) {
    x <- toupper(trimws(as.character(x)))
    out <- rep(NA, length(x))
    out[x %in% c("TRUE", "T", "1", "YES", "Y")] <- TRUE
    out[x %in% c("FALSE", "F", "0", "NO", "N")] <- FALSE
    out
  }
  
  summarise_maude_dir <- function(input_dir, suffix) {
    rds_dir <- file.path(input_dir, "rds")
    
    files <- list.files(
      rds_dir,
      pattern = "^MAUDE_guide_stats.*\\.rds$",
      full.names = TRUE
    )
    
    if (length(files) == 0) {
      stop("No MAUDE_guide_stats*.rds files found in: ", rds_dir)
    }
    
    df <- lapply(files, function(f) {
      x <- readRDS(f)
      
      required_cols <- c("sgRNA", "input", "lower", "upper")
      missing_cols <- setdiff(required_cols, colnames(x))
      
      if (length(missing_cols) > 0) {
        stop(
          "File is missing required columns: ",
          paste(missing_cols, collapse = ", "),
          "\nFile: ",
          f
        )
      }
      
      x %>%
        select(sgRNA,input,upper,lower) %>% 
        transmute(
          sgRNA = as.character(.data$sgRNA),
          input = suppressWarnings(as.numeric(as.character(.data$input))),
          lower = suppressWarnings(as.numeric(as.character(.data$lower))),
          upper = suppressWarnings(as.numeric(as.character(.data$upper)))
        )
    }) %>%
      bind_rows() %>%
      group_by(sgRNA) %>%
      summarise(
        input = mean_na(input),
        lower = mean_na(lower),
        upper = mean_na(upper),
        .groups = "drop"
      ) %>%
      rename(
        !!paste0("input_", suffix) := input,
        !!paste0("lower_", suffix) := lower,
        !!paste0("upper_", suffix) := upper
      )
    
    df
  }
  
  PA_df  <- summarise_maude_dir(PA_dir,  "PA")
  DCA_df <- summarise_maude_dir(DCA_dir, "DCA")
  GAL_df <- summarise_maude_dir(GAL_dir, "GAL")
  
  combined_df <- list(PA_df, DCA_df, GAL_df) %>%
    Reduce(function(x, y) full_join(x, y, by = "sgRNA"), .)
  
  sgrna_col <- if ("sgrna_id" %in% colnames(merged_sgRNA_df)) {
    "sgrna_id"
  } else if ("sgRNA" %in% colnames(merged_sgRNA_df)) {
    "sgRNA"
  } else {
    stop("merged_sgRNA_df must contain either `sgrna_id` or `sgRNA`.")
  }
  
  if (!"entrez" %in% colnames(merged_sgRNA_df)) {
    stop("merged_sgRNA_df must contain an `entrez` column.")
  }
  
  sgRNA_to_entrez <- merged_sgRNA_df %>%
    transmute(
      sgRNA = as.character(.data[[sgrna_col]]),
      entrez = as.character(.data$entrez)
    ) %>%
    distinct()
  
  overlap_entrez <- as.character(conjugates_overlap$entrez)
  
  combined_df <- combined_df %>%
    left_join(sgRNA_to_entrez, by = "sgRNA") %>%
    filter(!is.na(entrez), entrez %in% overlap_entrez)
  
  metric_suffixes <- c("PA", "DCA", "GAL")
  
  if (mode == "use_all") {
    combined_df <- combined_df %>%
      mutate(
        mean_input = row_mean_na(across(all_of(paste0("input_", metric_suffixes)))),
        mean_upper = row_mean_na(across(all_of(paste0("upper_", metric_suffixes)))),
        mean_lower = row_mean_na(across(all_of(paste0("lower_", metric_suffixes))))
      )
  }
  
  if (mode == "use_only_confirmed_hits") {
    required_overlap_cols <- c("entrez", "PA", "DCA", "GalNAc")
    missing_overlap_cols <- setdiff(required_overlap_cols, colnames(conjugates_overlap))
    
    if (length(missing_overlap_cols) > 0) {
      stop(
        "conjugates_overlap is missing required columns: ",
        paste(missing_overlap_cols, collapse = ", ")
      )
    }
    
    hit_lookup <- conjugates_overlap %>%
      transmute(
        entrez = as.character(.data$entrez),
        PA = to_logical_flag(.data$PA),
        DCA = to_logical_flag(.data$DCA),
        GalNAc = to_logical_flag(.data$GalNAc)
      ) %>%
      group_by(entrez) %>%
      summarise(
        PA = any(PA, na.rm = TRUE),
        DCA = any(DCA, na.rm = TRUE),
        GalNAc = any(GalNAc, na.rm = TRUE),
        .groups = "drop"
      )
    
    combined_df <- combined_df %>%
      left_join(hit_lookup, by = "entrez")
    
    confirmed_row_mean <- function(df, metric) {
      value_cols <- paste0(metric, "_", c("PA", "DCA", "GAL"))
      flag_cols <- c("PA", "DCA", "GalNAc")
      
      value_mat <- as.matrix(df[, value_cols])
      storage.mode(value_mat) <- "numeric"
      
      flag_mat <- as.matrix(df[, flag_cols])
      flag_mat[is.na(flag_mat)] <- FALSE
      
      value_mat[!flag_mat] <- NA_real_
      
      row_mean_na(value_mat)
    }
    
    combined_df <- combined_df %>%
      mutate(
        mean_input = confirmed_row_mean(pick(everything()), "input"),
        mean_upper = confirmed_row_mean(pick(everything()), "upper"),
        mean_lower = confirmed_row_mean(pick(everything()), "lower")
      )
  }
  combined_df <- combined_df %>%
    mutate(
      mean = row_mean_na(across(all_of(c("mean_input", "mean_upper", "mean_lower"))))
    )
  
  # Keep best sgRNA per entrez based on highest overall mean
  mean_df <- combined_df %>%
    arrange(entrez, desc(mean), sgRNA) %>%
    group_by(entrez) %>%
    slice(1) %>%
    ungroup()
  
  # Cut merged_sgRNA_df down to selected sgRNAs only
  merged_sgRNA_df_selected <- merged_sgRNA_df %>%
    mutate(
      sgrna_id = as.character(.data[[sgrna_col]]),
      entrez = as.character(.data$entrez)
    ) %>%
    filter(sgrna_id %in% mean_df$sgRNA)
  
  # Make lookup; should now be one seq/sgRNA per entrez
  seq_sgrna_lookup <- merged_sgRNA_df_selected %>%
    select(
      entrez,
      new_seq = seq,
      new_sgRNA = sgrna_id
    ) %>%
    distinct(entrez, .keep_all = TRUE)
  
  # Modify overlap/conjugates dataframe
  modified_overlap_df <- conjugates_overlap %>%
    mutate(entrez = as.character(entrez)) %>%
    left_join(seq_sgrna_lookup, by = "entrez") %>%
    mutate(
      seq = if_else(!is.na(new_seq), new_seq, seq),
      sgRNA = if_else(!is.na(new_sgRNA), new_sgRNA, sgRNA)
    ) %>%
    select(-new_seq, -new_sgRNA)
  
  return(
    list(
      modified_overlap_df = modified_overlap_df,
      mean_df = mean_df
    )
  )
}




plot_l2fc_scatter <- function(df,
                              baseline,
                              x_axis,
                              y_axis,
                              save_image_file_path = NULL,
                              genes = NULL,
                              mark_N_special = NULL,
                              pseudocount = 1e-9) {
  # Packages used: dplyr, ggplot2, ggrepel, rlang
  if (!is.data.frame(df)) stop("df must be a data.frame / tibble.")
  if (!"Gene" %in% names(df)) stop("df must contain a column named 'Gene'.")
  
  baseline_nm <- rlang::as_string(rlang::ensym(baseline))
  x_nm        <- rlang::as_string(rlang::ensym(x_axis))
  y_nm        <- rlang::as_string(rlang::ensym(y_axis))
  
  missing_cols <- base::setdiff(c(baseline_nm, x_nm, y_nm), names(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns in df: ", paste(missing_cols, collapse = ", "))
  }
  
  means <- df |>
    dplyr::group_by(.data$Gene) |>
    dplyr::summarise(
      baseline_mean = mean(.data[[baseline_nm]], na.rm = TRUE),
      x_mean        = mean(.data[[x_nm]],        na.rm = TRUE),
      y_mean        = mean(.data[[y_nm]],        na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      L2FC_x = log2((x_mean + pseudocount) / (baseline_mean + pseudocount)),
      L2FC_y = log2((y_mean + pseudocount) / (baseline_mean + pseudocount)),
      delta  = L2FC_y - L2FC_x
    )
  
  genes <- if (is.null(genes)) character(0) else unique(as.character(genes))
  
  top_genes <- bottom_genes <- character(0)
  if (!is.null(mark_N_special) && is.finite(mark_N_special) && mark_N_special > 0) {
    n <- as.integer(mark_N_special)
    tmp <- means |>
      dplyr::filter(!is.na(.data$delta))
    
    top_genes <- tmp |>
      dplyr::slice_max(.data$delta, n = n, with_ties = FALSE) |>
      dplyr::pull(.data$Gene)
    
    bottom_genes <- tmp |>
      dplyr::slice_min(.data$delta, n = n, with_ties = FALSE) |>
      dplyr::pull(.data$Gene)
  }
  
  plot_df <- means |>
    dplyr::mutate(
      category = dplyr::case_when(
        .data$Gene %in% genes       ~ "highlight",
        .data$Gene %in% top_genes   ~ "top",
        .data$Gene %in% bottom_genes~ "bottom",
        TRUE                        ~ "other"
      ),
      label = dplyr::if_else(.data$category == "other", NA_character_, .data$Gene)
    )
  
  print(paste("NAs in L2FC_x:", sum(is.na(means$L2FC_x))))
  print(paste("NAs in L2FC_y:", sum(is.na(means$L2FC_y))))
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$L2FC_x, y = .data$L2FC_y, color = .data$category)) +
    ggplot2::geom_point(alpha = 0.8, size = 1.7) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$label),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.2,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = c(
      other     = "grey80",
      top       = "red3",
      bottom    = "dodgerblue3",
      highlight = "darkgreen"
    )) +
    ggplot2::labs(
      title = paste0("L2FC vs ", baseline_nm),
      x = paste0("log2FC(", x_nm, " / ", baseline_nm, ")"),
      y = paste0("log2FC(", y_nm, " / ", baseline_nm, ")"),
      color = NULL
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    ggplot2::theme_classic()
  
  if (!is.null(save_image_file_path) && nzchar(save_image_file_path)) {
    ggplot2::ggsave(filename = save_image_file_path, plot = p, width = 7, height = 6, dpi = 300)
  }
  
  print(p)
  invisible(list(means = means, plot = p))
}

compare_bcwithqc_align_UMI_tools <- function(
    data_type,
    cfg,
    bcwithqc_function = process_bcwithqc_data,
    addage = "",
    bcwithqc_dir = cfg$paths$bcwithqc_output_folder,
    read_diff_histo_x_lim = c(-1000,1000),
    read_diff_histo_bin_width = 10,
    # read_diff_density_99_x_lim = c(-1000,1000),
    # read_diff_density_95_x_lim = c(-1000,1000),
    align_UMI_tools_umi_folder = cfg$paths$dedup_output_folder,
    align_UMI_tools_read_folder = cfg$paths$mapped_output_folder){
  
  cfg$paths$bcwithqc_output_folder <- bcwithqc_dir
  cfg$paths$dedup_output_folder <- align_UMI_tools_umi_folder
  cfg$paths$mapped_output_folder <- align_UMI_tools_read_folder
  cfg$counting$data_type <- data_type
  
  if (!(data_type %in% c("umis","reads"))){
    print('ERROR: data_type must be "umis" or "reads"')
    return(NULL)
  }  
  
  if (addage != ""){
    addage <- paste0("(",addage,")")
  }
  
  bcwithqc_df_long <- bcwithqc_function(
    cfg = cfg,
    data_type = data_type,
    skip_list = cfg$skip$files
  )
  if (data_type == "umis"){
    align_UMI_tools_df_long <- process_folder_files(
      cfg = cfg,
      folder_path = align_UMI_tools_umi_folder,
      skip_list = cfg$skip$files) #Add threshold df if thresholds should be applied
    data_name <- "UMI"
  }
  if (data_type == "reads"){
    align_UMI_tools_df_long <- process_folder_files(
      cfg = cfg,
      folder_path = align_UMI_tools_read_folder,
      skip_list = cfg$skip$files) #Add threshold df if thresholds should be applied
    data_name <- "Read"
  } 
  
  # Join the data frames on sgRNA and exp
  matched_df <- full_join(bcwithqc_df_long, align_UMI_tools_df_long, by = c("sgRNA", "condition","sublib","sample","group_category","exp"), suffix = c("_bcwithqc", "_align_UMI_tools")) %>%
    mutate(
      # Replace missing counts with 0
      count_bcwithqc = ifelse(is.na(count_bcwithqc), 0, count_bcwithqc),
      count_align_UMI_tools = ifelse(is.na(count_align_UMI_tools), 0, count_align_UMI_tools),
      
      # Fill in missing 'other' columns from available one
      # Example: if there’s another column like 'condition' or 'replicate'
      # replicate = coalesce(replicate_bcwithqc, replicate_align_UMI_tools)
      # (Uncomment and adapt as needed!)
    )
  
  # Compute the difference in count per matched entry
  matched_df <- matched_df %>%
    mutate(diff = count_bcwithqc - count_align_UMI_tools)
  
  # Total counts
  total_bcwithqc <- sum(matched_df$count_bcwithqc, na.rm = TRUE)
  total_align    <- sum(matched_df$count_align_UMI_tools, na.rm = TRUE)
  
  # Total difference
  total_diff <- total_bcwithqc - total_align
  
  # Percentage difference relative to align_UMI_tools
  perc_diff <- round(
    100 * total_diff / total_align,
    digits = 1
  )
  
  # Mean difference
  mean_diff <- mean(matched_df$diff, na.rm = TRUE)
  
  # Median difference
  median_diff <- median(matched_df$diff, na.rm = TRUE)
  
  # Standard deviation
  sd_diff <- sd(matched_df$diff, na.rm = TRUE)
  
  # Print results
  cat(sprintf("%-45s %15s\n", paste("Total", data_name, addage, "bcwithqc:"), total_bcwithqc))
  cat(sprintf("%-45s %15s\n", paste("Total", data_name, addage, "align_UMI_tools:"), total_align))
  cat(sprintf("%-45s %15s\n", paste("Total", data_name, addage, "Difference:"), total_diff))
  cat(sprintf("%-45s %14.1f%%\n", paste("Percentage", data_name, addage, "Difference:"), perc_diff))
  cat(sprintf("%-45s %15s\n", paste("Median", addage, "Difference:"), median_diff))
  cat(sprintf("%-45s %15.2f\n", paste("Mean", data_name, addage, "Difference:"), mean_diff))
  cat(sprintf("%-45s %15.2f\n", paste(data_name, addage, "Standard Deviation:"), sd_diff))
  
  total_df <- tibble(
    source = c("bcwithqc", "align_UMI_tools"),
    total_count = c(sum(matched_df$count_bcwithqc), sum(matched_df$count_align_UMI_tools))
  )
  
  p <- ggplot(total_df, aes(x = source, y = total_count, fill = source)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = total_count), vjust = -0.5) +
    labs(
      title = paste("Total",data_name, addage,"Counts"),
      x = "",
      y = paste("Total",data_name, addage,"Counts")
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  print(p)
  
  nonzero_df <- tibble(
    source = c("bcwithqc", "align_UMI_tools", "Max"),
    nonzero_count = c(sum(matched_df$count_bcwithqc > 0), sum(matched_df$count_align_UMI_tools > 0), nrow(matched_df))
  )
  
  p <- ggplot(nonzero_df, aes(x = source, y = nonzero_count, fill = source)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = nonzero_count), vjust = -0.5) +
    labs(
      title = paste("Coverage of",data_name, addage,"Counts"),
      x = "",
      y = paste("Coverage of",data_name, addage,"Counts")
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  print(p)
  
  # Assuming you've matched the two runs: count_bcwithqc vs count_align_UMI_tools
  p <- ggplot(matched_df, aes(x = count_bcwithqc + 1, y = count_align_UMI_tools + 1)) +
    geom_count(
      aes(
        size = after_stat(n),
        colour = after_stat(n)
      ),
      alpha = 0.9
    ) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_size_area(
      max_size = 4,
      trans = "log10",
      name = "Count overlap"
    ) +
    scale_colour_gradientn(
      colours = c("grey90", "skyblue2", "dodgerblue3", "navy"),
      trans = "log10",
      name = "Count overlap"
    ) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = paste("Comparison of sgRNA", data_name, addage, "Counts"),
      x = paste("bcwithqc", data_name, addage, "Counts"),
      y = paste("align_UMI_tools", data_name, addage, "Counts")
    ) +
    theme_minimal()
  
  print(p)
  # Histogram
  p <- ggplot(matched_df, aes(x = diff)) +
    geom_histogram(
      binwidth = read_diff_histo_bin_width,
      boundary = 0,
      fill = "steelblue",
      color = "white"
    ) +
    annotate(
      "label",
      x = read_diff_histo_x_lim[1],
      y = Inf,
      label = paste0("Binwidth: ", read_diff_histo_bin_width),
      hjust = 0,
      vjust = 1.1,
      fill = "white",
      color = "black",
      label.size = 0.3
    ) +
    labs(
      title = paste("Distribution of", data_name, addage, "Count Differences (bcwithqc - align_UMI_tools)"),
      x = paste(data_name, addage, "Count Difference"),
      y = "Frequency"
    ) +
    xlim(read_diff_histo_x_lim[1], read_diff_histo_x_lim[2]) +
    theme_minimal()
  print(p)
  # Compute the 2.5th and 97.5th percentiles (95% interval)
  x_limits <- quantile(matched_df$diff, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Plot with limited x-axis
  p <- ggplot(matched_df, aes(x = diff)) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    coord_cartesian(xlim = x_limits) +  # Zoom in without removing data
    labs(title = paste("Density of",data_name, addage, "Count Differences (95% Range)"),
         x = "bcwithqc - align_UMI_tools") +
    theme_minimal()
  print(p)
  x_limits <- quantile(matched_df$diff, probs = c(0.005, 0.995), na.rm = TRUE)
  
  # Plot with limited x-axis
  p <- ggplot(matched_df, aes(x = diff)) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    coord_cartesian(xlim = x_limits) +  # Zoom in without removing data
    labs(title = paste("Density of",data_name, addage, "Count Differences (99% Range)"),
         x = "bcwithqc - align_UMI_tools") +
    theme_minimal()
  print(p)
  matched_df <- matched_df %>%
    mutate(
      average = (count_bcwithqc + count_align_UMI_tools) / 2,
      diff = count_bcwithqc - count_align_UMI_tools
    )
  
  p <- ggplot(matched_df, aes(x = average, y = diff)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Bland–Altman Plot for",data_name, addage,"Counts"),
         x = paste("Average", data_name, addage, "Count"),
         y = "Difference (bcwithqc - align_UMI_tools)") +
    theme_minimal()
  print(p)
  return(matched_df)
}


plot_mapping_qc_boxes <- function(xlsx_paths, sample_names, output_dir = NULL) {
  
  stopifnot(length(xlsx_paths) == length(sample_names))
  
  cols_to_plot <- c(
    "correct_reads",
    "correct_perc",
    "wrong_reads",
    "wrong_perc",
    "coverage"
  )
  
  # Colours are assigned in pairs:
  # xlsx 1/2 = light/dark blue
  # xlsx 3/4 = light/dark red
  # xlsx 5/6 = light/dark green
  # etc.
  paired_colours <- c(
    "#9ecae1", "#08519c",  # light blue, dark blue
    "#fcae91", "#a50f15",  # light red, dark red
    "#a1d99b", "#006d2c",  # light green, dark green
    "#dadaeb", "#54278f",  # light purple, dark purple
    "#fdd0a2", "#a63603",  # light orange, dark orange
    "#bdbdbd", "#252525"   # light grey, dark grey
  )
  
  if (length(xlsx_paths) > length(paired_colours)) {
    stop(
      "Not enough predefined colours for ",
      length(xlsx_paths),
      " input files. Please add more colours to paired_colours."
    )
  }
  
  colour_df <- data.frame(
    dataset = sample_names,
    pair_colour = paired_colours[seq_along(sample_names)],
    stringsAsFactors = FALSE
  )
  
  read_one_file <- function(path, name) {
    df <- readxl::read_excel(path)
    
    colnames(df) <- colnames(df) %>%
      stringr::str_trim() %>%
      stringr::str_replace_all("\\s+", "_") %>%
      stringr::str_replace_all("-", "_") %>%
      tolower()
    
    missing_cols <- base::setdiff(cols_to_plot, colnames(df))
    
    if (length(missing_cols) > 0) {
      stop(
        "File ", path, " is missing columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    
    df %>%
      dplyr::select(dplyr::all_of(cols_to_plot)) %>%
      dplyr::mutate(
        dplyr::across(
          c(correct_perc, wrong_perc, coverage),
          ~ as.numeric(stringr::str_replace(as.character(.x), "%", ""))
        ),
        dplyr::across(
          c(correct_reads, wrong_reads),
          ~ as.numeric(.x)
        )
      ) %>%
      dplyr::mutate(dataset = name)
  }
  
  all_data <- purrr::map2_dfr(xlsx_paths, sample_names, read_one_file) %>%
    dplyr::left_join(colour_df, by = "dataset") %>%
    dplyr::mutate(dataset = factor(dataset, levels = sample_names))
  
  plot_data <- all_data %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(cols_to_plot),
      names_to = "metric",
      values_to = "value"
    )
  
  pretty_names <- c(
    correct_reads = "Correct reads",
    correct_perc  = "Correct reads [%]",
    wrong_reads   = "Wrong reads",
    wrong_perc    = "Wrong reads [%]",
    coverage      = "sgRNA Coverage [%]"
  )
  
  plots <- list()
  
  for (metric_name in cols_to_plot) {
    p <- plot_data %>%
      dplyr::filter(metric == metric_name) %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = dataset,
          y = value,
          fill = pair_colour
        )
      ) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_jitter(
        ggplot2::aes(color = pair_colour),
        width = 0.15,
        alpha = 0.5,
        size = 1
      ) +
      ggplot2::scale_fill_identity() +
      ggplot2::scale_color_identity() +
      ggplot2::labs(
        x = NULL,
        y = pretty_names[[metric_name]],
        title = paste("align_UMI_tools vs. bcwithqc -", pretty_names[[metric_name]])
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
    
    plots[[metric_name]] <- p
    
    if (!is.null(output_dir)) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      
      ggplot2::ggsave(
        filename = file.path(output_dir, paste0(metric_name, "_boxplot.png")),
        plot = p,
        width = 6,
        height = 4
      )
    }
  }
  
  for (p in plots) {
    print(p)
  }
  
  return(plots)
}



compare_hit_xlsx_pairs <- function(xlsx_paths,
                                   sample_names,
                                   output_dir = NULL,
                                   FDR_threshold = NULL,
                                   rank_method = "spearman",
                                   make_plots = TRUE) {
  
  stopifnot(length(xlsx_paths) == length(sample_names))
  
  if (length(xlsx_paths) %% 2 != 0) {
    stop("xlsx_paths must contain an even number of files because comparisons are done in pairs.")
  }
  
  required_cols <- c("entrez", "FDR", "significanceZ")
  
  paired_colours <- c(
    "#9ecae1", "#08519c",  # light blue, dark blue
    "#fcae91", "#a50f15",  # light red, dark red
    "#a1d99b", "#006d2c",  # light green, dark green
    "#dadaeb", "#54278f",  # light purple, dark purple
    "#fdd0a2", "#a63603",  # light orange, dark orange
    "#bdbdbd", "#252525"   # light grey, dark grey
  )
  
  if (length(xlsx_paths) > length(paired_colours)) {
    stop(
      "Not enough predefined colours for ",
      length(xlsx_paths),
      " input files. Please add more colours to paired_colours."
    )
  }
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  colour_df <- data.frame(
    dataset = sample_names,
    pair_colour = paired_colours[seq_along(sample_names)],
    stringsAsFactors = FALSE
  )
  
  read_one_file <- function(path, name) {
    
    df <- readxl::read_excel(path)
    
    colnames(df) <- colnames(df) %>%
      stringr::str_trim() %>%
      stringr::str_replace_all("\\s+", "_") %>%
      stringr::str_replace_all("-", "_") %>%
      tolower()
    
    # after lowercasing, required columns are:
    required_cols_lower <- tolower(required_cols)
    
    missing_cols <- base::setdiff(required_cols_lower, colnames(df))
    
    if (length(missing_cols) > 0) {
      stop(
        "File ", path, " is missing columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    
    df %>%
      dplyr::select(
        entrez,
        fdr,
        significancez
      ) %>%
      dplyr::mutate(
        entrez = as.character(entrez),
        fdr = as.numeric(fdr),
        significancez = as.numeric(significancez),
        dataset = name
      ) %>%
      dplyr::filter(!is.na(entrez), entrez != "")
  }
  
  all_data <- purrr::map2_dfr(xlsx_paths, sample_names, read_one_file) %>%
    dplyr::left_join(colour_df, by = "dataset") %>%
    dplyr::mutate(dataset = factor(dataset, levels = sample_names))
  
  if (!is.null(FDR_threshold)) {
    all_data_for_venn <- all_data %>%
      dplyr::filter(fdr <= FDR_threshold)
  } else {
    all_data_for_venn <- all_data
  }
  
  pair_df <- data.frame(
    pair_id = seq_len(length(sample_names) / 2),
    sample_1 = sample_names[seq(1, length(sample_names), by = 2)],
    sample_2 = sample_names[seq(2, length(sample_names), by = 2)],
    stringsAsFactors = FALSE
  )
  
  # ------------------------------------------------------------
  # 1. FDR boxplot
  # ------------------------------------------------------------
  
  fdr_boxplot <- ggplot2::ggplot(
    all_data,
    ggplot2::aes(
      x = dataset,
      y = fdr,
      fill = pair_colour
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      alpha = 0.35
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = pair_colour),
      width = 0.15,
      alpha = 0.5,
      size = 1
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      x = NULL,
      y = "FDR",
      title = "FDR distribution"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  if (make_plots) {
    print(fdr_boxplot)
  }
  
  if (!is.null(output_dir)) {
    ggplot2::ggsave(
      filename = file.path(output_dir, "FDR_boxplot.png"),
      plot = fdr_boxplot,
      width = 7,
      height = 4
    )
  }
  
  # Optional: often useful because FDR values can be strongly compressed near 0
  fdr_log_boxplot <- all_data %>%
    dplyr::mutate(
      neg_log10_fdr = -log10(pmax(fdr, .Machine$double.xmin))
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = dataset,
        y = neg_log10_fdr,
        fill = pair_colour
      )
    ) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = pair_colour),
      width = 0.15,
      alpha = 0.5,
      size = 1
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      x = NULL,
      y = "-log10(FDR)",
      title = "-log10(FDR) distribution"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  if (make_plots) {
    print(fdr_log_boxplot)
  }
  
  if (!is.null(output_dir)) {
    ggplot2::ggsave(
      filename = file.path(output_dir, "minus_log10_FDR_boxplot.png"),
      plot = fdr_log_boxplot,
      width = 7,
      height = 4
    )
  }
  
  # ------------------------------------------------------------
  # 2. Pairwise Venn diagrams
  # ------------------------------------------------------------
  
  make_pair_venn <- function(df1, df2, name1, name2, title = NULL) {
    
    entrez1 <- unique(df1$entrez)
    entrez2 <- unique(df2$entrez)
    
    area1 <- length(entrez1)
    area2 <- length(entrez2)
    cross_area <- length(intersect(entrez1, entrez2))
    
    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp)
    on.exit({
      grDevices::dev.off()
      unlink(tmp)
    }, add = TRUE)
    
    g_list <- VennDiagram::draw.pairwise.venn(
      area1 = area1,
      area2 = area2,
      cross.area = cross_area,
      category = c(name1, name2),
      fill = c("lightblue", "lightgreen"),
      alpha = c(0.6, 0.6),
      cat.cex = 0.9,
      cex = 1.3,
      col = c("darkblue", "darkgreen"),
      cat.fontface = "bold",
      cat.col = c("darkblue", "darkgreen"),
      label.col = c("darkblue", "black", "darkgreen"),
      cat.pos = c(-20, 20),
      cat.dist = c(0.03, 0.03),
      ind = TRUE
    )
    
    g <- grid::grobTree(gList = do.call(grid::gList, g_list))
    
    if (!is.null(title)) {
      title_grob <- grid::textGrob(
        title,
        gp = grid::gpar(fontsize = 14, fontface = "bold")
      )
      
      g <- gridExtra::arrangeGrob(
        title_grob,
        g,
        ncol = 1,
        heights = c(0.12, 0.88)
      )
    }
    
    return(g)
  }
  
  venns <- list()
  overlap_summary <- list()
  
  for (i in seq_len(nrow(pair_df))) {
    
    name1 <- pair_df$sample_1[i]
    name2 <- pair_df$sample_2[i]
    
    df1 <- all_data_for_venn %>% dplyr::filter(dataset == name1)
    df2 <- all_data_for_venn %>% dplyr::filter(dataset == name2)
    
    entrez1 <- unique(df1$entrez)
    entrez2 <- unique(df2$entrez)
    
    overlap_n <- length(intersect(entrez1, entrez2))
    
    overlap_summary[[i]] <- data.frame(
      pair = paste(name1, "vs", name2),
      dataset_1 = name1,
      dataset_2 = name2,
      n_dataset_1 = length(entrez1),
      n_dataset_2 = length(entrez2),
      n_overlap = overlap_n,
      perc_of_dataset_1 = round(100 * overlap_n / length(entrez1), 2),
      perc_of_dataset_2 = round(100 * overlap_n / length(entrez2), 2)
    )
    
    venn_title <- paste("Entrez overlap:", name1, "vs", name2)
    
    if (!is.null(FDR_threshold)) {
      venn_title <- paste0(
        venn_title,
        " | FDR <= ",
        FDR_threshold
      )
    }
    
    venn_grob <- make_pair_venn(
      df1 = df1,
      df2 = df2,
      name1 = name1,
      name2 = name2,
      title = venn_title
    )
    
    venns[[paste(name1, "vs", name2)]] <- venn_grob
    
    if (make_plots) {
      grid::grid.newpage()
      grid::grid.draw(venn_grob)
    }
    
    if (!is.null(output_dir)) {
      png_file <- file.path(
        output_dir,
        paste0(
          "venn_",
          gsub("[^A-Za-z0-9]+", "_", name1),
          "_vs_",
          gsub("[^A-Za-z0-9]+", "_", name2),
          ".png"
        )
      )
      
      grDevices::png(png_file, width = 1600, height = 1200, res = 200)
      grid::grid.draw(venn_grob)
      grDevices::dev.off()
    }
  }
  
  overlap_summary <- dplyr::bind_rows(overlap_summary)
  
  # ------------------------------------------------------------
  # 3. Rank-order similarity
  # ------------------------------------------------------------
  
  calculate_rank_similarity <- function(df1,
                                        df2,
                                        name1,
                                        name2,
                                        subset_name = "all",
                                        sign_filter = NULL) {
    
    if (!is.null(sign_filter)) {
      if (sign_filter == "positive") {
        df1 <- df1 %>% dplyr::filter(significancez >= 0)
        df2 <- df2 %>% dplyr::filter(significancez >= 0)
      }
      
      if (sign_filter == "negative") {
        df1 <- df1 %>% dplyr::filter(significancez < 0)
        df2 <- df2 %>% dplyr::filter(significancez < 0)
      }
    }
    
    df1_ranked <- df1 %>%
      dplyr::filter(!is.na(entrez), !is.na(fdr)) %>%
      dplyr::group_by(entrez) %>%
      dplyr::slice_min(order_by = fdr, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(fdr) %>%
      dplyr::mutate(rank_1 = dplyr::row_number()) %>%
      dplyr::select(entrez, fdr_1 = fdr, significancez_1 = significancez, rank_1)
    
    df2_ranked <- df2 %>%
      dplyr::filter(!is.na(entrez), !is.na(fdr)) %>%
      dplyr::group_by(entrez) %>%
      dplyr::slice_min(order_by = fdr, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(fdr) %>%
      dplyr::mutate(rank_2 = dplyr::row_number()) %>%
      dplyr::select(entrez, fdr_2 = fdr, significancez_2 = significancez, rank_2)
    
    paired <- dplyr::inner_join(df1_ranked, df2_ranked, by = "entrez")
    
    if (!is.null(FDR_threshold)) {
      paired <- paired %>%
        dplyr::filter(fdr_1 <= FDR_threshold | fdr_2 <= FDR_threshold)
    }
    
    if (nrow(paired) < 2) {
      rank_cor <- NA_real_
    } else {
      rank_cor <- cor(
        paired$rank_1,
        paired$rank_2,
        method = rank_method,
        use = "complete.obs"
      )
    }
    
    data.frame(
      pair = paste(name1, "vs", name2),
      dataset_1 = name1,
      dataset_2 = name2,
      subset = subset_name,
      n_shared_entrez = nrow(paired),
      rank_correlation = round(rank_cor, 3),
      method = rank_method,
      stringsAsFactors = FALSE
    )
  }
  
  rank_similarity <- list()
  
  for (i in seq_len(nrow(pair_df))) {
    
    name1 <- pair_df$sample_1[i]
    name2 <- pair_df$sample_2[i]
    
    df1 <- all_data %>% dplyr::filter(dataset == name1)
    df2 <- all_data %>% dplyr::filter(dataset == name2)
    
    rank_similarity[[paste0(i, "_all")]] <- calculate_rank_similarity(
      df1 = df1,
      df2 = df2,
      name1 = name1,
      name2 = name2,
      subset_name = "all",
      sign_filter = NULL
    )
    
    rank_similarity[[paste0(i, "_positive")]] <- calculate_rank_similarity(
      df1 = df1,
      df2 = df2,
      name1 = name1,
      name2 = name2,
      subset_name = "significanceZ >= 0",
      sign_filter = "positive"
    )
    
    rank_similarity[[paste0(i, "_negative")]] <- calculate_rank_similarity(
      df1 = df1,
      df2 = df2,
      name1 = name1,
      name2 = name2,
      subset_name = "significanceZ < 0",
      sign_filter = "negative"
    )
  }
  
  rank_similarity <- dplyr::bind_rows(rank_similarity)
  
  print(overlap_summary)
  print(rank_similarity)
  
  if (!is.null(output_dir)) {
    writexl::write_xlsx(
      list(
        entrez_overlap_summary = overlap_summary,
        FDR_rank_similarity_summary = rank_similarity
      ),
      path = file.path(output_dir, "hit_comparison_summary.xlsx")
    )
  }
  
  return(
    list(
      all_data = all_data,
      overlap_summary = overlap_summary,
      rank_similarity = rank_similarity,
      fdr_boxplot = fdr_boxplot,
      fdr_log_boxplot = fdr_log_boxplot,
      venns = venns
    )
  )
}




# Helper function to apply FDR filter
apply_FDR_filter <- function(df, FDR_threshold) {
  df <- df %>%
    filter(FDR <= FDR_threshold)
  return(df)
}

# Helper function to calculate number of overlapping entrez
calculate_overlap <- function(reference_df, target_df) {
  reference_entrez <- reference_df$entrez
  target_entrez <- target_df$entrez
  overlap_count <- sum(target_entrez %in% reference_entrez)
  return(overlap_count)
}
# Helper function to calculate the correlation between a refrence_df
# and a target_df
calculate_correlation <- function(reference_df,
                                  target_df,
                                  method = "spearman",
                                  use = "complete.obs",
                                  FDR_threshold = 0.05,
                                  significance_filter = NULL) {
  # keep only needed cols, drop duplicates, sort by FDR, assign ranks
  if (!is.null(significance_filter)) {
    reference_df <- reference_df %>% filter(significanceZ * significance_filter > 0)
    target_df <- target_df %>% filter(significanceZ * significance_filter > 0)
  }
  
  ref <- reference_df[, c("entrez", "FDR")]
  tgt <- target_df[, c("entrez", "FDR")]
  
  ref <- ref[!is.na(ref$entrez) & !is.na(ref$FDR), ]
  tgt <- tgt[!is.na(tgt$entrez) & !is.na(tgt$FDR), ]
  
  ref <- ref %>% rename(ref_FDR = FDR)
  tgt <- tgt %>% rename(tgt_FDR = FDR)
  
  # pair by entrez (inner join)
  paired <- merge(ref, tgt, by = "entrez") %>%
    filter(ref_FDR <= FDR_threshold | tgt_FDR <= FDR_threshold)
  # if not enough paired points, correlation is undefined
  if (nrow(paired) < 2) return(NA_real_)
  
  correlation <- cor(paired$ref_FDR, paired$tgt_FDR, method = method, use = use)
  correlation <- round(correlation, 3)
  return(correlation)
}
# Function to process each scenario
# generates a dataframe with hits, and overlap to reference for each df in df_list

make_compare_df <- function(df_list, name_list, FDR_threshold, significance_filter = NULL, top_N = NULL) {
  
  if (length(df_list) != length(name_list)){
    stop(paste("ERROR: in make_compare_df -> length of df_list:",length(df_list), " must be equal to lenght of name_list:",length(name_list)))
  }
  # First dataframe as the reference (100% by definition)
  reference_df <- df_list[[1]] %>%
    apply_FDR_filter(FDR_threshold)
  if (!is.null(significance_filter)) {
    reference_df <- reference_df %>% filter(significanceZ * significance_filter > 0)
  }
  
  # Initialize empty lists for storing results
  Overlap <- list()
  Overlap_percentage <- list()
  Total_Hits <- list()
  
  # Loop through the dataframes and calculate overlap for each subsample percentage
  for (i in seq_along(df_list)) {
    
    target_df <- df_list[[i]]
    
    # Apply FDR filter to the target dataframe
    target_df <- apply_FDR_filter(target_df, FDR_threshold)
    
    # Apply additional significance filter if provided
    if (!is.null(significance_filter)) {
      target_df <- target_df %>% filter(significanceZ * significance_filter > 0)
    }
    
    # Check overlap with the reference
    overlap_count <- calculate_overlap(reference_df, target_df)
    
    # Calculate percentage overlap
    overlap_percentage_value <- round((overlap_count / nrow(reference_df)) * 100,2)
    
    total_hits_value <- nrow(target_df)
    
    # Append results to the lists
    Overlap <- c(Overlap, overlap_count)
    Overlap_percentage <- c(Overlap_percentage, overlap_percentage_value)
    Total_Hits <- c(Total_Hits, total_hits_value)
  }
  # Combine the lists into a final dataframe
  final_results <- data.frame(
    Name = name_list,
    Overlap = unlist(Overlap),
    Overlap_percentage = unlist(Overlap_percentage),
    Total_Hits = unlist(Total_Hits)
  )
  
  # Return the final dataframe
  return(final_results)
}
make_venn_from_overlap_df <- function(df,
                                      ncol = 3,
                                      main_title = "Venn Diagramm",
                                      comparing = "") {
  stopifnot(all(c("Name", "Overlap", "Total_Hits") %in% colnames(df)))
  if (nrow(df) < 2) stop("df needs at least 2 rows (1 reference + >=1 target).")
  
  if (comparing != ""){
    if (comparing == "subsamples"){
      df <- df %>% 
        mutate(Name = paste0("sub_perc_",Name))
    } else {
      stop("currently only '' (sublib) and 'subsamples' are viable for comparing ")
    }
  }
  
  
  ref_name <- df$Name[1]
  ref_total <- df$Total_Hits[1]
  
  grobs <- list()
  #This is supposed to stop it drawing a venn diagramm with draw.pairwise.venn
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  
  for (i in 2:nrow(df)) {
    tgt_name  <- df$Name[i]
    tgt_total <- df$Total_Hits[i]
    overlap   <- df$Overlap[i]
    
    # sanity (avoid negative/invalid overlaps)
    overlap <- max(0, overlap)
    overlap <- min(overlap, ref_total, tgt_total)
    
    
    g_list <- VennDiagram::draw.pairwise.venn(
      area1 = ref_total,
      area2 = tgt_total,
      cross.area = overlap,
      category = c(ref_name, tgt_name),
      fill = c("lightgreen", "lightblue"),
      alpha = c(0.6, 0.6),
      cat.cex = 1.0,
      cex = 1.3,
      col = c("darkgreen", "darkblue"),
      cat.fontface = "bold",
      cat.col = c("darkgreen", "darkblue"),      # category label colours
      label.col = c("darkgreen", "black", "darkblue"),  # left / overlap / right numbers
      cat.pos  = c(-60, 60),
      cat.dist = c(0.2, 0.2),
      ind = TRUE
    )
    
    # wrap list-of-grobs into a single grobTree so gridExtra can place it
    g <- grid::grobTree(gList = do.call(grid::gList, g_list))
    g <- grid::grobTree(g, vp = grid::viewport(width = 0.5, height = 0.7))
    grobs[[tgt_name]] <- g
  }
  
  title_grob <- grid::textGrob(main_title, gp = grid::gpar(fontsize = 14, fontface = "bold"))
  
  return(
    gridExtra::arrangeGrob(
      title_grob,
      gridExtra::arrangeGrob(grobs = grobs, ncol = ncol),
      ncol = ncol,
      heights = c(0.08, 0.92)
    )
  )
}
# Make the Correlation Heatmap
make_correlation_heatmap <- function(df_list,
                                     name_list,
                                     comparing = "",
                                     plot_title = NULL,
                                     significance_filter = NULL,
                                     FDR_threshold = 0.05) {
  
  if (comparing == "subsamples"){
    name_list = paste0("sub_perc_",name_list)
  }
  
  n <- length(df_list)
  
  # correlation matrix
  cor_mat <- matrix(
    NA_real_,
    nrow = n, ncol = n,
    dimnames = list(name_list, name_list)
  )
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      cor_mat[i, j] <- calculate_correlation(
        df_list[[i]],
        df_list[[j]],
        FDR_threshold = FDR_threshold,
        significance_filter = significance_filter
      )
    }
  }
  
  cor_df <- as.data.frame(as.table(cor_mat))
  colnames(cor_df) <- c("Reference", "Target", "Correlation")
  
  p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Target, y = Reference, fill = Correlation)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(Correlation), "", sprintf("%.3f", Correlation))),
      size = 3
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = paste("Pairwise Spearman Correlation of FDR Ranks for",plot_title),
      x = NULL, y = NULL,
      fill = "Rank corr"
    )
  
  return(p)
}
plot_overlap <- function(df, title, comparing = "") {
  # Ensure the input dataframe contains the necessary columns
  if (!all(c("Name", "Overlap", "Overlap_percentage", "Total_Hits") %in% colnames(df))) {
    stop("The dataframe must contain 'Subsample_percentage', 'Overlap', 'Total_Hits' and 'Overlap_percentage' columns.")
  }
  
  # Create the plot
  if (comparing == ""){
    p <- ggplot(df, aes(x = Name, y = Overlap_percentage)) +
      geom_point(size = 3) +
      theme_bw() +
      labs(
        title = paste("Percentage of correct Hits for", title),
        x = "",
        y = "Overlap Percentage"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
      )  
  }
  if (comparing == "subsamples") {
    
    # Extract trailing number from tags like "A25" / "subX10" / "25"
    df$Number <- as.integer(sub("^\\D*([0-9]+)$", "\\1", df$Name))
    
    # Safety check (optional)
    if (any(is.na(df$Number))) {
      stop("Some 'Name' entries do not end in a number (cannot extract Number).")
    }
    
    # Sort by numeric value so the line goes left -> right correctly
    df <- df[order(df$Number), ]
    
    if (anyDuplicated(df$Number) == 0) {
      # mbers are unique → simple plot
      p <- ggplot(df, aes(x = Number)) +
        geom_line(aes(y = Overlap_percentage), size = 1, group = 1) +
        geom_point(aes(y = Overlap_percentage), size = 3) +
        labs(
          title = paste("Percentage of correct Hits for", title),
          x = "Subsample Percentage",
          y = "Overlap Percentage"
        ) +
        theme_bw()
      
    } else {
      # Numbers NOT unique → plot mean per Number + error bars
      
      # Summarise mean + SE (swap to sd if you prefer)
      tmp <- split(df$Overlap_percentage, df$Number)
      df_sum <- data.frame(
        Number = as.integer(names(tmp)),
        mean   = sapply(tmp, mean),
        sd     = sapply(tmp, sd),
        n      = sapply(tmp, length)
      )
      
      df_sum$se <- df_sum$sd / sqrt(df_sum$n)
      df_sum$se[is.na(df_sum$se)] <- 0  # happens when n=1 (sd=NA)
      
      df_sum <- df_sum[order(df_sum$Number), ]
      
      p <- ggplot(df_sum, aes(x = Number)) +
        geom_line(aes(y = mean), size = 1, group = 1) +
        geom_point(aes(y = mean), size = 3) +
        geom_errorbar(
          aes(ymin = mean - se, ymax = mean + se),
          width = 0.2
        ) +
        labs(
          title = paste("Percentage of correct Hits with SE for", title),
          x = "Subsample Percentage",
          y = "Overlap Percentage"
        ) +
        theme_bw()
      
      # Optional: show the individual points too (nice for debugging)
      # p <- p + geom_point(data = df, aes(x = Number, y = Overlap_percentage), alpha = 0.4)
    }
  }
  
  return(p)
}

# Function for a line plot comparing total hits
# Currently depricated and replaced by venn diagramms
plot_total_hits <- function(df, title, comparing = "") {
  # Ensure the input dataframe contains the necessary columns
  if (!all(c("Name", "Overlap", "Overlap_percentage", "Total_Hits") %in% colnames(df))) {
    stop("The dataframe must contain 'Name', 'Overlap', 'Total_Hits' and 'Overlap_percentage' columns.")
  }
  if (comparing == ""){
    x_lab_text = ""
  }
  if (comparing == "subsamples"){
    x_lab_text = "Subsample Percentage"
  }
  # Create the plot
  p <- ggplot(df, aes(x = Name)) +
    geom_line(aes(y = Total_Hits), size = 1) + 
    geom_point(aes(y = Total_Hits), size = 3) +
    labs(
      title = paste("Number of all Hits (not necessarily correct ones) for",title),
      x = x_lab_text,
      y = "Total Hits"
    ) +
    theme_bw()
  
  # Return the plot
  return(p)
}

compare_df_list <- function(df_list,
                            name_list,
                            FDR_threshold = 0.05,
                            top_N = NULL,
                            comparing = "",
                            correlation_heatmap = TRUE,
                            overlap = TRUE,
                            venn_diagram = TRUE) {
  
  
  all_hits_df <- make_compare_df(df_list, name_list, FDR_threshold)
  print(all_hits_df)
  sig_up_hits_df <- make_compare_df(df_list, name_list, FDR_threshold, significance_filter = 1)
  sig_down_hits_df <- make_compare_df(df_list, name_list, FDR_threshold, significance_filter = -1)
  results_list <- list(all_hits_df, sig_up_hits_df, sig_down_hits_df)
  names(results_list) <- c("All_Hits", "Up_Hits" ,"Down_Hits")
  names(results_list) <- c("All_Hits", "Up_Hits" ,"Down_Hits")
  signif_map <- list(
    All_Hits  = NULL,
    Up_Hits   = 1,
    Down_Hits = -1
  )
  
  # Loop over all dataframes in the results list and plot them
  for (scenario_name in names(results_list)) {
    signif_filter_for_scenario <- signif_map[[scenario_name]]
    if (correlation_heatmap == TRUE){
      q <- make_correlation_heatmap(df_list,
                                    name_list,
                                    comparing = comparing,
                                    plot_title = scenario_name,
                                    significance_filter = signif_filter_for_scenario,
                                    FDR_threshold = FDR_threshold)
      print(q)      
    }
    if (overlap == TRUE){
      p <- plot_overlap(results_list[[scenario_name]],
                        scenario_name,
                        comparing = comparing)
      print(p)       
    }
    if (venn_diagram == TRUE){
      r <- make_venn_from_overlap_df(results_list[[scenario_name]],
                                     ncol = 1,
                                     comparing = comparing,
                                     main_title = scenario_name)
      grid::grid.newpage()
      grid::grid.draw(r)      
    }
  }
  return(results_list)
}


move_subsamples <- function(subsample_master_dir,
                            subsample_prefix = "",
                            output_folder = get("rds_output_folder",
                                                envir = .GlobalEnv),
                            file_info_suffix = get("file_info_suffix",
                                                   envir = .GlobalEnv)) {
  
  # Basic checks
  if (subsample_prefix != ""){
    is_alpha_only <- function(x) grepl("^[A-Za-z]+$", x)
    if (!is_alpha_only(subsample_prefix)){
      stop("subsample prefix can only contain alphatbetic characters: ",subsample_prefix)
    }
  }
  if (!dir.exists(subsample_master_dir)) {
    stop("subsample_master_dir does not exist: ", subsample_master_dir)
  }
  if (!dir.exists(output_folder)) {
    stop("rds_output_folder does not exist: ", output_folder)
  }
  if (!is.character(file_info_suffix) || length(file_info_suffix) != 1) {
    stop("file_info_suffix must be a single character string.")
  }
  
  # Find subdirectories like subsample_X
  subsample_dirs <- list.dirs(subsample_master_dir, recursive = FALSE, full.names = TRUE)
  subsample_dirs <- subsample_dirs[grepl("^subsample_", basename(subsample_dirs))]
  
  if (length(subsample_dirs) == 0) {
    message("No subdirectories named 'subsample_#' found in: ", subsample_master_dir)
    return(invisible(NULL))
  }
  
  # Build regex to match files ending in <file_info_suffix>.rds
  # Example: ".*_fileinfo\\.rds$" depending on suffix
  pattern <- paste0(".*", file_info_suffix, "\\.rds$")
  
  copied <- character(0)
  skipped <- character(0)
  
  for (sd in subsample_dirs) {
    sub_name <- basename(sd)                 # "subsample_X"
    X <- sub("^subsample_", "", sub_name)    # keep whatever is after "subsample_"
    
    rds_dir <- file.path(sd, "rds")
    if (!dir.exists(rds_dir)) {
      message("Skipping, no rds dir in: ", sd)
      next
    }
    
    files <- list.files(rds_dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0){
      message("Skipping, no files matching patterin in dir: ", rds_dir)
      next
    } 
    for (file in files) {
      base <- basename(file) # e.g. "something_suffix.rds"
      new_base <- sub("\\.rds$", paste0("_subsample_", subsample_prefix, X, ".rds"), base)
      dest <- file.path(output_folder, new_base)
      
      ok <- file.copy(from = file, to = dest, overwrite = TRUE)
      if (ok) {
        copied <- c(copied, dest)
      } else {
        skipped <- c(skipped, file)
      }
    }
  }
  
  message("Copied ", length(copied), " file(s) to: ", output_folder)
  if (length(skipped) > 0) {
    message("WARNING: Failed to copy ", length(skipped), " file(s).")
  }
  
  invisible(list(copied = copied, failed = skipped))
}


automate_subsample_comparison <- function(file_info_prefix = get("file_info_suffix",
                                                                 envir = .GlobalEnv),
                                          FDR_threshold = 0.05,
                                          correlation_heatmap = TRUE,
                                          overlap = TRUE,
                                          venn_diagram = TRUE) {
  # Get the current settings from the global environment
  if(exists("Hits_current_settings")){
    Hits_current_settings <- get("Hits_current_settings", envir = .GlobalEnv)
  } else {
    Hits_current_settings <- add_info_wrapper(file_info_prefix)
  }
  
  rds_output_folder <- get("rds_output_folder", envir = .GlobalEnv)
  
  # 1) list files (now expecting subsample_<letters><digits>.rds)
  rds_files <- list.files(
    rds_output_folder,
    pattern = paste0(".*", file_info_prefix, "_subsample_[A-Za-z]+[0-9]+\\.rds$"),
    full.names = TRUE
  )
  
  # 2) extract the combined prefix+number part (e.g. "A25" from align_UMI_tools_reads_control_median_RI_subsample_A25.rds)
  subsample_tag <- sub(
    paste0(".*", file_info_prefix, "_subsample_([A-Za-z]+[0-9]+)\\.rds$"),
    "\\1",
    basename(rds_files)
  )
  subsample_tag_unique <- unique(subsample_tag)
  # 4) split into prefixes and numbers
  prefixes <- sub("^([A-Za-z]+)[0-9]+$", "\\1", subsample_tag_unique)
  numbers  <- as.integer(sub("^[A-Za-z]+([0-9]+)$", "\\1", subsample_tag_unique))
  
  # Generate the percentage list (starting from 100 and then descending from the available numbers)
  tag_df <- data.frame(prefix = prefixes,
                       number = numbers,
                       tag = subsample_tag_unique,
                       stringsAsFactors = FALSE)
  
  tag_df <- tag_df[order(tag_df$number, decreasing = TRUE), ]
  # Generate Hits_sub_x using add_info_wrapper for each percentage in percentage_list
  df_list <- lapply(tag_df$tag, function(tag) {
    add_info_wrapper(paste0(file_info_prefix, "_subsample_", tag))
  })
  print(length(df_list))
  # Call compare_df_list with the generated df_list and percentage_list
  result <- compare_df_list(df_list = c(list(Hits_current_settings), df_list),
                            name_list = c(100, tag_df$tag),
                            FDR_threshold = FDR_threshold,
                            comparing = "subsamples",
                            top_N = NULL,
                            correlation_heatmap = correlation_heatmap,
                            overlap = overlap,
                            venn_diagram = venn_diagram)
  
  return(result)
}
automate_sublibrary_comparison <- function(file_info_suffix = get("file_info_suffix", envir = .GlobalEnv),
                                           FDR_threshold = 0.05,
                                           correlation_heatmap = TRUE,
                                           overlap = TRUE,
                                           venn_diagram = TRUE) {
  # Get the current settings from the global environment
  Hits_all_sublibs <- add_info_wrapper(file_info_suffix)
  rds_output_folder <- get("rds_output_folder", envir = .GlobalEnv)
  
  # List all files in the rds output folder
  rds_files <- list.files(rds_output_folder, pattern = paste0(".*", file_info_suffix, "_no_sublib", "_L[0-9]+\\.rds$"), full.names = TRUE)
  
  # Extract everything between file_info_suffix + "_" and ".rds"
  # e.g. "...<file_info_suffix>_no_sublib_L1.rds" -> "no_sublib_L1"
  sublib_names <- sub(
    paste0(".*", file_info_suffix, "_(.+)\\.rds$"),
    "\\1",
    basename(rds_files)
  )
  sublib_names <- unique(sublib_names)
  sublib_suffixes <- paste0(file_info_suffix, "_", sublib_names)
  # Generate Hits_df_list using add_info_wrapper for each sublib_suffixes
  Hits_df_list <- lapply(sublib_names, function(sublib_names) {
    add_info_wrapper(paste0(file_info_suffix, "_", sublib_names))
  })
  
  # Call compare_sublibraries with the generated df_list and sublib_names
  result <- compare_df_list(df_list = c(list(Hits_all_sublibs), Hits_df_list),
                            name_list = c("all_sublibs", sublib_names),
                            FDR_threshold = FDR_threshold,
                            top_N = NULL,
                            correlation_heatmap = correlation_heatmap,
                            overlap = overlap,
                            venn_diagram = venn_diagram)
  
  return(result)
}
plot_three_mean_hit_heatmap <- function(
    PA_means,
    DCA_means,
    GALNAC_means,
    all_conjugates_overlap_classified,
    gene_col = "entrez",
    value_col = "significanceZ",
    fdr_col = "FDR",
    symbol_col = "symbol",
    show_significance_stars = TRUE,
    show_clustering_tree = TRUE,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_colnames = FALSE,
    rank_significanceZ = TRUE,
    n_clusters = NULL,
    row_labels = c("PA", "DCA", "GalNAc"),
    create_png_path = NULL
) {
  
  required_cols <- c(gene_col, value_col, fdr_col)
  
  check_df <- function(df, df_name) {
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
      stop(df_name, " is missing required columns: ",
           paste(missing_cols, collapse = ", "))
    }
  }
  
  check_df(PA_means, "PA_means")
  check_df(DCA_means, "DCA_means")
  check_df(GALNAC_means, "GALNAC_means")
  
  target_genes <- unique(as.character(all_conjugates_overlap_classified[[gene_col]]))
  
  prep_one <- function(df, dataset_name, rank_significanceZ = FALSE) {
    
    df <- df %>%
      dplyr::mutate(
        gene_id = as.character(.data[[gene_col]])
      ) %>%
      dplyr::filter(gene_id %in% target_genes) %>%
      dplyr::filter(!is.na(.data[[symbol_col]]))
    
    if (rank_significanceZ) {
      if (value_col != "significanceZ") {
        stop("rank_significanceZ = TRUE currently only makes sense when value_col = 'significanceZ'.")
      }
      
      df <- df %>%
        dplyr::mutate(
          value_rank = rank(abs(significanceZ), ties.method = "average", na.last = "keep"),
          value_rank = value_rank * sign(significanceZ)
        )
    }
    
    df %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        value = dplyr::first(if (rank_significanceZ) value_rank else .data[[value_col]]),
        FDR = dplyr::first(.data[[fdr_col]]),
        .groups = "drop"
      ) %>%
      dplyr::mutate(dataset = dataset_name)
  }
  
  pa_df  <- prep_one(PA_means, row_labels[1], rank_significanceZ)
  dca_df <- prep_one(DCA_means, row_labels[2], rank_significanceZ)
  gal_df <- prep_one(GALNAC_means, row_labels[3], rank_significanceZ)
  
  combined_df <- dplyr::bind_rows(pa_df, dca_df, gal_df)
  
  if (show_colnames) {
    required_symbol_col <- c(gene_col, symbol_col)
    
    df_list_for_labels <- list(
      PA_means = PA_means,
      DCA_means = DCA_means,
      GALNAC_means = GALNAC_means
    )
    
    for (nm in names(df_list_for_labels)) {
      df_obj <- df_list_for_labels[[nm]]
      missing_cols <- setdiff(required_symbol_col, colnames(df_obj))
      if (length(missing_cols) > 0) {
        stop(nm, " is missing required columns for labels: ",
             paste(missing_cols, collapse = ", "))
      }
    }
  }
  
  gene_order <- target_genes[target_genes %in% unique(combined_df$gene_id)]
  
  # value matrix
  value_wide <- combined_df %>%
    dplyr::select(dataset, gene_id, value) %>%
    tidyr::pivot_wider(names_from = gene_id, values_from = value) %>%
    as.data.frame()
  
  rownames(value_wide) <- value_wide$dataset
  value_wide$dataset <- NULL
  value_mat <- as.matrix(value_wide)
  
  value_mat <- value_mat[, gene_order, drop = FALSE]
  mode(value_mat) <- "numeric"
  
  # optional column labels based on symbol, with NA/blank symbols replaced by entrez
  if (show_colnames) {
    gene_labels_df <- dplyr::bind_rows(
      PA_means %>% dplyr::select(dplyr::all_of(c(gene_col, symbol_col))),
      DCA_means %>% dplyr::select(dplyr::all_of(c(gene_col, symbol_col))),
      GALNAC_means %>% dplyr::select(dplyr::all_of(c(gene_col, symbol_col)))
    ) %>%
      dplyr::mutate(
        gene_id = as.character(.data[[gene_col]]),
        symbol_filled = as.character(.data[[symbol_col]]),
        symbol_filled = ifelse(
          is.na(symbol_filled) | trimws(symbol_filled) == "",
          gene_id,
          symbol_filled
        )
      ) %>%
      dplyr::filter(gene_id %in% gene_order) %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        label = dplyr::first(symbol_filled),
        .groups = "drop"
      )
    
    gene_labels <- stats::setNames(gene_labels_df$label, gene_labels_df$gene_id)
    colnames(value_mat) <- gene_labels[colnames(value_mat)]
  }
  
  # FDR matrix
  fdr_wide <- combined_df %>%
    dplyr::select(dataset, gene_id, FDR) %>%
    tidyr::pivot_wider(names_from = gene_id, values_from = FDR) %>%
    as.data.frame()
  
  rownames(fdr_wide) <- fdr_wide$dataset
  fdr_wide$dataset <- NULL
  fdr_mat <- as.matrix(fdr_wide)
  
  fdr_mat <- fdr_mat[, gene_order, drop = FALSE]
  mode(fdr_mat) <- "numeric"
  colnames(fdr_mat) <- gene_order
  
  # Count for each gene in how many datasets FDR <= 0.05
  sig_counts <- colSums(fdr_mat <= 0.05, na.rm = TRUE)
  
  if (show_significance_stars) {
    star_vec <- c("", "*", "*\n*", "*\n*\n*")[sig_counts + 1]
    display_numbers <- matrix(
      rep(star_vec, each = nrow(value_mat)),
      nrow = nrow(value_mat),
      ncol = ncol(value_mat)
    )
    rownames(display_numbers) <- rownames(value_mat)
    colnames(display_numbers) <- colnames(value_mat)
  } else {
    display_numbers <- FALSE
  }
  
  # Use a version without NA for clustering distances
  clustering_mat <- value_mat
  clustering_mat[is.na(clustering_mat)] <- 0
  
  col_dist <- stats::dist(t(clustering_mat))
  row_dist <- stats::dist(clustering_mat)
  
  max_abs <- max(abs(value_mat), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) {
    max_abs <- 1
  }
  
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
  
  p <- pheatmap::pheatmap(
    mat = value_mat,
    color = heatmap_colors,
    breaks = breaks,
    cluster_cols = cluster_cols,
    cluster_rows = cluster_rows,
    clustering_distance_cols = col_dist,
    clustering_distance_rows = row_dist,
    cutree_cols = n_clusters, 
    display_numbers = display_numbers,
    number_color = "black",
    show_colnames = show_colnames,
    show_rownames = TRUE,
    border_color = NA,
    treeheight_col = if (show_clustering_tree && cluster_cols) 50 else 0,
    treeheight_row = if (show_clustering_tree && cluster_rows) 50 else 0,
    fontsize_row = 11,
    fontsize_col = 4,
    main = paste0(value_col, " heatmap for genes in all_conjugates_overlap_classified")
  )
  
  if (!is.null(create_png_path)) {
    dir.create(dirname(create_png_path), recursive = TRUE, showWarnings = FALSE)
    
    grDevices::png(
      filename = create_png_path,
      width = 16000,
      height = 1800,
      res = 300
    )
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
    grDevices::dev.off()
  }
  
  return(list(
    plot = p,
    value_matrix = value_mat,
    fdr_matrix = fdr_mat,
    significance_counts = sig_counts,
    filtered_PA = pa_df,
    filtered_DCA = dca_df,
    filtered_GALNAC = gal_df
  ))
}
plot_three_screen_significance_umap <- function(
    PA_means,
    DCA_means,
    GALNAC_means,
    all_conjugates_overlap_classified = NULL,
    gene_col = "entrez",
    value_col = "significanceZ",
    symbol_col = "symbol",
    use_only_overlap_genes = TRUE,
    fill_missing = 0,
    color_by = NULL,
    label_points = FALSE,
    point_size = 2,
    point_alpha = 0.8,
    seed = 123,
    n_neighbors = 15,
    min_dist = 0.2,
    metric = "euclidean",
    create_png_path = NULL
) {
  
  required_cols <- c(gene_col, value_col)
  
  check_df <- function(df, df_name) {
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
      stop(df_name, " is missing required columns: ",
           paste(missing_cols, collapse = ", "))
    }
  }
  
  check_df(PA_means, "PA_means")
  check_df(DCA_means, "DCA_means")
  check_df(GALNAC_means, "GALNAC_means")
  
  prep_df <- function(df, screen_name) {
    df %>%
      dplyr::mutate(
        gene_id = as.character(.data[[gene_col]])
      ) %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        value = dplyr::first(.data[[value_col]]),
        symbol = dplyr::first(if (symbol_col %in% colnames(df)) .data[[symbol_col]] else NA),
        .groups = "drop"
      ) %>%
      dplyr::rename(!!screen_name := value)
  }
  
  pa_df  <- prep_df(PA_means, "PA")
  dca_df <- prep_df(DCA_means, "DCA")
  gal_df <- prep_df(GALNAC_means, "GalNAc")
  
  merged_df <- pa_df %>%
    dplyr::full_join(dca_df, by = c("gene_id", "symbol")) %>%
    dplyr::full_join(gal_df, by = c("gene_id", "symbol"))
  
  # fill missing symbols with gene_id
  merged_df <- merged_df %>%
    dplyr::mutate(
      symbol = as.character(symbol),
      symbol = ifelse(is.na(symbol) | trimws(symbol) == "", gene_id, symbol)
    )
  
  if (use_only_overlap_genes) {
    if (is.null(all_conjugates_overlap_classified)) {
      stop("use_only_overlap_genes = TRUE requires all_conjugates_overlap_classified.")
    }
    
    target_genes <- unique(as.character(all_conjugates_overlap_classified[[gene_col]]))
    merged_df <- merged_df %>%
      dplyr::filter(gene_id %in% target_genes)
  }
  
  # optional annotation column for coloring
  if (!is.null(color_by)) {
    if (is.null(all_conjugates_overlap_classified)) {
      stop("color_by requires all_conjugates_overlap_classified to be provided.")
    }
    if (!(color_by %in% colnames(all_conjugates_overlap_classified))) {
      stop("Column '", color_by, "' not found in all_conjugates_overlap_classified.")
    }
    
    annot_df <- all_conjugates_overlap_classified %>%
      dplyr::mutate(gene_id = as.character(.data[[gene_col]])) %>%
      dplyr::group_by(gene_id) %>%
      dplyr::summarise(
        color_group = dplyr::first(.data[[color_by]]),
        .groups = "drop"
      )
    
    merged_df <- merged_df %>%
      dplyr::left_join(annot_df, by = "gene_id")
  }
  
  umap_input <- merged_df %>%
    dplyr::select(PA, DCA, GalNAc)
  
  if (!all(vapply(umap_input, is.numeric, logical(1)))) {
    stop("The selected value columns must be numeric.")
  }
  
  umap_input <- as.data.frame(umap_input)
  umap_input[is.na(umap_input)] <- fill_missing
  
  set.seed(seed)
  umap_layout <- uwot::umap(
    umap_input,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = metric,
    verbose = FALSE,
    ret_model = FALSE
  )
  
  plot_df <- merged_df %>%
    dplyr::mutate(
      UMAP1 = umap_layout[, 1],
      UMAP2 = umap_layout[, 2]
    )
  
  if (is.null(color_by)) {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP1, y = UMAP2)) +
      ggplot2::geom_point(size = point_size, alpha = point_alpha) +
      ggplot2::theme_classic() +
      ggplot2::labs(
        title = "UMAP of genes based on significanceZ across three screens",
        x = "UMAP1",
        y = "UMAP2"
      )
  } else {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP1, y = UMAP2, color = color_group)) +
      ggplot2::geom_point(size = point_size, alpha = point_alpha) +
      ggplot2::theme_classic() +
      ggplot2::labs(
        title = paste0("UMAP of genes based on significanceZ across three screens"),
        x = "UMAP1",
        y = "UMAP2",
        color = color_by
      )
  }
  
  if (label_points) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = symbol),
      size = 3,
      max.overlaps = 50
    )
  }
  
  if (!is.null(create_png_path)) {
    dir.create(dirname(create_png_path), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
      filename = create_png_path,
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
  
  return(list(
    plot = p,
    umap_df = plot_df,
    umap_input = umap_input
  ))
}

scatter_merged_vs_input <- function(
    merged_sgRNA_df,
    count_df_long,
    normalization = c("none", "CPM"),
    group_category = "all",                 # "all" or c("entry1","entry2")
    condition_value = "input",
    sgrna_col_long = "sgRNA",
    sgrna_col_merged = "sgrna_id",
    count_col_long = "count",
    count_col_merged = "count",
    group_col = "group_category",
    condition_col = "condition",
    sample_col_long = "sample",             # optional; improves CPM
    sublib_col = "sublib",
    sublib_chosen = "all",
    transform = c("none", "log10", "log1p"),
    add_lm = TRUE,
    annotate_cor = TRUE,
    mark_outlier = TRUE,
    n_outliers = 10,
    symbol_col_merged = "symbol",           # used for labels; fall back to sgrna_id if missing/NA
    point_alpha = 0.5,
    point_size = 1,
    mark_essential = TRUE,
    essential_col_merged = "is_essential",
    essential_color = "red",
    essential_alpha = 0.9,
    essential_size = NULL,
    remove_essential = FALSE
) {
  normalization <- match.arg(normalization)
  transform <- match.arg(transform)
  
  stopifnot(
    all(c(sgrna_col_long, count_col_long, group_col, condition_col) %in% colnames(count_df_long)),
    all(c(sgrna_col_merged, count_col_merged) %in% colnames(merged_sgRNA_df))
  )
  if (!is.null(sample_col_long) && !(sample_col_long %in% colnames(count_df_long))) {
    stop("sample_col_long was provided but is not a column in count_df_long.")
  }
  
  # 1) filter to input + group_category selection
  df_in <- count_df_long %>%
    filter(.data[[condition_col]] == .env$condition_value)
  
  n_before <- nrow(df_in)
  if (!(length(group_category) == 1 && group_category == "all")) {
    df_in <- df_in %>%
      filter(.data[[group_col]] %in% .env$group_category)
  }
  n_after <- nrow(df_in)
  message("Rows before/after group filter: ", n_before, " -> ", n_after)
  
  n_before <- nrow(df_in)
  if (!(length(sublib_chosen) == 1 && sublib_chosen == "all")) {
    df_in <- df_in %>%
      filter(.data[[sublib_col]] %in% .env$sublib_chosen)
  }
  n_after <- nrow(df_in)
  message("Rows before/after sublibrary filter: ", n_before, " -> ", n_after)
  
  
  
  essential_df <- readRDS(file.path(data_dir, "dependency_df.rds")) %>%
    mutate(entrez = as.character(entrez_id)) %>%
    select(-c(gene_symbol, entrez_id)) %>%
    mutate(
      is_essential = if_any(
        c(essential_HepG2, essential_general, essential_liver),
        ~ tidyr::replace_na(.x, FALSE)
      )
    )
  
  merged_sgRNA_df <- merged_sgRNA_df %>%
    mutate(entrez = as.character(entrez)) %>%   # ensure same type for join
    left_join(essential_df, by = "entrez")
  
  # 2) normalize (optional) then mean per sgRNA
  if (normalization == "CPM") {
    if (!is.null(sample_col_long)) {
      df_in <- df_in %>%
        group_by(.data[[sample_col_long]]) %>%
        mutate(.cpm = (.data[[count_col_long]] / sum(.data[[count_col_long]], na.rm = TRUE)) * 1e6) %>%
        ungroup()
      
      in_sum <- df_in %>%
        group_by(.data[[sgrna_col_long]]) %>%
        summarise(input_mean = mean(.cpm, na.rm = TRUE), .groups = "drop")
    } else {
      total <- sum(df_in[[count_col_long]], na.rm = TRUE)
      df_in <- df_in %>% mutate(.cpm = (.data[[count_col_long]] / total) * 1e6)
      
      in_sum <- df_in %>%
        group_by(.data[[sgrna_col_long]]) %>%
        summarise(input_mean = mean(.cpm, na.rm = TRUE), .groups = "drop")
    }
    
    merged_plot_df <- merged_sgRNA_df %>%
      mutate(merged_value = (.data[[count_col_merged]] / sum(.data[[count_col_merged]], na.rm = TRUE)) * 1e6)
  } else {
    in_sum <- df_in %>%
      group_by(.data[[sgrna_col_long]]) %>%
      summarise(input_mean = mean(.data[[count_col_long]], na.rm = TRUE), .groups = "drop")
    
    merged_plot_df <- merged_sgRNA_df %>%
      mutate(merged_value = .data[[count_col_merged]])
  }
  
  # Build label: prefer symbol, fall back to sgrna_id; if symbol col missing, use sgrna_id
  has_symbol <- !is.null(symbol_col_merged) && (symbol_col_merged %in% colnames(merged_plot_df))
  merged_plot_df <- merged_plot_df %>%
    mutate(.label = if (has_symbol) {
      dplyr::coalesce(as.character(.data[[symbol_col_merged]]),
                      as.character(.data[[sgrna_col_merged]]))
    } else {
      as.character(.data[[sgrna_col_merged]])
    })
  
  # 3) join (match sgRNA ids) — use all_of() in tidyselect contexts
  in_sum2 <- in_sum %>%
    rename(sgrna_id_join = all_of(sgrna_col_long))
  
  merged_min <- merged_plot_df %>%
    select(all_of(c(sgrna_col_merged, "merged_value", ".label")),
           any_of(essential_col_merged)) %>%                # safe if column missing
    rename(sgrna_id_join = all_of(sgrna_col_merged)) %>%
    mutate(.essential = dplyr::coalesce(.data[[essential_col_merged]], FALSE))
  
  
  joined <- in_sum2 %>% inner_join(merged_min, by = "sgrna_id_join")
  
  if (remove_essential){
    joined <- joined %>% filter(!is_essential)
  }
  
  if (nrow(joined) == 0) {
    stop("Join produced 0 rows. Check that sgRNA IDs match between the two data frames.")
  }
  
  # Fit LM once (needed for outlier distances and optionally for drawing the line)
  fit <- lm(input_mean ~ merged_value, data = joined)
  b <- unname(coef(fit)[["(Intercept)"]])
  m <- unname(coef(fit)[["merged_value"]])
  
  # Perpendicular distance from point to line y = m x + b:
  # distance = | -m*x + y - b | / sqrt(m^2 + 1)
  if (mark_outlier) {
    joined <- joined %>%
      mutate(.dist = abs((-m) * merged_value + input_mean - b) / sqrt(m^2 + 1))
    
    n_keep <- min(n_outliers, nrow(joined))
    
    outliers <- joined %>%
      arrange(desc(.dist)) %>%
      slice_head(n = n_keep)
  }
  
  # 4) plot
  p <- ggplot(joined, aes(x = merged_value, y = input_mean)) +
    geom_point(alpha = point_alpha, size = point_size) +
    labs(
      x = paste0("plasmid count", ifelse(normalization == "CPM", " (CPM)", "")),
      y = paste0("mean input ", count_col_long, ifelse(normalization == "CPM", " (CPM)", "")),
      title = paste(
        "Plasmid counts vs mean input counts for",
        paste(group_category, collapse = ", "),
        "Norm:", normalization
      )
    ) +
    theme_classic()
  
  
  if (transform == "log10") {
    p <- p + scale_x_log10() + scale_y_log10()
  } else if (transform == "log1p") {
    p <- p + scale_x_continuous(trans = "log1p") + scale_y_continuous(trans = "log1p")
  }
  
  if (annotate_cor) {
    r <- suppressWarnings(cor(joined$merged_value, joined$input_mean, use = "pairwise.complete.obs"))
    p <- p + annotate(
      "text",
      x = Inf, y = Inf, hjust = 1.1, vjust = 1.2,
      label = sprintf("Pearson r = %.3f\nn = %d", r, nrow(joined))
    )
  }
  
  # label top outliers
  if (mark_outlier && nrow(outliers) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = outliers,
        aes(label = .label),
        max.overlaps = Inf
      )
    } else {
      p <- p + geom_text(
        data = outliers,
        aes(label = .label),
        vjust = -0.5,
        check_overlap = TRUE
      )
    }
  }
  if (mark_essential) {
    if (is.null(essential_size)) essential_size <- point_size
    p <- p + geom_point(
      data = dplyr::filter(joined, .essential),
      aes(x = merged_value, y = input_mean),
      color = essential_color,
      alpha = essential_alpha,
      size = essential_size
    )
  }
  # draw LM line (using the same fit we used for distances)
  if (add_lm) {
    p <- p + geom_abline(intercept = b, slope = m, color = "blue")
  }
  
  print(p)
  list(plot = p, data = joined)
}

waterfall_roche_compare <- function(Hits_df,
                                    comp_df,
                                    signif_lines = TRUE,
                                    signif_level = 0.05,
                                    box_padding = 0.8,
                                    no_text = FALSE,
                                    mark_special = NULL,
                                    rho_neg_thresh = 0.0001,
                                    rho_pos_thresh = 0.0001,
                                    lfc_neg_thresh = 1.5,
                                    lfc_pos_thresh = 1.5) {
  
  # Create rank column
  Hits_df <- Hits_df %>%
    mutate(rank = rank(significanceZ, ties.method = "first"))
  signif_df <- Hits_df
  # Define the Limits for annotation boxes
  y_min   <- min(Hits_df$significanceZ, na.rm = TRUE)
  y_max   <- max(Hits_df$significanceZ, na.rm = TRUE)
  y_range <- y_max - y_min
  
  # Box geometry (relative to overall range)
  box_height <- 0.06 * y_range   # 8% of range high
  box_gap    <- 0.02 * y_range   # 3% of range between boxes
  
  # Centers of the two boxes (on y-axis)
  y_center1 <- y_min           # lower box
  y_center2 <- y_center1 + box_height + box_gap    # upper box
  
  # x range (can reuse for both boxes)
  x_min_box <- max(Hits_df$rank) - 3000
  x_max_box <- max(Hits_df$rank) + 500
  
  # Filter for Rho_neg and Rho_pos hits from comp_df
  rho_neg_hits_df <- comp_df %>%
    filter(Rho_neg <= rho_neg_thresh & 2^LFC_neg >= lfc_neg_thresh) %>%
    rename(symbol_rho = symbol)  # Rename the symbol to avoid conflict
  
  rho_pos_hits_df <- comp_df %>%
    filter(Rho_pos <= rho_pos_thresh & 2^LFC_pos >= lfc_pos_thresh) %>%
    rename(symbol_rho = symbol)  # Rename the symbol to avoid conflict
  
  # Filter only the rows in Hits_df that are in rho_neg_hits_df and rho_pos_hits_df
  rho_neg_hits_df <- inner_join(Hits_df, rho_neg_hits_df, by = "entrez")
  rho_pos_hits_df <- inner_join(Hits_df, rho_pos_hits_df, by = "entrez")
  
  # Start base plot
  p <- ggplot(Hits_df, aes(x = rank, y = significanceZ)) +
    geom_point(color = "lightgrey", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.3) +
    theme_bw() +
    labs(
      x = "Gene rank",
      y = "Significance Z-score",
      title = "Significance Z-score by gene rank"
    )
  
  # Mark Rho_neg hits
  if (nrow(rho_neg_hits_df) > 0) {
    p <- p + 
      geom_point(
        data = rho_neg_hits_df,
        aes(x = rank, y = significanceZ),
        color = "cyan", size = 2
      ) +
      annotate("rect",
               xmin = x_min_box,
               xmax = x_max_box,
               ymin = y_center2 - box_height / 2,
               ymax = y_center2 + box_height / 2,
               fill = "white", color = "black", size = 0.3) +
      annotate("point",
               x = x_min_box + 200,
               y = y_center2,
               color = "cyan", size = 2) +
      annotate("text",
               x = x_min_box + 400,
               y = y_center2,
               label = "Rho Neg hits", hjust = 0, size = 3)
  }
  
  # Mark Rho_pos hits
  if (nrow(rho_pos_hits_df) > 0) {
    p <- p + 
      geom_point(
        data = rho_pos_hits_df,
        aes(x = rank, y = significanceZ),
        color = "orange", size = 2
      ) +
      annotate("rect",
               xmin = x_min_box,
               xmax = x_max_box,
               ymin = y_center1 - box_height / 2,
               ymax = y_center1 + box_height / 2,
               fill = "white", color = "black", size = 0.3) +
      annotate("point",
               x = x_min_box + 200,
               y = y_center1,
               color = "orange", size = 2) +
      annotate("text",
               x = x_min_box + 400,
               y = y_center1,
               label = "Rho Pos hits", hjust = 0, size = 3)
  }
  
  # Add lines if requested (now based on input threshold for significance)
  
  if (isTRUE(signif_lines)) {
    
    ## Lowest positive significanceZ (closest to zero, > 0)
    pos_df <- signif_df %>% filter(significanceZ > 0)
    if (nrow(pos_df) > 0) {
      y_pos <- min(pos_df$significanceZ, na.rm = TRUE)
      p <- p +
        geom_hline(
          yintercept = y_pos,
          linetype = "dashed",
          color = "red",
          size = 0.3
        )
    }
    
    ## Negative significanceZ with smallest absolute value (closest to zero)
    neg_df <- signif_df %>% filter(significanceZ < 0)
    if (nrow(neg_df) > 0) {
      # largest negative value, e.g. -1 is "closer to zero" than -3
      y_neg <- max(neg_df$significanceZ, na.rm = TRUE)
      p <- p +
        geom_hline(
          yintercept = y_neg,
          linetype = "dashed",
          color = "red",
          size = 0.3
        )
    }
  }
  
  
  # Mark special symbols (can be multiple)
  if (!is.null(mark_special) && length(mark_special) > 0 && "symbol" %in% names(Hits_df)) {
    special_df <- Hits_df %>% filter(symbol %in% mark_special)
    
    if (nrow(special_df) > 0) {
      p <- p +
        geom_point(
          data = special_df,
          aes(x = rank, y = significanceZ),
          shape = 17,
          color = "orange",
          size = 4
        ) +
        geom_text_repel(
          data = special_df,
          aes(x = rank, y = significanceZ, label = symbol),
          color = "orange",
          size = 4,
          box.padding = 1,
          nudge_x = ifelse(special_df$significanceZ > 0, -500, 500),
          nudge_y = ifelse(special_df$significanceZ > 0, y_range * (0.05), y_range * (-0.05)),
          force = 4
        )
    }
  }
  
  return(p)
}

make_cor_heatmap <- function(df, name = ""){
  # Select only the numeric columns (exclude non-numeric columns like symbol, entrez, etc.)
  cor_df <- df %>% select(-c(symbol, entrez, numGuides, seq, sgRNA))
  
  # Compute the Pearson correlation matrix
  cor_matrix <- cor(cor_df, method = "pearson")
  
  # Reshape the correlation matrix for ggplot2
  melted_cor_matrix <- reshape::melt(cor_matrix)
  
  # Ensure the variable names are correctly assigned for ggplot
  colnames(melted_cor_matrix) <- c("Var1", "Var2", "value")  # Rename columns
  
  # Create a heatmap of the correlation matrix
  p <- ggplot(melted_cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", high = "green", mid = "white", midpoint = 0) +
    theme_bw() +
    labs(title = paste("Pearson Correlation Heatmap", name)) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)) 
  
  print(p)
  return(cor_matrix)
}




build_maude_guide_z_rank_df <- function(
    merged_sgRNA_df,
    DCA_dir,
    PA_dir,
    GAL_dir,
    gene_col = NULL
) {
  
  mean_na <- function(x) {
    x <- suppressWarnings(as.numeric(as.character(x)))
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }
  
  # Detect sgRNA ID column
  sgrna_col <- if ("sgrna_id" %in% colnames(merged_sgRNA_df)) {
    "sgrna_id"
  } else if ("sgRNA" %in% colnames(merged_sgRNA_df)) {
    "sgRNA"
  } else {
    stop("merged_sgRNA_df must contain either `sgrna_id` or `sgRNA`.")
  }
  
  if (!"entrez" %in% colnames(merged_sgRNA_df)) {
    stop("merged_sgRNA_df must contain an `entrez` column.")
  }
  
  # Use explicit gene_col if supplied, otherwise prefer gene/symbol, else entrez
  if (is.null(gene_col)) {
    candidate_gene_cols <- c("gene", "symbol", "gene_symbol", "target", "entrez")
    gene_col <- candidate_gene_cols[candidate_gene_cols %in% colnames(merged_sgRNA_df)][1]
  }
  
  if (is.na(gene_col) || !gene_col %in% colnames(merged_sgRNA_df)) {
    stop("Could not infer gene column. Please provide `gene_col = ...`.")
  }
  
  # Same control handling as before
  merged_sgRNA_df <- merged_sgRNA_df %>%
    mutate(
      sgrna_id_tmp = as.character(.data[[sgrna_col]]),
      entrez = if_else(
        is.na(.data$entrez),
        sub("^CONTROL_C_(.*)_\\d+$", "\\1", sgrna_id_tmp),
        as.character(.data$entrez)
      )
    )
  
  if (!"seq" %in% colnames(merged_sgRNA_df)) {
    stop("merged_sgRNA_df must contain a `seq` column.")
  }
  
  sgRNA_to_gene <- merged_sgRNA_df %>%
    transmute(
      sgRNA = as.character(.data[[sgrna_col]]),
      entrez = as.character(.data$entrez),
      gene = as.character(.data[[gene_col]]),
      seq = as.character(.data$seq)
    ) %>%
    mutate(
      gene = if_else(is.na(gene) | gene == "", entrez, gene)
    ) %>%
    distinct(sgRNA, .keep_all = TRUE)
  
  summarise_maude_z_dir <- function(input_dir, origin) {
    
    rds_dir <- file.path(input_dir, "rds")
    
    files <- list.files(
      rds_dir,
      pattern = "^MAUDE_guide_stats.*\\.rds$",
      full.names = TRUE
    )
    
    if (length(files) == 0) {
      stop("No MAUDE_guide_stats*.rds files found in: ", rds_dir)
    }
    
    lapply(files, function(f) {
      x <- readRDS(f)
      
      required_cols <- c("sgRNA", "Z")
      missing_cols <- setdiff(required_cols, colnames(x))
      
      if (length(missing_cols) > 0) {
        stop(
          "File is missing required columns: ",
          paste(missing_cols, collapse = ", "),
          "\nFile: ",
          f
        )
      }
      
      x %>%
        transmute(
          sgRNA = as.character(.data$sgRNA),
          Z = suppressWarnings(as.numeric(as.character(.data$Z))),
          source_file = basename(f)
        )
    }) %>%
      bind_rows() %>%
      group_by(sgRNA) %>%
      summarise(
        mean_Z = mean_na(Z),
        n_Z = sum(!is.na(Z)),
        n_files = n_distinct(source_file),
        .groups = "drop"
      ) %>%
      mutate(origin = origin, .before = 1)
  }
  
  PA_df  <- summarise_maude_z_dir(PA_dir,  "PA")
  DCA_df <- summarise_maude_z_dir(DCA_dir, "DCA")
  GAL_df <- summarise_maude_z_dir(GAL_dir, "GalNAc")
  
  ranked_df <- bind_rows(PA_df, DCA_df, GAL_df) %>%
    left_join(sgRNA_to_gene, by = "sgRNA") %>%
    mutate(
      abs_mean_Z = abs(mean_Z),
      Z_direction = case_when(
        mean_Z > 0 ~ "positive",
        mean_Z < 0 ~ "negative",
        mean_Z == 0 ~ "zero",
        TRUE ~ NA_character_
      )
    ) %>%
    group_by(origin, gene) %>%
    mutate(
      rank = dense_rank(desc(abs_mean_Z))
    ) %>%
    ungroup() %>%
    arrange(origin, gene, rank, sgRNA)
  
  ranked_df
}


qq_rectangle_plot <- function(
    align_rds_fpath,
    bcwithqc_rds_fpath,
    title = "QQ Rectangle Plot",
    FDR_thresh = 0.05,
    quantile_probs = c(0, 0.05, 0.5, 0.95, 1),
    alpha_both_sig = 0.85,
    alpha_other = 0.15,
    point_size = 1.2,
    return_data = FALSE,
    rank_gamma = 2,
    grey_frac = 0.08,
    save_fpath = NULL
) {
  
  prep_df <- function(fpath) {
    readRDS(fpath) %>%
      select(symbol, entrez, significanceZ, FDR) %>%
      filter(!is.na(entrez), !is.na(significanceZ)) %>%
      arrange(FDR, desc(abs(significanceZ))) %>%
      distinct(entrez, .keep_all = TRUE) %>%
      mutate(
        rank_sigZ = rank(significanceZ, ties.method = "average")
      )
  }
  
  get_fdr_cutoffs <- function(df) {
    sig_df <- df %>%
      filter(!is.na(FDR), FDR <= FDR_thresh)
    
    neg_ranks <- sig_df %>%
      filter(significanceZ < 0) %>%
      pull(rank_sigZ)
    
    pos_ranks <- sig_df %>%
      filter(significanceZ > 0) %>%
      pull(rank_sigZ)
    
    neg_cutoff <- if (length(neg_ranks) > 0) max(neg_ranks, na.rm = TRUE) else NA_real_
    pos_cutoff <- if (length(pos_ranks) > 0) min(pos_ranks, na.rm = TRUE) else NA_real_
    
    c(negative = neg_cutoff, positive = pos_cutoff)
  }
  
  make_quantile_axis <- function(n, probs) {
    breaks <- round(1 + probs * (n - 1))
    
    tibble(
      prob = probs,
      breaks = breaks,
      labels = paste0(probs * 100, "%")
    ) %>%
      distinct(breaks, .keep_all = TRUE)
  }
  
  centered_rank_transform <- function(rank, n, gamma = 2) {
    mid <- (n + 1) / 2
    sign(rank - mid) * abs(rank - mid)^gamma
  }
  
  align_df <- prep_df(align_rds_fpath) %>%
    rename(
      symbol_align = symbol,
      significanceZ_align = significanceZ,
      FDR_align = FDR,
      rank_align = rank_sigZ
    )
  
  bcwithqc_df <- prep_df(bcwithqc_rds_fpath) %>%
    rename(
      symbol_bcwithqc = symbol,
      significanceZ_bcwithqc = significanceZ,
      FDR_bcwithqc = FDR,
      rank_bcwithqc = rank_sigZ
    )
  
  align_cutoffs <- get_fdr_cutoffs(align_df %>%
                                     rename(
                                       symbol = symbol_align,
                                       significanceZ = significanceZ_align,
                                       FDR = FDR_align,
                                       rank_sigZ = rank_align
                                     ))
  
  bcwithqc_cutoffs <- get_fdr_cutoffs(bcwithqc_df %>%
                                        rename(
                                          symbol = symbol_bcwithqc,
                                          significanceZ = significanceZ_bcwithqc,
                                          FDR = FDR_bcwithqc,
                                          rank_sigZ = rank_bcwithqc
                                        ))
  
  n_align <- nrow(align_df)
  n_bcwithqc <- nrow(bcwithqc_df)
  
  plot_df <- inner_join(
    align_df,
    bcwithqc_df,
    by = "entrez"
  ) %>%
    mutate(
      symbol = coalesce(symbol_align, symbol_bcwithqc),
      delta_significanceZ = significanceZ_bcwithqc - significanceZ_align,
      both_FDR_sig = coalesce(FDR_align <= FDR_thresh, FALSE) |
        coalesce(FDR_bcwithqc <= FDR_thresh, FALSE),
      point_alpha = if_else(both_FDR_sig, alpha_both_sig, alpha_other),
      rank_align_trans = centered_rank_transform(rank_align, n_align, gamma = rank_gamma),
      rank_bcwithqc_trans = centered_rank_transform(rank_bcwithqc, n_bcwithqc, gamma = rank_gamma)
    )
  
  x_axis <- make_quantile_axis(n_align, quantile_probs)
  y_axis <- make_quantile_axis(n_bcwithqc, quantile_probs)
  
  x_axis <- x_axis %>%
    mutate(breaks_trans = centered_rank_transform(breaks, n_align, gamma = rank_gamma))
  
  y_axis <- y_axis %>%
    mutate(breaks_trans = centered_rank_transform(breaks, n_bcwithqc, gamma = rank_gamma))
  
  max_abs_delta <- max(abs(plot_df$delta_significanceZ), na.rm = TRUE)
  if (!is.finite(max_abs_delta) || max_abs_delta == 0) {
    max_abs_delta <- 1
  }
  
  x_lines <- align_cutoffs[!is.na(align_cutoffs)]
  y_lines <- bcwithqc_cutoffs[!is.na(bcwithqc_cutoffs)]
  
  x_lines_trans <- centered_rank_transform(x_lines, n_align, gamma = rank_gamma)
  y_lines_trans <- centered_rank_transform(y_lines, n_bcwithqc, gamma = rank_gamma)

  grey_half_width <- max_abs_delta * grey_frac
  
  p <- ggplot(
    plot_df,
    aes(
      x = rank_align_trans,
      y = rank_bcwithqc_trans,
      colour = delta_significanceZ,
      alpha = point_alpha
    )
  ) +
    geom_point(size = point_size) +
    
    geom_vline(
      xintercept = x_lines_trans,
      colour = "red",
      linetype = "dotted",
      linewidth = 0.7
    ) +
    geom_hline(
      yintercept = y_lines_trans,
      colour = "red",
      linetype = "dotted",
      linewidth = 0.7
    ) +
    scale_x_continuous(
      breaks = x_axis$breaks_trans,
      labels = x_axis$labels,
      expand = expansion(mult = 0.02)
    ) +
    scale_y_continuous(
      breaks = y_axis$breaks_trans,
      labels = y_axis$labels,
      expand = expansion(mult = 0.02)
    ) +
  
  scale_colour_gradientn(
    colours = c(
      "#2B6CB0",  # strong blue
      "#2B6CB0",  # keep blue until close to center
      "#48494B",  # narrow grey center
      "#E69F00",  # quickly shift to orange
      "#E69F00"   # strong orange
    ),
    values = scales::rescale(
      c(
        -max_abs_delta,
        -grey_half_width,
        0,
        grey_half_width,
        max_abs_delta
      ),
      from = c(-max_abs_delta, max_abs_delta)
    ),
    limits = c(-max_abs_delta, max_abs_delta),
    name = expression(Delta * " significanceZ\nbcwithqc - align")
  ) +
    scale_alpha_identity() +
    
    labs(
      title = title,
      x = "align rank by significanceZ [quantile / rank]",
      y = "bcwithqc rank by significanceZ [quantile / rank]"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9)
    )

  if (!is.null(save_fpath)) {
    dir.create(dirname(save_fpath), recursive = TRUE, showWarnings = FALSE)
    
    ggsave(
      filename = save_fpath,
      plot = p,
      width = 10,
      height = 6,
      units = "in",
      dpi = 300,
      bg = "transparent"
    )
  }
  return(p)
}

make_combined_delta_df <- function(
    align_rds_fpath,
    bcwithqc_rds_fpath,
    FDR_thresh = 0.05,
    only_one_screen_sig = TRUE
) {
  
  prep_df <- function(fpath) {
    readRDS(fpath) %>%
      select(symbol, entrez, significanceZ, FDR) %>%
      filter(!is.na(entrez), !is.na(significanceZ)) %>%
      arrange(FDR, desc(abs(significanceZ))) %>%
      distinct(entrez, .keep_all = TRUE) %>%
      mutate(
        rank_sigZ = rank(significanceZ, ties.method = "average")
      )
  }
  
  align_df <- prep_df(align_rds_fpath) %>%
    rename(
      symbol = symbol,
      significanceZ_align = significanceZ,
      FDR_align = FDR,
      rank_align = rank_sigZ
    )
  
  bcwithqc_df <- prep_df(bcwithqc_rds_fpath) %>%
    rename(
      significanceZ_bcwithqc = significanceZ,
      FDR_bcwithqc = FDR,
      rank_bcwithqc = rank_sigZ
    ) %>%
    select(-symbol)
  
  combined_delta <- inner_join(
    align_df,
    bcwithqc_df,
    by = "entrez"
  ) %>%
    mutate(
      deltaZ = significanceZ_bcwithqc - significanceZ_align,
      
      sig_align = coalesce(FDR_align <= FDR_thresh, FALSE),
      sig_bcwithqc = coalesce(FDR_bcwithqc <= FDR_thresh, FALSE),
      
      sig_only_one_screen = xor(sig_align, sig_bcwithqc),
      
      significant_screen = case_when(
        sig_align & !sig_bcwithqc ~ "align_only",
        !sig_align & sig_bcwithqc ~ "bcwithqc_only",
        sig_align & sig_bcwithqc ~ "both",
        TRUE ~ "neither"
      ),
      
      deltaZ_direction = case_when(
        deltaZ > 0 ~ "higher_in_bcwithqc",
        deltaZ < 0 ~ "higher_in_align",
        TRUE ~ "same"
      )
    ) %>%
    select(
      entrez,
      symbol,
      rank_align,
      rank_bcwithqc,
      significanceZ_align,
      significanceZ_bcwithqc,
      deltaZ,
      FDR_align,
      FDR_bcwithqc,
      sig_align,
      sig_bcwithqc,
      sig_only_one_screen,
      significant_screen,
      deltaZ_direction
    )
  
  if (only_one_screen_sig) {
    combined_delta <- combined_delta %>%
      filter(sig_only_one_screen)
  }
  
  return(combined_delta)
}
  
  
  
  
  
  