# One manifest row represents one FASTQ file.
# Single-end sample: one row with read = NA.
# Paired-end sample: two rows with the same pipeline_name and read = R1/R2.

validate_fastq_manifest_layout <- function(manifest) {
  required_cols <- c("pipeline_name", "read", "fastq_id")
  missing_cols <- setdiff(required_cols, colnames(manifest))
  
  if (length(missing_cols) > 0) {
    stop_log(
      "Manifest is missing columns required to validate FASTQ grouping: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  groups <- split(seq_len(nrow(manifest)), manifest$pipeline_name)
  invalid_groups <- character(0)
  
  for (pipeline_name in names(groups)) {
    idx <- groups[[pipeline_name]]
    reads <- manifest$read[idx]
    
    valid_single_end <-
      length(idx) == 1 &&
      (is.na(reads[1]) || !nzchar(reads[1]))
    
    valid_paired_end <-
      length(idx) == 2 &&
      !any(is.na(reads)) &&
      !anyDuplicated(reads) &&
      setequal(reads, c("R1", "R2"))
    
    if (!valid_single_end && !valid_paired_end) {
      read_description <- paste(
        ifelse(is.na(reads) | !nzchar(reads), "<single-end>", reads),
        collapse = ", "
      )
      
      invalid_groups <- c(
        invalid_groups,
        paste0(
          pipeline_name,
          ": ", length(idx), " file(s); read values = ", read_description
        )
      )
    }
  }
  
  if (length(invalid_groups) > 0) {
    stop_log(
      "Invalid FASTQ grouping in the manifest. Each pipeline sample must have either:\n",
      "  - exactly one row with an empty `read` value (single-end), or\n",
      "  - exactly two rows with `read` values R1 and R2 (paired-end).\n\n",
      paste0("  - ", invalid_groups, collapse = "\n")
    )
  }
  
  if (anyDuplicated(manifest$fastq_id)) {
    duplicated_ids <- unique(manifest$fastq_id[duplicated(manifest$fastq_id)])
    stop_log(
      "Manifest contains duplicate `fastq_id` values:\n",
      paste0("  - ", duplicated_ids, collapse = "\n")
    )
  }
  
  invisible(TRUE)
}


prepare_fastq_inputs <- function(
    fastq_dir,
    fastq_name_table_file_path = NULL,
    output_symlink_dir,
    manifest_output_path = NULL,
    strict_file_match = TRUE,
    overwrite_symlinks = TRUE
) {
  # -----------------------------
  # Helper functions
  # -----------------------------
  
  normalize_optional_path <- function(x) {
    if (is.null(x)) return(NULL)
    if (length(x) == 0) return(NULL)
    if (is.na(x)) return(NULL)
    
    x <- trimws(as.character(x))
    
    if (x %in% c("", "NULL", "null", "NA", "None", "none")) {
      return(NULL)
    }
    
    x
  }
  
  get_fastq_extension <- function(x) {
    # Supported:
    # .txt, .txt.gz, .txt.gzip
    # .fq, .fq.gz, .fq.gzip
    # .fastq, .fastq.gz, .fastq.gzip
    
    m <- regexpr(
      "\\.(txt|fq|fastq)(\\.gz|\\.gzip)?$",
      basename(x),
      ignore.case = TRUE
    )
    
    ifelse(m == -1, NA_character_, regmatches(basename(x), m))
  }
  
  strip_fastq_extension <- function(x) {
    sub(
      "\\.(txt|fq|fastq)(\\.gz|\\.gzip)?$",
      "",
      basename(x),
      ignore.case = TRUE
    )
  }
  
  normalize_read <- function(x) {
    x <- trimws(toupper(as.character(x)))
    x[x %in% c("", "NA", "N/A", "NULL", "NONE", "SE", "SINGLE")] <- NA_character_
    x[x %in% c("1", "READ1", "READ_1")] <- "R1"
    x[x %in% c("2", "READ2", "READ_2")] <- "R2"
    x
  }
  
  validate_bin <- function(x) {
    !is.na(x) & x %in% c("I", "L", "U")
  }
  
  validate_sublibrary <- function(x) {
    !is.na(x) & grepl("^L[0-9]+$", x)
  }
  
  validate_sample <- function(x) {
    !is.na(x) & nzchar(x) & !grepl("_", x)
  }
  
  validate_read <- function(x) {
    is.na(x) | x %in% c("R1", "R2")
  }
  
  # Accepted names without an xlsx parser:
  #   I_L1_sample.fastq.gz
  #   I_L1_sample_R1.fastq.gz
  #   I_L1_sample_R2.fastq.gz
  validate_input_name <- function(x) {
    grepl("^(I|L|U)_(L[0-9]+)_[^_]+(_R[12])?$", x)
  }
  
  parse_input_name <- function(input_name) {
    parsed <- lapply(strsplit(input_name, "_", fixed = TRUE), function(parts) {
      if (length(parts) == 3) {
        data.frame(
          bin = parts[1],
          sublibrary = parts[2],
          sample = parts[3],
          read = NA_character_,
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          bin = parts[1],
          sublibrary = parts[2],
          sample = parts[3],
          read = parts[4],
          stringsAsFactors = FALSE
        )
      }
    })
    
    do.call(rbind, parsed)
  }
  
  make_error_list <- function(x, bullet = "  - ") {
    paste0(bullet, x, collapse = "\n")
  }
  
  # -----------------------------
  # Normalize inputs
  # -----------------------------
  
  fastq_name_table_file_path <- normalize_optional_path(
    fastq_name_table_file_path
  )
  
  if (missing(fastq_dir) || is.null(fastq_dir) || !nzchar(fastq_dir)) {
    stop_log("`fastq_dir` must be provided.")
  }
  
  if (!dir.exists(fastq_dir)) {
    stop_log(
      "Input FASTQ directory does not exist:\n  ",
      fastq_dir
    )
  }
  
  if (missing(output_symlink_dir) || is.null(output_symlink_dir) || !nzchar(output_symlink_dir)) {
    stop_log("`output_symlink_dir` must be provided.")
  }
  
  # -----------------------------
  # Detect input files
  # -----------------------------
  
  all_files <- list.files(
    fastq_dir,
    full.names = TRUE,
    recursive = FALSE,
    include.dirs = FALSE
  )
  
  all_files <- all_files[file.info(all_files)$isdir == FALSE]
  
  supported_ext <- get_fastq_extension(all_files)
  input_files <- all_files[!is.na(supported_ext)]
  
  if (length(input_files) == 0) {
    stop_log(
      "No supported input files were found in:\n  ",
      fastq_dir,
      "\n\nSupported extensions are:\n",
      "  - .txt\n",
      "  - .txt.gz\n",
      "  - .txt.gzip\n",
      "  - .fq\n",
      "  - .fq.gz\n",
      "  - .fq.gzip\n",
      "  - .fastq\n",
      "  - .fastq.gz\n",
      "  - .fastq.gzip"
    )
  }
  
  input_file_basename <- basename(input_files)
  input_file_ext <- get_fastq_extension(input_files)
  input_file_stem <- strip_fastq_extension(input_files)
  
  if (anyDuplicated(input_file_basename)) {
    duplicated_files <- unique(input_file_basename[duplicated(input_file_basename)])
    stop_log(
      "Duplicate input file names were found. File basenames must be unique:\n",
      make_error_list(duplicated_files)
    )
  }
  
  # -----------------------------
  # Case 1: xlsx parser provided
  # -----------------------------
  log_info("xlxs parser file with: filename <-> pipline information association provided.\n",
           "Proceeding...")
  
  if (!is.null(fastq_name_table_file_path)) {
    if (!file.exists(fastq_name_table_file_path)) {
      stop_log(
        "The provided input FASTQ name parser file does not exist:\n  ",
        fastq_name_table_file_path
      )
    }
    
    if (!grepl("\\.xlsx$", fastq_name_table_file_path, ignore.case = TRUE)) {
      stop_log(
        "The input FASTQ name parser must be an .xlsx file:\n  ",
        fastq_name_table_file_path
      )
    }
    
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop_log(
        "Package `readxl` is required to read the input FASTQ name parser .xlsx file.\n",
        "Install it with:\n  install.packages(\"readxl\")"
      )
    }
    
    available_sheets <- readxl::excel_sheets(fastq_name_table_file_path)
    
    if (!"fastq_names" %in% available_sheets) {
      stop_log(
        "The input FASTQ name parser xlsx must contain a sheet named exactly `fastq_names`.\n\n",
        "Available sheets are:\n",
        make_error_list(available_sheets)
      )
    }
    
    fastq_names_sheet_index <- which(available_sheets == "fastq_names")
    
    if (fastq_names_sheet_index != 2) {
      stop_log(
        "The sheet `fastq_names` must be the second sheet in the xlsx file.\n\n",
        "Current sheet order is:\n",
        paste0(seq_along(available_sheets), ". ", available_sheets, collapse = "\n")
      )
    }
    
    parser <- readxl::read_excel(
      fastq_name_table_file_path,
      sheet = "fastq_names"
    )
    
    parser <- as.data.frame(parser, stringsAsFactors = FALSE)
    
    required_cols <- c("original_file", "bin", "sublibrary", "sample")
    missing_cols <- setdiff(required_cols, colnames(parser))
    
    if (length(missing_cols) > 0) {
      stop_log(
        "The parser xlsx is missing required columns:\n",
        make_error_list(missing_cols),
        "\n\nRequired columns are:\n",
        make_error_list(required_cols),
        "\n\nFor paired-end samples, also add a `read` column containing R1 or R2."
      )
    }
    
    # `read` is optional for backward-compatible single-end parser sheets.
    if (!"read" %in% colnames(parser)) {
      parser$read <- NA_character_
    }
    
    parser <- parser[, c(required_cols, "read"), drop = FALSE]
    
    parser$original_file <- trimws(as.character(parser$original_file))
    parser$bin <- trimws(as.character(parser$bin))
    parser$sublibrary <- trimws(as.character(parser$sublibrary))
    parser$sample <- trimws(as.character(parser$sample))
    parser$read <- normalize_read(parser$read)
    
    if (any(is.na(parser$original_file) | parser$original_file == "")) {
      stop_log("The parser xlsx contains empty values in `original_file`.")
    }
    
    if (anyDuplicated(parser$original_file)) {
      duplicated_original <- unique(parser$original_file[duplicated(parser$original_file)])
      stop_log(
        "The parser xlsx contains duplicate `original_file` entries:\n",
        make_error_list(duplicated_original)
      )
    }
    
    invalid_bin <- parser$original_file[!validate_bin(parser$bin)]
    invalid_sublibrary <- parser$original_file[!validate_sublibrary(parser$sublibrary)]
    invalid_sample <- parser$original_file[!validate_sample(parser$sample)]
    invalid_read <- parser$original_file[!validate_read(parser$read)]
    
    if (length(invalid_bin) > 0) {
      stop_log(
        "Invalid `bin` values in parser xlsx.\n\n",
        "Allowed bin values are: I, L, U\n\n",
        "Affected files:\n",
        make_error_list(invalid_bin)
      )
    }
    
    if (length(invalid_sublibrary) > 0) {
      stop_log(
        "Invalid `sublibrary` values in parser xlsx.\n\n",
        "Allowed sublibrary pattern is: L followed by digits, for example L1, L01, L002.\n\n",
        "Affected files:\n",
        make_error_list(invalid_sublibrary)
      )
    }
    
    if (length(invalid_sample) > 0) {
      stop_log(
        "Invalid `sample` values in parser xlsx.\n\n",
        "Sample names must be non-empty and must not contain underscores.\n\n",
        "Affected files:\n",
        make_error_list(invalid_sample)
      )
    }
    
    if (length(invalid_read) > 0) {
      stop_log(
        "Invalid `read` values in parser xlsx.\n\n",
        "Use an empty value for single-end data, or R1/R2 for paired-end data.\n\n",
        "Affected files:\n",
        make_error_list(invalid_read)
      )
    }
    
    parser$pipeline_name <- paste(
      parser$bin,
      parser$sublibrary,
      parser$sample,
      sep = "_"
    )
    
    # Match xlsx entries to files in the input directory.
    files_in_dir <- input_file_basename
    files_in_xlsx <- parser$original_file
    
    missing_from_dir <- setdiff(files_in_xlsx, files_in_dir)
    missing_from_xlsx <- setdiff(files_in_dir, files_in_xlsx)
    
    if (length(missing_from_dir) > 0) {
      stop_log(
        "Some files listed in the parser xlsx were not found in the input directory:\n",
        make_error_list(missing_from_dir)
      )
    }
    
    if (strict_file_match && length(missing_from_xlsx) > 0) {
      stop_log(
        "Some supported input files in the directory are missing from the parser xlsx:\n",
        make_error_list(missing_from_xlsx),
        "\n\nBecause `strict_file_match = TRUE`, every supported input file must be listed in the parser xlsx."
      )
    }
    
    matched_idx <- match(parser$original_file, input_file_basename)
    
    manifest <- data.frame(
      original_file = input_files[matched_idx],
      original_file_basename = input_file_basename[matched_idx],
      original_extension = input_file_ext[matched_idx],
      bin = parser$bin,
      sublibrary = parser$sublibrary,
      sample = parser$sample,
      read = parser$read,
      pipeline_name = parser$pipeline_name,
      stringsAsFactors = FALSE
    )
    
  } else {
    # -----------------------------
    # Case 2: no xlsx parser provided
    # -----------------------------
    log_info("No xlxs parser file with: filename <-> pipline information association provided.\n",
    "Proceeding by using filenames.")
    valid_name <- validate_input_name(input_file_stem)
    
    if (!all(valid_name)) {
      invalid_files <- input_file_basename[!valid_name]
      
      stop_log(
        "Some input files do not match the required naming pattern.\n\n",
        "Required single-end pattern:\n",
        "  BIN_SUBLIBRARY_SAMPLE.fastq.gz\n\n",
        "Required paired-end patterns:\n",
        "  BIN_SUBLIBRARY_SAMPLE_R1.fastq.gz\n",
        "  BIN_SUBLIBRARY_SAMPLE_R2.fastq.gz\n\n",
        "Where:\n",
        "  BIN        = I, L, or U\n",
        "  SUBLIBRARY = L followed by digits, e.g. L1, L01, L002\n",
        "  SAMPLE     = any string without underscores\n\n",
        "Invalid files:\n",
        make_error_list(invalid_files)
      )
    }
    
    parsed <- parse_input_name(input_file_stem)
    
    manifest <- data.frame(
      original_file = input_files,
      original_file_basename = input_file_basename,
      original_extension = input_file_ext,
      bin = parsed$bin,
      sublibrary = parsed$sublibrary,
      sample = parsed$sample,
      read = parsed$read,
      pipeline_name = paste(
        parsed$bin,
        parsed$sublibrary,
        parsed$sample,
        sep = "_"
      ),
      stringsAsFactors = FALSE
    )
  }
  
  # Unique identifier for an individual FASTQ/QC job.
  manifest$fastq_id <- ifelse(
    is.na(manifest$read) | !nzchar(manifest$read),
    manifest$pipeline_name,
    paste(manifest$pipeline_name, manifest$read, sep = "_")
  )
  
  validate_fastq_manifest_layout(manifest)
  
  # -----------------------------
  # Create symlink names
  # -----------------------------
  
  manifest$symlink_extension <- ifelse(
    grepl("\\.(gz|gzip)$", manifest$original_extension, ignore.case = TRUE),
    ".fastq.gz",
    ".fastq"
  )
  
  manifest$symlink_file_basename <- paste0(
    manifest$fastq_id,
    manifest$symlink_extension
  )
  
  if (anyDuplicated(manifest$symlink_file_basename)) {
    duplicated_symlinks <- unique(
      manifest$symlink_file_basename[duplicated(manifest$symlink_file_basename)]
    )
    
    stop_log(
      "Duplicate symlink names would be created:\n",
      make_error_list(duplicated_symlinks)
    )
  }
  
  manifest$symlink_file <- file.path(
    output_symlink_dir,
    manifest$symlink_file_basename
  )
  
  # -----------------------------
  # Prepare symlink directory
  # -----------------------------
  
  if (!dir.exists(output_symlink_dir)) {
    dir.create(output_symlink_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (!dir.exists(output_symlink_dir)) {
    stop_log(
      "Could not create output symlink directory:\n  ",
      output_symlink_dir
    )
  }
  
  existing_symlinks <- manifest$symlink_file[
    file.exists(manifest$symlink_file) |
      nzchar(Sys.readlink(manifest$symlink_file))
  ]
  
  if (length(existing_symlinks) > 0) {
    if (!overwrite_symlinks) {
      stop_log(
        "Some target symlink files already exist:\n",
        make_error_list(existing_symlinks),
        "\n\nSet `overwrite_symlinks = TRUE` to replace them."
      )
    } else {
      unlink(existing_symlinks)
    }
  }
  
  # -----------------------------
  # Create symlinks
  # -----------------------------
  
  symlink_ok <- file.symlink(
    from = normalizePath(manifest$original_file, mustWork = TRUE),
    to = manifest$symlink_file
  )
  
  if (!all(symlink_ok)) {
    failed <- manifest$symlink_file[!symlink_ok]
    
    stop_log(
      "Failed to create one or more symlinks:\n",
      make_error_list(failed),
      "\n\nThis can happen on some systems if symlink creation is not permitted."
    )
  }
  
  # -----------------------------
  # Write manifest
  # -----------------------------
  
  if (!is.null(manifest_output_path)) {
    manifest_dir <- dirname(manifest_output_path)
    
    if (!dir.exists(manifest_dir)) {
      dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    write.table(
      manifest,
      file = manifest_output_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = ""
    )
  }
  
  log_info(paste("Using standardized FASTQ folder: ", output_symlink_dir))
  
  if (!is.null(manifest_output_path)) {
    log_info(paste("FASTQ manifest written to: ", manifest_output_path))
  }
  
  return(manifest)
}


prepare_bcwithqc_inputs <- function(
    manifest,
    output_symlink_dir,
    overwrite_symlinks = TRUE,
    manifest_output_path = NULL
) {
  if (is.null(manifest) || !is.data.frame(manifest)) {
    stop_log("`manifest` must be a data.frame.")
  }
  
  required_cols <- c(
    "pipeline_name",
    "read",
    "fastq_id",
    "symlink_file_basename",
    "qc_filtered_paths"
  )
  
  missing_cols <- setdiff(required_cols, colnames(manifest))
  
  if (length(missing_cols) > 0) {
    stop_log(
      "Manifest is missing required columns for bcwithqc input creation: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  validate_fastq_manifest_layout(manifest)
  
  if (!dir.exists(output_symlink_dir)) {
    dir.create(output_symlink_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (!dir.exists(output_symlink_dir)) {
    stop_log("Could not create bcwithqc symlink folder: ", output_symlink_dir)
  }
  
  manifest$bcwithqc_input_dir <- file.path(
    output_symlink_dir,
    manifest$pipeline_name
  )
  
  manifest$bcwithqc_fastq_path <- file.path(
    manifest$bcwithqc_input_dir,
    manifest$symlink_file_basename
  )
  
  # Rebuild every sample directory once. This avoids stale files when a sample
  # changes from single-end to paired-end or vice versa.
  sample_dirs <- unique(manifest$bcwithqc_input_dir)
  
  for (sample_dir in sample_dirs) {
    sample_dir_exists <- dir.exists(sample_dir) || nzchar(Sys.readlink(sample_dir))
    
    if (sample_dir_exists) {
      if (!overwrite_symlinks) {
        stop_log(
          "bcwithqc sample directory already exists: ", sample_dir,
          "\nSet overwrite_symlinks = TRUE to rebuild it."
        )
      }
      
      unlink(sample_dir, recursive = TRUE, force = TRUE)
    }
    
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (!dir.exists(sample_dir)) {
      stop_log("Could not create bcwithqc sample directory: ", sample_dir)
    }
  }
  
  for (i in seq_len(nrow(manifest))) {
    target_file <- manifest$bcwithqc_fastq_path[i]
    source_file <- manifest$qc_filtered_paths[i]
    
    if (is.na(source_file) || !nzchar(source_file)) {
      stop_log(
        "Missing QC-filtered FASTQ path for FASTQ: ",
        manifest$fastq_id[i]
      )
    }
    
    if (!file.exists(source_file)) {
      stop_log(
        "QC-filtered FASTQ file does not exist for ",
        manifest$fastq_id[i],
        ":\n  ",
        source_file
      )
    }
    
    ok <- file.symlink(
      from = normalizePath(source_file, mustWork = TRUE),
      to = target_file
    )
    
    if (!ok) {
      stop_log(
        "Failed to create bcwithqc symlink:\n",
        "  from: ", source_file, "\n",
        "  to:   ", target_file
      )
    }
  }
  
  if (!is.null(manifest_output_path)) {
    manifest_dir <- dirname(manifest_output_path)
    
    if (!dir.exists(manifest_dir)) {
      dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    write.table(
      manifest,
      file = manifest_output_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = ""
    )
  }
  
  manifest
}
