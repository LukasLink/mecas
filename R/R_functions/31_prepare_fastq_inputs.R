prepare_fastq_inputs <- function(
    fastq_dir,
    fastq_name_table_file_path,
    bins,
    output_symlink_dir,
    manifest_output_path = NULL,
    overwrite_symlinks = TRUE
) {
  
  # ---------------------------------------------------------------------------
  # Helper functions
  # ---------------------------------------------------------------------------
  
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
    
    ifelse(
      m == -1,
      NA_character_,
      regmatches(basename(x), m)
    )
  }
  
  normalize_read <- function(x) {
    
    x <- trimws(toupper(as.character(x)))
    
    x[x %in% c("", "NA", "N/A", "NULL", "NONE", "SE", "SINGLE")] <- NA_character_
    x[x %in% c("1", "READ1", "READ_1")] <- "R1"
    x[x %in% c("2", "READ2", "READ_2")] <- "R2"
    
    x
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
  
  make_error_list <- function(x, bullet = "  - ") {
    paste0(bullet, x, collapse = "\n")
  }
  
  fastq_dir <- path.expand(trimws(fastq_dir))
  fastq_name_table_file_path <- path.expand(trimws(fastq_name_table_file_path))
  output_symlink_dir <- path.expand(trimws(output_symlink_dir))

  # ---------------------------------------------------------------------------
  # Detect input files
  # ---------------------------------------------------------------------------
  
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
  
  if (anyDuplicated(input_file_basename)) {
    duplicated_files <- unique(
      input_file_basename[duplicated(input_file_basename)]
    )
    
    stop_log(
      "Duplicate input file names were found. File basenames must be unique:\n",
      make_error_list(duplicated_files)
    )
  }
  
  # ---------------------------------------------------------------------------
  # Load fastq_names sheet
  # ---------------------------------------------------------------------------
  
  log_info(
    "Using XLSX FASTQ parser file:\n  ",
    fastq_name_table_file_path
  )
  
  parser <- tryCatch(
    readxl::read_excel(
      path = fastq_name_table_file_path,
      sheet = "fastq_names"
    ),
    error = function(e) {
      stop_log(
        "Failed to read the `fastq_names` sheet from:\n  ",
        fastq_name_table_file_path,
        "\n\nUnderlying error:\n  ",
        conditionMessage(e)
      )
    }
  )
  
  parser <- as.data.frame(
    parser,
    stringsAsFactors = FALSE
  )
  
  required_cols <- c(
    "original_file",
    "bin_name",
    "sublibrary",
    "sample",
    "read"
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(parser)
  )
  
  if (length(missing_cols) > 0) {
    stop_log(
      "The `fastq_names` sheet is missing required columns:\n",
      make_error_list(missing_cols),
      "\n\nRequired columns are:\n",
      make_error_list(required_cols)
    )
  }
  
  parser <- parser[, required_cols, drop = FALSE]
  
  parser$original_file <- trimws(as.character(parser$original_file))
  parser$bin_name <- trimws(as.character(parser$bin_name))
  parser$sublibrary <- trimws(as.character(parser$sublibrary))
  parser$sample <- trimws(as.character(parser$sample))
  parser$read <- normalize_read(parser$read)
  
  # ---------------------------------------------------------------------------
  # Validate parser values
  # ---------------------------------------------------------------------------
  
  if (any(is.na(parser$original_file) | !nzchar(parser$original_file))) {
    stop_log(
      "The `fastq_names` sheet contains empty values in `original_file`."
    )
  }
  
  if (anyDuplicated(parser$original_file)) {
    duplicated_original <- unique(
      parser$original_file[duplicated(parser$original_file)]
    )
    
    stop_log(
      "The `fastq_names` sheet contains duplicate `original_file` entries:\n",
      make_error_list(duplicated_original)
    )
  }
  
  invalid_bin_name <- parser$original_file[
    is.na(parser$bin_name) |
      !nzchar(parser$bin_name) |
      !(parser$bin_name %in% bins$bin_name)
  ]
  
  invalid_sublibrary <- parser$original_file[
    !validate_sublibrary(parser$sublibrary)
  ]
  
  invalid_sample <- parser$original_file[
    !validate_sample(parser$sample)
  ]
  
  invalid_read <- parser$original_file[
    !validate_read(parser$read)
  ]
  
  if (length(invalid_bin_name) > 0) {
    stop_log(
      "Invalid `bin_name` values were found in the `fastq_names` sheet.\n\n",
      "Every `bin_name` must be defined in the `bins` sheet.\n\n",
      "Affected files:\n",
      make_error_list(invalid_bin_name)
    )
  }
  
  if (length(invalid_sublibrary) > 0) {
    stop_log(
      "Invalid `sublibrary` values were found in the `fastq_names` sheet.\n\n",
      "Allowed sublibrary pattern: L followed by digits, for example L1, L01, or L002.\n\n",
      "Affected files:\n",
      make_error_list(invalid_sublibrary)
    )
  }
  
  if (length(invalid_sample) > 0) {
    stop_log(
      "Invalid `sample` values were found in the `fastq_names` sheet.\n\n",
      "Sample names must be non-empty and must not contain underscores.\n\n",
      "Affected files:\n",
      make_error_list(invalid_sample)
    )
  }
  
  if (length(invalid_read) > 0) {
    stop_log(
      "Invalid `read` values were found in the `fastq_names` sheet.\n\n",
      "Use an empty value for single-end data, or R1/R2 for paired-end data.\n\n",
      "Affected files:\n",
      make_error_list(invalid_read)
    )
  }
  
  # ---------------------------------------------------------------------------
  # Construct internal names
  # ---------------------------------------------------------------------------
  
  parser$pipeline_name <- paste(
    parser$bin_name,
    parser$sublibrary,
    parser$sample,
    sep = "_"
  )
  
  # ---------------------------------------------------------------------------
  # Match parser entries to files in input directory
  # ---------------------------------------------------------------------------
  
  files_in_dir <- input_file_basename
  files_in_xlsx <- parser$original_file
  
  missing_from_dir <- setdiff(files_in_xlsx, files_in_dir)
  
  missing_from_xlsx <- setdiff(files_in_dir, files_in_xlsx)
  
  if (length(missing_from_dir) > 0) {
    stop_log(
      "Some files listed in the `fastq_names` sheet were not found in the input directory:\n",
      make_error_list(missing_from_dir)
    )
  }
  
  if (length(missing_from_xlsx) > 0) {
    logger::log_warn(
      "Some supported input files in the directory are missing from the `fastq_names` sheet:\n",
      make_error_list(missing_from_xlsx),
      "\n\nThese files will be ignored because they are not listed in the XLSX manifest."
    )
  }
  
  matched_idx <- match(parser$original_file, input_file_basename)
  
  # ---------------------------------------------------------------------------
  # Construct manifest
  # ---------------------------------------------------------------------------
  
  manifest <- data.frame(
    original_file = input_files[matched_idx],
    original_file_basename = input_file_basename[matched_idx],
    original_extension = input_file_ext[matched_idx],
    bin_name = parser$bin_name,
    sublibrary = parser$sublibrary,
    sample = parser$sample,
    read = parser$read,
    pipeline_name = parser$pipeline_name,
    stringsAsFactors = FALSE
  )
  
  # Unique identifier for an individual FASTQ/QC job.
  manifest$fastq_id <- ifelse(
    is.na(manifest$read) | !nzchar(manifest$read),
    manifest$pipeline_name,
    paste(
      manifest$pipeline_name,
      manifest$read,
      sep = "_"
    )
  )
  
  validate_fastq_manifest(manifest)
  
  # ---------------------------------------------------------------------------
  # Construct standardized symlink names
  # ---------------------------------------------------------------------------
  
  manifest$symlink_extension <- ifelse(
    grepl(
      "\\.(gz|gzip)$",
      manifest$original_extension,
      ignore.case = TRUE
    ),
    ".fastq.gz",
    ".fastq"
  )
  
  manifest$symlink_file_basename <- paste0(
    manifest$fastq_id,
    manifest$symlink_extension
  )
  
  if (anyDuplicated(manifest$symlink_file_basename)) {
    duplicated_symlinks <- unique(
      manifest$symlink_file_basename[
        duplicated(manifest$symlink_file_basename)
      ]
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
  
  # ---------------------------------------------------------------------------
  # Prepare symlink directory
  # ---------------------------------------------------------------------------
  
  if (!dir.exists(output_symlink_dir)) {
    dir.create(
      output_symlink_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
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
    }
    
    unlink(existing_symlinks, force = TRUE)
  }
  
  # ---------------------------------------------------------------------------
  # Create symlinks
  # ---------------------------------------------------------------------------
  
  symlink_ok <- file.symlink(
    from = normalizePath(
      manifest$original_file,
      mustWork = TRUE
    ),
    to = manifest$symlink_file
  )
  
  if (!all(symlink_ok)) {
    failed <- manifest$symlink_file[!symlink_ok]
    
    stop_log(
      "Failed to create one or more symlinks:\n",
      make_error_list(failed),
      "\n\nThis can happen if symlink creation is not permitted on the current filesystem."
    )
  }
  
  # ---------------------------------------------------------------------------
  # Write manifest
  # ---------------------------------------------------------------------------
  
  if (!is.null(manifest_output_path)) {
    manifest_dir <- dirname(
      manifest_output_path
    )
    
    if (!dir.exists(manifest_dir)) {
      dir.create(
        manifest_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )
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
  
  log_info(
    "Using standardized FASTQ folder: ",
    output_symlink_dir
  )
  
  if (!is.null(manifest_output_path)) {
    log_info(
      "FASTQ manifest written to: ",
      manifest_output_path
    )
  }
  
  return(manifest)
}