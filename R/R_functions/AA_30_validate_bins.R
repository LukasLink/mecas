# R/R_functions/AA_30_validate_bins.R

validate_bins <- function(fastq_name_table_xlsx) {
  
  if (is.null(fastq_name_table_xlsx) ||
      length(fastq_name_table_xlsx) != 1 ||
      is.na(fastq_name_table_xlsx) ||
      !is.character(fastq_name_table_xlsx) ||
      !nzchar(trimws(fastq_name_table_xlsx))) {
    stop_log("`fastq_name_table_xlsx` must be a single non-empty file path.", call. = FALSE)
  }
  
  fastq_name_table_xlsx <- path.expand(trimws(fastq_name_table_xlsx))
  
  if (!file.exists(fastq_name_table_xlsx)) {
    stop_log(
      "The FASTQ name table does not exist:\n  ",
      fastq_name_table_xlsx,
      call. = FALSE
    )
  }
  
  if (!grepl("\\.xlsx$", fastq_name_table_xlsx, ignore.case = TRUE)) {
    stop_log(
      "The FASTQ name table must be an `.xlsx` file:\n  ",
      fastq_name_table_xlsx,
      call. = FALSE
    )
  }
  
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop_log(
      "Package `readxl` is required to read the FASTQ name table.\n",
      "Install it with:\n  install.packages(\"readxl\")",
      call. = FALSE
    )
  }
  
  available_sheets <- tryCatch(
    readxl::excel_sheets(fastq_name_table_xlsx),
    error = function(e) {
      stop_log(
        "Failed to open the FASTQ name table:\n  ",
        fastq_name_table_xlsx,
        "\n\nUnderlying error:\n  ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  required_sheets <- c("fastq_names", "bins")
  missing_sheets <- setdiff(required_sheets, available_sheets)
  
  if (length(missing_sheets) > 0) {
    stop_log(
      "The FASTQ name table is missing required sheet(s):\n",
      paste0("  - ", missing_sheets, collapse = "\n"),
      "\n\nRequired sheets are:\n",
      paste0("  - ", required_sheets, collapse = "\n"),
      "\n\nAvailable sheets are:\n",
      paste0("  - ", available_sheets, collapse = "\n"),
      call. = FALSE
    )
  }
  
  bins <- tryCatch(
    readxl::read_excel(
      path = fastq_name_table_xlsx,
      sheet = "bins"
    ),
    error = function(e) {
      stop_log(
        "Failed to read the `bins` sheet from:\n  ",
        fastq_name_table_xlsx,
        "\n\nUnderlying error:\n  ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  fastq_names <- tryCatch(
    readxl::read_excel(
      path = fastq_name_table_xlsx,
      sheet = "fastq_names"
    ),
    error = function(e) {
      stop_log(
        "Failed to read the `fastq_names` sheet from:\n  ",
        fastq_name_table_xlsx,
        "\n\nUnderlying error:\n  ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  bins <- as.data.frame(bins, stringsAsFactors = FALSE)
  fastq_names <- as.data.frame(fastq_names, stringsAsFactors = FALSE)
  
  required_bin_columns <- c(
    "bin_name",
    "sorted_or_unsorted",
    "bin_fraction_min",
    "bin_fraction_max"
  )
  
  missing_bin_columns <- setdiff(required_bin_columns, colnames(bins))
  
  if (length(missing_bin_columns) > 0) {
    stop_log(
      "The `bins` sheet is missing required columns:\n",
      paste0("  - ", missing_bin_columns, collapse = "\n"),
      "\n\nRequired columns are:\n",
      paste0("  - ", required_bin_columns, collapse = "\n"),
      call. = FALSE
    )
  }
  
  bins <- bins[, required_bin_columns, drop = FALSE]
  
  min_was_provided <- !is.na(bins$bin_fraction_min) & nzchar(trimws(as.character(bins$bin_fraction_min)))
  max_was_provided <- !is.na(bins$bin_fraction_max) & nzchar(trimws(as.character(bins$bin_fraction_max)))
  
  bins$bin_name <- trimws(as.character(bins$bin_name))
  bins$sorted_or_unsorted <- tolower(trimws(as.character(bins$sorted_or_unsorted)))
  bins$bin_fraction_min <- suppressWarnings(as.numeric(bins$bin_fraction_min))
  bins$bin_fraction_max <- suppressWarnings(as.numeric(bins$bin_fraction_max))
  
  if (any(is.na(bins$bin_name) | !nzchar(bins$bin_name))) {
    stop_log("The `bins` sheet contains empty values in `bin_name`.", call. = FALSE)
  }
  
  invalid_bin_name_format <- !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", bins$bin_name)
  
  if (any(invalid_bin_name_format)) {
    stop_log(
      "Invalid `bin_name` values were found.\n\n",
      "Bin names may contain letters, numbers, underscores, and hyphens, and must start with a letter or number.\n\n",
      "Invalid names:\n",
      paste0("  - ", bins$bin_name[invalid_bin_name_format], collapse = "\n"),
      call. = FALSE
    )
  }
  
  if (anyDuplicated(bins$bin_name)) {
    duplicated_bin_names <- unique(bins$bin_name[duplicated(bins$bin_name)])
    
    stop_log(
      "The `bins` sheet contains duplicate `bin_name` values:\n",
      paste0("  - ", duplicated_bin_names, collapse = "\n"),
      call. = FALSE
    )
  }
  
  allowed_bin_types <- c("sorted", "unsorted")
  invalid_bin_type_rows <- is.na(bins$sorted_or_unsorted) | !nzchar(bins$sorted_or_unsorted) | !bins$sorted_or_unsorted %in% allowed_bin_types
  
  if (any(invalid_bin_type_rows)) {
    invalid_bin_types <- unique(bins$sorted_or_unsorted[invalid_bin_type_rows])
    invalid_bin_types[is.na(invalid_bin_types) | !nzchar(invalid_bin_types)] <- "<empty>"
    
    stop_log(
      "Invalid values were found in `sorted_or_unsorted`:\n",
      paste0("  - ", invalid_bin_types, collapse = "\n"),
      "\n\nAllowed values are `sorted` and `unsorted`.",
      call. = FALSE
    )
  }
  
  unsorted_rows <- bins$sorted_or_unsorted == "unsorted"
  sorted_rows <- bins$sorted_or_unsorted == "sorted"
  
  if (sum(unsorted_rows) != 1) {
    stop_log(
      "Exactly one unsorted bin is required, but ",
      sum(unsorted_rows),
      " were found.",
      call. = FALSE
    )
  }
  
  if (!any(sorted_rows)) {
    log_warn("No sorted bins were specified. This will lead to a failure when running MAUDE.", call. = FALSE)
  }
  
  if (any(min_was_provided[unsorted_rows]) || any(max_was_provided[unsorted_rows])) {
    log_warn(
      "Values for `bin_fraction_min` or `bin_fraction_max` were found for the unsorted bin. ",
      "These values will be ignored because the unsorted bin has no fraction limits."
    )
  }
  
  bins$bin_fraction_min[unsorted_rows] <- NA_real_
  bins$bin_fraction_max[unsorted_rows] <- NA_real_
  
  invalid_sorted_min <- sorted_rows & (!min_was_provided | is.na(bins$bin_fraction_min))
  invalid_sorted_max <- sorted_rows & (!max_was_provided | is.na(bins$bin_fraction_max))
  
  if (any(invalid_sorted_min | invalid_sorted_max)) {
    invalid_bins <- bins$bin_name[invalid_sorted_min | invalid_sorted_max]
    
    stop_log(
      "All sorted bins must have numeric `bin_fraction_min` and `bin_fraction_max` values.\n\n",
      "Affected bins:\n",
      paste0("  - ", invalid_bins, collapse = "\n"),
      call. = FALSE
    )
  }
  
  outside_allowed_range <- sorted_rows & (
    bins$bin_fraction_min < 0 |
      bins$bin_fraction_min > 1 |
      bins$bin_fraction_max < 0 |
      bins$bin_fraction_max > 1
  )
  
  if (any(outside_allowed_range)) {
    stop_log(
      "Bin fractions must be between 0 and 1.\n\n",
      "Affected bins:\n",
      paste0("  - ", bins$bin_name[outside_allowed_range], collapse = "\n"),
      call. = FALSE
    )
  }
  
  invalid_bound_order <- sorted_rows & bins$bin_fraction_min >= bins$bin_fraction_max
  
  if (any(invalid_bound_order)) {
    stop_log(
      "`bin_fraction_min` must be smaller than `bin_fraction_max`.\n\n",
      "Affected bins:\n",
      paste0("  - ", bins$bin_name[invalid_bound_order], collapse = "\n"),
      call. = FALSE
    )
  }
  
  full_distribution_bins <- sorted_rows & bins$bin_fraction_min == 0 & bins$bin_fraction_max == 1
  
  if (any(full_distribution_bins)) {
    stop_log(
      "A sorted bin cannot span the complete 0–1 distribution.\n\n",
      "Affected bins:\n",
      paste0("  - ", bins$bin_name[full_distribution_bins], collapse = "\n"),
      call. = FALSE
    )
  }
  
  sorted_bins <- bins[sorted_rows, , drop = FALSE]
  sorted_bins <- sorted_bins[order(sorted_bins$bin_fraction_min, sorted_bins$bin_fraction_max), , drop = FALSE]
  
  if (nrow(sorted_bins) > 1) {
    overlapping_rows <- sorted_bins$bin_fraction_min[-1] < sorted_bins$bin_fraction_max[-nrow(sorted_bins)]
    
    if (any(overlapping_rows)) {
      overlap_indices <- which(overlapping_rows)
      
      overlap_descriptions <- vapply(
        overlap_indices,
        function(i) {
          paste0(
            sorted_bins$bin_name[i],
            " (", sorted_bins$bin_fraction_min[i], "–", sorted_bins$bin_fraction_max[i], ") overlaps ",
            sorted_bins$bin_name[i + 1],
            " (", sorted_bins$bin_fraction_min[i + 1], "–", sorted_bins$bin_fraction_max[i + 1], ")"
          )
        },
        character(1)
      )
      
      stop_log(
        "Overlapping sorted-bin ranges were found:\n",
        paste0("  - ", overlap_descriptions, collapse = "\n"),
        call. = FALSE
      )
    }
  }
  
  required_fastq_columns <- c(
    "original_file",
    "bin_name",
    "sublibrary",
    "sample",
    "read"
  )
  
  missing_fastq_columns <- setdiff(required_fastq_columns, colnames(fastq_names))
  
  if (length(missing_fastq_columns) > 0) {
    stop_log(
      "The `fastq_names` sheet is missing required columns:\n",
      paste0("  - ", missing_fastq_columns, collapse = "\n"),
      "\n\nRequired columns are:\n",
      paste0("  - ", required_fastq_columns, collapse = "\n"),
      call. = FALSE
    )
  }
  
  fastq_bin_names <- trimws(as.character(fastq_names$bin_name))
  
  if (any(is.na(fastq_bin_names) | !nzchar(fastq_bin_names))) {
    stop_log("The `fastq_names` sheet contains empty `bin_name` values.", call. = FALSE)
  }
  
  undefined_fastq_bins <- setdiff(unique(fastq_bin_names), bins$bin_name)
  
  if (length(undefined_fastq_bins) > 0) {
    stop_log(
      "The following bin names occur in `fastq_names` but are not defined in the `bins` sheet:\n",
      paste0("  - ", undefined_fastq_bins, collapse = "\n"),
      call. = FALSE
    )
  }
  
  unused_bin_definitions <- setdiff(bins$bin_name, unique(fastq_bin_names))
  
  if (length(unused_bin_definitions) > 0) {
    stop_log(
      "The following bins are defined in the `bins` sheet but do not occur in `fastq_names`:\n",
      paste0("  - ", unused_bin_definitions, collapse = "\n"),
      call. = FALSE
    )
  }
  
  bins <- dplyr::bind_rows(
    bins[bins$sorted_or_unsorted == "unsorted", , drop = FALSE],
    sorted_bins
  )
  
  reserved_bin_names <- c(
    "sgRNA",
    "sgrna_id",
    "count",
    "group_category",
    "bin_name",
    "sublib",
    "sublibrary",
    "sample",
    "exp",
    "isNontargeting"
  )
  
  invalid_reserved_bin_names <- intersect(
    bins$bin_name,
    reserved_bin_names
  )
  
  if (length(invalid_reserved_bin_names) > 0) {
    stop_log(
      "The following `bin_name` values are reserved pipeline column names:\n",
      paste0("  - ", invalid_reserved_bin_names, collapse = "\n")
    )
  }
  
  rownames(bins) <- NULL
  
  log_info(
    "Successfully validated ",
    nrow(bins),
    " bin definitions from: ",
    fastq_name_table_xlsx
  )
  
  return(bins)
}
