# 01a_check_config.R
#-------------------------------------------------------------------------------
# Config/input validation helpers
#-------------------------------------------------------------------------------

# This file provides one main function:
#   check_input(option, name)
#
# It validates a single config/setup option based on its name.
#
# Example:
#   check_input(opt$machine, "machine")
#   check_input(opt$slurm_wall_time, "slurm_wall_time")
#
# The function returns the input value invisibly if valid.
# It stops with an informative error if invalid.

#-------------------------------------------------------------------------------
# Small helpers
#-------------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.is_missing <- function(x) {
  is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))
}

.is_empty_string <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x) && !nzchar(x)
}

.is_null_or_empty <- function(x) {
  .is_missing(x) || .is_empty_string(x)
}

.fail_input <- function(name, option, message) {
  stop(
    "Invalid config option `", name, "`.\n",
    message, "\n",
    "Current value: ", paste(capture.output(str(option)), collapse = " "),
    call. = FALSE
  )
}

#-------------------------------------------------------------------------------
# Generic validators
#-------------------------------------------------------------------------------

.validate_required <- function(option, name) {
  if (.is_missing(option)) {
    .fail_input(name, option, "This option is required and cannot be NULL or NA.")
  }
  
  invisible(option)
}

.validate_optional <- function(option, name) {
  invisible(option)
}

.validate_logical <- function(option, name, required = TRUE) {
  if (!required && .is_missing(option)) {
    return(invisible(option))
  }
  
  .validate_required(option, name)
  
  if (!is.logical(option) || length(option) != 1 || is.na(option)) {
    .fail_input(name, option, "Expected a single TRUE or FALSE value.")
  }
  
  invisible(option)
}

.validate_string <- function(option,
                             name,
                             required = TRUE,
                             allow_empty = FALSE,
                             allow_null = FALSE) {
  if (!required && .is_missing(option)) {
    return(invisible(option))
  }
  
  if (allow_null && is.null(option)) {
    return(invisible(option))
  }
  
  .validate_required(option, name)
  
  if (!is.character(option) || length(option) != 1 || is.na(option)) {
    .fail_input(name, option, "Expected a single character string.")
  }
  
  if (!allow_empty && !nzchar(option)) {
    .fail_input(name, option, "Expected a non-empty character string.")
  }
  
  invisible(option)
}

.validate_choice <- function(choices,
                             required = TRUE,
                             allow_empty = FALSE) {
  force(choices)
  
  function(option, name) {
    if (!required && .is_missing(option)) {
      return(invisible(option))
    }
    
    .validate_string(
      option = option,
      name = name,
      required = required,
      allow_empty = allow_empty
    )
    
    if (!option %in% choices) {
      .fail_input(
        name,
        option,
        paste0(
          "Expected one of: ",
          paste(shQuote(choices), collapse = ", "),
          "."
        )
      )
    }
    
    invisible(option)
  }
}

.validate_integer <- function(option,
                              name,
                              required = TRUE,
                              min = NULL,
                              max = NULL,
                              allow_na = FALSE) {
  if (!required && .is_missing(option)) {
    return(invisible(option))
  }
  
  if (allow_na && length(option) == 1 && is.na(option)) {
    return(invisible(option))
  }
  
  .validate_required(option, name)
  
  numeric_option <- suppressWarnings(as.numeric(option))
  
  if (length(numeric_option) != 1 ||
      is.na(numeric_option) ||
      numeric_option != as.integer(numeric_option)) {
    .fail_input(name, option, "Expected a single integer value.")
  }
  
  if (!is.null(min) && numeric_option < min) {
    .fail_input(name, option, paste0("Expected an integer >= ", min, "."))
  }
  
  if (!is.null(max) && numeric_option > max) {
    .fail_input(name, option, paste0("Expected an integer <= ", max, "."))
  }
  
  invisible(option)
}

.validate_path_string <- function(option,
                                  name,
                                  required = TRUE,
                                  must_exist = FALSE,
                                  allow_empty = FALSE) {
  if (!required && .is_missing(option)) {
    return(invisible(option))
  }
  
  .validate_string(
    option = option,
    name = name,
    required = required,
    allow_empty = allow_empty
  )
  
  if (must_exist && nzchar(option) && !file.exists(option)) {
    .fail_input(name, option, "Path does not exist.")
  }
  
  invisible(option)
}

.validate_file_path <- function(option,
                                name,
                                required = TRUE,
                                must_exist = FALSE,
                                allow_empty = FALSE) {
  .validate_path_string(
    option = option,
    name = name,
    required = required,
    must_exist = must_exist,
    allow_empty = allow_empty
  )
  
  if (must_exist && nzchar(option) && !file.exists(option)) {
    .fail_input(name, option, "File does not exist.")
  }
  
  invisible(option)
}

.validate_dir_path <- function(option,
                               name,
                               required = TRUE,
                               must_exist = FALSE,
                               allow_empty = FALSE) {
  .validate_path_string(
    option = option,
    name = name,
    required = required,
    must_exist = must_exist,
    allow_empty = allow_empty
  )
  
  if (must_exist && nzchar(option) && !dir.exists(option)) {
    .fail_input(name, option, "Directory does not exist.")
  }
  
  invisible(option)
}

.validate_list_or_character <- function(option,
                                        name,
                                        required = FALSE) {
  if (!required && .is_missing(option)) {
    return(invisible(option))
  }
  
  if (!(is.character(option) || is.list(option))) {
    .fail_input(
      name,
      option,
      "Expected a character vector, comma-separated string, list, NULL, or empty value."
    )
  }
  
  invisible(option)
}

.validate_percentage <- function(option,
                                 name,
                                 required = TRUE,
                                 min = 0,
                                 max = 100) {
  .validate_numeric(
    option = option,
    name = name,
    required = required,
    min = min,
    max = max
  )
  
  invisible(option)
}

.validate_mem_string <- function(option,
                                 name,
                                 required = TRUE,
                                 allow_empty = FALSE) {
  .validate_string(
    option = option,
    name = name,
    required = required,
    allow_empty = allow_empty
  )
  
  if (allow_empty && .is_empty_string(option)) {
    return(invisible(option))
  }
  
  # Accept common SLURM memory strings, e.g. 4000M, 4G, 65g, 100gb
  valid <- grepl("^[0-9]+(K|M|G|T|KB|MB|GB|TB|k|m|g|t|kb|mb|gb|tb)?$", option)
  
  if (!valid) {
    .fail_input(
      name,
      option,
      "Expected a memory string such as '4000M', '4G', '65g', or '128GB'."
    )
  }
  
  invisible(option)
}

.validate_wall_time <- function(option,
                                name,
                                required = TRUE,
                                allow_empty = FALSE) {
  .validate_string(
    option = option,
    name = name,
    required = required,
    allow_empty = allow_empty
  )
  
  if (allow_empty && .is_empty_string(option)) {
    return(invisible(option))
  }
  
  # Accept HH:MM:SS, H:MM:SS, D-HH:MM:SS
  valid <- grepl("^([0-9]+-)?[0-9]+:[0-5][0-9]:[0-5][0-9]$", option)
  
  if (!valid) {
    .fail_input(
      name,
      option,
      "Expected SLURM wall time format 'HH:MM:SS' or 'D-HH:MM:SS', e.g. '48:00:00'."
    )
  }
  
  invisible(option)
}

.validate_email <- function(option,
                            name,
                            required = FALSE,
                            allow_empty = TRUE) {
  if (!required && .is_missing(option)) {
    return(invisible(option))
  }
  
  .validate_string(
    option = option,
    name = name,
    required = required,
    allow_empty = allow_empty
  )
  
  if (allow_empty && .is_empty_string(option)) {
    return(invisible(option))
  }
  
  valid <- grepl("^[^@\\s]+@[^@\\s]+\\.[^@\\s]+$", option)
  
  if (!valid) {
    .fail_input(name, option, "Expected a valid email address or empty/null.")
  }
  
  invisible(option)
}

.validate_regex_string <- function(option,
                                   name,
                                   required = FALSE,
                                   allow_empty = TRUE) {
  if (!required && .is_missing(option)) {
    return(invisible(option))
  }
  
  .validate_string(
    option = option,
    name = name,
    required = required,
    allow_empty = allow_empty
  )
  
  if (allow_empty && .is_empty_string(option)) {
    return(invisible(option))
  }
  
  # umi_tools uses Python-style regex syntax, including:
  #   (?P<umi_1>...)
  #   (?P<discard_1>...)
  #   fuzzy matching such as {e<=1}
  #
  # These are not valid base R/TRE regex patterns, so do NOT validate
  # with grepl().
  
  if (!grepl("\\(\\?P<", option, fixed = FALSE)) {
    .fail_input(
      name,
      option,
      "Expected a umi_tools-style regex with named groups such as (?P<umi_1>...) or (?P<discard_1>...)."
    )
  }
  
  if (!grepl("\\(\\?P<umi_[0-9]+>", option, fixed = FALSE)) {
    .fail_input(
      name,
      option,
      "Expected at least one umi_tools UMI group, e.g. (?P<umi_1>...)."
    )
  }
  
  open_parens <- gregexpr("\\(", option)[[1]]
  close_parens <- gregexpr("\\)", option)[[1]]
  
  n_open <- if (identical(open_parens, -1L)) 0L else length(open_parens)
  n_close <- if (identical(close_parens, -1L)) 0L else length(close_parens)
  
  if (n_open != n_close) {
    .fail_input(
      name,
      option,
      "Regex appears to have unbalanced parentheses."
    )
  }
  
  invisible(option)
}

.validate_executable <- function(option,
                                 name,
                                 required = FALSE,
                                 allow_empty = TRUE) {
  if (!required && .is_missing(option)) {
    return(invisible(option))
  }
  
  .validate_string(
    option = option,
    name = name,
    required = required,
    allow_empty = allow_empty
  )
  
  if (allow_empty && .is_empty_string(option)) {
    return(invisible(option))
  }
  
  # If option contains a slash, treat as file path.
  if (grepl("/", option, fixed = TRUE)) {
    if (!file.exists(option)) {
      .fail_input(name, option, "Executable path does not exist.")
    }
    
    if (file.access(option, mode = 1) != 0) {
      .fail_input(name, option, "File exists but is not executable.")
    }
    
    return(invisible(option))
  }
  
  # Otherwise treat as command name on PATH.
  found <- nzchar(Sys.which(option))
  
  if (!found) {
    .fail_input(
      name,
      option,
      paste0(
        "Executable was not found on PATH. ",
        "Provide an absolute path or make sure the command is available on PATH."
      )
    )
  }
  
  invisible(option)
}

.validate_numeric <- function(option,
                              name,
                              required = FALSE,
                              min = -Inf,
                              max = Inf,
                              allow_na = FALSE) {
  if (!required && .is_missing(option)) {
    return(invisible(option))
  }
  
  if (allow_na && length(option) == 1 && is.na(option)) {
    return(invisible(option))
  }
  
  if (!is.numeric(option) || length(option) != 1 || is.na(option)) {
    .fail_input(name, option, "Expected a single numeric value.")
  }
  
  if (option < min || option > max) {
    .fail_input(
      name,
      option,
      paste0("Expected a numeric value between ", min, " and ", max, ".")
    )
  }
  
  invisible(option)
}

.validate_numeric_vector <- function(option,
                                     name,
                                     required = FALSE,
                                     min = -Inf,
                                     max = Inf,
                                     allow_empty = TRUE) {
  if (!required && .is_missing(option)) {
    return(invisible(option))
  }
  
  if (allow_empty && length(option) == 0) {
    return(invisible(option))
  }
  
  if (!is.numeric(option)) {
    .fail_input(name, option, "Expected a numeric vector.")
  }
  
  if (any(is.na(option))) {
    .fail_input(name, option, "Numeric vector must not contain NA values.")
  }
  
  if (any(option < min | option > max)) {
    .fail_input(
      name,
      option,
      paste0("Expected all values to be between ", min, " and ", max, ".")
    )
  }
  
  invisible(option)
}
#-------------------------------------------------------------------------------
# Option-specific validator map
#-------------------------------------------------------------------------------
# Canonical validator names now mirror the nested cfg/YAML structure.
# Example:
#   check_input(cfg$filtering$drop_0s, "filtering.drop_0s")
#   check_input(cfg$paths$output_folder, "paths.output_folder")

.input_validators <- list(
  # Run behavior
  "run.first_time" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "run.start_with" = .validate_choice(
    c("beginning", "read_counting", "MAUDE_analysis", "generate_plots")
  ),
  
  "run.machine" = .validate_choice(c("local", "slurm")),
  
  # Paths
  "paths.output_folder" = function(option, name) {
    .validate_dir_path(option, name, required = TRUE, must_exist = FALSE)
  },
  
  "paths.input_folder" = function(option, name) {
    .validate_dir_path(option, name, required = TRUE, must_exist = TRUE)
  },
  
  "paths.library_path" = function(option, name) {
    .validate_file_path(option, name, required = TRUE, must_exist = TRUE)
  },
  
  "paths.star_index_folder" = function(option, name) {
    .validate_dir_path(option, name, required = TRUE, must_exist = TRUE)
  },
  
  "paths.fastq_name_table_xlsx" = function(option, name) {
    .validate_file_path(
      option,
      name,
      required = FALSE,
      must_exist = FALSE,
      allow_empty = TRUE
    )
  },
  
  "paths.strict_file_match" = function(option, name) {
    .validate_logical(option, name)
  },
  
  # Files
  "files.strict_file_match" = function(option, name) {
    .validate_logical(option, name)
  },
  
  # Modules
  "modules.use_modules" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "modules.seqtk" = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  "modules.star" = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  "modules.samtools" = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  "modules.umi_tools" = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  # bcwithqc options
  "bcwithqc.bcwithqc_bin" = function(option, name) {
    .validate_executable(
      option,
      name,
      required = FALSE,
      allow_empty = TRUE
    )
  },
  
  "bcwithqc.bcwithqc_dir" = function(option, name) {
    .validate_dir_path(
      option,
      name,
      required = FALSE,
      must_exist = FALSE,
      allow_empty = TRUE
    )
  },
  
  "bcwithqc.bcwithqc_config_path" = function(option, name) {
    .validate_file_path(
      option,
      name,
      required = FALSE,
      must_exist = FALSE,
      allow_empty = TRUE
    )
  },
  
  # align_UMI_tools options
  "align_UMI_tools.UMI_regex" = function(option, name) {
    .validate_regex_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  # QC filtering
  "qc_filtering.run" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "qc_filtering.min_qual" = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 0)
  },
  
  "qc_filtering.qual_offset" = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 0)
  },
  
  "qc_filtering.min_length" = function(option, name) {
    .validate_integer(option, name, required = FALSE, min = 1, allow_na = TRUE)
  },
  
  # Skip lists
  "skip.files" = function(option, name) {
    .validate_list_or_character(option, name, required = FALSE)
  },
  
  "skip.sublibraries" = function(option, name) {
    .validate_list_or_character(option, name, required = FALSE)
  },
  
  "skip.samples" = function(option, name) {
    .validate_list_or_character(option, name, required = FALSE)
  },
  
  # Controls
  "controls.include_controls" = function(option, name) {
    .validate_list_or_character(option, name, required = FALSE)
  },
  
  "controls.use_only_these_controls" = function(option, name) {
    .validate_list_or_character(option, name, required = FALSE)
  },
  
  "controls.same_controls_in_all_sublibraries" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "controls.subsample_controls" = function(option, name) {
    .validate_logical(option, name)
  },
  
  # Counting
  "counting.read_counting" = .validate_choice(c("bcwithqc", "align_UMI_tools")),
  
  "counting.data_type" = .validate_choice(c("reads", "umis")),
  
  # Replicates / grouping
  "replicates.method" = function(option, name) {
    .validate_string(option, name, required = TRUE, allow_empty = TRUE)
  },
  
  "replicates.combine_for_guide_stats" = function(option, name) {
    .validate_string(option, name, required = TRUE, allow_empty = TRUE)
  },
  
  "replicates.combine_for_gene_stats" = function(option, name) {
    .validate_string(option, name, required = TRUE, allow_empty = TRUE)
  },
  
  "replicates.auto_combine_replicates" = function(option, name) {
    .validate_logical(option, name)
  },
  
  # Normalization
  "normalization.norm_method" = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  "normalization.recover_input" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "normalization.use_custom_bins" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "normalization.upper_lower_percentage" = function(option, name) {
    .validate_percentage(option, name, required = TRUE, min = 0, max = 100)
  },
  
  # Filtering / MAUDE behavior
  "filtering.drop_0s" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "filtering.strict_mode" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "filtering.min_guides_per_gene" = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 0)
  },
  
  # Output naming
  "output.extra_suffix" = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  # SLURM
  "slurm.account" = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  "slurm.qos" = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  "slurm.cpus" = function(option, name) {
    .validate_integer(option, name, required = FALSE, min = 1)
  },
  
  "slurm.mem" = function(option, name) {
    .validate_mem_string(option, name, required = FALSE)
  },
  
  "slurm.wall_time" = function(option, name) {
    .validate_wall_time(option, name, required = FALSE)
  },
  
  "slurm.partition" = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  "slurm.array" = function(option, name) {
    .validate_integer(option, name, required = FALSE, min = 1)
  },
  
  "slurm.email" = function(option, name) {
    .validate_email(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  # Consensus calling
  "consensus.run" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "consensus.n_reps" = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 1)
  },
  
  "consensus.high_confidence.FDR_threshold" = function(option, name) {
    .validate_percentage(option, name, required = TRUE, min = 0, max = 1)
  },
  
  "consensus.high_confidence.hits_in_X_reps" = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 1)
  },
  
  "consensus.high_confidence.correlation_heatmap" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "consensus.high_confidence.overlap" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "consensus.high_confidence.venn_diagram" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "consensus.explorative.FDR_threshold" = function(option, name) {
    .validate_percentage(option, name, required = TRUE, min = 0, max = 1)
  },
  
  "consensus.explorative.hits_in_X_reps" = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 1)
  },
  
  "consensus.explorative.correlation_heatmap" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "consensus.explorative.overlap" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "consensus.explorative.venn_diagram" = function(option, name) {
    .validate_logical(option, name)
  },
  
  # General plots
  "plots.run_after_count" = function(option, name) {
    .validate_logical(option, name)
  },
  
  # Violin plot
  "plots.read_count_violin.violin_y_limit" = function(option, name) {
    .validate_integer(option, name, min = 1)
  },
  
  # Waterfall plot
  "plots.waterfall.mark_cntrl" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "plots.waterfall.mark_special" = function(option, name) {
    .validate_list_or_character(option, name, required = FALSE)
  },
  
  "plots.waterfall.mark_N_top_hits" = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 0)
  },
  
  "plots.waterfall.box_padding" = function(option, name) {
    .validate_numeric(option, name, required = TRUE, min = 0)
  },
  
  "plots.waterfall.no_text" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "plots.waterfall.signif_lines" = function(option, name) {
    .validate_logical(option, name)
  },
  
  "plots.waterfall.mark_all_signif_level" = function(option, name) {
    .validate_numeric(option, name, required = FALSE, min = 0, max = 1, allow_na = TRUE)
  },
  
  "plots.waterfall.break_in_plot" = function(option, name) {
    .validate_numeric_vector(option, name, required = FALSE, min = 0)
  },
  
  "plots.waterfall.top_padding" = function(option, name) {
    .validate_numeric(option, name, required = TRUE, min = 0)
  },
  
  "plots.waterfall.custom_title" = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE, allow_null = TRUE)
  },
  
  "plots.waterfall.width" = function(option, name) {
    .validate_numeric(option, name, required = TRUE, min = 1)
  },
  
  "plots.waterfall.height" = function(option, name) {
    .validate_numeric(option, name, required = TRUE, min = 1)
  },
  
  "plots.waterfall.file_format" = .validate_choice(c("png", "pdf", "svg"))
)

#-------------------------------------------------------------------------------
# Cross-option / dependency validation
#-------------------------------------------------------------------------------

.require_if <- function(condition,
                        option,
                        name,
                        reason,
                        validator = NULL) {
  if (!isTRUE(condition)) {
    return(invisible(option))
  }
  
  if (.is_null_or_empty(option)) {
    stop(
      "Missing required config option `", name, "`.\n",
      "Reason: ", reason,
      call. = FALSE
    )
  }
  
  if (!is.null(validator)) {
    validator(option, name)
  }
  
  invisible(option)
}

.check_allowed_when <- function(condition,
                                option,
                                name,
                                reason) {
  if (isTRUE(condition)) {
    return(invisible(option))
  }
  
  if (!.is_null_or_empty(option)) {
    warning(
      "Config option `", name, "` was provided but may be ignored.\n",
      "Reason: ", reason,
      call. = FALSE
    )
  }
  
  invisible(option)
}

check_config_dependencies <- function(cfg) {
  #-----------------------------------------------------------------------------
  # Basic required options for dependency checks
  #-----------------------------------------------------------------------------
  
  machine <- cfg$run$machine
  read_counting <- cfg$counting$read_counting
  
  check_input(machine, "run.machine")
  check_input(read_counting, "counting.read_counting")
  
  #-----------------------------------------------------------------------------
  # SLURM options are only required when machine == "slurm"
  #-----------------------------------------------------------------------------
  
  is_slurm <- identical(machine, "slurm")
  
  .require_if(
    condition = is_slurm,
    option = cfg$slurm$cpus,
    name = "slurm.cpus",
    reason = "`run.machine` is set to 'slurm'.",
    validator = function(option, name) {
      .validate_integer(option, name, required = TRUE, min = 1)
    }
  )
  
  .require_if(
    condition = is_slurm,
    option = cfg$slurm$mem,
    name = "slurm.mem",
    reason = "`run.machine` is set to 'slurm'.",
    validator = function(option, name) {
      .validate_mem_string(option, name, required = TRUE)
    }
  )
  
  .require_if(
    condition = is_slurm,
    option = cfg$slurm$wall_time,
    name = "slurm.wall_time",
    reason = "`run.machine` is set to 'slurm'.",
    validator = function(option, name) {
      .validate_wall_time(option, name, required = TRUE)
    }
  )
  
  .require_if(
    condition = is_slurm,
    option = cfg$slurm$array,
    name = "slurm.array",
    reason = "`run.machine` is set to 'slurm'.",
    validator = function(option, name) {
      .validate_integer(option, name, required = TRUE, min = 1)
    }
  )
  
  # Optional SLURM options: validate if provided
  if (!.is_null_or_empty(cfg$slurm$account)) {
    check_input(cfg$slurm$account, "slurm.account")
  }
  
  if (!.is_null_or_empty(cfg$slurm$qos)) {
    check_input(cfg$slurm$qos, "slurm.qos")
  }
  
  if (!.is_null_or_empty(cfg$slurm$partition)) {
    check_input(cfg$slurm$partition, "slurm.partition")
  }
  
  if (!.is_null_or_empty(cfg$slurm$email)) {
    check_input(cfg$slurm$email, "slurm.email")
  }
  
  .check_allowed_when(
    condition = is_slurm,
    option = cfg$slurm$account,
    name = "slurm.account",
    reason = "`run.machine` is not set to 'slurm'."
  )
  
  .check_allowed_when(
    condition = is_slurm,
    option = cfg$slurm$qos,
    name = "slurm.qos",
    reason = "`run.machine` is not set to 'slurm'."
  )
  
  .check_allowed_when(
    condition = is_slurm,
    option = cfg$slurm$partition,
    name = "slurm.partition",
    reason = "`run.machine` is not set to 'slurm'."
  )
  
  .check_allowed_when(
    condition = is_slurm,
    option = cfg$slurm$email,
    name = "slurm.email",
    reason = "`run.machine` is not set to 'slurm'."
  )
  
  #-----------------------------------------------------------------------------
  # bcwithqc-specific requirements
  #-----------------------------------------------------------------------------
  
  uses_bcwithqc <- identical(read_counting, "bcwithqc")
  
  if (uses_bcwithqc) {
    has_bcwithqc_dir <- !.is_null_or_empty(cfg$bcwithqc$bcwithqc_dir)
    has_bcwithqc_config <- !.is_null_or_empty(cfg$bcwithqc$bcwithqc_config_path)
    
    if (!has_bcwithqc_dir && !has_bcwithqc_config) {
      stop(
        "Missing required bcwithqc configuration.\n",
        "Reason: `counting.read_counting` is set to 'bcwithqc'.\n",
        "Provide at least one of:\n",
        "  - `bcwithqc.bcwithqc_dir`\n",
        "  - `bcwithqc.bcwithqc_config_path`",
        call. = FALSE
      )
    }
    
    .require_if(
      condition = TRUE,
      option = cfg$bcwithqc$bcwithqc_bin,
      name = "bcwithqc.bcwithqc_bin",
      reason = "`counting.read_counting` is set to 'bcwithqc'.",
      validator = function(option, name) {
        .validate_executable(
          option,
          name,
          required = TRUE,
          allow_empty = FALSE
        )
      }
    )
    
    if (has_bcwithqc_dir) {
      .validate_dir_path(
        cfg$bcwithqc$bcwithqc_dir,
        "bcwithqc.bcwithqc_dir",
        required = TRUE,
        must_exist = TRUE
      )
    }
    
    if (has_bcwithqc_config) {
      .validate_file_path(
        cfg$bcwithqc$bcwithqc_config_path,
        "bcwithqc.bcwithqc_config_path",
        required = TRUE,
        must_exist = TRUE
      )
    }
  }
  
  #-----------------------------------------------------------------------------
  # Module dependencies
  #-----------------------------------------------------------------------------
  
  if (isTRUE(cfg$modules$use_modules)) {
    if (isTRUE(cfg$qc_filtering$run)) {
      .require_if(
        condition = TRUE,
        option = cfg$modules$seqtk,
        name = "modules.seqtk",
        reason = "`modules.use_modules` is TRUE and `qc_filtering.run` is TRUE.",
        validator = function(option, name) {
          .validate_string(option, name, required = TRUE, allow_empty = FALSE)
        }
      )
    }
    
    if (identical(cfg$counting$read_counting, "align_UMI_tools")) {
      .require_if(
        condition = TRUE,
        option = cfg$modules$star,
        name = "modules.star",
        reason = "`modules.use_modules` is TRUE and `counting.read_counting` is 'align_UMI_tools'.",
        validator = function(option, name) {
          .validate_string(option, name, required = TRUE, allow_empty = FALSE)
        }
      )
      
      .require_if(
        condition = TRUE,
        option = cfg$modules$samtools,
        name = "modules.samtools",
        reason = "`modules.use_modules` is TRUE and `counting.read_counting` is 'align_UMI_tools'.",
        validator = function(option, name) {
          .validate_string(option, name, required = TRUE, allow_empty = FALSE)
        }
      )
      
      if (identical(cfg$counting$data_type, "umis")) {
        .require_if(
          condition = TRUE,
          option = cfg$modules$umi_tools,
          name = "modules.umi_tools",
          reason = "`modules.use_modules` is TRUE, `counting.read_counting` is 'align_UMI_tools', and `counting.data_type` is 'umis'.",
          validator = function(option, name) {
            .validate_string(option, name, required = TRUE, allow_empty = FALSE)
          }
        )
      }
    }
  }
  
  #-----------------------------------------------------------------------------
  # align_UMI_tools requirements
  #-----------------------------------------------------------------------------
  
  uses_umi_tools <- identical(cfg$counting$read_counting, "align_UMI_tools")
  
  .require_if(
    condition = uses_umi_tools,
    option = cfg$align_UMI_tools$UMI_regex,
    name = "align_UMI_tools.UMI_regex",
    reason = "`counting.read_counting` is set to 'align_UMI_tools'.",
    validator = function(option, name) {
      .validate_regex_string(option, name, required = TRUE, allow_empty = FALSE)
    }
  )
  
  #-----------------------------------------------------------------------------
  # consensus/start_with dependency
  #-----------------------------------------------------------------------------
  
  if (isTRUE(cfg$consensus$run) &&
      identical(cfg$run$start_with, "generate_plots")) {
    warning(
      "For `consensus.run`, `run` is set to TRUE but `run.start_with` is 'generate_plots'. ",
      "Consensus calling requires the pipeline to start with at least the MAUDE run.\n",
      " -> Consensus calling will be skipped!\n",
      "Set `run.start_with` to `MAUDE_analysis` or earlier to actually run consensus calling.",
      call. = FALSE
    )
  }
  
  #-----------------------------------------------------------------------------
  # consensus_call dependency
  #-----------------------------------------------------------------------------
  
  if (isTRUE(cfg$consensus$run)) {
    if (cfg$consensus$high_confidence$hits_in_X_reps >
        cfg$consensus$n_reps + 1) {
      stop_log(
        "`consensus.high_confidence.hits_in_X_reps` cannot be larger than ",
        "`consensus.n_reps + 1`.\n",
        "Current values:\n",
        "  consensus.n_reps: ",
        cfg$consensus$n_reps,
        "\n  consensus.high_confidence.hits_in_X_reps: ",
        cfg$consensus$high_confidence$hits_in_X_reps
      )
    }
    
    if (cfg$consensus$explorative$hits_in_X_reps >
        cfg$consensus$n_reps + 1) {
      stop_log(
        "`consensus.explorative.hits_in_X_reps` cannot be larger than ",
        "`consensus.n_reps + 1`.\n",
        "Current values:\n",
        "  consensus.n_reps: ",
        cfg$consensus$n_reps,
        "\n  consensus.explorative.hits_in_X_reps: ",
        cfg$consensus$explorative$hits_in_X_reps
      )
    }
  }
  
  invisible(TRUE)
}

#-------------------------------------------------------------------------------
# Main public function
#-------------------------------------------------------------------------------

check_input <- function(option, name) {
  if (!is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name)) {
    stop("`name` must be a single non-empty character string.", call. = FALSE)
  }
  
  validator <- .input_validators[[name]]
  
  if (is.null(validator)) {
    warning(
      "No validator registered for config option `", name, "`. ",
      "Input was not checked.",
      call. = FALSE
    )
    return(option)
  }
  
  validator(option, name)
  
  return(option)
}