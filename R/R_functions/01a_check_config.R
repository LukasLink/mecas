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

.validate_numeric <- function(option,
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
  
  if (length(numeric_option) != 1 || is.na(numeric_option)) {
    .fail_input(name, option, "Expected a single numeric value.")
  }
  
  if (!is.null(min) && numeric_option < min) {
    .fail_input(name, option, paste0("Expected a numeric value >= ", min, "."))
  }
  
  if (!is.null(max) && numeric_option > max) {
    .fail_input(name, option, paste0("Expected a numeric value <= ", max, "."))
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

.input_validators <- list(
  # Run behavior
  first_time = function(option, name) {
    .validate_logical(option, name)
  },
  
  start_with = .validate_choice(c("beginning", "read_counting", "MAUDE_analysis", "generate_plots")),
  
  machine = .validate_choice(c("local", "slurm")),
  
  # Paths
  output_folder = function(option, name) {
    .validate_dir_path(option, name, required = TRUE, must_exist = FALSE)
  },
  
  input_folder = function(option, name) {
    .validate_dir_path(option, name, required = TRUE, must_exist = TRUE)
  },
  
  library_path = function(option, name) {
    .validate_file_path(option, name, required = TRUE, must_exist = TRUE)
  },
  
  star_index_folder = function(option, name) {
    .validate_dir_path(option, name, required = TRUE, must_exist = TRUE)
  },
  
  fastq_name_table_xlsx = function(option, name) {
    .validate_file_path(option, name, required = FALSE, must_exist = FALSE, allow_empty = TRUE)
  },
  
  strict_file_match = function(option, name) {
    .validate_logical(option, name)
  },
  
  # Modules
  use_modules = function(option, name) {
    .validate_logical(option, name)
  },
  
  seqtk_module = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  star_module = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  samtools_module = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  umi_tools_module = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  # bcwithqc options
  bcwithqc_bin = function(option, name) {
    .validate_executable(
      option,
      name,
      required = FALSE,
      allow_empty = TRUE
    )
  },
  bcwithqc_dir = function(option, name) {
    .validate_dir_path(option, name, required = FALSE, must_exist = FALSE, allow_empty = TRUE)
  },
  
  bcwithqc_config_path = function(option, name) {
    .validate_file_path(option, name, required = FALSE, must_exist = FALSE, allow_empty = TRUE)
  },
  
  # align_UMI_tools options 
  UMI_regex = function(option, name) {
    .validate_regex_string(option, name, required = FALSE, allow_empty = TRUE)
  },  
  
  # QC filtering
  qc_filtering_run = function(option, name) {
    .validate_logical(option, name)
  },
  
  qc_min_qual = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 0)
  },
  
  qc_qual_offset = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 0)
  },
  
  qc_min_length = function(option, name) {
    .validate_integer(option, name, required = FALSE, min = 1, allow_na = TRUE)
  },
  
  # Lists
  skip_list = function(option, name) {
    .validate_list_or_character(option, name, required = FALSE)
  },
  
  skip_list_sublib = function(option, name) {
    .validate_list_or_character(option, name, required = FALSE)
  },
  
  skip_list_sample = function(option, name) {
    .validate_list_or_character(option, name, required = FALSE)
  },
  
  include_controls_list = function(option, name) {
    .validate_list_or_character(option, name, required = FALSE)
  },
  
  use_only_these_controls_list = function(option, name) {
    .validate_list_or_character(option, name, required = FALSE)
  },
  
  # Analysis options
  read_counting = .validate_choice(c("bcwithqc", "align_UMI_tools")),
  
  data_type = .validate_choice(c("reads", "umis")),
  
  method = function(option, name) {
    .validate_string(option, name, required = TRUE, allow_empty = TRUE)
  },
  
  norm_method = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  combine_for_guide_stats = function(option, name) {
    .validate_string(option, name, required = TRUE, allow_empty = TRUE)
  },
  
  combine_for_gene_stats = function(option, name) {
    .validate_string(option, name, required = TRUE, allow_empty = TRUE)
  },
  
  # Behavior options
  recover_input = function(option, name) {
    .validate_logical(option, name)
  },
  
  subsample_controls = function(option, name) {
    .validate_logical(option, name)
  },
  
  use_custom_bins = function(option, name) {
    .validate_logical(option, name)
  },
  
  same_controls_in_all_sublibraries = function(option, name) {
    .validate_logical(option, name)
  },
  
  # Suffix / thresholds
  extra_suffix = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  upper_lower_percentage = function(option, name) {
    .validate_percentage(option, name, required = TRUE, min = 0, max = 100)
  },
  
  # Filtering / MAUDE behavior
  drop_0s = function(option, name) {
    .validate_logical(option, name)
  },
  
  strict_mode = function(option, name) {
    .validate_logical(option, name)
  },
  
  min_guides_per_gene = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 0)
  },
  
  auto_combine_replicates = function(option, name) {
    .validate_logical(option, name)
  },
  
  # SLURM
  slurm_account = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  slurm_qos = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  slurm_cpus = function(option, name) {
    .validate_integer(option, name, required = FALSE, min = 1)
  },
  
  slurm_mem = function(option, name) {
    .validate_mem_string(option, name, required = FALSE)
  },
  
  slurm_wall_time = function(option, name) {
    .validate_wall_time(option, name, required = FALSE)
  },
  
  slurm_partition = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  slurm_array = function(option, name) {
    .validate_integer(option, name, required = FALSE, min = 1)
  },
  
  slurm_email = function(option, name) {
    .validate_email(option, name, required = FALSE, allow_empty = TRUE)
  },
  
  # Consensus calling
  consensus_run = function(option, name) {
    .validate_logical(option, name)
  },
  
  consensus_n_reps = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 1)
  },
  
  consensus_high_confidence_FDR_threshold = function(option, name) {
    .validate_percentage(option, name, required = TRUE, min = 0, max = 1)
  },
  
  consensus_high_confidence_hits_in_X_reps = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 1)
  },
  
  consensus_high_confidence_correlation_heatmap = function(option, name) {
    .validate_logical(option, name)
  },
  
  consensus_high_confidence_overlap = function(option, name) {
    .validate_logical(option, name)
  },
  
  consensus_high_confidence_venn_diagram = function(option, name) {
    .validate_logical(option, name)
  },
  
  consensus_explorative_FDR_threshold = function(option, name) {
    .validate_percentage(option, name, required = TRUE, min = 0, max = 1)
  },
  
  consensus_explorative_hits_in_X_reps = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 1)
  },
  
  consensus_explorative_correlation_heatmap = function(option, name) {
    .validate_logical(option, name)
  },
  
  consensus_explorative_overlap = function(option, name) {
    .validate_logical(option, name)
  },
  
  consensus_explorative_venn_diagram = function(option, name) {
    .validate_logical(option, name)
  },
  
  # Waterfall plot
  plots_waterfall_mark_cntrl = function(option, name) {
    .validate_logical(option, name)
  },
  
  plots_waterfall_mark_special = function(option, name) {
    .validate_list_or_character(option, name, required = FALSE)
  },
  
  plots_waterfall_mark_N_top_hits = function(option, name) {
    .validate_integer(option, name, required = TRUE, min = 0)
  },
  
  plots_waterfall_box_padding = function(option, name) {
    .validate_numeric(option, name, required = TRUE, min = 0)
  },
  
  plots_waterfall_no_text = function(option, name) {
    .validate_logical(option, name)
  },
  
  plots_waterfall_signif_lines = function(option, name) {
    .validate_logical(option, name)
  },
  
  plots_waterfall_mark_all_signif_level = function(option, name) {
    .validate_numeric(option, name, required = FALSE, min = 0, max = 1, allow_na = TRUE)
  },
  
  plots_waterfall_break_in_plot = function(option, name) {
    .validate_numeric_vector(option, name, required = FALSE, min = 0)
  },
  
  plots_waterfall_top_padding = function(option, name) {
    .validate_numeric(option, name, required = TRUE, min = 0)
  },
  
  plots_waterfall_custom_title = function(option, name) {
    .validate_string(option, name, required = FALSE, allow_empty = TRUE, allow_null = TRUE)
  },
  
  plots_waterfall_width = function(option, name) {
    .validate_numeric(option, name, required = TRUE, min = 1)
  },
  
  plots_waterfall_height = function(option, name) {
    .validate_numeric(option, name, required = TRUE, min = 1)
  },
  
  plots_waterfall_file_format = .validate_choice(c("png", "pdf", "svg"))
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

check_config_dependencies <- function(opt) {
  #-----------------------------------------------------------------------------
  # Basic required options for dependency checks
  #-----------------------------------------------------------------------------
  
  machine <- opt$machine
  read_counting <- opt$read_counting
  
  check_input(machine, "machine")
  check_input(read_counting, "read_counting")
  
  #-----------------------------------------------------------------------------
  # SLURM options are only required when machine == "slurm"
  #-----------------------------------------------------------------------------
  
  is_slurm <- identical(machine, "slurm")
  
  .require_if(
    condition = is_slurm,
    option = opt$slurm_cpus,
    name = "slurm_cpus",
    reason = "`run.machine` is set to 'slurm'.",
    validator = function(option, name) {
      .validate_integer(option, name, required = TRUE, min = 1)
    }
  )
  
  .require_if(
    condition = is_slurm,
    option = opt$slurm_mem,
    name = "slurm_mem",
    reason = "`run.machine` is set to 'slurm'.",
    validator = function(option, name) {
      .validate_mem_string(option, name, required = TRUE)
    }
  )
  
  .require_if(
    condition = is_slurm,
    option = opt$slurm_wall_time,
    name = "slurm_wall_time",
    reason = "`run.machine` is set to 'slurm'.",
    validator = function(option, name) {
      .validate_wall_time(option, name, required = TRUE)
    }
  )
  
  .require_if(
    condition = is_slurm,
    option = opt$slurm_array,
    name = "slurm_array",
    reason = "`run.machine` is set to 'slurm'.",
    validator = function(option, name) {
      .validate_integer(option, name, required = TRUE, min = 1)
    }
  )
  
  # Optional SLURM options: validate if provided
  if (!.is_null_or_empty(opt$slurm_account)) {
    check_input(opt$slurm_account, "slurm_account")
  }
  
  if (!.is_null_or_empty(opt$slurm_qos)) {
    check_input(opt$slurm_qos, "slurm_qos")
  }
  
  if (!.is_null_or_empty(opt$slurm_partition)) {
    check_input(opt$slurm_partition, "slurm_partition")
  }
  
  if (!.is_null_or_empty(opt$slurm_email)) {
    check_input(opt$slurm_email, "slurm_email")
  }
  
  .check_allowed_when(
    condition = is_slurm,
    option = opt$slurm_account,
    name = "slurm_account",
    reason = "`run.machine` is not set to 'slurm'."
  )
  
  .check_allowed_when(
    condition = is_slurm,
    option = opt$slurm_qos,
    name = "slurm_qos",
    reason = "`run.machine` is not set to 'slurm'."
  )
  
  .check_allowed_when(
    condition = is_slurm,
    option = opt$slurm_partition,
    name = "slurm_partition",
    reason = "`run.machine` is not set to 'slurm'."
  )
  
  .check_allowed_when(
    condition = is_slurm,
    option = opt$slurm_email,
    name = "slurm_email",
    reason = "`run.machine` is not set to 'slurm'."
  )
  
  #-----------------------------------------------------------------------------
  # bcwithqc-specific requirements
  #-----------------------------------------------------------------------------
  
  uses_bcwithqc <- identical(read_counting, "bcwithqc")
  
  if (uses_bcwithqc) {
    has_bcwithqc_dir <- !.is_null_or_empty(opt$bcwithqc_dir)
    has_bcwithqc_config <- !.is_null_or_empty(opt$bcwithqc_config_path)
    
    if (!has_bcwithqc_dir && !has_bcwithqc_config) {
      stop(
        "Missing required bcwithqc configuration.\n",
        "Reason: `read_counting` is set to 'bcwithqc'.\n",
        "Provide at least one of:\n",
        "  - `bcwithqc_dir`\n",
        "  - `bcwithqc_config_path`",
        call. = FALSE
      )
    }
    
    .require_if(
      condition = TRUE,
      option = opt$bcwithqc_bin,
      name = "bcwithqc_bin",
      reason = "`read_counting` is set to 'bcwithqc'.",
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
        opt$bcwithqc_dir,
        "bcwithqc_dir",
        required = TRUE,
        must_exist = TRUE
      )
    }
    
    if (has_bcwithqc_config) {
      .validate_file_path(
        opt$bcwithqc_config_path,
        "bcwithqc_config_path",
        required = TRUE,
        must_exist = TRUE
      )
    }
  }

  #-----------------------------------------------------------------------------
  # Module dependencies
  #-----------------------------------------------------------------------------
  
  if (isTRUE(opt$use_modules)) {
    
    if (isTRUE(opt$qc_filtering_run)) {
      .require_if(
        condition = TRUE,
        option = opt$seqtk_module,
        name = "seqtk_module",
        reason = "`modules.use_modules` is true and `qc_filtering.run` is true.",
        validator = function(option, name) {
          .validate_string(option, name, required = TRUE, allow_empty = FALSE)
        }
      )
    }
    
    if (identical(opt$read_counting, "align_UMI_tools")) {
      .require_if(
        condition = TRUE,
        option = opt$star_module,
        name = "star_module",
        reason = "`modules.use_modules` is true and `read_counting` is 'align_UMI_tools'.",
        validator = function(option, name) {
          .validate_string(option, name, required = TRUE, allow_empty = FALSE)
        }
      )
      
      .require_if(
        condition = TRUE,
        option = opt$samtools_module,
        name = "samtools_module",
        reason = "`modules.use_modules` is true and `read_counting` is 'align_UMI_tools'.",
        validator = function(option, name) {
          .validate_string(option, name, required = TRUE, allow_empty = FALSE)
        }
      )
      
      if (identical(opt$data_type, "umis")) {
        .require_if(
          condition = TRUE,
          option = opt$umi_tools_module,
          name = "umi_tools_module",
          reason = "`modules.use_modules` is true, `read_counting` is 'align_UMI_tools', and `data_type` is 'umis'.",
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
  
  uses_umi_tools <- identical(opt$read_counting, "align_UMI_tools")
  
  .require_if(
    condition = uses_umi_tools,
    option = opt$UMI_regex,
    name = "UMI_regex",
    reason = "`read_counting` is set to 'align_UMI_tools'.",
    validator = function(option, name) {
      .validate_regex_string(option, name, required = TRUE, allow_empty = FALSE)
    }
  )

  
  #-----------------------------------------------------------------------------
  # data_type dependency
  #-----------------------------------------------------------------------------
  
  if (isTRUE(opt$consensus_run) && identical(opt$start_with, "generate_plots")) {
    warning(
      "For `consensus` `run` is set to 'true' but `start_with` is 'generate_plots'. ",
      "Consensus calling requires the pipeline to start with at least the MAUDE run.\n",
      " -> Consesus Calling will be skipped! \nSet `start_with` to `MAUDE_analysis`",
      " or earlier to actually run Consensus Calling.",
      call. = FALSE
    )
  }
  
  #-----------------------------------------------------------------------------
  # consensus_call dependency
  #-----------------------------------------------------------------------------
  if (isTRUE(opt$consensus_run)) {
    
    if (opt$consensus_high_confidence_hits_in_X_reps >
        opt$consensus_n_reps + 1) {
      stop_log(
        "`consensus_high_confidence_hits_in_X_reps` cannot be larger than ",
        "`consensus_n_reps + 1`.\n",
        "Current values:\n",
        "  consensus_n_reps: ",
        opt$consensus_n_reps,
        "\n  consensus_high_confidence_hits_in_X_reps: ",
        opt$consensus_high_confidence_hits_in_X_reps
      )
    }
    
    if (opt$consensus_explorative_hits_in_X_reps >
        opt$consensus_n_reps + 1) {
      stop_log(
        "`consensus_explorative_hits_in_X_reps` cannot be larger than ",
        "`consensus_n_reps + 1`.\n",
        "Current values:\n",
        "  consensus_n_reps: ",
        opt$consensus_n_reps,
        "\n  consensus_explorative_hits_in_X_reps: ",
        opt$consensus_explorative_hits_in_X_reps
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