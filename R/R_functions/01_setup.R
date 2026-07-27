# project_setup.R

# ---- global knitr / options ----
knitr::opts_chunk$set(echo = FALSE)
options(bitmapType = "cairo")

# ---- packages ----

cran_packages <- c(
  "tidyverse",
  "Matrix",
  "conflicted",
  "ggplot2",
  "ggrepel",
  "writexl",
  "stringr",
  "ggbreak",
  "yaml",
  "tools",
  "logger"
)

bioc_packages <- c(
  "AnnotationDbi",
  "org.Hs.eg.db"
)

other_packages <- c(
  "MAUDE"
)

check_packages <- function(packages, source = "CRAN") {
  missing <- packages[
    !vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]
  
  if (length(missing) > 0) {
    message("Missing ", source, " package(s): ", paste(missing, collapse = ", "))
    
    if (source == "CRAN") {
      message("Install with:")
      message("install.packages(c(", paste(sprintf('"%s"', missing), collapse = ", "), "))")
    }
    
    if (source == "Bioconductor") {
      message("Install with:")
      message("BiocManager::install(c(", paste(sprintf('"%s"', missing), collapse = ", "), "))")
    }
    
    stop("Required package(s) missing.", call. = FALSE)
  }
}

check_packages(cran_packages, "CRAN")
check_packages(bioc_packages, "Bioconductor")
check_packages(other_packages, "other source")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(conflicted)
  library(MAUDE)
  library(ggplot2)
  library(ggrepel)
  library(writexl)
  library(stringr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggbreak)
  library(yaml)
  library(tools)
  library(logger)
})

suppressMessages({
  conflicted::conflicts_prefer(dplyr::rename)
  conflicted::conflicts_prefer(dplyr::filter)
  conflicted::conflicts_prefer(dplyr::select)
  conflicted::conflicts_prefer(dplyr::slice)
  conflicted::conflicts_prefer(dplyr::first)
  conflicted::conflicts_prefer(dplyr::desc)
  conflicted::conflicts_prefer(base::setdiff)
  conflicted::conflicts_prefer(base::intersect)
  conflicted::conflicts_prefer(base::unname)
  conflicted::conflicts_prefer(base::setequal)
})

#-------------------------------------------------------------------------------
# Small helpers
#-------------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

parse_comma_list <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(character())
  }
  
  x <- as.character(x)
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(character())
  }
  
  out <- unlist(strsplit(x, ",", fixed = TRUE))
  out <- trimws(out)
  out <- out[nzchar(out)]
  out
}

is_named_list <- function(x) {
  is.list(x) && !is.null(names(x))
}

deep_modify <- function(x, y) {
  if (is.null(y) || length(y) == 0) {
    return(x)
  }
  
  y <- y[!vapply(y, is.null, logical(1))]
  
  for (nm in names(y)) {
    if (is_named_list(x[[nm]]) && is_named_list(y[[nm]])) {
      x[[nm]] <- deep_modify(x[[nm]], y[[nm]])
    } else {
      x[[nm]] <- y[[nm]]
    }
  }
  
  x
}

params_to_overrides <- function(params) {
  if (is.null(params) || length(params) == 0) {
    return(list())
  }
  
  params <- as.list(params)
  
  # `config` is only used to locate the YAML file.
  params$config <- NULL
  
  # NULL means: do not override YAML.
  params <- params[!vapply(params, is.null, logical(1))]
  
  params
}

require_yaml_config <- function(x) {
  if (!is.list(x)) {
    stop("The parsed YAML config must be a list.", call. = FALSE)
  }
  
  required_top_level <- c("paths")
  
  missing <- required_top_level[!required_top_level %in% names(x)]
  
  if (length(missing) > 0) {
    stop(
      "Missing required top-level config section(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  
  invisible(x)
}

parse_qc_min_length <- function(option, name = "qc_filtering.min_length") {
  empty_result <- list(
    min_length_state = "infer",
    min_length_single = NA_integer_,
    min_length_R1 = NA_integer_,
    min_length_R2 = NA_integer_
  )
  
  # NULL or NA means infer from FASTQ files.
  if (is.null(option) || length(option) == 0L || (!is.list(option) && length(option) == 1L && is.na(option))) {
    return(empty_result)
  }
  
  # One scalar integer applies to all reads.
  if (!is.list(option) && length(option) == 1L) {
    .validate_integer(option = option, name = name, required = TRUE, min = 1)
    
    return(list(
      min_length_state = "provided_single",
      min_length_single = as.integer(option),
      min_length_R1 = NA_integer_,
      min_length_R2 = NA_integer_
    ))
  }
  
  # Separate R1/R2 configuration.
  if (!is.list(option)) {.fail_input(name, option, paste0(
        "Expected NULL, one integer, or a named mapping with ",
        "`R1` and `R2`."
      )
    )
  }
  
  option_names <- names(option)
  
  if (is.null(option_names) || any(is.na(option_names)) || any(!nzchar(option_names))) {
    .fail_input(name, option, "Separate paired-end lengths must be named `R1` and `R2`.")
  }
  
  allowed_names <- c("R1", "R2")
  unknown_names <- setdiff(option_names, allowed_names)
  
  if (length(unknown_names) > 0L) {
    .fail_input(
      name,
      option,
      paste0(
        "Unknown minimum-length field(s): ",
        paste(unknown_names, collapse = ", "),
        ". Expected only `R1` and `R2`."
      )
    )
  }
  
  missing_names <- setdiff(allowed_names, option_names)
  
  if (length(missing_names) > 0L) {
    .fail_input(
      name,
      option,
      paste0(
        "Separate paired-end minimum lengths require both `R1` and `R2`. ",
        "Missing: ",
        paste(missing_names, collapse = ", "),
        "."
      )
    )
  }
  
  .validate_integer(
    option = option$R1,
    name = paste0(name, ".R1"),
    required = TRUE,
    min = 1
  )
  
  .validate_integer(
    option = option$R2,
    name = paste0(name, ".R2"),
    required = TRUE,
    min = 1
  )
  
  list(
    min_length_state = "provided_R1R2",
    min_length_single = NA_integer_,
    min_length_R1 = as.integer(option$R1),
    min_length_R2 = as.integer(option$R2)
  )
}
#-------------------------------------------------------------------------------
# Logging
#-------------------------------------------------------------------------------

log_project_setup <- function(cfg) {
  logger::log_info("==========================================================")
  logger::log_info("Project setup options")
  logger::log_info("==========================================================")

  logger::log_info("setup_mode:                        {cfg$run$setup_mode}")
  logger::log_info("----------------------------------------------------------")
  
  logger::log_info("output_folder:                     {cfg$paths$output_folder}")
  logger::log_info("input_folder:                      {cfg$paths$input_folder}")
  logger::log_info("library_path:                      {cfg$paths$library_path}")
  logger::log_info("bcwithqc_config_path:              {cfg$paths$bcwithqc_config_path}")
  logger::log_info("fastq_name_table_xlsx:             {cfg$paths$fastq_name_table_xlsx}")
  logger::log_info("strict_file_match:                 {cfg$paths$strict_file_match}")
  logger::log_info("----------------------------------------------------------")

  logger::log_info("skip_list:                         {paste(cfg$skip$files, collapse = ', ')}")
  logger::log_info("skip_list_sublib:                  {paste(cfg$skip$sublibraries, collapse = ', ')}")
  logger::log_info("skip_list_sample:                  {paste(cfg$skip$samples, collapse = ', ')}")
  logger::log_info("include_controls_list:             {paste(cfg$controls$include_controls, collapse = ', ')}")
  logger::log_info(
    "use_only_these_controls_list:      {x}",
    x = if (length(cfg$controls$use_only_these_controls) > 0) {
      paste(cfg$controls$use_only_these_controls, collapse = ", ")
    } else {
      "<not set>"
    }
  )
  logger::log_info("----------------------------------------------------------")
  
  logger::log_info("data_type:                         {cfg$counting$data_type}")
  logger::log_info("method:                            {cfg$replicates$method}")
  logger::log_info("norm_method:                       {cfg$normalization$norm_method}")
  logger::log_info("combine_for_guide_stats:           {cfg$replicates$combine_for_guide_stats}")
  logger::log_info("combine_for_gene_stats:            {cfg$replicates$combine_for_gene_stats}")
  logger::log_info("----------------------------------------------------------")
  
  logger::log_info("recover_input:                     {cfg$normalization$recover_input}")
  logger::log_info("subsample_controls:                {cfg$controls$subsample_controls}")
  logger::log_info("use_custom_bins:                   {cfg$normalization$use_custom_bins}")
  logger::log_info("same_controls_in_all_sublibraries: {cfg$controls$same_controls_in_all_sublibraries}")
  logger::log_info("----------------------------------------------------------")
  
  logger::log_info("drop_0s:                           {cfg$filtering$drop_0s}")
  logger::log_info("strict_mode:                       {cfg$filtering$strict_mode}")
  logger::log_info("min_guides_per_gene:               {cfg$filtering$min_guides_per_gene}")
  logger::log_info("auto_combine_replicates:           {cfg$replicates$auto_combine_replicates}")
  logger::log_info("----------------------------------------------------------")
  
  logger::log_info("file_suffix:                       {cfg$suffix$file_suffix}")
  logger::log_info("file_info_suffix:                  {cfg$suffix$file_info_suffix}")
  logger::log_info("log_file:                          {cfg$paths$log_file}")
  
  invisible(cfg)
}

#-------------------------------------------------------------------------------
# Project setup
#-------------------------------------------------------------------------------

project_setup <- function(project_root_dir,
                          config = NULL,
                          config_path = NULL,
                          params = NULL,
                          setup_mode = "count",
                          overrides = list(),
                          envir = .GlobalEnv,
                          only_one_logger = TRUE) {
  stopifnot(is.character(project_root_dir), length(project_root_dir) == 1)
  
  allowed_setup_modes <- c(
    "setup",
    "QC_filtering",
    "count",
    "MAUDE",
    "plot"
  )
  
  if (!(setup_mode %in% allowed_setup_modes)) {
    stop(
      paste0(
        "Unsupported setup_mode: ",
        setup_mode,
        ". This is a programming error; the user cannot fix this."
      ),
      call. = FALSE
    )
  }
  
  if (is.null(config_path) && !is.null(params) && !is.null(params$config)) {
    config_path <- params$config
  }
  
  if (!is.null(config_path)) {
    config_path <- normalizePath(config_path, mustWork = TRUE)
    message("Using config file: ", config_path)
    config <- yaml::read_yaml(config_path)
  }
  
  if (is.null(config)) {
    stop("Provide either `config` or `config_path`.", call. = FALSE)
  }
  
  require_yaml_config(config)
  
  # Rmd params first, then explicit overrides. Explicit overrides win.
  # These should be nested lists, e.g.
  # overrides = list(output = list(extra_suffix = "run_1"))
  config <- deep_modify(config, params_to_overrides(params))
  config <- deep_modify(config, overrides)
  
  #=============================================================================
  # Validate / normalize config into cfg
  #=============================================================================
  
  cfg <- list()
  cfg$run <- list(
    setup_mode = setup_mode,
    run_id = format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  )
  
  cfg$paths <- list(
    project_root_dir = project_root_dir,
    
    output_folder = check_input(config$paths$output_folder, "paths.output_folder"),
    input_folder = check_input(config$paths$input_folder %||% "", "paths.input_folder"),
    library_path = check_input(config$paths$library_path %||% "", "paths.library_path"),
    fastq_name_table_xlsx = check_input(
      config$paths$fastq_name_table_xlsx %||% NULL,
      "paths.fastq_name_table_xlsx"
    ),
    bcwithqc_config_path = check_input(
      config$paths$bcwithqc_config_path %||% "",
      "paths.bcwithqc_config_path"
    ),
    strict_file_match = check_input(
      config$paths$strict_file_match %||% TRUE,
      "paths.strict_file_match"
    )
  )
  
  if (is.null(cfg$paths$output_folder) || !nzchar(cfg$paths$output_folder)) {
    stop("`paths.output_folder` must be provided in the YAML config.", call. = FALSE)
  }

  cfg$skip <- list(
    files = check_input(
      parse_comma_list(config$skip$files %||% c()),
      "skip.files"
    ),
    sublibraries = check_input(
      parse_comma_list(config$skip$sublibraries %||% c()),
      "skip.sublibraries"
    ),
    samples = check_input(
      parse_comma_list(config$skip$samples %||% c()),
      "skip.samples"
    )
  )
  
  cfg$controls <- list(
    include_controls = check_input(
      parse_comma_list(config$controls$include_controls %||% c()),
      "controls.include_controls"
    ),
    use_only_these_controls = check_input(
      parse_comma_list(config$controls$use_only_these_controls %||% c()),
      "controls.use_only_these_controls"
    ),
    same_controls_in_all_sublibraries = check_input(
      config$controls$same_controls_in_all_sublibraries %||% TRUE,
      "controls.same_controls_in_all_sublibraries"
    ),
    subsample_controls = check_input(
      config$controls$subsample_controls %||% FALSE,
      "controls.subsample_controls"
    )
  )
  
  cfg$counting <- list(
    data_type = check_input(
      config$counting$data_type %||% "reads",
      "counting.data_type"
    )
  )
  
  qc_min_length <- parse_qc_min_length(
    config$qc_filtering$min_length %||% NA
  )
  
  cfg$qc_filtering <- list(
    run = check_input(config$qc_filtering$run %||% TRUE, "qc_filtering.run"),
    min_qual = check_input(config$qc_filtering$min_qual %||% 20, "qc_filtering.min_qual"),
    min_length = check_input(config$qc_filtering$min_length %||% NA, "qc_filtering.min_length"),
    min_length_state = qc_min_length$min_length_state,
    min_length_single = qc_min_length$min_length_single,
    min_length_R1 = qc_min_length$min_length_R1,
    min_length_R2 = qc_min_length$min_length_R2
  )
  
  cfg$replicates <- list(
    method = check_input(config$replicates$method %||% "", "replicates.method"),
    combine_for_guide_stats = check_input(
      config$replicates$combine_for_guide_stats %||% "sample",
      "replicates.combine_for_guide_stats"
    ),
    combine_for_gene_stats = check_input(
      config$replicates$combine_for_gene_stats %||% "none",
      "replicates.combine_for_gene_stats"
    ),
    auto_combine_replicates = check_input(
      config$replicates$auto_combine_replicates %||% FALSE,
      "replicates.auto_combine_replicates"
    )
  )
  
  cfg$normalization <- list(
    norm_method = check_input(
      config$normalization$norm_method %||% "control_median",
      "normalization.norm_method"
    ),
    recover_input = check_input(
      config$normalization$recover_input %||% TRUE,
      "normalization.recover_input"
    ),
    use_custom_bins = check_input(
      config$normalization$use_custom_bins %||% FALSE,
      "normalization.use_custom_bins"
    ),
    upper_lower_percentage = check_input(
      config$normalization$upper_lower_percentage %||% 0.10,
      "normalization.upper_lower_percentage"
    )
  )
  
  cfg$filtering <- list(
    drop_0s = check_input(config$filtering$drop_0s %||% FALSE, "filtering.drop_0s"),
    strict_mode = check_input(config$filtering$strict_mode %||% FALSE, "filtering.strict_mode"),
    min_guides_per_gene = check_input(
      config$filtering$min_guides_per_gene %||% 0,
      "filtering.min_guides_per_gene"
    )
  )
  
  cfg$consensus <- list(
    run = check_input(config$consensus$run %||% TRUE, "consensus.run"),
    n_reps = check_input(config$consensus$n_reps %||% 19, "consensus.n_reps"),
    
    high_confidence = list(
      FDR_threshold = check_input(
        config$consensus$high_confidence$FDR_threshold %||% 0.05,
        "consensus.high_confidence.FDR_threshold"
      ),
      hits_in_X_reps = check_input(
        config$consensus$high_confidence$hits_in_X_reps %||% 16,
        "consensus.high_confidence.hits_in_X_reps"
      ),
      correlation_heatmap = check_input(
        config$consensus$high_confidence$correlation_heatmap %||% FALSE,
        "consensus.high_confidence.correlation_heatmap"
      ),
      overlap = check_input(
        config$consensus$high_confidence$overlap %||% FALSE,
        "consensus.high_confidence.overlap"
      ),
      venn_diagram = check_input(
        config$consensus$high_confidence$venn_diagram %||% FALSE,
        "consensus.high_confidence.venn_diagram"
      )
    ),
    
    explorative = list(
      FDR_threshold = check_input(
        config$consensus$explorative$FDR_threshold %||% 0.05,
        "consensus.explorative.FDR_threshold"
      ),
      hits_in_X_reps = check_input(
        config$consensus$explorative$hits_in_X_reps %||% 10,
        "consensus.explorative.hits_in_X_reps"
      ),
      correlation_heatmap = check_input(
        config$consensus$explorative$correlation_heatmap %||% TRUE,
        "consensus.explorative.correlation_heatmap"
      ),
      overlap = check_input(
        config$consensus$explorative$overlap %||% FALSE,
        "consensus.explorative.overlap"
      ),
      venn_diagram = check_input(
        config$consensus$explorative$venn_diagram %||% FALSE,
        "consensus.explorative.venn_diagram"
      )
    )
  )
  
  cfg$output <- list(
    extra_suffix = check_input(config$output$extra_suffix %||% "", "output.extra_suffix")
  )
  
  cfg$plots <- list(
    run_after_count = check_input(
      config$plots$run_after_count %||% TRUE,
      "plots.run_after_count"
    ),
    
    read_count_violin = list(
      violin_y_limit = check_input(
        config$plots$read_count_violin$violin_y_limit %||% 8000,
        "plots.read_count_violin.violin_y_limit"
      )
    ),
    
    waterfall = list(
      mark_cntrl = check_input(
        config$plots$waterfall$mark_cntrl %||% TRUE,
        "plots.waterfall.mark_cntrl"
      ),
      mark_special = check_input(
        config$plots$waterfall$mark_special %||% NULL,
        "plots.waterfall.mark_special"
      ),
      mark_N_top_hits = check_input(
        config$plots$waterfall$mark_N_top_hits %||% 3,
        "plots.waterfall.mark_N_top_hits"
      ),
      box_padding = check_input(
        config$plots$waterfall$box_padding %||% 0.8,
        "plots.waterfall.box_padding"
      ),
      no_text = check_input(
        config$plots$waterfall$no_text %||% FALSE,
        "plots.waterfall.no_text"
      ),
      signif_lines = check_input(
        config$plots$waterfall$signif_lines %||% TRUE,
        "plots.waterfall.signif_lines"
      ),
      mark_all_signif_level = check_input(
        config$plots$waterfall$mark_all_signif_level %||% 0.05,
        "plots.waterfall.mark_all_signif_level"
      ),
      break_in_plot = check_input(
        config$plots$waterfall$break_in_plot %||% c(),
        "plots.waterfall.break_in_plot"
      ),
      top_padding = check_input(
        config$plots$waterfall$top_padding %||% 0,
        "plots.waterfall.top_padding"
      ),
      custom_title = check_input(
        config$plots$waterfall$custom_title %||% NULL,
        "plots.waterfall.custom_title"
      ),
      width = check_input(
        config$plots$waterfall$width %||% 8,
        "plots.waterfall.width"
      ),
      height = check_input(
        config$plots$waterfall$height %||% 6,
        "plots.waterfall.height"
      ),
      file_format = check_input(
        config$plots$waterfall$file_format %||% "png",
        "plots.waterfall.file_format"
      )
    )
  )
  
  #=============================================================================
  # Construct suffixes
  #=============================================================================
  
  cfg$suffix <- list(
    recover_input_suffix = if (cfg$normalization$recover_input) "RI" else "",
    subsample_controls_suffix = if (cfg$controls$subsample_controls) "ss_cntrl" else "",
    custom_bins_suffix = if (cfg$normalization$use_custom_bins) "custom_bins" else "",
    drop_0s_suffix = if (cfg$filtering$drop_0s) "D0" else "",
    strict_mode_suffix = if (cfg$filtering$strict_mode) "strict" else "",
    auto_combine_replicates_suffix = if (cfg$replicates$auto_combine_replicates) "acr" else "",
    min_guides_per_gene_suffix = if (cfg$filtering$min_guides_per_gene > 0) {
      paste0("min_guides_", cfg$filtering$min_guides_per_gene)
    } else {
      ""
    },
    combine_for_guide_stats_suffix = if (cfg$replicates$combine_for_guide_stats == "") {
      ""
    } else {
      paste0("comb_", cfg$replicates$combine_for_guide_stats)
    },
    combine_for_gene_stats_suffix = if (cfg$replicates$combine_for_gene_stats == "") {
      ""
    } else {
      paste0("comb_", cfg$replicates$combine_for_gene_stats)
    }
  )
  
  skip_list_and_suffix <- create_skip_list_and_suffix(
    cfg$skip$files,
    cfg$skip$sublibraries,
    cfg$skip$samples
  )
  
  cfg$skip$files <- skip_list_and_suffix[[1]]
  cfg$suffix$skip_suffix <- skip_list_and_suffix[[2]]
  

  fs_parts <- c(
    cfg$counting$data_type,
    cfg$replicates$method,
    cfg$normalization$norm_method,
    cfg$suffix$recover_input_suffix,
    cfg$suffix$subsample_controls_suffix,
    cfg$suffix$drop_0s_suffix,
    cfg$suffix$strict_mode_suffix,
    cfg$suffix$custom_bins_suffix,
    cfg$suffix$min_guides_per_gene_suffix,
    cfg$suffix$combine_for_guide_stats_suffix,
    cfg$suffix$combine_for_gene_stats_suffix,
    cfg$suffix$auto_combine_replicates_suffix,
    cfg$suffix$skip_suffix,
    cfg$output$extra_suffix
  )
 
  
  fs_parts <- fs_parts[!is.na(fs_parts) & fs_parts != ""]
  
  cfg$suffix$file_suffix <- paste0("_", paste(fs_parts, collapse = "_"), ".rds")
  cfg$suffix$file_info_suffix <- paste(fs_parts, collapse = "_")
  
  #=============================================================================
  # Construct paths
  #=============================================================================
  
  cfg$paths$data_dir <- get_file_path(project_root_dir, "data")
  
  cfg$paths$bcwithqc_output_folder <- make_clean_dir(cfg$paths$output_folder, "bcwithqc_output")
  cfg$paths$qc_filtered_folder <- make_clean_dir(cfg$paths$output_folder, "QC_filtered")
  cfg$paths$rds_output_folder <- make_clean_dir(cfg$paths$output_folder, "rds")
  cfg$paths$results_output_folder <- make_clean_dir(cfg$paths$output_folder, "results")
  
  cfg$paths$plots_output_folder <- make_clean_dir(cfg$paths$output_folder, "plots")
  
  cfg$paths$fastq_symlinks_folder <- make_clean_dir(cfg$paths$output_folder, "fastq_symlinks")
  cfg$paths$bcwithqc_symlinks_folder <- make_clean_dir(cfg$paths$output_folder, "bcwithqc_symlinks")
  
  cfg$paths$config_path <- config_path %||% ""
  
  cfg$paths$count_df_fpath <- file.path(cfg$paths$rds_output_folder, "count_df.rds")
  #=============================================================================
  # Snakemake handeling
  #=============================================================================
  cfg$paths$snake <- list(
    state_dir = make_clean_dir(cfg$paths$output_folder, ".pipeline_state")
    )
  
  cfg$paths$manifest <- file.path(cfg$paths$snake$state_dir,"fastq_manifest.tsv")
  cfg$paths$snake$resolved_config_rds <- file.path(cfg$paths$snake$state_dir, "resolved_config.rds")
  cfg$paths$snake$resolved_config_yaml <- file.path(cfg$paths$snake$state_dir, "resolved_config.yaml")
  cfg$paths$snake$QC_filtering_params_sh <- file.path(cfg$paths$snake$state_dir,"QC_filtering_params.sh")
  
  cfg$paths$snake$done <- list(
    setup = file.path(cfg$paths$snake$state_dir, "01_setup.done"),
    infer_QC_filter_params = file.path(cfg$paths$snake$state_dir, "02_infer_QC_filter_params.done"),
    count = file.path(cfg$paths$snake$state_dir, "06_count.done"),
    MAUDE = file.path(cfg$paths$snake$state_dir, "07_MAUDE.done"),
    plot = file.path(cfg$paths$snake$state_dir, "08_plot.done")
  )
  #=============================================================================
  # Logging 
  #=============================================================================
  cfg$paths$log_folder <- make_clean_dir(cfg$paths$output_folder, "logs")
  # This is for the old logging version, which has one log file for everything
  if (isTRUE(only_one_logger)){
    cfg$paths$log_file <- file.path(
      cfg$paths$log_folder,
      paste0(
        "log_",
        cfg$suffix$file_info_suffix,
        "_",
        format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),
        ".log"
      )
    )
    
    initialize_pipeline_logger(cfg$paths$log_file)
    
  } else {
    # This is for the new snakemake log version. 
    cfg$paths$log_files <- list(
      setup = file.path(cfg$paths$log_folder, paste0(cfg$run$run_id, "_01_setup.log")),
      infer_QC_filter_params = file.path(cfg$paths$log_folder, paste0(cfg$run$run_id, "_02_infer_QC_filter_params.log")),
      count = file.path(cfg$paths$log_folder, paste0(cfg$run$run_id, "_06_count.log")),
      MAUDE = file.path(cfg$paths$log_folder, paste0(cfg$run$run_id, "_07_MAUDE.log")),
      plot = file.path(cfg$paths$log_folder, paste0(cfg$run$run_id, "_08_plot.log"))
    )
    if (is.null(cfg$paths$log_files[[setup_mode]])) {
      stop("No log file configured for setup_mode: ", setup_mode, call. = FALSE)
    }
    initialize_pipeline_logger(cfg$paths$log_files[[setup_mode]])
  }
  
  #=============================================================================
  # Load library
  #=============================================================================
  
  cfg$merged_sgRNA_df <- read_library_file(
    library_path = cfg$paths$library_path
  )
  
  #=============================================================================
  # Print setup summary
  #=============================================================================
  
  log_project_setup(cfg)
  
  invisible(cfg)
}