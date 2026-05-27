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
  missing <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  
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
  conflicts_prefer(dplyr::rename)
  conflicts_prefer(dplyr::filter)
  conflicts_prefer(dplyr::select)
  conflicts_prefer(dplyr::slice)
  conflicts_prefer(dplyr::first)
  conflicts_prefer(base::setdiff)
  conflicts_prefer(base::intersect)
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
  
  # Handles both single comma-separated strings and direct vector values
  # from YAML, Rmd params, or explicit overrides.
  out <- unlist(strsplit(x, ",", fixed = TRUE))
  out <- trimws(out)
  out <- out[nzchar(out)]
  out
}

apply_overrides <- function(x, overrides = list()) {
  if (is.null(overrides) || length(overrides) == 0) {
    return(x)
  }
  
  overrides <- overrides[!vapply(overrides, is.null, logical(1))]
  
  if (length(overrides) == 0) {
    return(x)
  }
  
  for (nm in names(overrides)) {
    x[[nm]] <- overrides[[nm]]
  }
  
  x
}

params_to_overrides <- function(params) {
  if (is.null(params) || length(params) == 0) {
    return(list())
  }
  
  params <- as.list(params)
  
  # `config` is used to locate the YAML file and should not become a pipeline option.
  params$config <- NULL
  
  # Rmd params default to NULL. NULL means: do not override the YAML value.
  params <- params[!vapply(params, is.null, logical(1))]
  
  params
}

looks_like_yaml_config <- function(x) {
  is.list(x) && any(c("run", "paths", "counting", "resources", "normalization") %in% names(x))
}

looks_like_user_opts <- function(x) {
  is.list(x) && any(c("output_folder", "input_folder", "read_counting", "extra_suffix") %in% names(x))
}

#-------------------------------------------------------------------------------
# Convert config.yaml into the legacy user_opts object
#-------------------------------------------------------------------------------
# This is the single place where YAML names are mapped to the legacy/global names
# used by the rest of the Rmd. Add new YAML options here, not in the Rmd body.

config_to_user_opts <- function(config) {
  stopifnot(is.list(config))
  
  list(
    # Run behavior
    first_time = config$run$first_time %||% FALSE,
    machine = config$run$machine %||% "local",
    
    # Directories / paths
    output_folder = config$paths$output_folder,
    input_folder = config$paths$input_folder %||% "",
    library_path = config$paths$library_path %||% "",
    
    star_index_folder = config$paths$star_index_folder %||% "",
    
    # Optional FASTQ name parser table
    fastq_name_table_xlsx = config$paths$fastq_name_table_xlsx %||% "",
    strict_file_match = config$paths$strict_file_match %||% TRUE,
    
    # Optional bcwithqc paths
    bcwithqc_dir = config$bcwithqc$bcwithqc_dir %||% "",
    bcwithqc_config_path = config$bcwithqc$bcwithqc_config_path %||% "",
    
    # align_UMI_tools options
    UMI_regex = config$align_UMI_tools$umi_regex %||% "",
    
    # Skip lists
    skip_list = config$skip$files %||% c(),
    skip_list_sublib = config$skip$sublibraries %||% c(),
    skip_list_sample = config$skip$samples %||% c(),
    
    # Controls
    include_controls_list = config$controls$include_controls %||% c(),
    use_only_these_controls_list = config$controls$use_only_these_controls %||% c(),
    same_controls_in_all_sublibraries =
      config$controls$same_controls_in_all_sublibraries %||% TRUE,
    subsample_controls = config$controls$subsample_controls %||% FALSE,
    
    # Counting options
    read_counting = config$counting$read_counting %||% "align_UMI_tools",
    data_type = config$counting$data_type %||% "reads",
    
    # SLURM / resource settings
    slurm_account = config$slurm$account %||% "",
    slurm_qos = config$slurm$qos %||% "",
    slurm_cpus = config$slurm$cpus %||% 1,
    slurm_mem = config$slurm$mem %||% "4g",
    slurm_wall_time = config$slurm$wall_time %||% "01:00:00",
    slurm_partition = config$slurm$partition %||% "",
    slurm_array = config$slurm$array %||% 1,
    slurm_email = config$slurm$email %||% "",
    
    # Modules
    use_modules       = config$modules$use_modules %||% FALSE,
    seqtk_module      = config$modules$seqtk %||% "",
    star_module       = config$modules$star %||% "",
    samtools_module   = config$modules$samtools %||% "",
    umi_tools_module  = config$modules$umi_tools %||% "",
    
    # QC filtering
    qc_filtering_run  = config$qc_filtering$run %||% TRUE,
    qc_min_qual       = config$qc_filtering$min_qual %||% 20,
    qc_qual_offset    = config$qc_filtering$qual_offset %||% 20,
    qc_min_length     = config$qc_filtering$min_length %||% NA,
    
    # Replicates / grouping
    method = config$replicates$method %||% "",
    combine_for_guide_stats = config$replicates$combine_for_guide_stats %||% "sample",
    combine_for_gene_stats = config$replicates$combine_for_gene_stats %||% "none",
    auto_combine_replicates = config$replicates$auto_combine_replicates %||% FALSE,
    
    # Normalization / bins
    norm_method = config$normalization$norm_method %||% "control_median",
    recover_input = config$normalization$recover_input %||% TRUE,
    use_custom_bins = config$normalization$use_custom_bins %||% FALSE,
    upper_lower_percentage = config$normalization$upper_lower_percentage %||% 0.10,
    
    # Filtering
    drop_0s = config$filtering$drop_0s %||% FALSE,
    strict_mode = config$filtering$strict_mode %||% FALSE,
    min_guides_per_gene = config$filtering$min_guides_per_gene %||% 0,
    
    # Output naming
    extra_suffix = config$output$extra_suffix %||% ""
  )
}

#-------------------------------------------------------------------------------
# Project setup
#-------------------------------------------------------------------------------
# New recommended usage from the Rmd:
#   project_setup(project_root_dir = project_root_dir,
#                 config_path = config_path,
#                 params = params,
#                 use_old_suffix_construction = TRUE)
#
# For programmatic multiple runs from the same config:
#   project_setup(project_root_dir = project_root_dir,
#                 config_path = config_path,
#                 overrides = list(extra_suffix = "run_1"))
#
# Still supported legacy usage:
#   config <- yaml::read_yaml(config_path)
#   project_setup(project_root_dir, config, use_old_suffix_construction = TRUE)
#   project_setup(project_root_dir, user_opts, use_old_suffix_construction = TRUE)

project_setup <- function(project_root_dir,
                          config_or_user_opts = NULL,
                          config_path = NULL,
                          params = NULL,
                          setup_mode = "count",
                          overrides = list(),
                          envir = .GlobalEnv,
                          use_old_suffix_construction = FALSE) {
  stopifnot(is.character(project_root_dir), length(project_root_dir) == 1)
  
  if (!(setup_mode %in% c("count","make_reference"))){
    stop(paste0(
      "unsupported setup_mode: ",
      setup_mode,
      "this is a programm error, user can do nothing to fix this."),
      call. = FALSE)
  }
  
  # ---- accept either config object, config path, or legacy user_opts ----
  # If config_path was not supplied explicitly, allow params$config to provide it.
  if (is.null(config_path) && !is.null(params) && !is.null(params$config)) {
    config_path <- params$config
  }
  
  if (!is.null(config_path)) {
    config_path <- normalizePath(config_path, mustWork = TRUE)
    message("Using config file: ", config_path)
    config_or_user_opts <- yaml::read_yaml(config_path)
  }
  
  if (is.null(config_or_user_opts)) {
    stop("Provide either `config_or_user_opts` or `config_path`.", call. = FALSE)
  }
  
  if (looks_like_yaml_config(config_or_user_opts)) {
    user_opts <- config_to_user_opts(config_or_user_opts)
  } else if (looks_like_user_opts(config_or_user_opts)) {
    user_opts <- config_or_user_opts
  } else {
    stop(
      "`config_or_user_opts` must be either the parsed YAML config or the legacy user_opts list.",
      call. = FALSE
    )
  }
  
  # Apply Rmd params first, then explicit overrides. Explicit overrides win if both are supplied.
  user_opts <- apply_overrides(user_opts, params_to_overrides(params))
  user_opts <- apply_overrides(user_opts, overrides)
  
  # Backfill defaults, so older YAML/user_opts files do not fail when new fields are added.
  defaults <- config_to_user_opts(list())
  user_opts <- modifyList(defaults, user_opts)
  
  if (is.null(user_opts$output_folder) || !nzchar(user_opts$output_folder)) {
    stop("`output_folder` must be provided in config$paths$output_folder or user_opts$output_folder.", call. = FALSE)
  }
  
  # add if (setup_mode == ) check here. 
  # Use the fully resolved options from YAML + params + overrides.
  opt <- user_opts
  
  #=============================================================================
  # Normalize resolved options and comma-separated/vector list fields
  #=============================================================================
  first_time                        <- check_input(opt$first_time, "first_time")
  machine                           <- check_input(opt$machine, "machine")
  
  output_folder                     <- check_input(opt$output_folder, "output_folder")
  input_folder                      <- check_input(opt$input_folder, "input_folder")
  library_path                      <- check_input(opt$library_path, "library_path")
  star_index_folder                 <- check_input(opt$star_index_folder, "star_index_folder")
  fastq_name_table_xlsx             <- check_input(opt$fastq_name_table_xlsx, "fastq_name_table_xlsx")
  strict_file_match                 <- check_input(opt$strict_file_match, "strict_file_match")
  
  bcwithqc_dir                      <- check_input(opt$bcwithqc_dir, "bcwithqc_dir")
  bcwithqc_config_path              <- check_input(opt$bcwithqc_config_path, "bcwithqc_config_path")
  UMI_regex                         <- check_input(opt$UMI_regex, "UMI_regex")
  
  skip_list                         <- check_input(parse_comma_list(opt$skip_list), "skip_list")
  skip_list_sublib                  <- check_input(parse_comma_list(opt$skip_list_sublib), "skip_list_sublib")
  skip_list_sample                  <- check_input(parse_comma_list(opt$skip_list_sample), "skip_list_sample")
  include_controls_list             <- check_input(parse_comma_list(opt$include_controls_list), "include_controls_list")
  use_only_these_controls_list      <- check_input(parse_comma_list(opt$use_only_these_controls_list), "use_only_these_controls_list")
  
  read_counting                     <- check_input(opt$read_counting, "read_counting")
  data_type                         <- check_input(opt$data_type, "data_type")
  method                            <- check_input(opt$method, "method")
  norm_method                       <- check_input(opt$norm_method, "norm_method")
  combine_for_guide_stats           <- check_input(opt$combine_for_guide_stats, "combine_for_guide_stats")
  combine_for_gene_stats            <- check_input(opt$combine_for_gene_stats, "combine_for_gene_stats")
  
  # Modules
  use_modules                      <- check_input(opt$use_modules, "use_modules")
  seqtk_module                     <- check_input(opt$seqtk_module, "seqtk_module")
  star_module                      <- check_input(opt$star_module, "star_module")
  samtools_module                  <- check_input(opt$samtools_module, "samtools_module")
  umi_tools_module                 <- check_input(opt$umi_tools_module, "umi_tools_module")
  
  # QC filtering
  qc_filtering_run                 <- check_input(opt$qc_filtering_run, "qc_filtering_run")
  qc_min_qual                      <- check_input(opt$qc_min_qual, "qc_min_qual")
  qc_qual_offset                   <- check_input(opt$qc_qual_offset, "qc_qual_offset")
  qc_min_length                    <- check_input(opt$qc_min_length, "qc_min_length")
  
  recover_input                     <- check_input(opt$recover_input, "recover_input")
  subsample_controls                <- check_input(opt$subsample_controls, "subsample_controls")
  use_custom_bins                   <- check_input(opt$use_custom_bins, "use_custom_bins")
  same_controls_in_all_sublibraries <- check_input(opt$same_controls_in_all_sublibraries, "same_controls_in_all_sublibraries")
  
  extra_suffix                      <- check_input(opt$extra_suffix, "extra_suffix")
  upper_lower_percentage            <- check_input(opt$upper_lower_percentage, "upper_lower_percentage")
  
  drop_0s                           <- check_input(opt$drop_0s, "drop_0s")
  strict_mode                       <- check_input(opt$strict_mode, "strict_mode")
  min_guides_per_gene               <- check_input(opt$min_guides_per_gene, "min_guides_per_gene")
  auto_combine_replicates           <- check_input(opt$auto_combine_replicates, "auto_combine_replicates")
  
  slurm_account                     <- check_input(opt$slurm_account, "slurm_account")
  slurm_qos                         <- check_input(opt$slurm_qos, "slurm_qos")
  slurm_cpus                        <- check_input(opt$slurm_cpus, "slurm_cpus")
  slurm_mem                         <- check_input(opt$slurm_mem, "slurm_mem")
  slurm_wall_time                   <- check_input(opt$slurm_wall_time, "slurm_wall_time")
  slurm_partition                   <- check_input(opt$slurm_partition, "slurm_partition")
  slurm_array                       <- check_input(opt$slurm_array, "slurm_array")
  slurm_email                       <- check_input(opt$slurm_email, "slurm_email")
  
  check_config_dependencies(opt)
  #=============================================================================
  # Construct special suffixes
  #=============================================================================
  recover_input_suffix <- if (recover_input) "RI" else ""
  subsample_controls_suffix <- if (subsample_controls) "ss_cntrl" else ""
  custom_bins_suffix <- if (use_custom_bins) "custom_bins" else ""
  drop_0s_suffix <- if (drop_0s) "D0" else ""
  strict_mode_suffix <- if (strict_mode) "strict" else ""
  auto_combine_replicates_suffix <- if (auto_combine_replicates) "acr" else ""
  min_guides_per_gene_suffix <- if (min_guides_per_gene > 0) paste0("min_guides_", min_guides_per_gene) else ""
  combine_for_guide_stats_suffix <- if (combine_for_guide_stats == "") "" else paste0("comb_", combine_for_guide_stats)
  combine_for_gene_stats_suffix <- if (combine_for_gene_stats == "") "" else paste0("comb_", combine_for_gene_stats)
  
  #=============================================================================
  # Construct the skip list
  #=============================================================================
  skip_list_and_suffix <- create_skip_list_and_suffix(
    skip_list,
    skip_list_sublib,
    skip_list_sample
  )
  skip_list <- skip_list_and_suffix[[1]]
  skip_suffix <- skip_list_and_suffix[[2]]
  
  #=============================================================================
  # Construct the file suffix
  #=============================================================================
  if (use_old_suffix_construction) {
    if (read_counting == "align_UMI_tools") {
      read_counting_old_name <- "lukas"
    } else if (read_counting == "bcwithqc") {
      read_counting_old_name <- "john"
    } else {
      stop(
        read_counting,
        " is not a valid read_counting value. Please enter either 'align_UMI_tools' or 'bcwithqc'.",
        call. = FALSE
      )
    }
    
    fs_parts <- c(
      read_counting_old_name,
      data_type,
      method,
      norm_method,
      recover_input_suffix,
      drop_0s_suffix,
      strict_mode_suffix,
      min_guides_per_gene_suffix,
      combine_for_guide_stats_suffix,
      auto_combine_replicates_suffix,
      skip_suffix,
      extra_suffix
    )
  } else {
    fs_parts <- c(
      read_counting,
      data_type,
      method,
      norm_method,
      recover_input_suffix,
      subsample_controls_suffix,
      drop_0s_suffix,
      strict_mode_suffix,
      custom_bins_suffix,
      min_guides_per_gene_suffix,
      combine_for_guide_stats_suffix,
      combine_for_gene_stats_suffix,
      auto_combine_replicates_suffix,
      skip_suffix,
      extra_suffix
    )
  }
  
  fs_parts <- fs_parts[!is.na(fs_parts) & fs_parts != ""]
  
  file_suffix <- paste0("_", paste(fs_parts, collapse = "_"), ".rds")
  file_info_suffix <- paste(fs_parts, collapse = "_")
  
  
  
  #=============================================================================
  # Construct File Paths
  #=============================================================================
  data_dir <- get_file_path(project_root_dir, "data")
  
  # Keep subfolder names relative. Avoid leading slashes here.
  genome_output_folder <- make_clean_dir(output_folder, "ref")
  dedup_output_folder <- make_clean_dir(output_folder, "dedup")
  bcwithqc_output_folder <- make_clean_dir(output_folder, "bcwithqc_output")
  mapped_output_folder <- make_clean_dir(output_folder, "mapped")
  qc_filtered_folder <- make_clean_dir(output_folder, "QC_filtered")
  rds_output_folder <- make_clean_dir(output_folder, "rds")
  results_output_folder <- make_clean_dir(output_folder, "results")
  
  symlinks_folder <- get_file_path(output_folder, "symlinks")
  fastq_symlinks_folder <- make_clean_dir(output_folder, "fastq_symlinks")
  bcwithqc_symlinks_folder <- make_clean_dir(output_folder, "bcwithqc_symlinks")
  
  log_folder <- make_clean_dir(output_folder, "logs")
  log_file <- file.path(
    log_folder,
    paste0("log_", file_info_suffix, "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".log")
  )
  log_appender(appender_tee(log_file))
  logger::log_layout(
    logger::layout_glue_generator(
      format = "{format(time, '%Y-%m-%d %H:%M:%S')} [{level}] {msg}"
    )
  )
  log_info(paste0("Logger initialized. Find log file at: ", log_file))
  

  merged_sgRNA_df <- read_library_file(
    library_path = library_path
  )
  
  # Return everything the Rmd uses.
  return_list <- list(
    # =========================
    # Options
    # =========================
    
    opt                               = opt,
    first_time                        = first_time,
    machine                           = machine,
    
    output_folder                     = output_folder,
    input_folder                      = input_folder,
    library_path                      = library_path,
    star_index_folder                 = star_index_folder,
    fastq_name_table_xlsx             = fastq_name_table_xlsx,
    strict_file_match                 = strict_file_match,
    
    bcwithqc_dir                      = bcwithqc_dir,
    bcwithqc_config_path              = bcwithqc_config_path,
    UMI_regex                         = UMI_regex,
    
    skip_list                         = skip_list,
    skip_list_sublib                  = skip_list_sublib,
    skip_list_sample                  = skip_list_sample,
    include_controls_list             = include_controls_list,
    use_only_these_controls_list      = use_only_these_controls_list,
    
    read_counting                     = read_counting,
    data_type                         = data_type,
    method                            = method,
    norm_method                       = norm_method,
    combine_for_guide_stats           = combine_for_guide_stats,
    combine_for_gene_stats            = combine_for_gene_stats,
    
    recover_input                     = recover_input,
    subsample_controls                = subsample_controls,
    use_custom_bins                   = use_custom_bins,
    same_controls_in_all_sublibraries = same_controls_in_all_sublibraries,
    
    extra_suffix                      = extra_suffix,
    upper_lower_percentage            = upper_lower_percentage,
    
    drop_0s                           = drop_0s,
    strict_mode                       = strict_mode,
    min_guides_per_gene               = min_guides_per_gene,
    auto_combine_replicates           = auto_combine_replicates,
    
    slurm_account                     = slurm_account,
    slurm_qos                         = slurm_qos,
    slurm_cpus                        = slurm_cpus,
    slurm_mem                         = slurm_mem,
    slurm_wall_time                   = slurm_wall_time,
    slurm_partition                   = slurm_partition,
    slurm_array                       = slurm_array,
    slurm_email                       = slurm_email,
    
    # =========================
    # Derived values and paths
    # =========================
    skip_suffix       = skip_suffix,
    file_suffix       = file_suffix,
    file_info_suffix  = file_info_suffix,
    
    data_dir                  = data_dir,
    genome_output_folder      = genome_output_folder,
    qc_filtered_folder        = qc_filtered_folder,
    dedup_output_folder       = dedup_output_folder,
    bcwithqc_output_folder    = bcwithqc_output_folder,
    mapped_output_folder      = mapped_output_folder,
    rds_output_folder         = rds_output_folder,
    results_output_folder     = results_output_folder,
    fastq_symlinks_folder     = fastq_symlinks_folder,
    bcwithqc_symlinks_folder  = bcwithqc_symlinks_folder,
    log_folder                = log_folder, 
    
    log_file = log_file,
    
    merged_sgRNA_df = merged_sgRNA_df
  )
  
  list2env(return_list, envir = envir)
  
  # Preserve old downstream behavior: if no explicit control subset is supplied,
  # remove the global object instead of leaving it as character(0).
  if (length(use_only_these_controls_list) == 0 &&
      exists("use_only_these_controls_list", envir = envir, inherits = FALSE)) {
    rm(use_only_these_controls_list, envir = envir)
  }
  #=============================================================================
  # Print the options for verification
  #=============================================================================
  logger::log_info("==========================================================")
  logger::log_info("Returning project setup options for verification:")
  logger::log_info("==========================================================")
  
  logger::log_info("first_time:                        {first_time}")
  logger::log_info("machine:                           {machine}")
  logger::log_info("----------------------------------------------------------")

  logger::log_info("output_folder:                     {output_folder}")
  logger::log_info("input_folder:                      {input_folder}")
  logger::log_info("library_path:                      {library_path}")
  logger::log_info("star_index_folder:                 {star_index_folder}")
  logger::log_info("fastq_name_table_xlsx:             {fastq_name_table_xlsx}")
  logger::log_info("strict_file_match:                 {strict_file_match}")
  logger::log_info("----------------------------------------------------------")
  
  logger::log_info("bcwithqc_dir:                      {bcwithqc_dir}")
  logger::log_info("bcwithqc_config_path:              {bcwithqc_config_path}")
  logger::log_info("UMI_regex:                         {UMI_regex}")
  logger::log_info("----------------------------------------------------------")
  
  logger::log_info("skip_list:                         {paste(skip_list, collapse = ', ')}")
  logger::log_info("skip_list_sublib:                  {paste(skip_list_sublib, collapse = ', ')}")
  logger::log_info("skip_list_sample:                  {paste(skip_list_sample, collapse = ', ')}")
  logger::log_info("include_controls_list:             {paste(include_controls_list, collapse = ', ')}")
  logger::log_info(
    "use_only_these_controls_list:      {x}",
    x = if (exists("use_only_these_controls_list")) {
      paste(use_only_these_controls_list, collapse = ", ")
    } else {
      "<not set>"
    }
  )
  logger::log_info("----------------------------------------------------------")
  
  logger::log_info("read_counting:                     {read_counting}")
  logger::log_info("data_type:                         {data_type}")
  logger::log_info("method:                            {method}")
  logger::log_info("norm_method:                       {norm_method}")
  logger::log_info("combine_for_guide_stats:           {combine_for_guide_stats}")
  logger::log_info("combine_for_gene_stats:            {combine_for_gene_stats}")
  logger::log_info("----------------------------------------------------------")
  
  logger::log_info("recover_input:                     {recover_input}")
  logger::log_info("subsample_controls:                {subsample_controls}")
  logger::log_info("use_custom_bins:                   {use_custom_bins}")
  logger::log_info("same_controls_in_all_sublibraries: {same_controls_in_all_sublibraries}")
  logger::log_info("----------------------------------------------------------")
  
  logger::log_info("extra_suffix:                      {extra_suffix}")
  logger::log_info("upper_lower_percentage:            {upper_lower_percentage}")
  logger::log_info("----------------------------------------------------------")
  
  logger::log_info("drop_0s:                           {drop_0s}")
  logger::log_info("strict_mode:                       {strict_mode}")
  logger::log_info("min_guides_per_gene:               {min_guides_per_gene}")
  logger::log_info("auto_combine_replicates:           {auto_combine_replicates}")
  logger::log_info("----------------------------------------------------------")
  
  logger::log_info("slurm_account:                     {slurm_account}")
  logger::log_info("slurm_qos:                         {slurm_qos}")
  logger::log_info("slurm_cpus:                        {slurm_cpus}")
  logger::log_info("slurm_mem:                         {slurm_mem}")
  logger::log_info("slurm_wall_time:                   {slurm_wall_time}")
  logger::log_info("slurm_partition:                   {slurm_partition}")
  logger::log_info("slurm_array:                       {slurm_array}")
  logger::log_info("slurm_email:                       {slurm_email}")
  logger::log_info("----------------------------------------------------------")
  
  logger::log_info("File Suffix: {file_suffix}")
  
  invisible(return_list)
}
