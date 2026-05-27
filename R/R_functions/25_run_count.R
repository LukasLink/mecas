# run_count.R

run_count <- function(config_path, project_root_dir, cli_args){
  #-----------------------------------------------------------------------------
  # Run setup
  #-----------------------------------------------------------------------------
  logger::log_info("Begining Project setup...")
  
  config <- yaml::read_yaml(config_path)
  # Priority explicit overrides > Rmd params > config.yaml
  project_setup(
    project_root_dir = project_root_dir,
    config_path = config_path,
    setup_mode = "count",
    use_old_suffix_construction = FALSE
  )
  logger::log_info("Finished Project setup.")
  #-----------------------------------------------------------------------------
  # Optional: handle use_only_these_controls
  #-----------------------------------------------------------------------------
  
  # If config$controls$use_only_these_controls contains entries, use them.
  # Otherwise, remove the object to preserve old downstream behavior.
  
  if (exists("use_only_these_controls_list") &&
      length(use_only_these_controls_list) == 0) {
    rm(use_only_these_controls_list)
  }
  #-----------------------------------------------------------------------------
  # Create Symlinks
  #-----------------------------------------------------------------------------
  logger::log_info("Creating symlinks of input files...")
  manifest_output_path <- get_file_path(
    rds_output_folder,
    paste0("standardized_fastq_manifest_", file_info_suffix, ".tsv")
  )
  manifest <- prepare_fastq_inputs(
    fastq_dir = input_folder,
    fastq_name_table_file_path = fastq_name_table_xlsx,
    output_symlink_dir = fastq_symlinks_folder,
    manifest_output_path = manifest_output_path,
    strict_file_match = strict_file_match
  )
  log_info(paste("Using standardized FASTQ folder: ", fastq_symlinks_folder))
  log_info(paste("FASTQ manifest written to: ", manifest_output_path))
  logger::log_info("Finished Creating symlinks of input files.")
  #-----------------------------------------------------------------------------
  # Make slurm settings 
  #-----------------------------------------------------------------------------
  slurm_settings <- list(
    account    = slurm_account,
    qos        = slurm_qos,
    cpus       = slurm_cpus,
    mem        = slurm_mem,
    wall_time  = slurm_wall_time,
    partition  = slurm_partition,
    array      = slurm_array,
    email      = slurm_email
  )
  # slurm_settings <- list(
  #   account    = slurm_account,
  #   qos        = slurm_qos,
  #   cpus       = 1,
  #   mem        = "1g",
  #   wall_time  = "01:01:00",
  #   partition  = slurm_partition,
  #   array      = 1,
  #   email      = slurm_email
  # )
  #-----------------------------------------------------------------------------
  # Run QC filtering
  #-----------------------------------------------------------------------------
  if (opt$qc_filtering_run) {
    logger::log_info("Begining QC filtering of reads...")
    
    if (is.na(opt$qc_min_length) || is.null(opt$qc_min_length)) {
      opt$qc_min_length <- infer_qc_min_length(
        manifest = manifest,
        fastq_col = "symlink_file",
        n_lines = 10000
      )
    }
    
    manifest <- manifest %>%
      mutate(
        qc_filtered_paths = file.path(
          qc_filtered_folder,
          ensure_gz_suffix(symlink_file_basename)
        )
      )
    
    run_shell_step(
      step_name = "QC_filtering",
      script_path = file.path(project_root_dir, "shell", "QC_filtering.sh"),
      args = c(
        "--manifest", manifest_output_path,
        "--output-dir", qc_filtered_folder,
        "--min-qual", as.character(opt$qc_min_qual),
        "--qual-offset", as.character(opt$qc_qual_offset),
        "--min-length", as.character(opt$qc_min_length),
        "--use-modules", if (opt$use_modules) "true" else "false",
        "--seqtk-module", opt$seqtk_module
      ),
      slurm_settings = slurm_settings,
      machine = opt$machine,
      log_dir = log_folder
    )
    
    logger::log_info("Finished QC filtering of reads.")
  } else {
    logger::log_info("Skipping QC filtering of reads, as QC_filtering.run in config yaml file is set to false.")
    manifest <- manifest %>%
      mutate(
        qc_filtered_paths = symlink_file
      )
  }
  #-----------------------------------------------------------------------------
  # Optional bcwithqc symlink creation goes here
  #-----------------------------------------------------------------------------
  
  if (read_counting == "bcwithqc") {
    manifest <- prepare_bcwithqc_inputs(
      manifest = manifest,
      output_symlink_dir = bcwithqc_symlinks_folder,
      overwrite_symlinks = TRUE,
      manifest_output_path = file.path(bcwithqc_symlinks_folder,
                                       "bcwithqc_symlink_manifest.tsv")
    )
    
    log_info(paste("Using bcwithqc FASTQ folder: ", bcwithqc_symlinks_folder))
  }
  #-----------------------------------------------------------------------------
  # Read counting
  #-----------------------------------------------------------------------------
  # Make sure the manifest on disk contains qc_filtered_paths
  write.table(
    manifest,
    file = manifest_output_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  # Figure out how much memory STAR is allowed to use
  limit_bam_sort_ram <- floor(slurm_mem_to_bytes(opt$slurm_mem) * 0.90)
  
  # For align_UMI_tools
  if (identical(opt$read_counting, "align_UMI_tools")) {
    logger::log_info("Launching read counting for align_UMI_tools...")
    
    run_shell_step(
      step_name = "align_UMI_tools",
      script_path = file.path(project_root_dir, "shell", "align_UMI_tools.sh"),
      args = c(
        "--manifest", manifest_output_path,
        "--output-dir", opt$output_folder,
        "--star-index-folder", opt$star_index_folder,
        "--data-type", opt$data_type,
        "--umi-regex", opt$UMI_regex,
        "--threads", as.character(opt$slurm_cpus),
        "--limit-bam-sort-ram", as.character(limit_bam_sort_ram),
        "--use-modules", if (isTRUE(opt$use_modules)) "true" else "false",
        "--star-module", opt$star_module,
        "--samtools-module", opt$samtools_module,
        "--umi-tools-module", opt$umi_tools_module
      ),
      slurm_settings = slurm_settings,
      machine = opt$machine,
      log_dir = log_folder
    )
  }
  # For align_UMI_tools
  if (identical(opt$read_counting, "bcwithqc")) {
    logger::log_info("Launching read counting for bcwithqc...")
    
    run_shell_step(
      step_name = "bcwithqc",
      script_path = file.path(project_root_dir, "shell", "run_bcwithqc.sh"),
      args = c(
        "--manifest", manifest_output_path,
        "--output-dir", opt$output_folder,
        "--bcwithqc-config", opt$bcwithqc_config_path,
        "--star-index-folder", opt$star_index_folder,
        "--bcwithqc-bin", file.path(opt$bcwithqc_dir, "bcwithqc"),
        "--threads", as.character(opt$slurm_cpus),
        "--keep-intermediary", "false",
        "--existing-results-mode", "override",
        "--verbosity", "-vv",
        "--use-modules", if (isTRUE(opt$use_modules)) "true" else "false",
        "--star-module", opt$star_module
      ),
      slurm_settings = slurm_settings,
      machine = opt$machine,
      log_dir = log_folder,
      extra_slurm_sbatch_lines = c("#SBATCH --constraint=avx512")
    )
  }
  logger::log_info("Finished Read Counting.")
  #-----------------------------------------------------------------------------
  # Get count_df_long 
  #-----------------------------------------------------------------------------
  logger::log_info("Reading in read/UMI counts...")
  
  if (identical(opt$read_counting, "bcwithqc")){
    count_df_long <- process_bcwithqc_data(
      data_type = data_type,
      skip_list = skip_list)
  } else if (identical(opt$read_counting, "align_UMI_tools")) {
    if (identical(opt$data_type, "umis")){
      count_df_long <- process_folder_files(dedup_output_folder,
                                            skip_list = skip_list) #Add threshold df if thresholds should be applied 
    }
    if (identical(opt$data_type, "reads")){
      count_df_long <- process_folder_files(mapped_output_folder,
                                            skip_list = skip_list) #Add threshold df if thresholds should be applied 
    }
  } else {
    stop(opt$read_counting, " is not a valid read_counting value. Exiting. read_counting must be 'bcwithqc' or 'align_UMI_tools'.")
  }
  logger::log_info("Finished Reading in read/UMI counts.")
  logger::log_info("sgRNAs aligned to the wrong sublibrary were excluded from the analysis.")
  
  #-----------------------------------------------------------------------------
  # Optional count_df_long include_controls_list if use_only_these_controls_list is given
  #-----------------------------------------------------------------------------
  if (exists("use_only_these_controls_list")) {
    if (length(use_only_these_controls_list) > 0){
      # list of sgRNAs that were not targeting (before) and not in allowed controls
      excluded_controls <- count_df_long %>%
        filter(group_category != "targeting", !sgRNA %in% use_only_these_controls_list) %>%
        distinct(sgRNA) %>%
        pull(sgRNA)
      
      include_controls_list <- c(include_controls_list, excluded_controls)
    }
  }
  if (length(include_controls_list) > 0){
    for (control_gene in include_controls_list){
      count_df_long$group_category[grepl(control_gene, count_df_long$sgRNA)] <- "targeting"
    }
  }
  
  #-----------------------------------------------------------------------------
  # Optional Combine either samples or sublibraries of the same condition. 
  #-----------------------------------------------------------------------------
  
  if (combine_for_guide_stats != ""){
    if (combine_for_guide_stats == "sample"){
      count_df_long <- count_df_long %>%
        group_by(sgRNA, sublib, condition) %>%
        summarise(
          # set sample name for the combined rows
          sample = "sample_1",
          count = sum(count, na.rm = TRUE),
          exp = dplyr::first(exp),
          group_category = dplyr::first(group_category),
          
          # keep one sublib/sgRNA (also fine even though they're grouping keys)
          .groups = "drop"
        )
    }
    if (combine_for_guide_stats == "sublib"){
      count_df_long <- count_df_long %>%
        group_by(sgRNA, sample, condition) %>%
        summarise(
          sublib = "sublib_1",
          count = sum(count, na.rm = TRUE),
          exp = dplyr::first(exp),
          group_category = dplyr::first(group_category),
          
          # keep one sublib/sgRNA (also fine even though they're grouping keys)
          .groups = "drop"
        )    
    }
  }
  #-----------------------------------------------------------------------------
  # Optional Violin Plots
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Make coverage file
  #-----------------------------------------------------------------------------
  logger::log_info("Creating Coverage file for Read/UMI counts...")
  df_cov <- parse_coverage_file(log_file)
  
  mapping_results_df <- add_star_log_stats(df_cov)
  
  rm(df_cov)
  
  overall_targeting <- count_df_long %>%
    summarise(
      total_counts     = sum(count, na.rm = TRUE),
      targeting_counts = sum(count[group_category == "targeting"], na.rm = TRUE),
      targeting_perc   = 100 * targeting_counts / total_counts
    ) %>%
    mutate(targeting_perc = sprintf("%.2f%%", targeting_perc))
  
  targeting_by_group <- count_df_long %>%
    group_by(condition, sublib, sample) %>%
    summarise(
      total_counts     = sum(count, na.rm = TRUE),
      targeting_counts = sum(count[group_category == "targeting"], na.rm = TRUE),
      targeting_perc   = 100 * targeting_counts / total_counts,
      .groups = "drop"
    ) %>%
    mutate(targeting_perc = sprintf("%.2f%%", targeting_perc))
  
  overall_targeting_merged <- merged_sgRNA_df %>%
    summarise(
      total_counts     = sum(count, na.rm = TRUE),
      targeting_counts = sum(count[!is.na(entrez)], na.rm = TRUE),
      targeting_perc   = 100 * targeting_counts / total_counts
    ) %>%
    mutate(targeting_perc = sprintf("%.2f%%", targeting_perc))
  
  
  # ---- Write everything into ONE excel file as separate sheets ----
  write_xlsx(
    list(
      mapping_results          = mapping_results_df,
      overall_targeting        = overall_targeting,
      targeting_by_group       = targeting_by_group,
      overall_targeting_reference = overall_targeting_merged
    ),
    get_file_path(results_output_folder,
                  paste0(file_info_suffix, "_mapping_results.xlsx"))
  )
  
  logger::log_info("Finished Creating Coverage file for Read/UMI counts.")
  logger::log_info(paste("The coverage file is:",
                         get_file_path(results_output_folder,
                                       paste0(file_info_suffix, "_mapping_results.xlsx"))))
  #-----------------------------------------------------------------------------
  # prepare data for MAUDE
  #-----------------------------------------------------------------------------
  logger::log_info("Preparing Data for MAUDE...")
  if (subsample_controls == TRUE){
    # count_df_long_old <- subsample_controls_func_old(count_df_long, merged_sgRNA_df)
    count_df_long <- subsample_controls_func(count_df_long, merged_sgRNA_df)
  }
  
  if (!(norm_method %in% c("","control_median"))){
    stop("Error: 'norm_method' must be one of '', 'control_median'. The script will now stop.")
  }
  if (norm_method == "control_median"){
    count_df_long <- normalize_count_df_long(count_df_long,
                                             norm_method = norm_method)
  }
  
  
  
  maude_counts_df <- count_df_long_to_wide(count_df_long = count_df_long,
                                           print = FALSE,
                                           drop_0s = drop_0s,
                                           recover_input = recover_input)
  if (!(method %in% c("","rep","sum",'rep_sample','rep_sublib'))){
    stop("Error: 'method' must be one of '', 'rep', 'rep_sample', 'rep_sublib', or 'sum'. The script will now stop.")
  }
  
  if (method == ""){
    maude_counts_df <- maude_counts_df %>% 
      mutate(exp = "rep1")
  }
  if (method == "rep"){
    
  }
  if (method == "rep_sample"){
    maude_counts_df <- maude_counts_df %>%
      mutate(exp = sample)
  }
  if (method == "rep_sublib"){
    maude_counts_df <- maude_counts_df %>%
      mutate(exp = sublib)
  }
  if (method == "sum"){
    stop("Method 'sum' is deprecated, use the option combine_for_guide_stats instead")
    maude_counts_df <- count_df_long_to_wide(count_df_long = count_df_long,
                                             print = FALSE,
                                             drop_0s = drop_0s,
                                             recover_input = TRUE,
                                             for_sum = TRUE)
    # Group by sgRNA and summarize the required columns
    maude_counts_df <- maude_counts_df %>%
      group_by(sgRNA, sublib) %>%
      summarize(
        input = pmax(sum(input, na.rm = TRUE), 0),  # Sum and ensure minimum is 0
        upper = pmax(sum(upper, na.rm = TRUE), 0),  # Sum and ensure minimum is 0
        lower = pmax(sum(lower, na.rm = TRUE), 0),  # Sum and ensure minimum is 0
        isNontargeting = dplyr::first(isNontargeting),  # Take the first value of isNontargeting (same for all in the group)
        .groups = 'drop'  # Drop the group structure after summarizing
      ) %>%
      mutate(
        exp = "rep1",
        input = input + 1,
        upper = upper + 1,
        lower = lower + 1
      )
  }
  
  if (strict_mode){
    if (pseudocount_added){
      umi_threshold <- 2
    } else {
      umi_threshold <- 1
    } 
    cat("Strict Mode enabled\n")
    cat("Rows before strict mode: \t", nrow(maude_counts_df),"\n")
    maude_counts_df <- maude_counts_df %>%
      filter(if_any(c(input, upper, lower), ~ . <= umi_threshold))
    cat("Rows after strict mode: \t", nrow(maude_counts_df),"\n")
  }
  if (length(include_controls_list) > 0) {
    for (control_gene in include_controls_list) {
      maude_counts_df$isNontargeting[grepl(control_gene, maude_counts_df$sgRNA)] <- FALSE
    }
  }
  if (exists("use_only_these_controls_list")) {
    if (length(use_only_these_controls_list) > 0){
      maude_counts_df$isNontargeting[ !(maude_counts_df$sgRNA %in% use_only_these_controls_list) ] <- FALSE
    }
  }
  #-----------------------------------------------------------------------------
  # Pre MAUDE plots
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # Run MAUDE
  #-----------------------------------------------------------------------------
  logger::log_info("Starting initial MAUDE run ...")
  
  unique_exp <- unique(maude_counts_df$exp)
  
  # Define bin stats with 10% for lower/upper each
  lower_bin_end = upper_lower_percentage
  upper_bin_start = 1 - upper_lower_percentage
  
  maude_bins <- tibble(Bin = rep(c('upper', 'lower'), length(unique_exp)),  # Repeat 'upper' and 'lower' for each exp
                       exp = rep(unique_exp, each = 2),  # Repeat each exp value twice for 'upper' and 'lower'
                       binStartQ = ifelse(rep(c('upper', 'lower'), length(unique_exp)) == 'lower', 0.001, upper_bin_start),
                       binEndQ = ifelse(rep(c('upper', 'lower'), length(unique_exp)) == 'lower', lower_bin_end, 0.999),
                       fraction = binEndQ - binStartQ,
                       binStartZ = qnorm(binStartQ),
                       binEndZ = qnorm(binEndQ)) %>%
    select(Bin, binStartQ, binEndQ, fraction, binStartZ, binEndZ, exp) %>%
    as.data.frame()
  
  
  
  if (first_time == TRUE){
    ## The input dataframe needs to have the lower, upper and input columns.
    ## use maude to calculate guide level statistics.
    maude_guide_stats <- findGuideHitsAllScreens(
      experiments = unique(maude_counts_df['exp']),
      countDataFrame = maude_counts_df,
      binStats = maude_bins,
      sortBins = c('lower', 'upper'),
      unsortedBin = 'input',
      negativeControl = 'isNontargeting'
    )
    
    saveRDS(maude_guide_stats,file.path(rds_output_folder,
                                        paste0("MAUDE_guide_stats", file_suffix)))
  } else {
    maude_guide_stats <- readRDS(file.path(rds_output_folder,
                                           paste0("MAUDE_guide_stats", file_suffix)))
  }
  
  
  if (first_time == TRUE){
    
    maude_guide_stats <- maude_guide_stats %>%
      left_join(
        merged_sgRNA_df %>%
          select(sgrna_id, entrez) %>%
          distinct(sgrna_id, .keep_all = TRUE),
        by = c("sgRNA" = "sgrna_id")
      ) %>%
      mutate(entrez = coalesce(as.character(entrez), sgRNA))
    
    # any entries from include_controls_list are manually turned into genes
    if (length(include_controls_list) > 0) {
      for (control_gene in include_controls_list) {
        # here we remove everything after the last _, so stuff like AAVS1_9 and
        # AAVS1_13 are both treated as AAVS1
        control_gene <- sub("_[^_]*$", "", control_gene)
        maude_guide_stats$entrez[grepl(control_gene, maude_guide_stats$sgRNA)] <- control_gene
      }
    }
    
    if (combine_for_gene_stats != "none"){
      if (!(combine_for_gene_stats %in% c("all","sublib","sample"))){
        stop("combine_for_gene_stats must be one of: 'all','none','sublib','sample'")
      }
      if (combine_for_gene_stats == "all"){
        maude_guide_stats$exp <- "rep1"
      }
      if (combine_for_gene_stats == "sublib"){
        maude_guide_stats <- maude_guide_stats %>% 
          mutate(exp = sample)
      }
      if (combine_for_gene_stats == "sample"){
        maude_guide_stats <- maude_guide_stats %>% 
          mutate(exp = sublib)
      }    
    }
    
    ## calculate gene-level summarized scores
    maude_gene_stats <- getElementwiseStats(
      experiments = unique(maude_guide_stats['exp']),
      normNBSummaries = maude_guide_stats,
      negativeControl = 'isNontargeting',
      elementIDs = 'entrez'
    )
    
    # Filter out all genes with not enough guides pointing to them
    maude_gene_stats <- maude_gene_stats %>%
      filter(numGuides >= min_guides_per_gene)
    
    saveRDS(maude_gene_stats,file.path(rds_output_folder,
                                       paste0("MAUDE_gene_stats", file_suffix)))
  } else {
    maude_gene_stats <- readRDS(file.path(rds_output_folder,
                                          paste0("MAUDE_gene_stats", file_suffix)))
  }
  
  maude_guide_stats <- readRDS(file.path(rds_output_folder,
                                         paste0("MAUDE_guide_stats", file_suffix)))
  maude_gene_stats <- readRDS(file.path(rds_output_folder,
                                        paste0("MAUDE_gene_stats", file_suffix)))
  
  
  
  #-----------------------------------------------------------------------------
  # post MAUDE plots
  #-----------------------------------------------------------------------------
  
  
  #-----------------------------------------------------------------------------
  # Finding Consensus Hits: MAUDE
  #-----------------------------------------------------------------------------
  
  
  
  
  
  
  
  
  
  
  
  run_shell_step(
    step_name = "test_dummy",
    script_path = file.path(project_root_dir, "shell", "test_dummy.sh"),
    slurm_settings = slurm_settings,
    machine = "slurm", # LOREM -> change this to machine = machine in final setup
    log_dir = log_folder,
    extra_slurm_sbatch_lines = c("#SBATCH --constraint=avx512")
  )
}


