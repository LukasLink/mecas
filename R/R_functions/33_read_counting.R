run_read_counting <- function(manifest,
                              manifest_output_path,
                              run_read_counting_stage,
                              opt,
                              slurm_settings,
                              project_root_dir,
                              log_folder) {
  
  if (isTRUE(run_read_counting_stage)) {
    
    #-----------------------------------------------------------------------------
    # Run QC filtering
    #-----------------------------------------------------------------------------
    
    if (opt$qc_filtering_run) {
      logger::log_info("Beginning QC filtering of reads...")
      
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
      
      write.table(
        manifest,
        file = manifest_output_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
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
      logger::log_info("Skipping QC filtering of reads.")
      
      manifest <- manifest %>%
        mutate(qc_filtered_paths = symlink_file)
      
      write.table(
        manifest,
        file = manifest_output_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
    
    #-----------------------------------------------------------------------------
    # Optional bcwithqc symlink creation
    #-----------------------------------------------------------------------------
    
    if (identical(opt$read_counting, "bcwithqc")) {
      manifest <- prepare_bcwithqc_inputs(
        manifest = manifest,
        output_symlink_dir = bcwithqc_symlinks_folder,
        overwrite_symlinks = TRUE,
        manifest_output_path = file.path(
          bcwithqc_symlinks_folder,
          "bcwithqc_symlink_manifest.tsv"
        )
      )
      
      logger::log_info("Using bcwithqc FASTQ folder: {bcwithqc_symlinks_folder}")
    }
    
    #-----------------------------------------------------------------------------
    # Read counting shell jobs
    #-----------------------------------------------------------------------------
    
    write.table(
      manifest,
      file = manifest_output_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    limit_bam_sort_ram <- floor(slurm_mem_to_bytes(opt$slurm_mem) * 0.90)
    
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
    
    if (identical(opt$read_counting, "bcwithqc")) {
      logger::log_info("Launching read counting for bcwithqc...")
      
      run_shell_step(
        step_name = "bcwithqc",
        script_path = file.path(project_root_dir, "shell", "run_bcwithqc.sh"),
        args = c(
          "--bcwithqc-bin", opt$bcwithqc_bin,
          "--manifest", manifest_output_path,
          "--output-dir", opt$output_folder,
          "--bcwithqc-config", opt$bcwithqc_config_path,
          "--star-index-folder", opt$star_index_folder,
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
    
    logger::log_info("Finished read counting.")
    
  } else {
    
    logger::log_info(
      "Skipping QC/read-counting stage because start_with is: {opt$start_with}"
    )
  }
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
  
  return(count_df_long)
}