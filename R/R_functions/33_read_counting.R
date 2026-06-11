run_read_counting <- function(manifest,
                              manifest_output_path,
                              run_read_counting_stage,
                              cfg,
                              project_root_dir) {
  
  log_dir <- cfg$paths$log_folder %||% cfg$paths$logs_folder
  # This is hopefully now a Chesterton’s fence, if I correctly renamed everything. 
  
  if (isTRUE(run_read_counting_stage)) {
    
    #-----------------------------------------------------------------------------
    # Run QC filtering
    #-----------------------------------------------------------------------------
    
    if (isTRUE(cfg$qc_filtering$run)) {
      logger::log_info("Beginning QC filtering of reads...")
      
      if (is.na(cfg$qc_filtering$min_length) || is.null(cfg$qc_filtering$min_length)) {
        cfg$qc_filtering$min_length <- infer_qc_min_length(
          manifest = manifest,
          fastq_col = "symlink_file",
          n_lines = 10000
        )
      }
      
      manifest <- manifest %>%
        dplyr::mutate(
          qc_filtered_paths = file.path(
            cfg$paths$qc_filtered_folder,
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
        cfg = cfg,
        args = c(
          "--manifest", manifest_output_path,
          "--output-dir", cfg$paths$qc_filtered_folder,
          "--min-qual", as.character(cfg$qc_filtering$min_qual),
          "--qual-offset", as.character(cfg$qc_filtering$qual_offset),
          "--min-length", as.character(cfg$qc_filtering$min_length),
          "--use-modules", if (isTRUE(cfg$modules$use_modules)) "true" else "false",
          "--seqtk-module", cfg$modules$seqtk
        ),
        log_dir = log_dir
      )
      
      logger::log_info("Finished QC filtering of reads.")
      
    } else {
      logger::log_info("Skipping QC filtering of reads.")
      
      manifest <- manifest %>%
        dplyr::mutate(qc_filtered_paths = symlink_file)
      
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
    
    if (identical(cfg$counting$read_counting, "bcwithqc")) {
      manifest <- prepare_bcwithqc_inputs(
        manifest = manifest,
        output_symlink_dir = cfg$paths$bcwithqc_symlinks_folder,
        overwrite_symlinks = TRUE,
        manifest_output_path = file.path(
          cfg$paths$bcwithqc_symlinks_folder,
          "bcwithqc_symlink_manifest.tsv"
        )
      )
      
      logger::log_info(
        "Using bcwithqc FASTQ folder: {cfg$paths$bcwithqc_symlinks_folder}"
      )
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
    
    limit_bam_sort_ram <- floor(slurm_mem_to_bytes(cfg$slurm$mem) * 0.90)
    
    if (identical(cfg$counting$read_counting, "align_UMI_tools")) {
      logger::log_info("Launching read counting for align_UMI_tools...")
      
      run_shell_step(
        step_name = "align_UMI_tools",
        script_path = file.path(project_root_dir, "shell", "align_UMI_tools.sh"),
        cfg = cfg,
        args = c(
          "--manifest", manifest_output_path,
          "--output-dir", cfg$paths$output_folder,
          "--star-index-folder", cfg$paths$star_index_folder,
          "--data-type", cfg$counting$data_type,
          "--umi-regex", cfg$align_UMI_tools$UMI_regex,
          "--threads", as.character(cfg$slurm$cpus),
          "--limit-bam-sort-ram", as.character(limit_bam_sort_ram),
          "--use-modules", if (isTRUE(cfg$modules$use_modules)) "true" else "false",
          "--star-module", cfg$modules$star,
          "--samtools-module", cfg$modules$samtools,
          "--umi-tools-module", cfg$modules$umi_tools
        ),
        log_dir = log_dir
      )
    }
    
    if (identical(cfg$counting$read_counting, "bcwithqc")) {
      logger::log_info("Launching read counting for bcwithqc...")
      
      run_shell_step(
        step_name = "bcwithqc",
        script_path = file.path(project_root_dir, "shell", "run_bcwithqc.sh"),
        cfg = cfg,
        args = c(
          "--bcwithqc-bin", cfg$bcwithqc$bcwithqc_bin,
          "--manifest", manifest_output_path,
          "--output-dir", cfg$paths$output_folder,
          "--bcwithqc-config", cfg$bcwithqc$bcwithqc_config_path,
          "--star-index-folder", cfg$paths$star_index_folder,
          "--threads", as.character(cfg$slurm$cpus),
          "--keep-intermediary", "false",
          "--existing-results-mode", "override",
          "--verbosity", "-vv",
          "--use-modules", if (isTRUE(cfg$modules$use_modules)) "true" else "false",
          "--star-module", cfg$modules$star
        ),
        log_dir = log_dir,
        extra_slurm_sbatch_lines = c("#SBATCH --constraint=avx512")
      )
    }
    
    logger::log_info("Finished read counting.")
    
  } else {
    
    logger::log_info(
      "Skipping QC/read-counting stage because start_with is: {cfg$run$start_with}"
    )
  }
  
  #-----------------------------------------------------------------------------
  # Read counts back into R
  #-----------------------------------------------------------------------------
  
  logger::log_info("Reading in read/UMI counts...")
  
  if (identical(cfg$counting$read_counting, "bcwithqc")) {
    
    count_df_long <- process_bcwithqc_data(
      cfg = cfg,
      data_type = cfg$counting$data_type,
      skip_list = cfg$skip$files
    )
    
  } else if (identical(cfg$counting$read_counting, "align_UMI_tools")) {
    
    if (identical(cfg$counting$data_type, "umis")) {
      count_df_long <- process_folder_files(
        cfg = cfg,
        folder_path = cfg$paths$dedup_output_folder,
        skip_list = cfg$skip$files
      )
    }
    
    if (identical(cfg$counting$data_type, "reads")) {
      count_df_long <- process_folder_files(
        cfg = cfg,
        folder_path = cfg$paths$mapped_output_folder,
        skip_list = cfg$skip$files
      )
    }
    
  } else {
    stop(
      cfg$counting$read_counting,
      " is not a valid read_counting value. ",
      "read_counting must be 'bcwithqc' or 'align_UMI_tools'.",
      call. = FALSE
    )
  }
  
  logger::log_info("Finished reading in read/UMI counts.")
  logger::log_info("sgRNAs aligned to the wrong sublibrary were excluded from the analysis.")
  
  #-----------------------------------------------------------------------------
  # Optional: use_only_these_controls / include_controls
  #-----------------------------------------------------------------------------
  
  use_only_these_controls_list <- cfg$controls$use_only_these_controls %||% character()
  include_controls_list <- cfg$controls$include_controls %||% character()
  
  if (length(use_only_these_controls_list) > 0) {
    excluded_controls <- count_df_long %>%
      dplyr::filter(
        group_category != "targeting",
        !sgRNA %in% use_only_these_controls_list
      ) %>%
      dplyr::distinct(sgRNA) %>%
      dplyr::pull(sgRNA)
    
    include_controls_list <- c(include_controls_list, excluded_controls)
  }
  
  if (length(include_controls_list) > 0) {
    for (control_gene in include_controls_list) {
      count_df_long$group_category[
        grepl(control_gene, count_df_long$sgRNA)
      ] <- "targeting"
    }
  }
  
  #-----------------------------------------------------------------------------
  # Optional: combine samples or sublibraries of the same condition
  #-----------------------------------------------------------------------------
  
  combine_for_guide_stats <- cfg$counting$combine_for_guide_stats %||% ""
  
  if (!identical(combine_for_guide_stats, "")) {
    
    if (identical(combine_for_guide_stats, "sample")) {
      count_df_long <- count_df_long %>%
        dplyr::group_by(sgRNA, sublib, condition) %>%
        dplyr::summarise(
          sample = "sample_1",
          count = sum(count, na.rm = TRUE),
          exp = dplyr::first(exp),
          group_category = dplyr::first(group_category),
          .groups = "drop"
        )
    }
    
    if (identical(combine_for_guide_stats, "sublib")) {
      count_df_long <- count_df_long %>%
        dplyr::group_by(sgRNA, sample, condition) %>%
        dplyr::summarise(
          sublib = "sublib_1",
          count = sum(count, na.rm = TRUE),
          exp = dplyr::first(exp),
          group_category = dplyr::first(group_category),
          .groups = "drop"
        )
    }
  }
  
  #-----------------------------------------------------------------------------
  # Make coverage file
  #-----------------------------------------------------------------------------
  
  logger::log_info("Creating coverage file for read/UMI counts...")
  
  df_cov <- parse_coverage_file(cfg$paths$log_file)
  
  mapping_results_df <- add_star_log_stats(
    df = df_cov,
    cfg = cfg
  )
  
  rm(df_cov)
  
  overall_targeting <- count_df_long %>%
    dplyr::summarise(
      total_counts = sum(count, na.rm = TRUE),
      targeting_counts = sum(count[group_category == "targeting"], na.rm = TRUE),
      targeting_perc = 100 * targeting_counts / total_counts
    ) %>%
    dplyr::mutate(targeting_perc = sprintf("%.2f%%", targeting_perc))
  
  targeting_by_group <- count_df_long %>%
    dplyr::group_by(condition, sublib, sample) %>%
    dplyr::summarise(
      total_counts = sum(count, na.rm = TRUE),
      targeting_counts = sum(count[group_category == "targeting"], na.rm = TRUE),
      targeting_perc = 100 * targeting_counts / total_counts,
      .groups = "drop"
    ) %>%
    dplyr::mutate(targeting_perc = sprintf("%.2f%%", targeting_perc))
  
  overall_targeting_merged <- cfg$merged_sgRNA_df %>%
    dplyr::summarise(
      total_counts = sum(count, na.rm = TRUE),
      targeting_counts = sum(count[!is.na(entrez)], na.rm = TRUE),
      targeting_perc = 100 * targeting_counts / total_counts
    ) %>%
    dplyr::mutate(targeting_perc = sprintf("%.2f%%", targeting_perc))
  
  writexl::write_xlsx(
    list(
      mapping_results = mapping_results_df,
      overall_targeting = overall_targeting,
      targeting_by_group = targeting_by_group,
      overall_targeting_reference = overall_targeting_merged
    ),
    get_file_path(
      cfg$paths$results_output_folder,
      paste0(cfg$suffix$file_info_suffix, "_mapping_results.xlsx")
    )
  )
  
  logger::log_info("Finished creating coverage file for read/UMI counts.")
  logger::log_info(
    paste(
      "The coverage file is:",
      get_file_path(
        cfg$paths$results_output_folder,
        paste0(cfg$suffix$file_info_suffix, "_mapping_results.xlsx")
      )
    )
  )
  
  return(count_df_long)
}