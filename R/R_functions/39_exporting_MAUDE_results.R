handle_auto_combine_replicates_and_export <- function(file_info_suffix,
                                                      maude_results,
                                                      cfg,
                                                      file_suffix = cfg$suffix$file_suffix) {
  
  auto_combine_replicates <- cfg$replicates$auto_combine_replicates
  
  if (isTRUE(auto_combine_replicates)) {
    logger::log_info("Combining replicates as auto_combine_replicates is set to true.")
    
    if (length(unique(maude_results$maude_gene_stats$exp)) > 1) {
      
      genes_with_rep <- add_info_wrapper(
        suffix = file_info_suffix,
        cfg = cfg
      )
      
      export_df <- genes_with_rep %>%
        dplyr::group_by(symbol, entrez) %>%
        dplyr::summarise(
          numGuides = sum(numGuides, na.rm = TRUE),
          seq = dplyr::first(seq),
          sgRNA = dplyr::first(sgRNA),
          meanZ = mean(meanZ, na.rm = TRUE),
          
          stoufferZ = sum(significanceZ, na.rm = TRUE) /
            sqrt(sum(!is.na(significanceZ))),
          significanceZ = stoufferZ,
          
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          p.value = 2 * stats::pnorm(abs(stoufferZ), lower.tail = FALSE),
          FDR = stats::p.adjust(p.value, method = "BH")
        ) %>%
        dplyr::select(
          symbol,
          entrez,
          numGuides,
          stoufferZ,
          meanZ,
          significanceZ,
          p.value,
          FDR,
          seq,
          sgRNA
        )
      
    } else {
      message("No replicates present, skipping replicate combination")
      
      export_df <- add_info_wrapper(
        suffix = file_info_suffix,
        cfg = cfg
      )
    }
    
  } else {
    export_df <- add_info_wrapper(
      suffix = file_info_suffix,
      cfg = cfg
    )
  }
  
  include_controls_list <- cfg$controls$include_controls %||% character()
  
  if (length(include_controls_list) > 0) {
    for (control_gene in include_controls_list) {
      export_df$symbol[export_df$entrez == control_gene] <- control_gene
    }
  }
  
  logger::log_info("Exporting results of initial MAUDE run...")
  
  csv_file_path <- sub(
    "\\.rds$",
    ".csv",
    file.path(
      cfg$paths$results_output_folder,
      paste0("MAUDE_Hits", file_suffix)
    )
  )
  
  readr::write_csv(export_df, csv_file_path)
  logger::log_info("Results of initial MAUDE run exported to: {csv_file_path}")
  
  excel_file_path <- sub(
    "\\.rds$",
    ".xlsx",
    file.path(
      cfg$paths$results_output_folder,
      paste0("MAUDE_Hits", file_suffix)
    )
  )
  
  writexl::write_xlsx(export_df, excel_file_path)
  logger::log_info("Results of initial MAUDE run exported to: {excel_file_path}")
  
  return(export_df)
}