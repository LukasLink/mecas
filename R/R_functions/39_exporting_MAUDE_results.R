handle_auto_combine_replicates_and_export <- function(file_info_suffix, opt){
  if (auto_combine_replicates){
    log_info("Combining replicates as auto_combine_replicates is set to true.")
    # if we have replicates, handle them now
    if (length(unique(maude_gene_stats$exp)) > 1){
      genes_with_rep <- add_info_wrapper(file_info_suffix)
      
      export_df <- genes_with_rep %>%
        group_by(symbol, entrez) %>%
        summarise(
          numGuides = sum(numGuides, na.rm = TRUE),
          seq       = first(seq),
          sgRNA     = first(sgRNA),
          meanZ     = mean(meanZ, na.rm = TRUE),
          
          # Stouffer’s method (equal weights); use count of non-NA Zs in the denominator
          stoufferZ = sum(significanceZ, na.rm = TRUE) /
            sqrt(sum(!is.na(significanceZ))),
          # For clarity: final significance Z = stoufferZ
          significanceZ = stoufferZ,
          
          .groups = "drop"
        ) %>%
        mutate(
          # two-sided p from Stouffer Z
          p.value = 2 * pnorm(abs(stoufferZ), lower.tail = FALSE),
          FDR     = p.adjust(p.value, method = "BH")
        )
      export_df <- export_df %>%
        select(symbol, entrez, numGuides, stoufferZ, meanZ, significanceZ, p.value, FDR, seq, sgRNA)
    } else {
      print("No replicates present, skipping replicate combination")
      export_df <- add_info_wrapper(file_info_suffix)
    }
  } else {
    export_df <- add_info_wrapper(file_info_suffix)
  }
  
  if (length(include_controls_list) > 0) {
    for (control_gene in include_controls_list) {
      export_df$symbol[export_df$entrez == control_gene] <- control_gene
    }
  }
  
  
  log_info("Exporting results of initial MAUDE run...")
  csv_file_path <- sub(".rds",".csv",file.path(results_output_folder,
                                               paste0("MAUDE_Hits", file_suffix)))
  write_csv(export_df, csv_file_path)
  log_info("Results of initial MAUDE run exported to: {csv_file_path}")
  excel_file_path <- sub("\\.rds$", ".xlsx", file.path(results_output_folder,
                                                       paste0("MAUDE_Hits", file_suffix)))
  write_xlsx(export_df, excel_file_path)
  log_info("Results of initial MAUDE run exported to: {excel_file_path}")
  
  return(export_df)
}

