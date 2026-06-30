add_info_to_gene_stats <- function(maude_guide_stats,
                                   maude_gene_stats,
                                   cfg,
                                   merged_sgRNA_df = cfg$merged_sgRNA_df,
                                   data_dir = cfg$paths$data_dir) {
  
  maude_guide_stats <- maude_guide_stats %>%
    dplyr::left_join(
      merged_sgRNA_df %>%
        dplyr::mutate(entrez = as.character(entrez)) %>% 
        dplyr::select(sgrna_id, entrez, seq, symbol) %>%
        dplyr::distinct(sgrna_id, .keep_all = TRUE),
      by = c("sgRNA" = "sgrna_id")
    ) %>%
    dplyr::mutate(entrez = dplyr::coalesce(as.character(entrez), sgRNA))
  
  maude_guide_stats <- maude_guide_stats %>%
    dplyr::mutate(abs_meanZ = abs(mean)) %>%
    dplyr::group_by(entrez) %>%
    dplyr::slice_max(order_by = abs_meanZ, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  maude_gene_stats <- maude_gene_stats %>%
    dplyr::left_join(
      maude_guide_stats %>%
        dplyr::select(entrez, seq, sgRNA, symbol),
      by = "entrez"
    )
  
  export_df <- maude_gene_stats %>% 
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
    ) %>% 
    dplyr::arrange(significanceZ)
  
  # Add Liangfu's expression data
  HepG2_cpm_Xie <- readRDS(file.path(data_dir, "Xie_hepato_df_cpm_RNAseq.rds"))
  
  export_df <- export_df %>%
    dplyr::left_join(
      HepG2_cpm_Xie %>%
        dplyr::filter(!is.na(entrez_id)) %>%
        dplyr::distinct(entrez_id, .keep_all = TRUE) %>% 
        dplyr::select(
          entrez_id,
          HepG2_CPM_Xie = HepG2,
          primary_hepato_CPM_Xie = primary
        ),
      by = c("entrez" = "entrez_id")
    )
  
  HepG2_tpm_Xie <- readRDS(file.path(data_dir, "Xie_hepato_df_tpm_RNAseq.rds"))
  
  export_df <- export_df %>%
    dplyr::left_join(
      HepG2_tpm_Xie %>%
        dplyr::filter(!is.na(entrez_id)) %>%
        dplyr::distinct(entrez_id, .keep_all = TRUE) %>% 
        dplyr::select(
          entrez_id,
          HepG2_TPM_Xie = HepG2,
          primary_hepato_TPM_Xie = primary
        ),
      by = c("entrez" = "entrez_id")
    )
  
  # Add Nico's HepG2 RNAseq expression data
  HepG2_tpm_Battisti <- readRDS(file.path(data_dir, "Battisti_HepG2_dual_rep_RNAseq.rds"))
  
  export_df <- export_df %>%
    dplyr::left_join(
      HepG2_tpm_Battisti %>%
        dplyr::filter(!is.na(entrez)) %>%
        dplyr::distinct(entrez, .keep_all = TRUE) %>% 
        dplyr::select(
          entrez,
          HepG2_CPM_Battisti = CPM_WT,
          HepG2_14D9_CPM_Battisti = CPM_dual_rep_14D9,
          HepG2_14B11_CPM_Battisti = CPM_dual_rep_14B11
        ),
      by = "entrez"
    )
  
  # Add information about essential genes
  essential_df <- readRDS(file.path(data_dir, "dependency_df.rds")) %>%
    dplyr::mutate(entrez = as.character(entrez_id)) %>% 
    dplyr::select(-c(gene_symbol, entrez_id))
  
  export_df <- export_df %>%
    dplyr::left_join(essential_df, by = "entrez")
  
  cardio_cpm_Xie <- readRDS(file.path(data_dir, "Xie_cardio_df_cpm_RNAseq.rds"))
  
  export_df <- export_df %>%
    dplyr::left_join(
      cardio_cpm_Xie %>%
        dplyr::filter(!is.na(entrez_id)) %>%
        dplyr::distinct(entrez_id, .keep_all = TRUE) %>% 
        dplyr::select(
          entrez_id,
          cardio_IPSC_CPM = IPS_d,
          cardio_primary_CPM = primary
        ),
      by = c("entrez" = "entrez_id")
    )
  
  cardio_tpm_Xie <- readRDS(file.path(data_dir, "Xie_cardio_df_tpm_RNAseq.rds"))
  
  export_df <- export_df %>%
    dplyr::left_join(
      cardio_tpm_Xie %>%
        dplyr::filter(!is.na(entrez_id)) %>%
        dplyr::distinct(entrez_id, .keep_all = TRUE) %>% 
        dplyr::select(
          entrez_id,
          cardio_IPSC_TPM = IPS_d,
          cardio_primary_TPM = primary
        ),
      by = c("entrez" = "entrez_id")
    )
  
  cardio_cpm_Xie_Liangfu_own_data <- readRDS(
    file.path(data_dir, "Xie_cardio_df_cpm_Liangfus_own_RNAseq.rds")
  )
  
  export_df <- export_df %>%
    dplyr::left_join(
      cardio_cpm_Xie_Liangfu_own_data %>%
        dplyr::filter(!is.na(entrez_id)) %>%
        dplyr::distinct(entrez_id, .keep_all = TRUE) %>%
        dplyr::mutate(
          both = rowMeans(dplyr::across(c(Xie_CM_Rep_1, Xie_CM_Rep_2)), na.rm = TRUE)
        ) %>% 
        dplyr::select(
          entrez_id,
          cardio_IPSC_CPM_Xie = both
        ),
      by = c("entrez" = "entrez_id")
    )
  
  cardio_tpm_Xie_Liangfus_own_data <- readRDS(
    file.path(data_dir, "Xie_cardio_df_tpm_Liangfus_own_RNAseq.rds")
  )
  
  export_df <- export_df %>%
    dplyr::left_join(
      cardio_tpm_Xie_Liangfus_own_data %>%
        dplyr::filter(!is.na(entrez_id)) %>%
        dplyr::distinct(entrez_id, .keep_all = TRUE) %>%
        dplyr::mutate(
          both = rowMeans(dplyr::across(c(Xie_CM_Rep_1, Xie_CM_Rep_2)), na.rm = TRUE)
        ) %>% 
        dplyr::select(
          entrez_id,
          cardio_IPSC_TPM_Xie = both
        ),
      by = c("entrez" = "entrez_id")
    )
  
  # Add GO term information
  go_map <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = as.character(export_df$entrez),
      keytype = "ENTREZID",
      columns = c("GO", "ONTOLOGY")
    )
  )
  
  go_collapsed <- go_map %>%
    dplyr::mutate(ENTREZID = as.character(ENTREZID)) %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::summarise(
      GO_terms = paste(unique(stats::na.omit(GO)), collapse = ";"),
      GO_ontology = paste(unique(stats::na.omit(ONTOLOGY)), collapse = ";"),
      go_count = dplyr::n_distinct(stats::na.omit(GO)),
      .groups = "drop"
    )
  
  export_df <- export_df %>%
    dplyr::left_join(go_collapsed, by = c("entrez" = "ENTREZID"))
  
  return(export_df)
}


add_info_wrapper <- function(
    cfg,
    file_suffix = cfg$suffix$file_suffix,
    folder = cfg$paths$rds_output_folder,
    prefix_gene = "MAUDE_gene_stats",
    prefix_guide = "MAUDE_guide_stats") {
  
  gene <- readRDS(file.path(folder, paste0(prefix_gene, file_suffix)))
  guide <- readRDS(file.path(folder, paste0(prefix_guide, file_suffix)))
  
  result <- add_info_to_gene_stats(
    maude_guide_stats = guide,
    maude_gene_stats = gene,
    cfg = cfg
  )
  
  return(result)
}

check_UMI_counts <- function(guide_df, reference_df) {
  guide_df <- guide_df %>%
    dplyr::left_join(
      reference_df %>%
        dplyr::select(sgrna_id, entrez), 
      by = c("sgRNA" = "sgrna_id")
    )
  
  result_df <- guide_df %>%
    dplyr::group_by(sublib, sample) %>%
    dplyr::summarise(
      input_sum = sum(input, na.rm = TRUE),
      lower_sum = sum(lower, na.rm = TRUE),
      upper_sum = sum(upper, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(result_df)
}