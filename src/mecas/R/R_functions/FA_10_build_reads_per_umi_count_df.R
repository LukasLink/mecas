# R/R_functions/FA_10_build_reads_per_umi_count_df.R

build_reads_per_umi_count_df <- function(
    bam_fpath,
    library,
    set_bin_name_to,
    set_sample_to,
    chunk_size = 10000000L,
    consolidate_every = 10L
) {
  
  sgrna_list <- library$seq
  n_sgrnas <- length(sgrna_list)
  
  bam <- Rsamtools::BamFile(bam_fpath, yieldSize = chunk_size)
  
  open(bam)
  on.exit(close(bam), add = TRUE)
  
  chunk_results <- list()
  chunks_since_consolidation <- 0L
  chunk_number <- 0L
  reads_processed <- 0L
  
  memory_used_gb <- function() {
    gc()[2, 2] / 1024
  }
  
  memory_rss_gb <- function() {
    status <- readLines("/proc/self/status")
    rss_line <- status[startsWith(status, "VmRSS:")]
    as.numeric(sub("VmRSS:\\s+([0-9]+) kB.*", "\\1", rss_line)) / 1024^2
  }
  
  repeat {
    
    chunk_start <- Sys.time()
    # Only the CB and UB optional BAM tags are loaded for chunk_size reads:
    #   CB: corrected sgRNA sequence
    #   UB: corrected the UMI sequence
    reads <- Rsamtools::scanBam(
      bam,
      param = Rsamtools::ScanBamParam(
        flag = Rsamtools::scanBamFlag(
          isSecondMateRead = FALSE
        ),
        what = character(),
        tag = c("CB", "UB")
      )
    )[[1]]
    
    if (length(reads$tag$CB) == 0L) { break } # End of BAm file
    
    #Progress counters
    chunk_number <- chunk_number + 1L
    chunks_since_consolidation <- chunks_since_consolidation + 1L
    reads_processed <- reads_processed + length(reads$tag$CB)
    
    sgRNA <- reads$tag$CB
    UMI <- reads$tag$UB
    
    keep <- !is.na(sgRNA) & !is.na(UMI)
    
    if (any(keep)) {
      
      # Remove records missing either of the required tags.
      sgRNA <- sgRNA[keep]
      UMI <- UMI[keep]
      
      # Convert each sgRNA sequence into its 1-based position in the sgrna_list
      sgRNA_id <- match(sgRNA, sgrna_list)
      
      keep_library_sgrna <- !is.na(sgRNA_id)
      
      if (any(keep_library_sgrna)) {
        
        sgRNA_id <- sgRNA_id[keep_library_sgrna]
        UMI <- UMI[keep_library_sgrna]
        
        # Construct a chunk-local lookup table for UMI sequences.
        # umi_values contains each UMI observed in this chunk once, while
        # UMI_id gives the 1-based position of every read's UMI in the table
        umi_values <- unique(UMI)
        UMI_id <- match(UMI, umi_values)
        
        # Encode each (sgRNA_id, UMI_id) pair as one numeric key:
        #
        #   pair_key = sgRNA_id + n_sgrnas * (UMI_id - 1)
        #
        # Because sgRNA_id ranges from 1 to n_sgrnas, every possible
        # (sgRNA_id, UMI_id) combination receives a unique key. This lets
        # the code count pairs using numeric sorting and run-length encoding
        # instead of grouping the original character vectors.
        #
        # Doubles are used to avoid integer overflow when many distinct UMIs
        # produce large encoded keys.
        pair_key <- as.double(sgRNA_id) + as.double(n_sgrnas) * (as.double(UMI_id) - 1)
        
        # Sorting places all identical (sgRNA, UMI) keys next to each other.
        ordered_key <- sort(pair_key)
        
        # Run-length encoding then returns:
        #   runs$values  = the unique encoded pair keys
        #   runs$lengths = the number of reads carrying each pair
        runs <- rle(ordered_key)
        
        unique_keys <- runs$values
        
        # Decode each pair key back into its original sgRNA and UMI IDs.
        #
        # Modulo recovers the position within an n_sgrnas-sized block,
        # which is the sgRNA ID.
        sgRNA_id <- ((unique_keys - 1) %% n_sgrnas) + 1L
        # Integer division by n_sgrnas identifies the block containing the
        # key, which is the chunk-local UMI ID.
        umi_id <- floor((unique_keys - 1) / n_sgrnas) + 1L
        
        # Store one row per distinct (sgRNA, UMI) pair observed in this
        # chunk. `count` is the number of reads carrying that exact pair.
        #
        # Identical pairs found in later chunks remain separate here and
        # are combined during periodic or final consolidation.
        chunk_results[[length(chunk_results) + 1L]] <- data.frame(
          sgRNA = sgrna_list[sgRNA_id],
          UMI = umi_values[umi_id],
          count = runs$lengths,
          stringsAsFactors = FALSE
        )
      }
    }
    
    if (chunks_since_consolidation >= consolidate_every) {
      
      if (length(chunk_results) > 0L) {
        chunk_results <- list(
          dplyr::bind_rows(chunk_results) %>%
            dplyr::group_by(sgRNA, UMI) %>%
            dplyr::summarise(count = sum(count), .groups = "drop")
        )
      }
      
      chunks_since_consolidation <- 0L
    }
    
    chunk_time <- as.numeric(
      difftime(Sys.time(), chunk_start, units = "secs")
    )
    
    log_info(sprintf(
      "Finished chunk %d | reads: %s | time: %.2f s | R memory: %.2f GB | RSS: %.2f GB",
      chunk_number,
      format(reads_processed, big.mark = ","),
      chunk_time,
      memory_used_gb(),
      memory_rss_gb()
    ))
  }
  
  
  if (length(chunk_results) == 0L) {
    result <- data.frame(
      sgRNA = character(),
      UMI = character(),
      count = integer(),
      stringsAsFactors = FALSE
    )
  } else {
    result <- dplyr::bind_rows(chunk_results) %>%
      dplyr::group_by(sgRNA, UMI) %>%
      dplyr::summarise(count = sum(count), .groups = "drop")
  }
  
  # Modify it to be in line with count_df_long formating
  return(
    result %>%
      left_join(
        library %>%
          select(seq, sgrna_id, type),
        by = c("sgRNA" = "seq")
      ) %>%
      rename(
        sgrna_seq = sgRNA,
        sgRNA = sgrna_id,
        group_category = type
      ) %>%
      mutate(
        bin_name = set_bin_name_to,
        sample = set_sample_to,
        sublib = UMI,
        exp = paste(sublib, sample, sep = "_")
      ) %>%
      select(c(sgRNA, sublib, bin_name, count, group_category, sample, exp))
  )
  
}