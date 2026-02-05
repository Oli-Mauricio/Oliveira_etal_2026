suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(rstudioapi)
  library(tidyverse)
})

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#----------------------------
# Helpers
#----------------------------

# Parse sample metadata from column name:
# e.g. "Camk2a15min_TRAP_Shock4"
parse_sample_info <- function(sample_name) {
  cell_type  <- str_extract(sample_name, "^(Camk2a|Pvalb|Sst)")
  timepoint  <- str_extract(sample_name, "(15min|60min)")
  condition  <- str_extract(sample_name, "(HC|Box|Shock)")
  tech       <- str_extract(sample_name, "_(TRAP|TL)_") %>% str_replace_all("_", "")
  
  tibble(
    sample = sample_name,
    cell_type = cell_type,
    timepoint = timepoint,
    condition = condition,
    tech = tech
  )
}

# Run DESeq2 for one (cell_type, timepoint) subset and return results + PCA plot
run_deseq_for_subset <- function(count_df, cell_type, timepoint) {
  # 1) choose TRAP columns for this subset
  sample_cols <- colnames(count_df) %>%
    setdiff("Geneid") %>%
    keep(~ str_detect(.x, "^" %||% "") ) # placeholder for safety
  
  # Explicit filter using patterns in your naming scheme
  sample_cols <- colnames(count_df) %>%
    setdiff("Geneid") %>%
    keep(~ str_detect(.x, paste0("^", cell_type, timepoint, "_TRAP_")))
  
  if (length(sample_cols) < 2) {
    return(list(
      dds = NULL,
      results = NULL,
      pca = NULL,
      message = paste("Not enough samples for", cell_type, timepoint)
    ))
  }
  
  # 2) build colData from column names
  meta <- bind_rows(lapply(sample_cols, parse_sample_info)) %>%
    filter(tech == "TRAP") %>%                        # safety; should already be TRAP
    mutate(condition = factor(condition, levels = c("HC", "Box", "Shock"))) %>%
    mutate(condition_number = case_when(condition == "HC" ~ 0, condition == "Box" ~ 1, condition == "Shock" ~ 2)) %>%
    column_to_rownames("sample")
  
  # 3) build count matrix
  mat <- count_df %>%
    select(Geneid, all_of(rownames(meta))) %>%
    tibble::column_to_rownames("Geneid") %>%
    as.matrix()
  
  # Ensure integer-ish counts
  storage.mode(mat) <- "integer"
  
  #INSERT HERE COMBAT-SEQ WITH CORRECT BATCH VALUES
  #CREATE ONE OBJECT CONTAINING THE BATCHES PER DATASET
  #SELECT THE CORRECT ONE ITERATIVELY, MATCHING THE RUN IN THE FOR LOOP.
  
  # 4) DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData   = meta,
    design    = ~ condition
  )
  
  # Basic filtering (optional but common)
  keep <- rowSums(counts(dds) >= 20) >= 5
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  
  # 5) results: Shock vs HC, Box vs HC
  res_shock <- results(dds, contrast = c("condition", "Shock", "HC"))
  res_box   <- results(dds, contrast = c("condition", "Box",   "HC"))
  
  # Convert to data.frames with Geneid column
  res_shock_df <- as.data.frame(res_shock) %>%
    tibble::rownames_to_column("Geneid") %>%
    arrange(padj)
  
  res_box_df <- as.data.frame(res_box) %>%
    tibble::rownames_to_column("Geneid") %>%
    arrange(padj)
  
  # 6) PCA plot using DESeq2::plotPCA (colored by condition)
  vsd <- vst(dds, blind = FALSE)
  pca_df <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pca_df, "percentVar"))
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 5) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(paste(cell_type, timepoint, "TRAP")) +
    theme_bw(base_size = 20)
  
  list(
    dds = dds,
    results = list(
      Shock_vs_HC = res_shock_df,
      Box_vs_HC   = res_box_df
    ),
    pca = p_pca,
    meta = meta
  )
}

#----------------------------
# Main driver
#----------------------------

load("All_SampleCounts.RData")

# Identify which (cell_type, timepoint) combos exist in your TRAP columns
trap_cols <- colnames(counts_wide) %>%
  setdiff("Geneid") %>%
  keep(~ str_detect(.x, "_TRAP_")) %>%
  keep(~ !str_detect(.x, "_TL_"))  # extra safety; TL shouldn't match anyway

trap_meta <- bind_rows(lapply(trap_cols, parse_sample_info)) %>%
  filter(tech == "TRAP") %>%
  distinct(cell_type, timepoint) %>%
  arrange(cell_type, timepoint)

# Run all datasets, return a nested list:
# out[["Camk2a"]][["15min"]] etc.
out <- list()
for (i in seq_len(nrow(trap_meta))) {
  ct <- trap_meta$cell_type[i]
  tp <- trap_meta$timepoint[i]
  
  if (is.null(out[[ct]])) out[[ct]] <- list()
  
  out[[ct]][[tp]] <- run_deseq_for_subset(counts_wide, ct, tp)
}

# out is your final list:
# out$Camk2a$`15min`$results$Shock_vs_HC
# out$Camk2a$`15min`$pca  (a ggplot object you can print)
out


# Keep ONLY DE results (and optionally meta / pca)
extract_results_only <- function(out, keep_meta = TRUE, keep_pca = FALSE) {
  res_only <- list()
  
  for (ct in names(out)) {
    res_only[[ct]] <- list()
    for (tp in names(out[[ct]])) {
      
      x <- out[[ct]][[tp]]
      
      res_only[[ct]][[tp]] <- list(
        results = x$results
      )
      
      if (keep_meta && !is.null(x$meta)) {
        # meta can include rownames; store as plain column
        res_only[[ct]][[tp]]$meta <- tibble::rownames_to_column(
          as.data.frame(x$meta), var = "sample"
        )
      }
      
      if (keep_pca && !is.null(x$pca)) {
        res_only[[ct]][[tp]]$pca <- x$pca
      }
    }
  }
  
  res_only
}

out_results <- extract_results_only(out, keep_meta = TRUE, keep_pca = FALSE)

# Save compactly (RDS compresses well)
save(out_results, file = "DESeq2_results.RData")
