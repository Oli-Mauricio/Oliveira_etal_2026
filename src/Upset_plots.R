suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(UpSetR)
  library(rstudioapi)
})

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load("DESeq2_results.RData")

# Build binary DE membership table (Geneid x {Box/Shock}x{15/60}) for one cell type
make_binary_de_table <- function(out_results,
                                 cell_type,
                                 padj_thresh = 0.05,
                                 timepoints = c("15min", "60min")) {
  
  stopifnot(cell_type %in% names(out_results))
  
  # Helper: extract one contrast table and return Geneid + binary
  extract_binary <- function(res_df, col_name, padj_thresh) {
    res_df %>%
      transmute(
        Geneid = Geneid,
        !!col_name := as.integer(!is.na(padj) & padj < padj_thresh)
      )
  }
  
  pieces <- list()
  
  for (tp in timepoints) {
    obj <- out_results[[cell_type]][[tp]]
    
    # safety checks
    if (is.null(obj) || is.null(obj$results)) next
    if (is.null(obj$results$Box_vs_HC) || is.null(obj$results$Shock_vs_HC)) next
    
    pieces[[paste0("Box_", tp)]] <- extract_binary(obj$results$Box_vs_HC, paste0("Box_", tp), padj_thresh)
    pieces[[paste0("Shock_", tp)]] <- extract_binary(obj$results$Shock_vs_HC, paste0("Shock_", tp), padj_thresh)
  }
  
  if (length(pieces) == 0) {
    return(NULL)
  }
  
  # Full join across columns to keep union of genes; missing => 0
  bin_df <- Reduce(function(x, y) full_join(x, y, by = "Geneid"), pieces) %>%
    mutate(across(-Geneid, ~ tidyr::replace_na(.x, 0L))) %>%
    distinct(Geneid, .keep_all = TRUE)
  
  bin_df
}

# Make UpSet plot for one cell type (draws plot; also returns data)
plot_upset_for_cell_type <- function(out_results,
                                     cell_type,
                                     padj_thresh = 0.05,
                                     timepoints = c("15min", "60min"),
                                     set_order = c("Box_15min", "Shock_15min",
                                                   "Box_60min", "Shock_60min"),
                                     nsets = 4,
                                     point.size = 5,
                                     line.size = 2,
                                     text.scale = 1.6,
                                     sets.bar.color = c("#008631","#c23b22","#9DD6AD","#ffecec"),
                                     matrix.color = "#008631",
                                     shade.color = "grey",
                                     main.bar.color = "#008631",
                                     mainbar.y.label = "Number of DEGs") {
  
  bin_df <- make_binary_de_table(out_results, cell_type, padj_thresh, timepoints)
  if (is.null(bin_df)) stop("No usable results found for cell type: ", cell_type)
  
  cols <- intersect(set_order, colnames(bin_df))
  upset_input <- bin_df %>% dplyr::select(dplyr::all_of(cols))
  
  # --- ensure sets.bar.color length matches number of sets exactly ---
  n_sets_used <- length(cols)
  if (length(sets.bar.color) == 1) {
    sets.bar.color <- rep(sets.bar.color, n_sets_used)
  } else if (length(sets.bar.color) != n_sets_used) {
    # repeat or truncate safely
    sets.bar.color <- rep(sets.bar.color, length.out = n_sets_used)
  }
  
  print(UpSetR::upset(
    upset_input,
    sets = cols,
    nsets = min(nsets, n_sets_used),
    order.by = "freq",
    point.size = point.size,
    line.size = line.size,
    text.scale = text.scale,
    sets.bar.color = sets.bar.color,
    matrix.color = matrix.color,
    shade.color = shade.color,
    main.bar.color = main.bar.color,
    mainbar.y.label = mainbar.y.label
  ))
  
  invisible(list(binary_table = bin_df, upset_input = upset_input, sets_used = cols))
}


# Run for all cell types; returns list of binary tables (plots are drawn when called)
make_all_binary_tables <- function(out_results,
                                   padj_thresh = 0.05,
                                   timepoints = c("15min", "60min")) {
  out <- list()
  for (ct in names(out_results)) {
    out[[ct]] <- make_binary_de_table(out_results, ct, padj_thresh, timepoints)
  }
  out
}


binary_tables <- make_all_binary_tables(out_results, padj_thresh = 0.05)

plot_upset_for_cell_type(
  out_results,
  cell_type = "Sst",
  sets.bar.color = c("#008631", "#c23b22", "#9DD6AD", "#ffecec"),
  matrix.color = "#008631",
  main.bar.color = "#008631"
)
