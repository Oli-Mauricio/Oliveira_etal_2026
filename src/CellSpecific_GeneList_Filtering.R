#

# This is the filtering routine to obtain a list of common genes identified in all sequencing experiments (home cage).
# It converts raw FeatureCount files into RPKM values, then  

# Load all data sets, converts raw values into RPKM, then produces one data frame of TRAP / Total Lysate relative vals
# Filters samples that were not enriched to signature genes Camk2a, Pvalb and Sst
# Filters dataset to contain only genes with 1.2-fold enrichment (~20%) over total lysate as cell type-specific

### NOTE:
# The high incidence of samples with low purity in interneurons is expected, due to low input.

#

suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(edgeR)
})

Raw_to_RPKM <- function(df_interest, length_df) {
  
  df_interest <- as.data.frame(df_interest)
  
  gene.length.temp <- length_df %>% 
    filter(Geneid %in% rownames(df_interest))
  
  df_interest <- df_interest %>% 
    filter(rownames(df_interest) %in% gene.length.temp$Geneid)
  
  match_idx <- match(rownames(df_interest), gene.length.temp$Geneid)
  df_interest <- df_interest[match_idx,]
  
  
  dge.list <- DGEList(df_interest)
  
  rpkm.values <- as.data.frame(rpkm(dge.list, gene.length = gene.length.temp$Length))
  
  return(rpkm.values)
  
}


clean_sample_name <- function(x) {
  b <- basename(x)
  
  # Drop junk prefixes up to last dot before cell-type token
  b <- sub("^.*\\.(?=(Camk2a|Pvalb|Sst))", "", b, perl = TRUE)
  
  # Capture "<CellTime>_(TL|TRAP)_<Group>(_<Extra>)?" before "_S<digits>"
  m <- str_match(b, "((?:Camk2a|Pvalb|Sst)[^_]*_(?:TL|TRAP)_[^_]+(?:_[^_]+)?)_S\\d+")
  if (!is.na(m[, 2])) return(m[, 2])
  
  # Fallback: trim typical suffixes
  b <- sub("Aligned.*$", "", b)
  b <- sub("\\.bam$", "", b)
  b <- sub("_L\\d+.*$", "", b)
  b
}

read_fc_select <- function(fp, col7 = 7) {
  # fread is much more forgiving than read.delim/read.table
  dt <- read.delim(fp, comment.char = "#")
  
  if (!("Geneid" %in% names(dt))) {
    stop(sprintf("File %s does not contain a 'Geneid' column.", fp))
  }
  if (ncol(dt) < col7) {
    stop(sprintf("File %s has %d columns; cannot extract column %d.", fp, ncol(dt), col7))
  }
  
  sample_nm <- clean_sample_name(names(dt)[col7])
  
  dt %>%
    dplyr::select(Geneid, !!sample_nm := all_of(names(dt)[col7]))
}

# ----------------------------------------------------------------------------------
# MAIN: READ FILES, CONVERT THEM INTO A SINGLE DATA FRAME, THEN CALCULATE ENRICHMENT
# ----------------------------------------------------------------------------------

fc_dir   <- "C:\\Users\\mauri\\Dropbox\\TRAP_Alt_Analysis\\FeatCounts"
fc_files <- sort(list.files(fc_dir, pattern = "\\.txt$|\\.tsv$|\\.tab$|L[0-9]+$", full.names = TRUE))

# Read all files, stacking them on a list. Then transform them into a data frame using reduce
count_list <- lapply(fc_files, function(f) {
  tryCatch(
    read_fc_select(f, col7 = 7),
    error = function(e) {
      message("\nFAILED on file:\n  ", f, "\nReason:\n  ", conditionMessage(e), "\n")
      stop(e)
    }
  )
})

counts_wide <- reduce(count_list, inner_join, by = "Geneid")
save(counts_wide, file = "All_counts.RData")

# Build length data frame to use in the conversion to RPKM.
last_fp <- tail(fc_files, 1)
last_df <- read.delim(
  last_fp, 
  comment.char = "#"
)

if (!all(c("Geneid", "Length") %in% names(last_df))) {
  stop(sprintf("Last file (%s) does not contain both 'Geneid' and 'Length' columns.", last_fp))
}

length_df <- last_df %>% dplyr::select(Geneid, Length)

#Conversion to RPKM
rownames(counts_wide) <- counts_wide$Geneid
RPKM.df <- Raw_to_RPKM(counts_wide[,-1], length_df)

#Enrichment data frame. Divide the counts in TRAP fraction by the counts in TL fraction
RPKM.div.df <- data.frame(
  "HC1_Camk2a" = RPKM.df$Camk2a15min_TRAP_HC1 / RPKM.df$Camk2a15min_TL_HC1,
  "HC2_Camk2a" = RPKM.df$Camk2a15min_TRAP_HC2 / RPKM.df$Camk2a15min_TL_HC2,
  "HC3_Camk2a" = RPKM.df$Camk2a15min_TRAP_HC3 / RPKM.df$Camk2a15min_TL_HC3,
  "HC4_Camk2a" = RPKM.df$Camk2a60min_TRAP_HC1 / RPKM.df$Camk2a60min_TL_HC1,
  "HC5_Camk2a" = RPKM.df$Camk2a60min_TRAP_HC2 / RPKM.df$Camk2a60min_TL_HC2,
  "HC6_Camk2a" = RPKM.df$Camk2a60min_TRAP_HC3 / RPKM.df$Camk2a60min_TL_HC3,
  "HC7_Camk2a" = RPKM.df$Camk2a60min_TRAP_HC4 / RPKM.df$Camk2a60min_TL_HC4,
  "HC1_Pvalb" = RPKM.df$Pvalb15min_TRAP_HC1 / RPKM.df$Pvalb15min_TL_HC1,
  #HC2 included due to comparable marker expression, but not included for cell type-specific gene filtering.
  #"HC2_Pvalb" = RPKM.df$Pvalb15min_TRAP_HC2 / 1, 
  "HC3_Pvalb" = RPKM.df$Pvalb15min_TRAP_HC3 / RPKM.df$Pvalb15min_TL_HC3,
  "HC4_Pvalb" = RPKM.df$Pvalb15min_TRAP_HC4 / RPKM.df$Pvalb15min_TL_HC4,
  "HC5_Pvalb" = RPKM.df$Pvalb60min_TRAP_HC1 / RPKM.df$Pvalb60min_TL_HC1,
  "HC6_Pvalb" = RPKM.df$Pvalb60min_TRAP_HC2 / RPKM.df$Pvalb60min_TL_HC2,
  "HC7_Pvalb" = RPKM.df$Pvalb60min_TRAP_HC3 / RPKM.df$Pvalb60min_TL_HC3,
  "HC8_Pvalb" = RPKM.df$Pvalb60min_TRAP_HC4 / RPKM.df$Pvalb60min_TL_HC4,
  "HC9_Pvalb" = RPKM.df$Pvalb60min_TRAP_HC5 / RPKM.df$Pvalb60min_TL_HC5,
  "HC10_Pvalb" = RPKM.df$Pvalb60min_TRAP_HC6 / RPKM.df$Pvalb60min_TL_HC6,
  "HC1_Sst" = RPKM.df$Sst15min_TRAP_HC1 / RPKM.df$Sst15min_TL_HC1,
  "HC2_Sst" = RPKM.df$Sst15min_TRAP_HC2 / RPKM.df$Sst15min_TL_HC2,
  "HC3_Sst" = RPKM.df$Sst15min_TRAP_HC3 / RPKM.df$Sst15min_TL_HC3,
  "HC4_Sst" = RPKM.df$Sst15min_TRAP_HC4 / RPKM.df$Sst15min_TL_HC4,
  "HC5_Sst" = RPKM.df$Sst60min_TRAP_HC1 / RPKM.df$Sst60min_TL_HC1,
  "HC6_Sst" = RPKM.df$Sst60min_TRAP_HC2 / RPKM.df$Sst60min_TL_HC2,
  "HC7_Sst" = RPKM.df$Sst60min_TRAP_HC3 / RPKM.df$Sst60min_TL_HC3,
  "HC8_Sst" = RPKM.df$Sst60min_TRAP_HC4 / RPKM.df$Sst60min_TL_HC4,
  "HC9_Sst" = RPKM.df$Sst60min_TRAP_HC5 / RPKM.df$Sst60min_TL_HC5,
  row.names = rownames(RPKM.df)
)

# Repairs raw-to-rpkm conversion errors (e.g., rawCounts=0 in total lysate - the denominator - generates Inf)
RPKM.div.df[is.na(RPKM.div.df)] <- 0
RPKM.div.df[RPKM.div.df == "Inf"] <- 0
RPKM.div.df[RPKM.div.df == "-Inf"] <- 0

#
## Filter out non-enriched samples
gene.markers <- c("ENSMUSG00000024617.16" ,"ENSMUSG00000005716.16", "ENSMUSG00000004366.4")

markers.df <- RPKM.div.df %>%
  filter(rownames(RPKM.div.df) %in% gene.markers)

#Split the main dataset to filter gene names that are enriched in each cell type

Camk2a.df <- RPKM.div.df %>%
  dplyr::select(contains("Camk2a")) %>%
  filter(rowSums(. >= 1.2) >= 4)

Pvalb.df <- RPKM.div.df %>%
  dplyr::select(contains("Pvalb")) %>%
  filter(rowSums(. >= 1.2) >= 6)

Sst.df <- RPKM.div.df %>%
  dplyr::select(contains("Sst")) %>%
  filter(rowSums(. >= 1.2) >= 6)

save(Camk2a.df, Pvalb.df, Sst.df, file = "Filtered_gene_cell_Specific.RData")
