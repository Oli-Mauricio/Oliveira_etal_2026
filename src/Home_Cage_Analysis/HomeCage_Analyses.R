suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(rtracklayer)
  library(ggrepel)
  library(ggpubr)
  library(tibble)
})

#FIGURE 1 - CELL TYPE-SPECIFIC TRANSLATOME#####
#Fig 1E
#-------------------------------------------------------------------------------
##PCA of all genes enriched in neurons
#first get a full list of all genes identified as enriched in all cell types
load("Filtered_gene_cell_Specific.RData")
all_genes <- c(rownames(Camk2a.df), rownames(Pvalb.df), rownames(Sst.df))
all_genes <- unique(all_genes)

#Prepare counts data frame only with Home Cage samples
load("All_counts.RData")
All.HC.samples <- counts_wide %>% 
  dplyr::select(Geneid, contains("TRAP")) %>%
  dplyr::select(Geneid, contains("HC")) %>%
  filter(Geneid %in% all_genes) %>%
  dplyr::distinct(Geneid, .keep_all = TRUE)
rownames(All.HC.samples) <- All.HC.samples$Geneid
All.HC.samples <- dplyr::select(All.HC.samples, -Geneid)

#Prepare Metadata for DESeq2 run
Metadata_df <- data.frame(
  samples = names(All.HC.samples),
  Neuron_type = str_replace(names(All.HC.samples), "[0-9]+min_TRAP_HC[0-9]$", ""),
  row.names = names(All.HC.samples)
)  

my_database_TRAP_all <- DESeqDataSetFromMatrix(countData = as.matrix(All.HC.samples),
                                               colData = Metadata_df,
                                               design = ~ Neuron_type)

keep_TRAP <- rowSums(counts(my_database_TRAP_all) >= 20) >= 5
my_database_TRAP_all <- my_database_TRAP_all[keep_TRAP, ]
my_database_TRAP_all <- estimateSizeFactors(my_database_TRAP_all)
normalized_db_TRAP_all <- vst(my_database_TRAP_all, blind = TRUE)

#PCA all
plotPCA(normalized_db_TRAP_all, intgroup = c("Neuron_type")) + 
  geom_point(size = 7) +
  scale_color_manual(values = c("#b1c6d9","#008631","#ffecec")) +
  #geom_text(aes(label = name)) +
  xlim(-40,40) + 
  ylim(-40, 40) + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        legend.position = "bottom",
        panel.grid = element_blank()) 


# -------------------------------------------------------------------------
# CELL MARKERS PLOT
# -------------------------------------------------------------------------

gtf_path <- "C:\\Users\\mauri\\Dropbox\\rMATS_AS_translatome project\\gencode.vM25.annotation.gtf"

# 1) Build mapping table from GTF
gtf <- rtracklayer::import(gtf_path)

gene_map <- as.data.frame(gtf) %>%
  filter(type == "gene") %>%
  transmute(
    Geneid    = as.character(gene_id),
    gene_name = as.character(gene_name)
  ) %>%
  distinct() %>%
  filter(!is.na(Geneid), !is.na(gene_name))

# If your rownames look like ENSMUSG... .xx, strip version so it matches GTF gene_id
strip_ensembl_version <- function(x) sub("\\.[0-9]+$", "", x)

marker_genes <- c(
  "Camk2a", "Slc17a7",
  "Gad1", "Gad2", "Syt2",
  "Pvalb", "Sst", "Aldh1l1",
  "Gfap", "Cnp", "Cx3cr1"
)

# 2) Convert rownames (Ensembl) -> symbols, then plot as before
markers.df <- RPKM.div.df %>%
  rownames_to_column("Geneid") %>%
  left_join(gene_map, by = "Geneid") %>%
  filter(gene_name %in% marker_genes)

# 3) Long format (tidyr instead of gather)
gather.df <- markers.df %>%
  pivot_longer(
    cols = HC1_Camk2a:HC9_Sst,
    names_to = "sample",
    values_to = "normalized_RPKM"
  ) %>%
  mutate(
    geneID = factor(gene_name, levels = marker_genes),
    `Cell type` = str_replace(sample, "^.*_", "")
  )

ggplot(gather.df, aes(`Cell type`, log2(normalized_RPKM), fill = `Cell type`)) + 
  geom_boxplot(linewidth = 0.8) + 
  geom_point() +
  geom_hline(yintercept = 0, linewidth = 1) +
  scale_fill_manual(values = c("#b1c6d9","#008631","#ffecec")) + 
  ylab("log2(RPKM TRAP /\n RPKM Total)") +
  scale_y_continuous(breaks = seq(-8, 8, by = 2)) +
  facet_wrap(~ geneID, nrow = 2) + 
  theme_bw() + 
  theme(
    axis.title.y = element_text(size = 24),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=24),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 24),
    strip.text = element_text(size = 22),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )



# -----------------------------------------------------------------------
# PRODUCTION OF HOME CAGE DESEQ2 OBJECT
# -----------------------------------------------------------------------

run_pairwise_neuron_type_results <- function(dds,
                                             neuron_col = "Neuron_type",
                                             alpha = 0.05,
                                             cooksCutoff = FALSE,
                                             independentFiltering = FALSE) {
  # Make sure the factor exists and is a factor
  stopifnot(neuron_col %in% colnames(colData(dds)))
  colData(dds)[[neuron_col]] <- factor(colData(dds)[[neuron_col]])
  
  # Fit model once
  dds <- DESeq(dds)
  
  # Levels present (e.g., Camk2a, Pvalb, Sst)
  lvls <- levels(colData(dds)[[neuron_col]])
  lvls <- lvls[lvls %in% unique(as.character(colData(dds)[[neuron_col]]))]  # keep used levels
  
  # All pairwise combinations
  pairs <- combn(lvls, 2, simplify = FALSE)
  
  # Store results here
  res_list <- list()
  
  for (p in pairs) {
    a <- p[1]
    b <- p[2]
    
    # Name like "Camk2a_vs_Pvalb"
    nm <- paste0(a, "_vs_", b)
    
    res <- results(
      dds,
      contrast = c(neuron_col, a, b),
      alpha = alpha,
      cooksCutoff = cooksCutoff,
      independentFiltering = independentFiltering
    )
    
    # Full, unfiltered table; keep Geneid as a column for easy joins later
    res_df <- as.data.frame(res) %>%
      rownames_to_column("Geneid") %>%
      arrange(padj)
    
    res_list[[nm]] <- res_df
  }
  
  list(dds = dds, results = res_list, neuron_levels = lvls)
}


out_HC <- run_pairwise_neuron_type_results(
  dds = my_database_TRAP_all,
  neuron_col = "Neuron_type",
  alpha = 0.05,
  cooksCutoff = FALSE,
  independentFiltering = FALSE
)

# Your results list:
# out_HC$results$Camk2a_vs_Pvalb
# out_HC$results$Camk2a_vs_Sst
# out_HC$results$Pvalb_vs_Sst


#-------------------------------------------------------------------------------
# VENN DIAGRAM OF HOME CAGE SAMPLES
#-------------------------------------------------------------------------------


Venn_df <- list(Camk2a = rownames(Camk2a.df),
           Pvalb = rownames(Pvalb.df),
           Sst = rownames(Sst.df))


attributes(Venn_df) <- list(names = names(Venn_df), 
                       row.names = 1:7000, #Random number much larger than the total number of a cell type-specific gene list 
                       class = 'data.frame')

#Venn diagram
ggVennDiagram(Venn_df, label_alpha = 0, 
              set_color = "midnightblue", 
              label_size = 8,
              label_percent_digit = 1,
              set_size = 8) + 
  scale_fill_gradient(low = "#F4FAFE", 
                      high = "#4981BF")

# ------------------------------------------
# COMPARISON OF CAMK2A-TO-INTERNEURON DEGs 
# ------------------------------------------

common_genes <- Reduce(intersect, list(rownames(Camk2a.df),
                                       rownames(Pvalb.df),
                                       rownames(Sst.df)))

res_Camk2a_vs_Pvalb <- out_HC$results$Camk2a_vs_Pvalb
res_Camk2a_vs_Sst   <- out_HC$results$Camk2a_vs_Sst
res_Pvalb_vs_Sst    <- out_HC$results$Pvalb_vs_Sst

df.final <- tibble(Geneid = common_genes) %>%
  left_join(
    res_Camk2a_vs_Pvalb %>% select(Geneid, Camk2avsPV = log2FoldChange),
    by = "Geneid"
  ) %>%
  left_join(
    res_Camk2a_vs_Sst %>% select(Geneid, Camk2avsSst = log2FoldChange),
    by = "Geneid"
  ) %>%
  left_join(
    res_Pvalb_vs_Sst %>% select(Geneid, padj_PvalbvsSst = padj),
    by = "Geneid"
  ) %>%
  mutate(
    sig_PvalbvsSst = if_else(!is.na(padj_PvalbvsSst) & padj_PvalbvsSst < 0.05, "yes", "no")
  )

gene.markers <- c("ENSMUSG00000024617.16" ,"ENSMUSG00000005716.16", "ENSMUSG00000004366.4")

ggplot(df.final, aes(Camk2avsPV, Camk2avsSst)) + 
  geom_point(data = df.final, aes(fill = sig_PvalbvsSst), size = 3, shape = 21) + 
  scale_fill_manual(values = c("yes" = "darkblue", "no" = "gray")) +
  geom_smooth(method = 'lm', colour = "black", linewidth = 1, alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  xlim(c(-3,3)) + ylim (c(-3,3)) +
  xlab(label = "log2FC - Camk2a/Pvalb") + ylab(label = "log2FC - Camk2a/Sst") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.position = "bottom") +
  stat_cor(method = "spearman", size = 7, label.x = -3, label.y = 2.5, cor.coef.name = "rho") +
  geom_label_repel(aes(label = ifelse(Geneid %in% gene.markers, row.names(df.final), '')),
                   max.overlaps = Inf, box.padding = 1
  )
