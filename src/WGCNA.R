# ============================================================
# WGCNA from DESeq2-VST + Cytoscape export (clean functional form)
# Keeps your original parameter values the same.
# ============================================================

suppressPackageStartupMessages({
  library("sva")
  library("ggplot2")
  library("WGCNA")
  library("tidyverse")
  library("clusterProfiler")
  library(rstudioapi)
  library(caret)
  library(mice)
  library(DESeq2)
  library(ggpubr)
  library(anRichment) #INSTALL THIS PACKAGE AGAIN
})

# WGCNA threads (optional; keep if you like)
WGCNA::allowWGCNAThreads()

# ---- small helper: capture base plots as objects (recordedplot) ----
.capture_plot <- function(expr) {
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 1400, height = 900, res = 150)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)
  
  force(expr)
  grDevices::recordPlot()
}

# ---- small helper: coerce vst input into WGCNA "datExpr" (samples x genes) ----
.as_wgcna_datExpr <- function(vst_df) {
  x <- as.data.frame(vst_df)
  
  # If it looks like genes x samples, transpose
  # (typical VST assay has genes in rows, samples in columns)
  if (nrow(x) > ncol(x)) {
    datExpr <- t(x)
  } else {
    # could already be samples x genes; keep as-is
    datExpr <- as.matrix(x)
  }
  
  datExpr <- as.data.frame(datExpr)
  datExpr
}

# ============================================================
# Main function
# ============================================================
run_wgcna_from_vst <- function(
    vst_df,
    metadata_df,
    neuronal_genes = NULL,                 # vector of gene symbols to keep (your all_neuronal_genes)
    trait_col = "condition",               # used to build Dummy vars + trait of interest
    trait_level_a = "Shock",               # "SHOCK"
    trait_level_b = "HC",                  # "HC"
    modules_for_cytoscape = c("black", "blue", "brown", "green", "turquoise", "red", "yellow"),
    
    # --- keep your original values as defaults ---
    power_soft = 20,
    powers_test = c(1:10, seq(from = 12, to = 20, by = 2)),
    TOMType = "signed",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = TRUE,
    saveTOMFileBase = "TRAP_invivo_TOM",
    maxBlockSize = 15000,
    cytoscape_threshold = 0.1
) {
  # ---- Input checks ----
  if (is.null(colnames(vst_df))) stop("vst_df must have column names (sample IDs).")
  if (is.null(rownames(vst_df))) stop("vst_df must have row names (gene IDs / symbols).")
  if (is.null(rownames(metadata_df))) {
    stop("metadata_df must have rownames matching sample IDs (same as colnames(vst_df) if genes x samples).")
  }
  
  # ---- Build datExpr (samples x genes) ----
  datExpr <- .as_wgcna_datExpr(vst_df)
  
  # Ensure sample alignment: datExpr rownames should be sample IDs
  # If vst_df was genes x samples, after transpose rownames become sample IDs.
  # If vst_df was samples x genes, rownames must already be sample IDs.
  if (!all(rownames(datExpr) %in% rownames(metadata_df))) {
    stop("Sample IDs in vst_df/datExpr do not match rownames(metadata_df).")
  }
  
  # Reorder metadata to datExpr sample order
  metadata_df <- metadata_df[rownames(datExpr), , drop = FALSE]
  
  # ---- Optional neuronal filtering (your all_neuronal_genes) ----
  if (!is.null(neuronal_genes)) {
    keep <- colnames(datExpr) %in% neuronal_genes
    datExpr <- datExpr[, keep, drop = FALSE]
  }
  
  # ---- Quality check: goodSamplesGenes ----
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
      printFlush(paste("Removing genes:", paste(colnames(datExpr)[!gsg$goodGenes], collapse = ", ")))
    }
    if (sum(!gsg$goodSamples) > 0) {
      printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))
    }
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
    metadata_df <- metadata_df[rownames(datExpr), , drop = FALSE]
  }
  
  # ---- Soft threshold diagnostic (same as your script) ----
  sft <- pickSoftThreshold(datExpr, powerVector = powers_test, verbose = 5)
  
  plot_soft_threshold <- .capture_plot({
    par(mfrow = c(1, 2))
    cex1 <- 0.9
    
    plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
         xlab = "Soft Threshold (power)",
         ylab = "Scale Free Topology Model Fit, signed R^2",
         type = "n", main = "Scale independence"
    )
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
         labels = powers_test, cex = cex1, col = "red"
    )
    abline(h = 0.90, col = "red")
    
    plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
         xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
         type = "n", main = "Mean connectivity"
    )
    text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers_test, cex = cex1, col = "red")
  })
  
  # ---- Create dummy trait matrix (same approach as your script) ----
  out <- binarizeCategoricalColumns(metadata_df,
                                    includePairwise = TRUE,
                                    includeLevelVsAll = FALSE
  )
  Dummy.vars <- data.frame(out, row.names = rownames(metadata_df))
  
  # Plot sample dendrogram + trait heatmap (your QC)
  sampleTree2 <- hclust(dist(datExpr), method = "average")
  traitColors <- numbers2colors(Dummy.vars, signed = TRUE)
  
  plot_sample_dendro_traits <- .capture_plot({
    par(mar = c(20, 6.5, 3, 3))
    plotDendroAndColors(sampleTree2, traitColors,
                        groupLabels = names(Dummy.vars),
                        main = "Sample dendrogram and trait heatmap"
    )
  })
  
  # ---- Build network/modules (same parameters as your script) ----
  cor <- WGCNA::cor
  
  net <- blockwiseModules(datExpr,
                          power = power_soft,
                          TOMType = TOMType,
                          minModuleSize = minModuleSize,
                          reassignThreshold = reassignThreshold,
                          mergeCutHeight = mergeCutHeight,
                          numericLabels = numericLabels,
                          pamRespectsDendro = pamRespectsDendro,
                          saveTOMs = saveTOMs,
                          saveTOMFileBase = saveTOMFileBase,
                          maxBlockSize = maxBlockSize,
                          verbose = 3
  )
  
  moduleLabels <- net$colors
  moduleColors <- labels2colors(net$colors)
  MEs <- net$MEs
  geneTree <- net$dendrograms[[1]]
  
  plot_dendro_modules <- .capture_plot({
    mergedColors <- labels2colors(net$colors)
    plotDendroAndColors(net$dendrograms[[1]],
                        mergedColors[net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05
    )
  })
  
  # ---- Trait of interest: treatment.SHOCK.vs.HC (same idea, but built robustly) ----
  target_pairwise_name <- paste0(trait_col, ".", trait_level_a, ".vs.", trait_level_b)
  if (!target_pairwise_name %in% colnames(Dummy.vars)) {
    stop(paste0(
      "Could not find dummy trait column '", target_pairwise_name, "'.\n",
      "Available Dummy.vars columns include:\n  ",
      paste(head(colnames(Dummy.vars), 30), collapse = ", "),
      if (ncol(Dummy.vars) > 30) " ... (truncated)" else ""
    ))
  }
  
  trait.interest <- as.data.frame(Dummy.vars[[target_pairwise_name]])
  names(trait.interest) <- paste0(trait_level_a, "vs", trait_level_b)
  
  # ---- Gene MM/GS tables (same math as your script) ----
  modNames <- substring(names(MEs), 3)
  nSamples <- nrow(datExpr)
  
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste("MM", modNames, sep = "")
  names(MMPvalue) <- paste("p.MM", modNames, sep = "")
  
  geneTraitSignificance <- as.data.frame(cor(datExpr, trait.interest, use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) <- paste("GS.", names(trait.interest), sep = "")
  names(GSPvalue) <- paste("p.GS.", names(trait.interest), sep = "")
  
  df.geneMMvsGS <- data.frame(
    geneModuleMembership,
    mergedColors = moduleColors,
    geneTraitSignificance,
    row.names = rownames(geneModuleMembership)
  )
  
  # ---- GO enrichment (same anRichment flow as your script) ----
  # convert2entrez requires organism name; your script used "mouse"
  symbol <- rownames(df.geneMMvsGS)
  entrez <- convert2entrez(organism = "mouse", symbol = symbol)
  df.geneMMvsGS$entrez <- entrez
  
  GOcollection <- buildGOcollection(organism = "mouse")
  GOenrichment <- enrichmentAnalysis(
    classLabels = df.geneMMvsGS$mergedColors,
    identifiers = df.geneMMvsGS$entrez,
    refCollection = GOcollection,
    threshold = 0.05,
    thresholdType = "Bonferroni",
    getOverlapEntrez = TRUE,
    getOverlapSymbols = TRUE,
    ignoreLabels = "grey"
  )
  GO_table <- GOenrichment$enrichmentTable
  
  # ---- Eigengene network + module-trait correlation heatmap (same as your script) ----
  MEs2 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  MET <- orderMEs(MEs2)
  
  moduleTraitCor <- cor(MEs2, Dummy.vars, use = "p")
  moduleTraitPValue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  plot_eigengene_network <- .capture_plot({
    par(cex = 0.9)
    plotEigengeneNetworks(MET, "",
                          marDendro = c(0, 4, 1, 2),
                          marHeatmap = c(3, 4, 1, 2),
                          cex.lab = 0.8,
                          xLabelsAngle = 90
    )
  })
  
  plot_module_trait_heatmap <- .capture_plot({
    textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPValue, 1), ")", sep = ""
    )
    dim(textMatrix) <- dim(moduleTraitCor)
    
    par(mar = c(10, 6.5, 3, 3))
    labeledHeatmap(
      Matrix = moduleTraitCor,
      xLabels = names(Dummy.vars),
      yLabels = names(MEs2),
      ySymbols = names(MEs2),
      colorLabels = FALSE,
      colors = greenWhiteRed(50),
      textMatrix = textMatrix,
      setStdMargins = FALSE,
      cex.text = 0.5,
      zlim = c(-1, 1)
    )
  })
  
  # ---- TOM plot (your heatmap) ----
  TOM <- TOMsimilarityFromExpr(datExpr, power = power_soft)
  dissTOM <- 1 - TOM
  plotTOM <- dissTOM^10
  diag(plotTOM) <- NA
  
  plot_tom_heatmap <- .capture_plot({
    TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
  })
  
  # ---- Cytoscape export (same module list + threshold) ----
  probes <- colnames(datExpr)
  
  inModule <- is.finite(match(moduleColors, modules_for_cytoscape))
  modProbes <- probes[inModule]
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(modProbes, modProbes)
  
  cyt <- exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste("CytoscapeInput-edges.txt", sep = ""),
    nodeFile = paste("CytoscapeInput-nodes.txt", sep = ""),
    weighted = TRUE,
    threshold = cytoscape_threshold,
    nodeNames = modProbes,
    nodeAttr = moduleColors[inModule]
  )
  
  # ---- Return a single “analysis bundle” list ----
  list(
    # main outputs
    datExpr = datExpr,
    metadata = metadata_df,
    Dummy.vars = Dummy.vars,
    
    net = net,
    moduleLabels = moduleLabels,
    moduleColors = moduleColors,
    MEs = MEs2,
    geneTree = geneTree,
    
    # tables
    df.geneMMvsGS = df.geneMMvsGS,
    moduleTraitCor = moduleTraitCor,
    moduleTraitPValue = moduleTraitPValue,
    GO_enrichment_table = GO_table,
    
    # cytoscape object (what you need to plot in Cytoscape)
    cytoscape = cyt,
    
    # QC plots stored as objects
    plots = list(
      soft_threshold = plot_soft_threshold,
      sample_dendrogram_traits = plot_sample_dendro_traits,
      dendrogram_modules = plot_dendro_modules,
      eigengene_network = plot_eigengene_network,
      module_trait_heatmap = plot_module_trait_heatmap,
      tom_heatmap = plot_tom_heatmap
    )
  )
}

# ============================================================
# MM vs GS plotting function (takes the required df)
# ============================================================
plot_mm_vs_gs <- function(
    df_geneMMvsGS,
    module_color,
    trait_name = NULL,
    point_cex = 1.2
) {
  stopifnot("mergedColors" %in% colnames(df_geneMMvsGS))
  
  # Pick a trait column automatically if not provided
  if (is.null(trait_name)) {
    trait_candidates <- grep("^GS\\.", colnames(df_geneMMvsGS), value = TRUE)
    if (length(trait_candidates) == 0) stop("No GS.* column found in df_geneMMvsGS.")
    trait_name <- trait_candidates[[1]]
  } else {
    # accept "ShockvsHC" or "GS.ShockvsHC"
    if (!startsWith(trait_name, "GS.")) trait_name <- paste0("GS.", trait_name)
    if (!trait_name %in% colnames(df_geneMMvsGS)) stop(paste0("Trait column not found: ", trait_name))
  }
  
  mm_col <- paste0("MM", module_color)
  if (!mm_col %in% colnames(df_geneMMvsGS)) {
    stop(paste0(
      "Module membership column not found: ", mm_col, "\n",
      "Available MM columns: ", paste(grep("^MM", colnames(df_geneMMvsGS), value = TRUE), collapse = ", ")
    ))
  }
  
  # subset module
  inMod <- df_geneMMvsGS$mergedColors == module_color
  dfMod <- df_geneMMvsGS[inMod, , drop = FALSE]
  
  # Return a plot object (recorded base plot) so it can be stored in a list
  .capture_plot({
    verboseScatterplot(
      abs(dfMod[[mm_col]]),
      abs(dfMod[[trait_name]]),
      xlab = paste("Module Membership in", module_color, "module"),
      ylab = paste("Gene significance for", sub("^GS\\.", "", trait_name)),
      main = paste("Module membership vs. gene significance\n", module_color, "module"),
      cex = point_cex
    )
  })
}


## ---------------------------------------------------------
# MAIN
## ---------------------------------------------------------
set.seed(123)

load("C:\\Users\\mauri\\Dropbox\\TRAP_Alt_Analysis\\WGCNA\\All_counts.RData")

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

trap_cols <- colnames(counts_wide) %>%
  setdiff("Geneid") %>%
  keep(~ str_detect(.x, "_TRAP_")) %>%
  keep(~ !str_detect(.x, "_TL_"))  # extra safety; TL shouldn't match anyway

WGCNA.df <- dplyr::select(counts_wide, contains("TRAP"))

trap_meta <- bind_rows(lapply(trap_cols, parse_sample_info)) %>%
  filter(tech == "TRAP") %>%
  as.data.frame()

rownames(trap_meta) <- trap_meta$sample

#Build the input_file

WGCNA.DEseq <- DESeqDataSetFromMatrix(countData = as.matrix(WGCNA.df),
                                      colData = trap_meta,
                                      design = ~ condition + cell_type)

keep_TRAP <- rowSums(counts(WGCNA.DEseq) >= 20) >= 3
WGCNA.DEseq <- WGCNA.DEseq[keep_TRAP, ]
WGCNA.DEseq <- estimateSizeFactors(WGCNA.DEseq)
WGCNA.DEnorm <- vst(WGCNA.DEseq, blind = TRUE)

#Keep rownames
WGCNA.df.filt <- WGCNA.df[keep_TRAP,]
WGCNA.final <- as.data.frame(assay(WGCNA.DEnorm))
rownames(WGCNA.final) <- rownames(WGCNA.df.filt)


#PCA
plotPCA(WGCNA.DEnorm, intgroup = c("condition", "cell_type")) + 
  geom_point(size = 4) +
  scale_color_manual(values = c("blue", "green", "red", "lightblue", "lightgreen", "pink",
                                "midnightblue", "forestgreen", "darkred")) +
  #geom_text(aes(label = name)) +
  #xlim(-10,12) + 
  #ylim(-8, 10) + 
  theme_bw() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20))

#LOAD THE CELL TYPE-SPECIFIC GENE LISTS, FIND ALL NEURONAL GENES


all_neuronal_genes <- unique(c(rownames(Camk2a.df), rownames(Pvalb.df), rownames(Sst.df)))

# vst_mat: genes x samples (rownames = genes, colnames = samples)
# metadata: rownames = samples
res <- run_wgcna_from_vst(
  vst_df = WGCNA.final,
  metadata_df = trap_meta,
  neuronal_genes = all_neuronal_genes
)

# Cytoscape tables:
edges <- res$cytoscape$edgeData
nodes <- res$cytoscape$nodeData

# Print any stored QC plot:
res$plots$module_trait_heatmap

# MM vs GS plot for (e.g.) turquoise:
p_mm_gs <- plot_mm_vs_gs(res$df.geneMMvsGS, module_color = "turquoise", trait_name = "ShockvsHC")
p_mm_gs
