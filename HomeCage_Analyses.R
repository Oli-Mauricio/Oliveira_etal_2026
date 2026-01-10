library(tidyverse)
library(DESeq2)

#FIGURE 1 - CELL TYPE-SPECIFIC TRANSLATOME#####
#Fig 1E
#-------------------------------------------------------------------------------
##PCA of all genes enriched in neurons
#first get a full list of all genes identified as enriched in all cell types
load("Filtered_gene_cell_Specific.RData")
all_genes <- c(rownames(Filt.Camk2a), rownames(Filt.Pvalb), rownames(Filt.Sst))
all_genes <- unique(all_genes)

load("All_samples_HC_TRAP.RData")
All.HC.samples <- All.HC.samples %>% 
  filter( rownames(All.HC.samples) %in% all_genes)

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

#FIG 1D
##CELL MARKERS PLOT
markers.df <- RPKM.div.df %>%
  filter(rownames(RPKM.div.df) %in% c("Camk2a", "Slc17a7",
                                      "Gad1", "Gad2", "Syt2",
                                      "Pvalb", "Sst", "Aldh1l1",
                                      "Gfap", "Cnp", "Cx3cr1"))
markers.df$geneID <- rownames(markers.df)
gather.df <- gather(markers.df, key = "sample", value = "normalized RPKM", HC1_Camk2a:HC11_Sst)
gather.df$geneID <- factor(gather.df$geneID, levels = c("Camk2a", "Slc17a7",
                                                        "Gad1", "Gad2", "Syt2",
                                                        "Pvalb", "Sst", "Aldh1l1",
                                                        "Gfap", "Cnp", "Cx3cr1"))
gather.df$`Cell type` <- str_replace(gather.df$sample, "^.*_", "")

ggplot(
  gather.df, aes(`Cell type`, log2(`normalized RPKM`), 
                 fill = `Cell type`)) + 
  geom_boxplot(linewidth = 0.8) + 
  geom_point() +
  geom_hline(yintercept = 0, size = 1) +
  scale_fill_manual(values = c("#b1c6d9","#008631","#ffecec")) + 
  ylab(label="log2(RPKM TRAP /\n RPKM Total")+
  scale_y_continuous(breaks = seq(-8, 8, by = 2)) +
  facet_wrap(.~ geneID, nrow = 2) + 
  theme_bw() + 
  theme(axis.title.y = element_text(size = 24),
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

##GENERATION OF DESEQ2 OBJECTS AND DIFFERENTIAL EXPRESSION
#Camk2a vs PV
my_database_TRAP_all <- DESeq(my_database_TRAP_all)
Results_Camk2avsPvalb <- results(my_database_TRAP_all, 
                                 contrast = c("Neuron_type", "Camk2a", "Pvalb"), 
                                 alpha = 0.05, 
                                 cooksCutoff = F, 
                                 independentFiltering = F)

My_results_all_TRAP_annotated.Camk2avsPV <- data.frame(Results_Camk2avsPvalb) 
Camk2a_Pvalb_sig <- subset(My_results_all_TRAP_annotated.Camk2avsPV, padj < 0.05)
Camk2a_Pvalb_sig <- Camk2a_Pvalb_sig %>%
  arrange(padj)

#Camk2a vs Sst
my_database_TRAP_all <- DESeq(my_database_TRAP_all)
Results_Camk2avsSst <- results(my_database_TRAP_all, 
                               contrast = c("Neuron_type", "Camk2a", "Sst"), 
                               alpha = 0.05, 
                               cooksCutoff = F, 
                               independentFiltering = F)

My_results_all_TRAP_annotated.Camk2avsSst <- data.frame(Results_Camk2avsSst) 
Camk2a_Sst_sig <- subset(My_results_all_TRAP_annotated.Camk2avsSst, padj < 0.05)
Camk2a_Sst_sig <- Camk2a_Sst_sig %>%
  arrange(padj)

#PV vs Sst
my_database_TRAP_all <- DESeq(my_database_TRAP_all)
Results_PvalbvsSst <- results(my_database_TRAP_all, 
                              contrast = c("Neuron_type", "Pvalb", "Sst"), 
                              alpha = 0.05, 
                              cooksCutoff = F, 
                              independentFiltering = F)

My_results_all_TRAP_annotated.PV.Sst <- data.frame(Results_PvalbvsSst) 
Pvalb_Sst_sig <- subset(My_results_all_TRAP_annotated.PV.Sst, padj < 0.05)
Pvalb_Sst_sig <- Pvalb_Sst_sig %>%
  arrange(padj)




#Fig 1F
#-------------------------------------------------------------------------------
#

load("Filtered_gene_cell_Specific.RData")

df <- list(Camk2a = rownames(Filt.Camk2a),
           Pvalb = rownames(Filt.Pvalb),
           Sst = rownames(Filt.Sst))


attributes(df) <- list(names = names(df), 
                       row.names = 1:6724, 
                       class = 'data.frame')

#Venn diagram
ggVennDiagram(df, label_alpha = 0, 
              set_color = "midnightblue", 
              label_size = 8,
              label_percent_digit = 1,
              set_size = 8) + 
  scale_fill_gradient(low = "#F4FAFE", 
                      high = "#4981BF")

#Fig 1G
#-------------------------------------------------------------------------------
#
#Find common genes for DEG analysis
common_genes <- Reduce(intersect, list(rownames(Filt.Camk2a),
                                       rownames(Filt.Pvalb),
                                       rownames(Filt.Sst)))

common.Camk2a.Pvalb <- My_results_all_TRAP_annotated.Camk2avsPV %>% 
  filter(rownames(My_results_all_TRAP_annotated.Camk2avsPV) %in% common_genes)

common.Camk2a.Sst <- My_results_all_TRAP_annotated.Camk2avsSst %>% 
  filter(rownames(My_results_all_TRAP_annotated.Camk2avsSst) %in% common_genes)

common.Pvalb.Sst <- My_results_all_TRAP_annotated.PV.Sst %>% 
  filter(rownames(My_results_all_TRAP_annotated.PV.Sst) %in% common_genes)


df.final <- data.frame("Camk2avsPV" = common.Camk2a.Pvalb$log2FoldChange, 
                       "Camk2avsSst" = common.Camk2a.Sst$log2FoldChange, 
                       row.names = rownames(common.Camk2a.Pvalb))

df.final <- df.final %>% 
  mutate(sig_PvalbvsSst = case_when(rownames(df.final) %in% rownames(Pvalb_Sst_sig) ~ "yes", 
                                    TRUE ~ "no"))



ggplot(df.final, aes(Camk2avsPV, Camk2avsSst)) + 
  geom_point(data = df.final, aes(fill = sig_PvalbvsSst), size = 3, shape = 21) + 
  scale_fill_manual(values = c("yes" = "darkblue", "no" = "gray")) +
  geom_smooth(method = 'lm', colour = "black", size = 1, alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  xlim(c(-15,15)) + ylim (c(-15,15)) +
  xlab(label = "log2FC - Camk2a/Pvalb") + ylab(label = "log2FC - Camk2a/Sst") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.position = "bottom") +
  stat_cor(method = "spearman", size = 6, label.x = -15, label.y = 12, cor.coef.name = "rho") +
  geom_label_repel(aes(label = ifelse(row.names(df.final) %in% c("Sst", "Pvalb"), row.names(df.final), '')),
                   max.overlaps = Inf, box.padding = 1
  )