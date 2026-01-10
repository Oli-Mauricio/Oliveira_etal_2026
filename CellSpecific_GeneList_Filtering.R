#

# This is the filtering routine to obtain a list of common genes identified in all sequencing experiments (home cage).
# It converts raw FeatureCount files into RPKM values, then  

# Load all data sets, converts raw values into RPKM, then produces one data frame of TRAP / Total Lysate relative vals
# Filters samples that were not enriched to signature genes Camk2a, Pvalb and Sst
# Filters dataset to contain only genes with 1.2-fold enrichment (~20%) over total lysate as cell type-specific

### NOTE:
# The high incidence of samples with low purity in interneurons is expected, due to low input.

#

library(edgeR)
library(tidyverse)

Raw_to_RPKM <- function(df_interest) {
  
  df_interest <- as.data.frame(df_interest)
  
  gene.length.temp <- gene.length %>% 
    filter(Geneid %in% rownames(df_interest))
  
  df_interest <- df_interest %>% 
    filter(rownames(df_interest) %in% gene.length.temp$Geneid)
  
  match_idx <- match(rownames(df_interest), gene.length.temp$Geneid)
  df_interest <- df_interest[match_idx,]
  
  
  dge.list <- DGEList(df_interest)
  
  rpkm.values <- as.data.frame(rpkm(dge.list, gene.length = gene.length.temp$Length))
  
  return(rpkm.values)
  
}

#
# INSERT HERE ROUTINE TO CONVERT ALL SAMPLES INTO RPKM
#
#

RPKM.Camk2a_60min <- RPKM.Camk2a_60min %>%
  filter(rownames(RPKM.Camk2a_60min) %in% rownames(RPKM.Camk2a_15min)) %>%
  mutate(geneID = rownames(RPKM.Camk2a_60min))
RPKM.Pvalb_15min <- RPKM.Pvalb_15min %>%
  filter(rownames(RPKM.Pvalb_15min) %in% rownames(RPKM.Camk2a_15min)) %>%
  mutate(geneID = rownames(RPKM.Pvalb_15min))
RPKM.Pvalb_60min <- RPKM.Pvalb_60min %>%
  filter(rownames(RPKM.Pvalb_60min) %in% rownames(RPKM.Camk2a_15min)) %>%
  mutate(geneID = rownames(RPKM.Pvalb_60min))
RPKM.Sst_15min <- RPKM.Sst_15min %>%
  filter(rownames(RPKM.Sst_15min) %in% rownames(RPKM.Camk2a_15min)) %>%
  mutate(geneID = rownames(RPKM.Sst_15min))
RPKM.Sst_60min <- RPKM.Sst_60min %>%
  filter(rownames(RPKM.Sst_60min) %in% rownames(RPKM.Camk2a_60min)) %>%
  mutate(geneID = rownames(RPKM.Sst_60min)) 

RPKM.Camk2a_15min$geneID <- rownames(RPKM.Camk2a_15min)


RPKM.df <- RPKM.Camk2a_15min %>%
  inner_join(RPKM.Camk2a_60min, by = "geneID") %>%
  inner_join(RPKM.Pvalb_15min, by = "geneID") %>%
  inner_join(RPKM.Pvalb_60min, by = "geneID") %>%
  inner_join(RPKM.Sst_15min, by = "geneID") %>%
  inner_join(RPKM.Sst_60min, by = "geneID")

rownames(RPKM.df) <- RPKM.df$geneID
RPKM.df <- dplyr::select(RPKM.df, -geneID)

RPKM.div.df <- data.frame(
  "HC1_Camk2a" = RPKM.df$HC_TRAP_Female2_C2 / RPKM.df$HC_TotLys_Female2_C2,
  "HC2_Camk2a" = RPKM.df$HC_TRAP_Male1_C2 / RPKM.df$HC_TotLys_Male1_C2,
  "HC3_Camk2a" = RPKM.df$HC_TRAP_Male2_C2 / RPKM.df$HC_TotLys_Male2_C2,
  "HC4_Camk2a" = RPKM.df$FeatCount_HC2_TRAP_Camk2a1h_S13_L008 / RPKM.df$FeatCount_HC2_TL_Camk2a1h_S8_L008,
  "HC5_Camk2a" = RPKM.df$FeatCount_HC3_TRAP_Camk2a1h_S14_L008 / RPKM.df$FeatCount_HC3_TL_Camk2a1h_S9_L008,
  "HC6_Camk2a" = RPKM.df$FeatCount_HC4_TRAP_Camk2a1h_S18_L008 / RPKM.df$FeatCount_HC4_TL_Camk2a1h_S27_L008,
  "HC7_Camk2a" = RPKM.df$FeatCount_HC5_TRAP_Camk2a1h_S19_L008 / RPKM.df$FeatCount_HC5_TL_Camk2a1h_S28_L008,
  "HC8_Camk2a" = RPKM.df$FeatCount_HC6_TRAP_Camk2a1h_S20_L008 / RPKM.df$FeatCount_HC6_TL_Camk2a1h_S29_L008,
  "HC1_Pvalb" = RPKM.df$HC_5_TRAP / RPKM.df$HC_5_TOTLYS,
  "HC2_Pvalb" = RPKM.df$HC1_TRAP / RPKM.df$HC1_TOTLYS,
  "HC3_Pvalb" = RPKM.df$HC3_TRAP / RPKM.df$HC3_TOTLYS,
  "HC4_Pvalb" = RPKM.df$HC4_TRAP / RPKM.df$HC4_TOTLYS,
  "HC5_Pvalb" = RPKM.df$HC1_TRAP_1h / RPKM.df$HC1_TOTLYS_1h,
  "HC6_Pvalb" = RPKM.df$HC2_TRAP_1h / RPKM.df$HC2_TOTLYS_1h,
  "HC7_Pvalb" = RPKM.df$HC3_TRAP_1h / RPKM.df$HC3_TOTLYS_1h,
  "HC8_Pvalb" = RPKM.df$HC4_TRAP_1h / RPKM.df$HC4_TOTLYS_1h,
  "HC9_Pvalb" = RPKM.df$HC5_TRAP_1h / RPKM.df$HC5_TOTLYS_1h,
  "HC10_Pvalb" = RPKM.df$HC6_TRAP_1h / RPKM.df$HC6_TOTLYS_1h,
  "HC11_Pvalb" = RPKM.df$HC7_TRAP_1h / RPKM.df$HC7_TOTLYS_1h,
  "HC12_Pvalb" = RPKM.df$HC8_TRAP_1h / RPKM.df$HC8_TOTLYS_1h,
  "HC1_Sst" = RPKM.df$FeatCount_HC1_TRAP_SST_S33_L003 / RPKM.df$FeatCount_HC1_TOTLYS_SST_S53_L003,
  "HC2_Sst" = RPKM.df$FeatCount_HC2_TRAP_SST_S34_L003 / RPKM.df$FeatCount_HC2_TOTLYS_SST_S54_L003,
  "HC3_Sst" = RPKM.df$FeatCount_HC3_TRAP_SST_S35_L003 / RPKM.df$FeatCount_HC3_TOTLYS_SST_S55_L003,
  "HC4_Sst" = RPKM.df$FeatCount_HC4_TRAP_SST_S36_L003 / RPKM.df$FeatCount_HC4_TOTLYS_SST_S56_L003,
  "HC5_Sst" = RPKM.df$FeatCount_HC5_TRAP_SST_S37_L003 / RPKM.df$FeatCount_HC5_TOTLYS_SST_S57_L003,
  "HC6_Sst" = RPKM.df$FeatCount_HC6_TRAP_SST_S38_L003 / RPKM.df$FeatCount_HC6_TOTLYS_SST_S58_L003,
  "HC7_Sst" = RPKM.df$FeatCount_HC7_TRAP_SST_S39_L003 / RPKM.df$FeatCount_HC7_TOTLYS_SST_S59_L003,
  "HC8_Sst" = RPKM.df$`FeatCount_TRAP-HC2C1_S16_L006` / RPKM.df$`FeatCount_TL-HC2C1_S2_L006`,
  "HC9_Sst" = RPKM.df$`FeatCount_TRAP-HC3C2_S23_L006` / RPKM.df$`FeatCount_TL-HC3C2_S9_L006`,
  "HC10_Sst" = RPKM.df$HC1C3_TRAP / RPKM.df$HC1C3_TL,
  "HC11_Sst" = RPKM.df$HC2C3_TRAP / RPKM.df$HC2C3_TL,
  row.names = rownames(RPKM.df)
)

# Repairs raw-to-rpkm conversion errors (e.g., rawCounts=0 in total lysate - the denominator - generates Inf)
RPKM.div.df[is.na(RPKM.div.df)] <- 0
RPKM.div.df[RPKM.div.df == "Inf"] <- 0
RPKM.div.df[RPKM.div.df == "-Inf"] <- 0

#
## Filter out non-enriched samples
gene.markers <- c("Camk2a" ,"Pvalb", "Sst")

#
# INSERT HERE ROUTINE THAT FINDS A COLUMN THAT HAS NO ENRICHMENT AND ELIMINATES IT.
#
#

RPKM.div.df <- RPKM.div.df[,-c(3,22,25,26,29)] #HC2,5, 6 and 9 in Sst and HC3 in Camk2a

Filt.Camk2a <- RPKM.div.df %>% 
  filter(rowSums(RPKM.div.df[,c(1:7)] >= 1.2) >= 3) # Camk2a samples

Filt.Pvalb <- RPKM.div.df %>%
  filter(rowSums(RPKM.div.df[,c(8:19)] >= 1.2) >= 6) # Pvalb samples

Filt.Sst <- RPKM.div.df %>%
  filter(rowSums(RPKM.div.df[,c(20:27)] >= 1.2) >= 4) # Sst samples

save(Filt.Camk2a, Filt.Pvalb, Filt.Sst, file = "Filtered_gene_cell_Specific.RData")
