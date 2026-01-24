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


#Iterate over all files, select count column and rename it. Then append to main data frame
filelist <- list.files("D:\\Dropbox\\TRAP_Alt_Analysis\\FeatCounts", full.names = T)
filelist <- filelist[!grepl("\\.summary$", filelist)]


all_data <- data.frame()
for(file in filelist) {
  
  temp_file <- read.delim(file, comment.char = "#")
  
  temp_file <- temp_file[,c(1,7)]
  
  names(temp_file) <- c("gene_id", str_replace(names(temp_file[2]), "X.scratch.mm10349.files_salmon.", ""))
  names(temp_file) <- c("gene_id", str_replace(names(temp_file[2]), "_S[0-9]+_L[0-9]+Aligned.sortedByCoord.out.bam", ""))
  
  if(length(all_data) == 0) {
    
    all_data <- temp_file
    
  } else {
        
        all_data <- all_data %>%
          left_join(temp_file, by = "gene_id")
        
    }
      
}

#Retrieve gene lengths, then convert to RPKM. Right now we will use all samples to identify purification efficiency
FeatCount_Camk2a15min_TRAP_Box1_S1_L001 <- read.delim("D:/Dropbox/TRAP_Alt_Analysis/FeatCounts/FeatCount_Camk2a15min_TRAP_Box1_S1_L001", comment.char="#")
gene.length <- dplyr::select(FeatCount_Camk2a15min_TRAP_Box1_S1_L001, Geneid, Length)
rownames(all_data) <- all_data$gene_id

all_data_rpkm <- Raw_to_RPKM(all_data)

stopifnot(sum(is.na(all_data_rpkm)) == 0)

hc.rpkm <- dplyr::select(all_data_rpkm, contains("HC"))

# CURRENTLY THERE ARE ONLY TRAP SAMPLES, RETRIEVE THE OTHER ONES FROM BIGPURPLE INTO GREENE, RUN PIPELINE AND COLLECT THEM AS WELL

#Create a data frame with the division between TRAP and TL for each sample in the dataset.
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
