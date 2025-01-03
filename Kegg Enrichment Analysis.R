library(KEGG.db)
library(enrichplot)
library(pathview)
library(clusterProfiler)
library(ggnewscale)
library(ggridges)
library(R.utils)
library(ggplot2)
library(ggupset)
library(forcats)
library(tidyverse)
library(dplyr)

# Set the option for downloading method in clusterProfiler
R.utils::setOption("clusterProfiler.download.method", 'auto')

# Set working directory
setwd("D:/230408 Metagenomic analysis/240311 Metagenome/240322ClusterprofilerKEGG")

# List all CSV files in the directory
FileName = list.files(getwd(), pattern = '*.csv')

# Create output names based on the file names
Outputname = lapply(strsplit(FileName, "_regulated"), "[", 1)

# Iterate over each file for analysis
for (i in 1:13) {
  
  # Read the data from the CSV file
  ko_deseq <- read.csv(FileName[i], header = TRUE)
  
  # Filter data based on adjusted p-value (padj)
  booleanSig <- ko_deseq$padj < 0.05
  length(booleanSig)
  
  # Keep only the significant rows
  geneSig <- ko_deseq[booleanSig,]
  dim(geneSig)
  
  # Further filter by absolute log2 fold change > 2
  booleanSigSize <- abs(geneSig$log2FoldChange) > 2
  length(booleanSigSize)
  
  # Keep rows with both significant p-value and large fold change
  geneSigSize <- geneSig[booleanSigSize,]
  dim(geneSigSize) # Size after filtering
  
  # Extract the log2 fold change values for further analysis
  geneSigSize_fc <- geneSigSize$log2FoldChange
  names(geneSigSize_fc) <- geneSigSize$X  # Assign KO numbers as names to the vector
  
  # Sort the fold change vector in decreasing order
  geneSigSize_fc = sort(geneSigSize_fc, decreasing = TRUE)
  geneSigSize_ko <- geneSigSize$X
  
  # Set the desired organism for KEGG enrichment analysis (KO is used for KEGG pathways)
  organism = "ko"
  
  # Perform KEGG enrichment analysis for down-regulated genes
  kk2_down <- enrichKEGG(
    gene = geneSigSizeKoDown,
    organism = "ko",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,  # New addition for q-value cutoff
    keyType = "kegg"
  )
  
  # Perform KEGG enrichment analysis for up-regulated genes
  kk2_up <- enrichKEGG(
    gene = geneSigSizeKoUp,
    organism = "ko",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,  # New addition for q-value cutoff
    keyType = "kegg"
  )
  
  # Generate dot plots for the results (showing top 20 categories)
  dotplot(kk2_down, showCategory = 20, label = 60)
  dotplot(kk2_up, showCategory = 20, label = 60)
  
  # Save the RData file with results
  save.image(paste0(file = Outputname[[i]], '.RData'))
}
