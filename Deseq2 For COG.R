# Install and load required libraries
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("DESeq2")

library(DESeq2)
library(SQMtools)
library(readxl)

# Set working directory
setwd("D:/230408 Metagenomic analysis/240311 Metagenome/240311Deseq2/cogInoculum/")

# Load SQM data and metadata
wormPlastic <- loadSQM("D:/230408 Metagenomic analysis/wormPlastic2/")
metadata <- read_excel("D:/230408 Metagenomic analysis/240311 Metagenome/metaForMetagenomeUpdated.xlsx")

# Load functional abundance data
cogMetaseq <- wormPlastic[["functions"]][["COG"]][["abund"]]
cogMetaseq <- t(cogMetaseq)

# Load comparison list
ComparedList <- read.csv("D:/230408 Metagenomic analysis/240311 Metagenome/240311Deseq2/cogInoculum/ComparedList/ComparedList.csv", header = TRUE)

# Loop over comparisons
for (i in 1:3) {
  
  # Merge functional abundance data with metadata
  cogMetaseq0 <- merge(cogMetaseq, metadata, by.x = "row.names", by.y = "SampleName")
  cogMetaseq0 <- subset(cogMetaseq0, (Position == ComparedList$a1[i] & (Type == ComparedList$a2[i] | Type == ComparedList$a3[i])))
  Type <- ComparedList$a4[i]

  # Select relevant columns for analysis
  a <- c(1, 33537:33542)
  group <- cogMetaseq0[, a]
  cogMetaseq0 <- cogMetaseq0[, 1:33536]
  row.names(cogMetaseq0) <- cogMetaseq0[, 1]
  cogMetaseq0 <- cogMetaseq0[, -1]
  cogMetaseq01 <- t(cogMetaseq0)

  # Define output prefix
  prefix <- paste0(Type, '_', "Shiftafterfeeding", '_')

  # Prepare DESeq2 dataset
  group$Stage <- factor(group$Position)
  FullCountTable <- DESeqDataSetFromMatrix(
    countData = cogMetaseq01,
    colData = group,
    design = ~Type
  )

  # Filter low-count data
  deseq2Data <- FullCountTable[rowSums(counts(FullCountTable)) > 100, ]
  ds_filtered <- DESeq(deseq2Data)

  # Perform differential analysis
  re_filtered <- results(ds_filtered)
  summary(re_filtered)
  
  # Diagnostic plots
  plotMA(re_filtered, ylim = c(-2, 2))
  plotDispEsts(ds_filtered, ylim = c(1e-6, 1e1))
  hist(re_filtered$pvalue, breaks = 20, col = "grey")

  # Save all results to CSV
  csvname1 <- paste0(prefix, "re_filtered_deseq.csv")
  write.csv(as.data.frame(re_filtered), csvname1)

  # Volcano plot
  with(re_filtered, plot(log2FoldChange, -log10(pvalue), pch = 20, main = "Volcano plot", xlim = c(-8, 8)))
  with(subset(re_filtered, padj < 0.01), points(log2FoldChange, -log10(pvalue), pch = 20, col = "blue"))
  with(subset(re_filtered, padj < 0.01 & abs(log2FoldChange) > .58), points(log2FoldChange, -log10(pvalue), pch = 20, col = "red"))
  with(subset(re_filtered, padj < 0.01 & abs(log2FoldChange) > 1), points(log2FoldChange, -log10(pvalue), pch = 20, col = "green"))

  # Data transformation and filtering
  vsdata <- varianceStabilizingTransformation(ds_filtered, blind = FALSE)
  reSig <- re_filtered[which(re_filtered$padj < 0.05), ]
  reSigup <- reSig[which(reSig$log2FoldChange > 1), ]
  reSigdown <- reSig[which(reSig$log2FoldChange < -1), ]

  # Save filtered results to CSV
  csvname2 <- paste0(prefix, "Upregulated.csv")
  write.csv(as.data.frame(reSigup), csvname2)
  csvname3 <- paste0(prefix, "Downregulated.csv")
  write.csv(as.data.frame(reSigdown), csvname3)
  csvname4 <- paste0(prefix, "Regulated.csv")
  write.csv(as.data.frame(reSig), csvname4)
  csvname5 <- paste0(prefix, "gene_counts.csv")
  write.csv(as.data.frame(cogMetaseq01), csvname5)
  csvname6 <- paste0(prefix, "metadata.csv")
  write.csv(as.data.frame(group), csvname6)
}
