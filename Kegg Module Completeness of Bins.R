# Load necessary libraries
library(reshape2)
library(tidyr)
library(SQMtools)
library(tidyverse)
library(ComplexHeatmap)
library(readxl)
library(eoffice)
library(circlize)
library(pheatmap)
library(dendextend)

# Load the SQM object for plastic degradation analysis
wormPlastic <- loadSQM('D:/230408 Metagenomic analysis/wormPlastic2/', tax_mode = "nofilter", trusted_functions_only = FALSE, engine = 'data.table')

# Load the list of abundant bins
BinAb <- read.csv("D:/230408 Metagenomic analysis/240311 Metagenome/240423 binning analysis/AbundantBins.csv", header = TRUE, row.names = 1)[,1]

# Simplify bin names in the list
temp1 <- sapply(BinAb, function(x) strsplit(x, "\\.")[[1]][1])
temp2 <- sapply(BinAb, function(x) strsplit(x, "\\.")[[1]][2])           
BinAb <- paste0('Bin.', temp1, '.', temp2)
BinAb <- gsub('concoct', 'c', BinAb)
BinAb <- gsub('metabat2', 'm', BinAb)

# Load contig-bin mapping and simplify bin names
ContigTb <- as.data.frame(wormPlastic[["contigs"]][["table"]])
colnames(ContigTb)[5] <- "BinID"
temp1 <- sapply(ContigTb[, 5], function(x) strsplit(x, "\\.")[[1]][1])
temp2 <- sapply(ContigTb[, 5], function(x) strsplit(x, "\\.")[[1]][2])           
ContigTb[, 5] <- paste0('Bin.', temp1, '.', temp2)
ContigTb[, 5] <- gsub('concoct', 'c', ContigTb[, 5])
ContigTb[, 5] <- gsub('metabat2', 'm', ContigTb[, 5])

# Filter for high-abundance bins and extract corresponding contig list
ContigTbAb <- subset(ContigTb, BinID %in% BinAb)
ContigList <- rownames(ContigTbAb)

# Load ORF table and filter out rows with no KEGG ID
OrfTb <- as.data.frame(wormPlastic[["orfs"]][["table"]])
colnames(OrfTb)[8] <- "KEGGID"
colnames(OrfTb)[1] <- "ContigID"
OrfTb <- subset(OrfTb, KEGGID != "")
OrfTb <- subset(OrfTb, ContigID %in% ContigList)

# Map contig IDs to bin IDs
ContigCol <- OrfTb$ContigID
BinCol <- sapply(ContigCol, function(x) ContigTbAb[x, 5])
OrfTb <- cbind(OrfTb, BinCol)

# Extract KEGG IDs associated with each bin
KeggBin <- select(OrfTb, c('BinCol', "KEGGID"))
KeggList <- sapply(unique(BinCol), function(BinTarget) {
  KEGGTemp <- subset(KeggBin, BinCol == BinTarget)
  paste0(KEGGTemp$KEGGID, collapse = ' ')
})

# Prepare result table with bin and KEGG IDs
Result <- as.data.frame(cbind(unique(BinCol), KeggList))
Result$KeggList <- gsub('\\*', '', Result$KeggList)
Result$KeggList <- gsub(';', ' ', Result$KeggList)
Result$KeggList <- gsub(' ', ',', Result$KeggList)
write.csv(Result, 'KEGGBIN.csv', quote = FALSE)

# Save KEGG list for each bin as individual files
for (i in 1:nrow(Result)) {
  write.table(Result[i, 2], Result[i, 1], quote = FALSE, row.names = FALSE, col.names = FALSE)  
}

# Export the total KEGG list
KeggBin$KEGGID <- gsub('\\*', '', KeggBin$KEGGID)
a <- paste0(KeggBin$KEGGID, collapse = ',')
write.csv(a, 'TotalKEGG.csv')

# Integrate the percentage of different bins
Total <- read.csv("D:/230408 Metagenomic analysis/240311 Metagenome/240514KeggCompleteness/result/Total.summary.kegg_pathways.tsv", sep = "\t")[,2:3]
Total <- column_to_rownames(Total, 'pathway_name')
setwd("D:/230408 Metagenomic analysis/240311 Metagenome/240514KeggCompleteness/result/result_selected")
file_names <- list.files()

for (i in 1:length(file_names)) {
  temp <- read.csv(file_names[i], sep = "\t")[,2:3]
  temp <- column_to_rownames(temp, 'pathway_name')
  Total <- merge(Total, temp, by.x = 0, by.y = 0, all.x = TRUE)
  Total <- column_to_rownames(Total, 'Row.names')
}

# Rename columns based on file names
setwd('D:/Liu Yi-Nan/Desktop/')
temp1 <- sapply(file_names, function(x) strsplit(x, "\\.")[[1]][1])
temp2 <- sapply(file_names, function(x) strsplit(x, "\\.")[[1]][2])
temp3 <- sapply(file_names, function(x) strsplit(x, "\\.")[[1]][3])
BinName <- paste0(temp1, '.', temp2, '.', temp3)
Total <- Total[, -1]
colnames(Total) <- BinName
Total <- Total %>% replace(is.na(.), 0)

# Load pathway annotations and filter by selected bins
KeggPreList <- read_excel("D:/Liu Yi-Nan/Desktop/240731SelectedKEGGBin.xlsx", col_names = FALSE)
KeggPreList <- KeggPreList[[1]]
TotalSig <- Total[KeggPreList,]
TotalSig[TotalSig < 90] <- 0
TotalSig1 <- TotalSig
TotalSig1 <- TotalSig1[which(rowSums(TotalSig1) != 0),]
TotalSig2 <- TotalSig1
TotalSig2[TotalSig2 > 0] <- 100

# Generate heatmap of pathway completeness
Heatmap(TotalSig2)

# Transpose the data for heatmap visualization
TotalSig2 <- t(TotalSig2)
col_fun3 <- colorRamp2(c(0, 100), c("#def1f0", "#035a8d"))
AnnoPathway <- read.csv('D:/230408 Metagenomic analysis/240311 Metagenome/240514KeggCompleteness/result/240523TotalCompleteness.csv', header = FALSE)
Anno <- sapply(AnnoPathway$V4, function(x) strsplit(x, "; ")[[1]][2])
names(Anno) <- AnnoPathway$V3

# Read bin abundance data
BinAbList <- rev(read.csv("D:/230408 Metagenomic analysis/240311 Metagenome/240526BinAbList.csv")$x)
BinAbPhylum <- rev(read.csv("D:/230408 Metagenomic analysis/240311 Metagenome/240526BinAbOrder.csv")$x)
TotalSig2 <- TotalSig2[BinAbList,]

# Set pathway types for annotations
PathwayType <- Anno[colnames(TotalSig2)]
PathwayType <- gsub(' ', '_', PathwayType)
PathwayType <- sort(PathwayType)
TotalSig2 <- TotalSig2[, names(PathwayType)]

# Annotate heatmap
TopAnno <- columnAnnotation(
  gap = unit(0.05, "cm"),
  height = unit(1, "cm"),
  annotation_height = unit(0.1, "cm"),
  Pathway_type = PathwayType,
  col = list(Pathway_type = c(Lipid_metabolism = "#7fc97f", Xenobiotics_biodegradation = "#beaed4",
                              Amino_acid_metabolism = "#fdc086", Carbohydrate_metabolism = "#f4cae4",
                              Nucleotide_metabolism = "#ffff99", Energy_metabolism = "#386cb0",
                              Metabolism_of_cofactors_and_vitamins = "#f0027f", Biosynthesis_of_terpenoids_and_polyketides = '#c15f16',
                              Glycan_metabolism = '#666666'))
)

BinAbPhylum <- gsub('assified ', '.', BinAbPhylum)
RightAnno2 <- rowAnnotation(
  width = unit(1.5, "cm"),
  annotation_height = unit(0.1, "cm"),
  Order = BinAbPhylum,
  col = list(Order = c(Bacteroidales = '#a6cee3', Burkholderiales = '#1f78b4', Corynebacteriales = '#b2df8a', Enterobacterales = '#33a02c', 
                       Flavobacteriales = '#fb9a99', Fusobacteriales = '#e31a1c', Hyphomicrobiales = '#fdbf6f', Lactobacillales = '#ff7f00', 
                       Moraxellales = '#cab2d6', Pseudomonadales = '#6a3d9a', Uncl.Tenericutes = '#ffff99', Xanthomonadales = '#b15928'))
)

# Load necessary libraries
library(ComplexHeatmap)
library(ggplot2)
library(readxl)
library(circlize)

# Load pathway abundance data
Pazy <- read.csv('D:/230408 Metagenomic analysis/240311 Metagenome/240717 代谢通路分析绘图与网络分析/AbPazyNum.csv', row.names = 1)
Pazy <- Pazy[rownames(TotalSig2),]

# Hierarchical clustering of pathways based on Euclidean distance
dend <- as.dendrogram(hclust(dist(t(Total)[BinAbList,], method = 'euclidean')), method = "average")

# Pathway names for annotation
PathwayName <- colnames(TotalSig2)
PathwayName1 <- PathwayName
column_names <- anno_text(PathwayName1, rot = 90, location = unit(1, "npc"), just = "right")

# Create the heatmap for pathway completeness and enzyme numbers
a <- Heatmap(TotalSig2,
             column_labels = PathwayName1,
             border = 'black',
             column_gap = unit(2, "mm"),
             row_gap = unit(2, "mm"),
             name = 'Pathway completeness (%)',
             rect_gp = gpar(col = "white", lty = 1.5, lwd = 4),
             cluster_columns = FALSE,
             col = col_fun3,
             row_names_side = 'right',
             cluster_rows = dend,
             top_annotation = TopAnno,
             right_annotation = RightAnno2,
             column_title = NULL,
             row_title = NULL,
             row_dend_gp = gpar(lwd = 2),
             heatmap_legend_param = list(title_position = "lefttop-rot"),
             row_dend_width = unit(1.5, "cm")) +
  Heatmap(Pazy,
          border = 'black',
          column_gap = unit(4, "mm"),
          name = 'Enzyme number',
          cluster_columns = FALSE,
          rect_gp = gpar(col = "white", lty = 1.5, lwd = 2),
          col = col_fun2,
          column_names_rot = 90,
          column_title = NULL,
          row_title = NULL,
          row_dend_gp = gpar(lwd = 2),
          heatmap_legend_param = list(title_position = "lefttop-rot"),
          row_dend_width = unit(1.5, "cm"))

# Save the heatmap as a PowerPoint file
topptx(a, '240802KEGGAbCompleteness1.pptx', width = 18, height = 8, append = TRUE)

# Save TotalSig2 data as a CSV file
write.csv(TotalSig2, '240802KEGGBinSelectByHand.csv')

# Calculate statistics for pathway completeness (mean, max, min, SD)
TotalMean <- apply(Total, 1, mean)
TotalMax <- apply(Total, 1, max)
NMax <- which(TotalMax > 70)
TotalMin <- apply(Total, 1, min)
NMin <- which(TotalMin < 70)
TotalSd <- apply(Total, 1, sd)
NMean <- which((TotalMean > 4.5) & (TotalMean < 46.9))
NSd <- which(TotalSd > 21.2)
NMeanSd <- intersect(NMean, NSd)
NMeanSdMax <- intersect(NMeanSd, NMax)
NMeanSdMaxMin <- intersect(NMeanSdMax, NMin)
TotalMeanSdMaxMin <- Total[NMeanSdMaxMin, ]
TTotalMeanSdMaxMin <- t(TotalMeanSdMaxMin)

# Define color functions for the heatmap
col_fun3 <- colorRamp2(c(0, 100), c("#def1f0", "#035a8d"))
col_fun2 <- colorRamp2(c(0, 8), c("white", "#83062c"))
col_fun1 <- colorRamp2(c(0, 100), c("white", "black"))

# Load and process pathway annotation data
AnnoPathway <- read.csv('D:/230408 Metagenomic analysis/240311 Metagenome/240514KeggCompleteness/result/240523TotalCompleteness.csv', header = FALSE)
Anno <- sapply(AnnoPathway$V4, function(x) strsplit(x, "; ")[[1]][2])
names(Anno) <- AnnoPathway$V3

# Load bin abundance and phylum data
BinAbList <- rev(read.csv("D:/230408 Metagenomic analysis/240311 Metagenome/240526BinAbList.csv")$x)
BinAbPhylum <- rev(read.csv("D:/230408 Metagenomic analysis/240311 Metagenome/240526BinAbOrder.csv")$x)

# Filter TTotalMeanSdMaxMin based on bin list
TTotalMeanSdMaxMin <- TTotalMeanSdMaxMin[BinAbList, ]

# Map pathway types based on annotations
PathwayType <- Anno[colnames(TTotalMeanSdMaxMin)]
PathwayType <- gsub(' ', '_', PathwayType)

# Remove "Gene_set" pathway type and save selected pathways
a <- which(PathwayType == "Gene_set")
TTotalMeanSdMaxMin <- TTotalMeanSdMaxMin[, -a]
PathwayType <- Anno[colnames(TTotalMeanSdMaxMin)]
PathwayType <- gsub(' ', '_', PathwayType)
write.csv(colnames(TTotalMeanSdMaxMin), '240715KeggPathwaySelected.csv')

# Create top annotation for heatmap
TopAnno <- columnAnnotation(
  gap = unit(0.05, "cm"),
  height = unit(1, "cm"),
  annotation_height = unit(0.1, "cm"),
  Pathway_type = PathwayType,
  col = list(Pathway_type = c(
    Lipid_metabolism = "#7fc97f", Xenobiotics_biodegradation = "#beaed4",
    Amino_acid_metabolism = "#fdc086", Carbohydrate_metabolism = "#f4cae4",
    Nucleotide_metabolism = "#ffff99", Energy_metabolism = "#386cb0",
    Metabolism_of_cofactors_and_vitamins = "#f0027f", Biosynthesis_of_terpenoids_and_polyketides = '#c15f16',
    Glycan_metabolism = '#666666'
  ))
)

# Clean up phylum names and create right annotation for heatmap
BinAbPhylum <- gsub('assified ', '.', BinAbPhylum)
RightAnno2 <- rowAnnotation(
  width = unit(1.5, "cm"),
  annotation_height = unit(0.1, "cm"),
  Order = BinAbPhylum,
  col = list(Order = c(
    Bacteroidales = '#a6cee3', Burkholderiales = '#1f78b4', Corynebacteriales = '#b2df8a', 
    Enterobacterales = '#33a02c', Flavobacteriales = '#fb9a99', Fusobacteriales = '#e31a1c', 
    Hyphomicrobiales = '#fdbf6f', Lactobacillales = '#ff7f00', Moraxellales = '#cab2d6', 
    Pseudomonadales = '#6a3d9a', Uncl.Tenericutes = '#ffff99', Xanthomonadales = '#b15928'
  ))
)

# Create dendrogram and heatmap for filtered data
dend <- as.dendrogram(hclust(dist(TTotalMeanSdMaxMin, method = 'euclidean')), method = "average")
PathwayName <- colnames(TTotalMeanSdMaxMin)
PathwayName1 <- sapply(PathwayName, function(x) strsplit(x, ",")[[1]][1])

a <- Heatmap(TTotalMeanSdMaxMin,
             column_labels = PathwayName1,
             border = 'black',
             column_gap = unit(2, "mm"),
             row_gap = unit(2, "mm"),
             name = 'Pathway completeness (%)',
             rect_gp = gpar(col = "white", lty = 1.5, lwd = 2),
             cluster_columns = TRUE,
             column_km = 7,
             col = col_fun3,
             row_names_side = 'right',
             cluster_rows = dend,
             top_annotation = TopAnno,
             right_annotation = RightAnno2,
             column_title = NULL,
             row_title = NULL,
             row_split = 5,
             row_dend_gp = gpar(lwd = 2),
             heatmap_legend_param = list(title_position = "lefttop-rot"),
             row_dend_width = unit(1.5, "cm")) +
  Heatmap(Pazy,
          border = 'black',
          column_gap = unit(4, "mm"),
          name = 'Enzyme number',
          cluster_columns = FALSE,
          rect_gp = gpar(col = "white", lty = 1.5, lwd = 2),
 

# Load necessary data
KeggPreList = read_excel("D:\\230408 Metagenomic analysis\\240311 Metagenome\\KeggFilListPre.xlsx", col_names = FALSE)
KeggPreList = as.vector(KeggPreList)[[1]]
TotalFil = Total[KeggPreList, ]
write.csv(TotalFil, 'TotalPercentagePathwayBins.csv')

# Load pathway annotation data
AnnoPathway = read.csv('D:\\230408 Metagenomic analysis\\240311 Metagenome\\240514KeggCompleteness\\result\\240523TotalCompleteness.csv', header = FALSE)
Anno = sapply(AnnoPathway$V4, function(x) strsplit(x, "; ")[[1]][2])
names(Anno) = AnnoPathway$V3

# Prepare bin annotation data
BinAbList = rev(read.csv("D:\\230408 Metagenomic analysis\\240311 Metagenome\\240526BinAbList.csv")$x)
BinAbPhylum = rev(read.csv("D:\\230408 Metagenomic analysis\\240311 Metagenome\\240526BinAbOrder.csv")$x)

# Subset data by bin list and transpose
TotalFil = t(TotalFil)[BinAbList, ]
PathwayType = Anno[colnames(TotalFil)]
PathwayType = gsub(' ', '_', PathwayType)

# Arrange data by pathway type
TotalFil1 = as.data.frame(t(TotalFil))
TotalFil1$type = PathwayType
TotalFil1 = arrange(TotalFil1, type)

# Remove the 33rd column and transpose again
TotalFil1 = TotalFil1[, -33]
TotalFil2 = t(TotalFil1)

# Define pathway annotation
TopAnno = columnAnnotation(
  gap = unit(0.05, "cm"),
  height = unit(1, "cm"),
  annotation_height = unit(0.1, "cm"),
  Pathway_type = PathwayType,
  col = list(
    Pathway_type = c(
      Lipid_metabolism = "#7fc97f",
      Xenobiotics_biodegradation = "#beaed4",
      Amino_acid_metabolism = "#fdc086",
      Carbohydrate_metabolism = "#f4cae4",
      Nucleotide_metabolism = "#ffff99",
      Energy_metabolism = "#386cb0",
      Metabolism_of_cofactors_and_vitamins = "#f0027f"
    )
  )
)

# Clean up bin phylogeny annotations
BinAbPhylum = gsub('assified ', '.', BinAbPhylum)

# Define right-side annotation with order of bin phyla
RightAnno2 = rowAnnotation(
  width = unit(1.5, "cm"),
  annotation_height = unit(0.1, "cm"),
  Order = BinAbPhylum,
  col = list(
    Order = c(
      Bacteroidales = '#a6cee3',
      Burkholderiales = '#1f78b4',
      Corynebacteriales = '#b2df8a',
      Enterobacterales = '#33a02c',
      Flavobacteriales = '#fb9a99',
      Fusobacteriales = '#e31a1c',
      Hyphomicrobiales = '#fdbf6f',
      Lactobacillales = '#ff7f00',
      Moraxellales = '#cab2d6',
      Pseudomonadales = '#6a3d9a',
      Uncl_Tenericutes = '#ffff99',
      Xanthomonadales = '#b15928'
    )
  )
)

# Define color functions for heatmaps
col_fun2 = colorRamp2(c(0, 100), c("white", "black"))
col_fun3 = colorRamp2(c(0, 100), c("white", "#f73133"))

# Create heatmap for pathway completeness
a = Heatmap(
  TotalFil2,
  column_gap = unit(4, "mm"),
  name = 'Pathway completeness',
  rect_gp = gpar(type = "none"),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  col = col_fun2,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.circle(
      x = x, y = y, r = 0.32 * min(unit.c(width, height)), 
      gp = gpar(fill = col_fun2(TotalFil2[i, j]), col = 'black', lwd = 1.5)
    )
  },
  top_annotation = TopAnno,
  left_annotation = RightAnno2,
  column_names_rot = 90,
  column_title = NULL,
  row_title = NULL,
  row_dend_gp = gpar(lwd = 2),
  heatmap_legend_param = list(title_position = "lefttop-rot"),
  row_dend_width = unit(1.5, "cm")
) + rowAnnotation(Bin_name = anno_text(rownames(TotalFil)))

# Save the heatmap to a PowerPoint file
topptx(a, '240629KEGGAbCompleteness.pptx', width = 11, height = 7.5)
