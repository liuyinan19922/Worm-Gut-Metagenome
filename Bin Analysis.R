library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggprism)
library(ggalluvial)
library(dplyr)
library(SQMtools)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(dendextend)

# Load the SQM data for worm plastic analysis
wormPlastic = loadSQM("D:/230408 Metagenomic analysis/wormPlastic2/")
metadata = read_excel("D:/230408 Metagenomic analysis/240311 Metagenome/metaForMetagenomeUpdated.xlsx")
metadata = as.data.frame(metadata)

# Extract bins and filter based on completeness and contamination
BinTable = wormPlastic[["bins"]][["table"]]
n1 = which(BinTable$Completeness >= 90) # Filter bins with completeness >= 90
n2 = which(BinTable$Contamination < 5)  # Filter bins with contamination < 5
n3 = intersect(n1, n2)  # Get bins satisfying both conditions
BinTable1 = BinTable[n3,]  # Subset of high-quality bins

# Copy the selected high-quality bins to a new folder
FileNameFrom = paste0("D:/230408 Metagenomic analysis/240311 Metagenome/240423 binning analysis/bins/", rownames(BinTable1), '.fa')
FileNameTo = paste0("D:/230408 Metagenomic analysis/240311 Metagenome/240423 binning analysis/BinsHig/", rownames(BinTable1), '.fa')

for (i in 1:length(FileNameFrom)) {
  file.copy(FileNameFrom[i], FileNameTo[i])  # Copy each bin file
}

# Filter abundance data based on CPM and calculate sums
data = as.data.frame(wormPlastic[["bins"]][["cpm"]][n3,])
data = data[, metadata$SampleName]
data$rSum = rowSums(data)
data1 = dplyr::arrange(data, desc(rSum))  # Sort bins by abundance
data1 = mutate(data1, mean = rSum / 33)   # Calculate mean abundance per sample
n4 = rownames(data1)[1:32]  # Select bins with the highest abundance (top 32)

# Filter out columns 34 and 35 (metadata-related columns)
data2 = data1[, -c(34, 35)]
data2t = t(data2)
data2t = as.data.frame(data2t)

# Aggregate data by sample group
data2_group = aggregate(t(data2), by = list(metadata$Full), mean)
data2_group = column_to_rownames(data2_group, "Group.1")
data2_groupt = t(data2_group)
data2_groupt = data2_groupt[, unique(metadata$Full)]

# Calculate standard deviation for each group
data2_group_sd = aggregate(t(data2), by = list(metadata$Full), sd)
data2_group_sd = column_to_rownames(data2_group_sd, "Group.1")
data2_groupt_sd = t(data2_group_sd)
data2_groupt_sd = data2_groupt_sd[, unique(metadata$Full)]

# Round the mean and standard deviation tables for output
data2_groupt_2 = round(data2_groupt, 2)
data2_groupt_sd_2 = round(data2_groupt_sd, 2)
data2_in_table = paste0(data2_groupt_2, ' ± ', data2_groupt_sd_2)

# Create a matrix of mean and standard deviation values
data2_in_table1 = matrix(data2_in_table, nrow = nrow(data2_groupt_2), ncol = ncol(data2_groupt_2), byrow = FALSE)
row.names(data2_in_table1) = row.names(data2_groupt_2)
colnames(data2_in_table1) = colnames(data2_groupt_2)

# Save the results to CSV
write.csv(data2_in_table1, '240721binAbundanceMeanSd.csv')

# Prepare data for heatmap visualization
data2_groupt_tr = t(scale(t(data2_groupt)))  # Normalize data for heatmap

# Extract metadata for sample types and positions
Type = metadata$Type
dim(Type) = c(3, 11)
Type = Type[1,]
Position = metadata$Position
dim(Position) = c(3, 11)
Position = Position[1,]

# Define heatmap annotations for sample types and plastic types
TopAnno = HeatmapAnnotation(
  gap = unit(0.05, "cm"),
  height = unit(1, "cm"),
  annotation_height = unit(0.1, "cm"),
  Plastic_type = Type,
  Sample_type = Position,
  col = list(
    Plastic_type = c(Control = '#000000', HDPE = "#fc6e63", PP = "#01bd0c", PS = "#4e9eff"),
    Sample_type = c(Inoculum = "#444444", Biofilm = "#888888", Planktonic = "#e7e7e7")
  )  # Define colors for each type
)

RightAnno = HeatmapAnnotation(
  Abundance = anno_boxplot(data2t, add_points = TRUE, annotation_label = 'Abundance (cpm)',
                           gp = gpar(col = '#00bd0e', fill = '#e2f0d9', lwd = 1), size = unit(1, "pt")),
  which = 'row', width = unit(3, "cm")
)
RightAnno@anno_list[["Abundance"]]@label = 'Abundance (cpm)'

# Define color scale for the heatmap
col_fun2 = colorRamp2(c(-1.5, 0, 1.5), c("#3894bf", "white", "#f48c8f"))
x = apply(data2_groupt, 1, mean)

# Define dendrogram for clustering rows
dend = as.dendrogram(hclust(dist(data2_groupt_tr)))
dend = color_branches(dend, k = 2, col = c("#3894bf", "#f48c8f"))

# Define column index for heatmap
ColIndex = c(1, 2, 2, 2, 3, 4, 4, 5, 5, 6, 6)

# Plot the heatmap
Heatmap(data2_groupt_tr,
        column_gap = unit(4, "mm"),
        name = 'Normalized_abundance',
        rect_gp = gpar(col = "white", lty = 1.5, lwd = 2),
        col = col_fun2,
        cluster_rows = dend,
        column_split = ColIndex,
        top_annotation = TopAnno,
        right_annotation = RightAnno,
        column_title = NULL,
        row_title = NULL,
        row_dend_gp = gpar(lwd = 2),
        row_split = 2,
        heatmap_legend_param = list(title_position = "lefttop-rot"),
        row_dend_width = unit(1.5, "cm")
)

# Load required libraries
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggprism)
library(ggalluvial)
library(dplyr)
library(SQMtools)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(dendextend)

# Load dataset
wormPlastic = loadSQM("D:/230408 Metagenomic analysis/wormPlastic2/")
metadata <- read_excel("D:/230408 Metagenomic analysis/240311 Metagenome/metaForMetagenomeUpdated.xlsx")
metadata = as.data.frame(metadata)

# Filter bin data based on completeness and contamination
BinTable = wormPlastic[["bins"]][["table"]]
n1 = which(BinTable$Completeness >= 90)  # Select bins with completeness >= 90%
n2 = which(BinTable$Contamination < 5)   # Select bins with contamination < 5%
n3 = intersect(n1, n2)  # Get bins satisfying both conditions
BinTable1 = BinTable[n3,]

# Copy high-quality bins to a new folder
FileNameFrom = paste0("D:/230408 Metagenomic analysis/240311 Metagenome/240423 binning analysis/bins/", rownames(BinTable1), '.fa')
FileNameTo = paste0("D:/230408 Metagenomic analysis/240311 Metagenome/240423 binning analysis/BinsHig/", rownames(BinTable1), '.fa')

for (i in 1:length(FileNameFrom)) {
  file.copy(FileNameFrom[i], FileNameTo[i])
}

# Subset data based on abundance
data = as.data.frame(wormPlastic[["bins"]][["cpm"]][n3,])
data = data[, metadata$SampleName]
data$rSum = rowSums(data)
data1 = dplyr::arrange(data, desc(rSum))
data1 = mutate(data1, mean = rSum / 33)
n4 = rownames(data1)[1:32]  # Select bins with cpm > 0.1

# Data transformation and aggregation
data2 = data1[1:32, ]
data2 = data2[, c(-34, -35)]  # Remove unnecessary columns
data2t = t(data2)
data2t = as.data.frame(data2t)
data2_group = aggregate(t(data2), by = list(metadata$Full), mean)
data2_group = column_to_rownames(data2_group, "Group.1")
data2_groupt = t(data2_group)
data2_groupt = data2_groupt[, unique(metadata$Full)]

# Calculate mean and standard deviation for abundance data
data2_group_sd = aggregate(t(data2), by = list(metadata$Full), sd)
data2_group_sd = column_to_rownames(data2_group_sd, "Group.1")
data2_groupt_sd = t(data2_group_sd)
data2_groupt_sd = data2_groupt_sd[, unique(metadata$Full)]

# Prepare final table for export
data2_groupt_2 = round(data2_groupt, 2)
data2_groupt_sd_2 = round(data2_groupt_sd, 2)
data2_in_table = paste0(data2_groupt_2, ' ± ', data2_groupt_sd_2)
data2_in_table1 = matrix(data2_in_table, nrow = nrow(data2_groupt_2), ncol = ncol(data2_groupt_2), byrow = FALSE)
row.names(data2_in_table1) = row.names(data2_groupt_2)
colnames(data2_in_table1) = colnames(data2_groupt_2)
write.csv(data2_in_table1, '240721binAbundanceMeanSd.csv')

# Generate heatmap data
data2_groupt_tr = t(scale(t(data2_groupt)))

Type = metadata$Type
dim(Type) = c(3, 11)
Type = Type[1,]

Position = metadata$Position
dim(Position) = c(3, 11)
Position = Position[1,]

# Heatmap annotation setup
TopAnno = HeatmapAnnotation(
  gap = unit(0.05, "cm"),
  height = unit(1, "cm"),
  annotation_height = unit(0.1, "cm"),
  Plastic_type = Type,
  Sample_type = Position,
  col = list(
    Plastic_type = c(Control = '#000000', HDPE = "#fc6e63", PP = "#01bd0c", PS = "#4e9eff"),
    Sample_type = c(Inoculum = "#444444", Biofilm = "#888888", Planktonic = "#e7e7e7"))
)

RightAnno = HeatmapAnnotation(Abundance = anno_boxplot(data2t, add_points = TRUE, annotation_label = 'Abundance (cpm)', 
                                                        gp = gpar(col = '#00bd0e', fill = '#e2f0d9', lwd = 1), size = unit(1, "pt")),
                              which = 'row', width = unit(3, "cm"))
RightAnno@anno_list[["Abundance"]]@label = 'Abundance (cpm)'

# Set color scheme for heatmap
col_fun2 = colorRamp2(c(-1.5, 0, 1.5), c("#3894bf", "white", "#f48c8f"))
x = apply(data2_groupt, 1, mean)

# Clustered heatmap with annotations
dend = as.dendrogram(hclust(dist(data2_groupt_tr)))
dend = color_branches(dend, k = 2, col = c("#3894bf", "#f48c8f"))

Heatmap(data2_groupt_tr,
        column_gap = unit(4, "mm"),
        name = 'Normalized_abundance',
        rect_gp = gpar(col = "white", lty = 1.5, lwd = 2),
        col = col_fun2,
        cluster_rows = dend,
        column_split = ColIndex,
        top_annotation = TopAnno,
        right_annotation = RightAnno,
        column_title = NULL,
        row_title = NULL,
        row_dend_gp = gpar(lwd = 2),
        row_split = 2,
        heatmap_legend_param = list(title_position = "lefttop-rot"),
        row_dend_width = unit(1.5, "cm"))

# Further data preparation for analysis
annotation = read.table("D:\\230408 Metagenomic analysis\\240311 Metagenome\\BinAnnotationFull.tsv")[1:70, 1:2]
annotation = column_to_rownames(annotation, 'V1')

tax_annotation = sapply(t(annotation), function(x) strsplit(x, ":")[[1]])[3,]
tax_id = as.numeric(sapply(t(annotation), function(x) strsplit(x, ":")[[1]])[4,])
tax_binid = row.names(annotation)

Anno = gsub('k__|\\|p__|\\|c__|\\|o__|\\|f__|\\|g__|\\|s__|\\|t__', 'SEP', tax_annotation) %>%
  sapply(function(x) strsplit(x, c("SEP")))
Anno = as.data.frame(t(as.data.frame(Anno)))
Anno = Anno[, -1]
Anno$tax_binid = tax_binid
Anno$iden = tax_id

# Rename columns and classify SGB
colnames(Anno) = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'SGBNum', 'BinID', "Identity")
Anno = mutate(Anno, SGB = case_when(Anno$Identity <= 0.05 ~ "kSGB", Anno$Identity > 0.05 ~ 'uSGB'))
Anno = Anno[, -c(8, 10)]
row.names(Anno) = NULL
Anno = column_to_rownames(Anno, 'BinID')

# Annotate and update taxonomy
rawtax = as.data.frame(wormPlastic[['bins']][["tax"]])[rownames(BinTable1),]
rawtax$SGB = rep('Novel', 99)
colnames(rawtax) = colnames(Anno)

for (i in 1:nrow(rawtax)) {
  if (rownames(rawtax)[i] %in% rownames(Anno)) {
    if (Anno[rownames(rawtax)[i], 'SGB'] == 'kSGB') {
      rawtax[i, ] = Anno[rownames(rawtax)[i], ]
    }
  }
}

# Modify bin names
temp1 = sapply(rownames(rawtax), function(x) strsplit(x, "\\.")[[1]][1])
temp2 = sapply(rownames(rawtax), function(x) strsplit(x, "\\.")[[1]][2])
rawtax$NewName = paste0('Bin.', temp1, '.', temp2)
rawtax$NewName = gsub('concoct', 'c', rawtax$NewName)
rawtax$NewName = gsub('metabat2', 'm', rawtax$NewName)

# Save the annotated table
rawtax2 = rawtax
BinTable1 = as.data.frame(BinTable1)
rawtax2 = cbind(rawtax2, BinTable1)
write.csv(rawtax2, '240603BinMeta.csv')

# Prepare table for itol visualization
BinTable2 = BinTable1[, 11:76]
num = rep(1:33)
num = (num * 2) - 1
BinTable3 = BinTable2[, num]
colnames(BinTable3) = gsub('Coverage ', '', colnames(BinTable3))
BinTable4 = BinTable3[, metadata$SampleName]
library(clipr)
write_clip(BinTable4)

# iTol visualization setup (Tree and annotations)
library(itol.toolkit)
library(data.table)
library(dplyr)
library(ape)
library(stringr)
library(tidyr)
library(RColorBrewer)

TreeFile = read.tree("D:/230408 Metagenomic analysis/240311 Metagenome/240423 binning analysis/Bins_output_FullMedQual/RAxML_bestTree.Bins_Medium_refined.tre")
rawtax2 = rownames_to_column(rawtax2, 'RawName')
DataFile = rawtax2

# Create hub and units for iTol visualization
hub = create_hub(TreeFile)
Unit1 = create_unit(data = rawtax2 %>% select(RawName, NewName), key = 'Rename', type = 'LABELS', tree = TreeFile)
Unit2 = create_unit(data = rawtax2 %>% select(RawName, Phylum), key = 'Phylum', type = 'DATASET_COLORSTRIP', tree = TreeFile)
Unit3 = create_unit(data = rawtax2 %>% select(RawName, SGB), key = "SGB", type = "DATASET_COLORSTRIP", tree = TreeFile, color = 'aaas')

hub3 = hub + Unit1 + Unit2 + Unit3
write_hub(hub3, getwd())

# Export the iTol units
unit7 = create_unit(data = rawtax2 %>% select(RawName, Phylum), key = "Phylum", type = "TREE_COLORS", subtype = "range", tree = TreeFile, color = "Pastel2")
write_unit(unit7)

unit8 = create_unit(data = rawtax2 %>% select(RawName, Species), key = 'Species', type = "DATASET_TEXT", size_factor = 1, rotation = 0, position = -1, color = "#000000", tree = TreeFile)
write_unit(unit8)

unit9 = create_unit(data = rawtax2 %>% select(RawName, Order), key = "Order", type = "TREE_COLORS", subtype = "range", tree = TreeFile, color = "Set3")
setwd("d:/Liu Yi-Nan/Desktop/")
write_unit(unit9)

# Analyze abundant bins
rawtaxAb = as.data.frame(wormPlastic[['bins']][["tax"]])[n4,]
rawtaxAb$SGB = rep('Novel', 32)
colnames(rawtaxAb) = colnames(Anno)

for (i in 1:nrow(rawtaxAb)) {
  if (rownames(rawtaxAb)[i] %in% rownames(Anno)) {
    if (Anno[rownames(rawtaxAb)[i], 'SGB'] == 'kSGB') {
      rawtaxAb[i,] = Anno[rownames(rawtaxAb)[i],]
    }
  }
}

# Modify bin names for abundant bins
temp1 = sapply(rownames(rawtaxAb), function(x) strsplit(x, "\\.")[[1]][1])
temp2 = sapply(rownames(rawtaxAb), function(x) strsplit(x, "\\.")[[1]][2])
rawtaxAb$NewName = paste0('Bin.', temp1, '.', temp2)
rawtaxAb$NewName = gsub('concoct', 'c', rawtaxAb$NewName)
rawtaxAb$NewName = gsub('metabat2', 'm', rawtaxAb$NewName)

rawtax2Ab = rawtaxAb
rawtax2Ab = rownames_to_column(rawtax2Ab, 'RawName')

# Load tree file
TreeFile = read.tree("D:\\【【工作 from 191029】】\\230408 Metagenomic analysis\\240311 Metagenome\\240423 binning analysis\\Bins_output_FullMedQual\\RAxML_bestTree.Bins_Medium_refined.tre")

# Initialize a temporary vector for column name mapping
temp = vector()   

# Map raw names to columns in TotalFil
for (j in 1:ncol(TotalFil)) {
  temp[j] = rawtax2Ab$RawName[which(colnames(TotalFil)[j] == rawtax2Ab$NewName)] 
}
colnames(TotalFil) = temp

# Transpose TotalFil and convert to a data frame
TotalFil = t(TotalFil)
TotalFil = as.data.frame(TotalFil)
TotalFil = rownames_to_column(TotalFil, 'BinRawName')

# Create hub and heatmap visualization for KEGG bins
hub = create_hub(TreeFile)
Unit1 = create_unit(data = TotalFil, key = 'KEGGBinsAb', type = 'DATASET_HEATMAP', tree = TreeFile)
Unit1@specific_themes$heatmap$tree$tree_display <- 0
Unit1@specific_themes$heatmap$color$min <- 'white'
Unit1@specific_themes$heatmap$color$max <- "#484cab"
Unit1@specific_themes$heatmap$use_mid <- 0
hub1 = hub + Unit1
write_hub(hub1, getwd())

# Map raw names to columns in TotalDirect
temp = vector()   
for (j in 1:ncol(TotalDirect)) {
  temp[j] = rawtax2Ab$RawName[which(colnames(TotalDirect)[j] == rawtax2Ab$NewName)] 
}
colnames(TotalDirect) = temp

# Transpose TotalDirect and convert to a data frame
TotalDirect = t(TotalDirect)
TotalDirect = as.data.frame(TotalDirect)
TotalDirect = rownames_to_column(TotalDirect, 'BinRawName')

# Create hub and heatmap visualization for direct KEGG bins
hub = create_hub(TreeFile)
Unit1 = create_unit(data = TotalDirect, key = 'KEGGBinsAbDirect', type = 'DATASET_HEATMAP', tree = TreeFile)
Unit1@specific_themes$heatmap$tree$tree_display <- 0
Unit1@specific_themes$heatmap$color$min <- 'white'
Unit1@specific_themes$heatmap$color$max <- "#484cab"
Unit1@specific_themes$heatmap$use_mid <- 0
hub1 = hub + Unit1
write_hub(hub1, getwd())

# Create hub for additional tree customization
hub = create_hub(TreeFile)

# Unit for pruning by abundance
Unit0 = create_unit(data = rawtax2Ab$RawName, key = "PruneByAbundance", type = "PRUNE", tree = TreeFile)

# Unit for renaming bins
Unit1 = create_unit(data = rawtax2Ab %>% select(RawName, NewName), key = 'Rename', type = 'LABELS', tree = TreeFile)

# Unit for displaying species as text labels
Unit2 = create_unit(data = rawtax2Ab %>% select(RawName, Species), key = 'Species', type = "DATASET_TEXT", size_factor = 1,
                    rotation = 0, position = -1, color = "#000000", tree = TreeFile)

# Unit for genus color strip visualization
Unit3 = create_unit(data = rawtax2Ab %>% select(RawName, Genus), color = 'Pastel1', key = "Genus", type = "TREE_COLORS", subtype = "range", tree = TreeFile)

# Unit for SGB color strip visualization
Unit4 = create_unit(data = rawtax2Ab %>% select(RawName, SGB), key = "SGB", type = "DATASET_COLORSTRIP", tree = TreeFile, color = 'aaas')

# Unit for clade genus visualization
Unit5 = create_unit(data = rawtax2Ab %>% select(RawName, Genus), key = "Clade_genus", type = "TREE_COLORS", subtype = "clade", 
                    size_factor = 4, tree = TreeFile)

# Combine all units into a hub
hub_total = hub + Unit0 + Unit1 + Unit2 + Unit3 + Unit4 + Unit5
write_hub(hub_total, getwd())

# Write unit for displaying species labels
write_hub(hub + Unit2)

# Unit for displaying order as a color range in the tree
unit9 = create_unit(data = rawtax2Ab %>% select(RawName, Order),
                    key = "Order", 
                    type = "TREE_COLORS", 
                    subtype = "range", 
                    tree = TreeFile,
                    color = "Paired")

write_unit(unit9)
