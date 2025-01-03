library(tidyverse)
library(ggprism)
library(stringr)
library(readxl)
library('SQMtools')
library(pheatmap)
library(Hmisc)
library("agricolae")
library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(eoffice)

# Data import
wormPlastic = loadSQM('D:/Path_to_data/wormPlastic2/', tax_mode = "nofilter", trusted_functions_only = FALSE, engine = "data.table")
metadata <- read_excel("D:/Path_to_data/metaForMetagenomeUpdated.xlsx")

PlasticDB = read.csv("D:/Path_to_data/240413PlasticDBAlign.csv", header = TRUE)
PlasticDB = mutate(PlasticDB, Recovery = (send - sstart + 1) / slen)
PlasticDB = subset(PlasticDB, pident > 50)
PlasticDB = subset(PlasticDB, Recovery > 0.50)
PlasticDB = subset(PlasticDB, evalue < 1e-10)

simplifyPlasticDB = PlasticDB

# Process the best annotations
simplifyPlasticDB = transform(simplifyPlasticDB, Index = sapply(simplifyPlasticDB$sseqid, function(x) strsplit(x, "\\|\\|")[[1]][1]))
simplifyPlasticDB = transform(simplifyPlasticDB, EnzymeName = sapply(simplifyPlasticDB$sseqid, function(x) strsplit(x, "\\|\\|")[[1]][2]))
simplifyPlasticDB = transform(simplifyPlasticDB, PlasitcType = sapply(simplifyPlasticDB$sseqid, function(x) strsplit(x, "\\|\\|")[[1]][4]))
simplifyPlasticDB = transform(simplifyPlasticDB, NameType = paste0(EnzymeName, ' (', PlasitcType, ") ", Index))

# Remove unrelated enzymes
simplifyPlasticDB = subset(simplifyPlasticDB, !(Index %in% c('00016','00169','00170','00071','00172','00114')))

# Include signal peptide information
Sp = read.csv("D:/Path_to_data/prediction_results.csv", header = TRUE)
Sp[,1] = sapply(Sp[,1], function(x) strsplit(x, " _")[[1]][1])

# Extract metagenomic data abundance based on best annotations
GoiList = simplifyPlasticDB$qseqid
TpmTable = wormPlastic[["orfs"]][["tpm"]]
TpmGoi = subset(TpmTable, rownames(TpmTable) %in% GoiList)
simplifyPlasticDB1 = merge(TpmGoi, simplifyPlasticDB, by.x = 0, by.y = "qseqid")
simplifyPlasticDB1 = merge(simplifyPlasticDB1, Sp, by.x = "Row.names", by.y = "X..ID")

# Annotate enzymes with signal peptide prediction
temp2 = simplifyPlasticDB1[, c("NameType", 'Prediction')]
temp2 = subset(temp2, Prediction != 'OTHER')
temp3 = paste0(temp2[, 1], '&', temp2[, 2])
temp4 = unique(temp3)
Position = sapply(temp4, function(x) strsplit(x, "&")[[1]][1])
S = sapply(temp4, function(x) strsplit(x, "&")[[1]][2])
GeneSp = data.frame(Position = Position, Signal_peptide = S)

# Add enzyme type annotation
simplifyPlasticDB1 = simplifyPlasticDB1 %>% mutate(Enzyme_type = case_when(
  Index %in% c('00174','00111','00209','00093','00180','00061','00104','00035') ~ 'Oxidase',  
  !(Index %in% c('00174','00111','00209','00093','00180','00061','00104','00035')) ~ 'Hydrolase'
))
temp5 = simplifyPlasticDB1[, c("NameType", 'Enzyme_type')]
temp6 = paste0(temp5[, 1], '&', temp5[, 2])
temp7 = unique(temp6)
Position1 = sapply(temp7, function(x) strsplit(x, "&")[[1]][1])
Enzyme_type = sapply(temp7, function(x) strsplit(x, "&")[[1]][2])
names(Enzyme_type) = Position1
Enzyme_type1 = Enzyme_type    

# Write results to CSV file
# write.csv(simplifyPlasticDB1, "240414plasticDBTpmTable.csv")

# Aggregate data by enzyme type
simplifyPlasticDB2 = aggregate(simplifyPlasticDB1[2:34], by = list(simplifyPlasticDB1$NameType), sum)

# Create heatmap
simplifyPlasticDB2heat = column_to_rownames(simplifyPlasticDB2, 'Group.1')
simplifyPlasticDB2heat = simplifyPlasticDB2heat[, metadata$SampleName]

# Prepare the data for analysis
data3 = simplifyPlasticDB1
data2 = simplifyPlasticDB2heat

# Aggregate data by group
data2_group = aggregate(t(data2), by = list(metadata$Full), mean)
data2_group_sd = aggregate(t(data2), by = list(metadata$Full), sd)
data2_group = column_to_rownames(data2_group, "Group.1")
data2_group_sd = column_to_rownames(data2_group_sd, "Group.1")

# Transpose the data for easier processing
data2_groupt = t(data2_group)
data2_groupt = data2_groupt[, unique(metadata$Full)]
data2_groupt_sd = t(data2_group_sd)
data2_groupt_sd = data2_groupt_sd[, unique(metadata$Full)]

# Prepare the table with the mean and standard deviation values
data2_groupt_table = paste0(round(data2_groupt, 2), ' ± ', round(data2_groupt_sd, 2))
data2_groupt_table = matrix(data2_groupt_table, nrow = nrow(data2_groupt), ncol = ncol(data2_groupt_sd), byrow = FALSE)
row.names(data2_groupt_table) = row.names(data2_groupt)
colnames(data2_groupt_table) = colnames(data2_groupt)
write.table(data2_groupt_table, "clipboard", sep = ",", quote = F)

# Prepare data for percentage calculation of different enzyme types
RedoxPlasticDB = subset(simplifyPlasticDB, Index %in% c('00174', '00111', '00209', '00093', '00180', '00061', '00104', '00035'))
HydrolysisPlasticDB = subset(simplifyPlasticDB, !(Index %in% c('00174', '00111', '00209', '00093', '00180', '00061', '00104', '00035')))
RedoxName = RedoxPlasticDB$NameType
HydrolysisName = HydrolysisPlasticDB$NameType

# Select the corresponding rows based on enzyme names
n1 = which(row.names(data2_groupt) %in% RedoxName)
n2 = which(row.names(data2_groupt) %in% HydrolysisName)
Redox = data2_groupt[n1, ]
Redox_sd = data2_groupt_sd[n1, ]
RedoxSum = colSums(Redox)

Hydrolysis = data2_groupt[n2, ]
HydrolysisSum = colSums(Hydrolysis)
Hydrolysis_sd = data2_groupt_sd[n2, ]

# Calculate percentages and standard deviations for Redox and Hydrolysis
RedoxPercentage = Redox
Redox_sdPercentage = Redox_sd
HydrolysisPercentage = Hydrolysis
Hydrolysis_sdPercentage = Hydrolysis_sd

for (j in 1:ncol(Redox)) {
  RedoxPercentage[, j] = round(Redox[, j] / RedoxSum[j], 4)
  Redox_sdPercentage[, j] = round(Redox_sd[, j] / RedoxSum[j], 4)
  HydrolysisPercentage[, j] = round(Hydrolysis[, j] / HydrolysisSum[j], 4)
  Hydrolysis_sdPercentage[, j] = round(Hydrolysis_sd[, j] / HydrolysisSum[j], 4)
}

# Narrow down to columns of interest
RedoxPercentage = RedoxPercentage[, 6:11]
Redox_sdPercentage = Redox_sdPercentage[, 6:11]
HydrolysisPercentage = HydrolysisPercentage[, 6:11]
Hydrolysis_sdPercentage = Hydrolysis_sdPercentage[, 6:11]

# Prepare the table for Redox and Hydrolysis
Redox_in_table = paste0(RedoxPercentage, ' ± ', Redox_sdPercentage)
Redox_in_table1 = matrix(Redox_in_table, nrow = nrow(RedoxPercentage), ncol = ncol(RedoxPercentage), byrow = FALSE)
row.names(Redox_in_table1) = row.names(Redox)
colnames(Redox_in_table1) = colnames(RedoxPercentage)
write.csv(Redox_in_table1, '240803RedoxPlasticDBMeanSd.csv')

Hydrolysis_in_table = paste0(HydrolysisPercentage, ' ± ', Hydrolysis_sdPercentage)
Hydrolysis_in_table1 = matrix(Hydrolysis_in_table, nrow = nrow(HydrolysisPercentage), ncol = ncol(HydrolysisPercentage), byrow = FALSE)
row.names(Hydrolysis_in_table1) = row.names(Hydrolysis)
colnames(Hydrolysis_in_table1) = colnames(HydrolysisPercentage)
write.csv(Hydrolysis_in_table1, '240803HydrolysisPlasticDBMeanSd.csv')

# Prepare data for visualization
data = RedoxPercentage
row.names(data) = gsub("Unclassified", 'Uncl.', row.names(data))
row.names(data) = gsub("Uncl. Bacteria", 'Unclassified Bacteria', row.names(data))

# Set metadata for the plot
metadata = c("HDPE", 'HDPE', 'PP', "PP", "PS", 'PS')
names(metadata) = colnames(data)
metadata = as.data.frame(metadata)

# Create the taxonomic stack plot
p = tax_stackplot(data, metadata, groupID = "metadata", style = "sample", topN = 9)
p
OTUdata1 = p[["data"]]

# Define custom colors for the plot
Yanse = c('#00B8DF', '#FF6C65', '#00B800', '#D39100', '#619BFA', '#DB72F6', '#9ed2ff', '#feb1ad', '#a0c3a0', '#dac597', '#a6b7fb', '#e8b2f7', '#c6c6c6', '#156f82', '#a2524f', '#177817', '#765714', '#4e6da1', '#9459a3')
Yanse1 = c(Yanse[2], Yanse[1], Yanse[5], Yanse[7], Yanse[8], Yanse[3], Yanse[9], Yanse[4], Yanse[12], Yanse[6], Yanse[13])

# Create the stacked bar plot
p1 = ggplot(OTUdata1, aes(x = variable, y = value, fill = Taxonomy)) +
  geom_bar(colour = 'white', alpha = 0.85, position = "fill", stat = "identity", width = 0.4, size = 0.4) +
  facet_grid(.~group, space = 'free', scale = 'free') +
  scale_fill_manual(values = Yanse) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Plastic type', y = 'Abundance') +
  guides(fill = guide_legend(reverse = F)) +
  theme_prism(border = TRUE) +
  theme(panel.spacing.x = unit(0.15, "cm"), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title = element_text(size = 18))

p1

# Export the plot to a PowerPoint file
library(eoffice)
topptx(p1, "240803RedoxContributionEnzyme.pptx", width = 6.5, height = 3.75)

# Prepare the data for heatmap annotation
data2_groupt_tr = t(scale(t(data2_groupt)))
data2t = as.data.frame(t(data2))

# Add signal peptide data
Position2 = setdiff(row.names(data2_groupt_tr), GeneSp$Position)
S2 = rep('N.D.', length(Position2))
GeneSp2 = data.frame(Position = Position2, Signal_peptide = S2)
GeneSp = rbind(GeneSp, GeneSp2)
GeneSp = GeneSp[-4, ]
GeneSp3 = GeneSp[[2]]
names(GeneSp3) = GeneSp[[1]]

# Prepare the heatmap annotations
Type = metadata$Type
dim(Type) = c(3, 11)
Type = Type[1, ]
Position = metadata$Position
dim(Position) = c(3, 11)
Position = Position[1, ]

TopAnno = HeatmapAnnotation(
  gap = unit(0.05, "cm"),
  height = unit(1, "cm"),
  annotation_height = unit(0.1, "cm"),
  Plastic_type = Type,
  Sample_type = Position,
  col = list(
    Plastic_type = c(Control = '#000000', HDPE = "#fc6e63", PP = "#01bd0c", PS = "#4e9eff"),
    Sample_type = c(Inoculum = "#444444", Biofilm = "#888888", Planktonic = "#e7e7e7")
  )
)

# Add additional annotations to the heatmap
RightAnno1 = HeatmapAnnotation(Count = anno_barplot(as.vector(table(data3$NameType)), gp = gpar(col = '#00bd0e', fill = '#e2f0d9', lwd = 1)), which = 'row', width = unit(1.5, "cm"))
RightAnno2 = HeatmapAnnotation(Abundance = anno_boxplot(data2t, add_points = T, annotation_label = 'Abundance (%)', gp = gpar(col = '#fc6e63', fill = '#f9b6b2', lwd = 1), size = unit(1, "pt")), which = 'row', width = unit(1.5, "cm"))
RightAnno2@anno_list[["Abundance"]]@label = 'Abundance (tpm)'
RightAnno3 = rowAnnotation(Signal_peptide = GeneSp3[row.names(data2_groupt_tr)], col = list(Signal_peptide = c(SP = "#6baed6", TAT = "#0053a1", LIPO = '#92cce4', N.D. = "#e7e7e7")))
RightAnno4 = rowAnnotation(Enzyme_type = Enzyme_type[row.names(data2_groupt_tr)], col = list(Enzyme_type = c(Oxidase = "#33a02b", Hydrolase = "#e0f3db")))

# Define color palette
col_fun2 <- colorRamp2(c(-1.5, 0, 1.5), c("#3894bf", "white", "#f48c8f"))

# Calculate row means
x = apply(data2_groupt, 1, mean)

# Define column indices for grouping
ColIndex = c(1, 2, 2, 2, 3, 4, 4, 5, 5, 6, 6)

# Create a dendrogram using hierarchical clustering
dend = as.dendrogram(hclust(dist(data2_groupt_tr)))
dend = color_branches(dend, k = 2, col = c("#f48c8f", "#3894bf"))

# Load required library
library('eoffice')

# Create heatmap with dendrogram and annotations
Fig = Heatmap(data2_groupt_tr,
              border = 'black',
              column_gap = unit(4, "mm"),
              name = 'Normalized_abundance',
              rect_gp = gpar(col = "white", lty = 1.5, lwd = 4),
              cluster_columns = F,
              col = col_fun2,
              row_names_side = 'right',
              cluster_rows = dend,
              column_split = ColIndex,
              top_annotation = TopAnno,
              column_title = NULL,
              row_title = NULL,
              row_split = 2,
              row_dend_gp = gpar(lwd = 2.5),
              heatmap_legend_param = list(title_position = "lefttop-rot"),
              row_dend_width = unit(1.5, "cm")) +
      RightAnno1 + RightAnno2 + RightAnno3 + RightAnno4 + 
      rowAnnotation(rn = anno_text(rownames(data2_groupt_tr)), annotation_name_gp = gpar(fontsize = 60))

# Save figure to PowerPoint
topptx(Fig, '240507PlasticDBHeatMapRmMeaninglessSP.pptx', width = 8, height = 7, append = TRUE)

# Statistical tests for enzyme types
# Perform pairwise t-test
TestData = data.frame(t(simplifyPlasticDB2heat))
TestData$Type = metadata$Full
Res.ttest = pairwise.t.test(TestData[, 1], metadata$Full, p.adjust.method = "BH")

# LSD test for comparisons between groups
TestRes = vector()
for (i in 1:(ncol(TestData) - 1)) {
    model <- aov(TestData[, i] ~ Type, data = TestData)
    out <- LSD.test(model, "Type", p.adj = "BH")
    TestRes1 = out[["groups"]][unique(metadata$Full), ]
    TestRes1 = TestRes1$groups
    TestRes = rbind(TestRes, TestRes1)
}
colnames(TestRes) = unique(metadata$Full)
rownames(TestRes) = colnames(TestData)[1:33]
write.csv(TestRes, "240715PlasticDBEnzymeEachConditionBH.csv")

# Test for group-wise comparison using LSD
TestData$Group = metadata$Group
TestRes = vector()
for (i in 1:(ncol(TestData) - 2)) {
    model <- aov(TestData[, i] ~ Group, data = TestData)
    out <- LSD.test(model, "Group", p.adj = "bonferroni")
    TestRes1 = out[["groups"]][unique(metadata$Group), ]
    TestRes1 = TestRes1$groups
    TestRes = rbind(TestRes, TestRes1)
}
colnames(TestRes) = unique(metadata$Group)
rownames(TestRes) = colnames(TestData[, 1:33])
write.csv(TestRes, "240715PlasticDBEnzymeGroup.csv")

# Modify Group variable based on metadata
TypePlus = vector()
for (i in 1:33) {
    if (metadata$Stage[i] == "Inoculum") {
        TypePlus[i] = metadata$Full[i]
    }
    if (metadata$Stage[i] == "Late") {
        TypePlus[i] = metadata$Type[i]
    }
}
TestData$Group = TypePlus
metadata$Group = TypePlus

# Test again with modified group
TestRes = vector()
for (i in 1:(ncol(TestData) - 2)) {
    model <- aov(TestData[, i] ~ Group, data = TestData)
    out <- LSD.test(model, "Group", p.adj = "bonferroni")
    TestRes1 = out[["groups"]][unique(metadata$Group), ]
    TestRes1 = TestRes1$groups
    TestRes = rbind(TestRes, TestRes1)
}
colnames(TestRes) = unique(metadata$Group)
rownames(TestRes) = colnames(TestData[, 1:33])
write.csv(TestRes, "240715PlasticDBEnzymeTypePlus.csv")

# Analyze hit samples
data1 = simplifyPlasticDB1[2:34][metadata$SampleName]
hit_sample = colSums(data1 != 0)
hit_sample_mean = aggregate(hit_sample, list(metadata$Full), mean)
hit_sample_mean = column_to_rownames(hit_sample_mean, 'Group.1')
hit_sample_mean = hit_sample_mean[unique(metadata$Full), ]
hit_sample_sd = aggregate(hit_sample, list(metadata$Full), sd)
hit_sample_sd = column_to_rownames(hit_sample_sd, 'Group.1')
hit_sample_sd = hit_sample_sd[unique(metadata$Full), ]
sample_info = unique(metadata$Full)

# Calculate sum and mean for samples
sum_sample = colSums(data1)
sum_sample_mean = aggregate(sum_sample, list(metadata$Full), mean)
sum_sample_mean = column_to_rownames(sum_sample_mean, 'Group.1')
sum_sample_mean = sum_sample_mean[unique(metadata$Full), ]

sum_sample_sd = aggregate(sum_sample, list(metadata$Full), sd)
sum_sample_sd = column_to_rownames(sum_sample_sd, 'Group.1')
sum_sample_sd = sum_sample_sd[unique(metadata$Full), ]

# Subset and calculate sum for hydrolase enzyme type
Enzyme_type1 = Enzyme_type1[row.names(data2)]
data2_hydro = subset(data2, Enzyme_type1 == "Hydrolase")
sum_sample_hydro = colSums(data2_hydro)

sum_sample_hydro_mean = aggregate(sum_sample_hydro, list(metadata$Full), mean)
sum_sample_hydro_mean = column_to_rownames(sum_sample_hydro_mean, 'Group.1')
sum_sample_hydro_mean = sum_sample_hydro_mean[unique(metadata$Full), ]

sum_sample_hydro_sd = aggregate(sum_sample_hydro, list(metadata$Full), sd)
sum_sample_hydro_sd = column_to_rownames(sum_sample_hydro_sd, 'Group.1')
sum_sample_hydro_sd = sum_sample_hydro_sd[unique(metadata$Full), ]

Enzyme_type = rep('Hydrolase', 11)
sample_info = unique(metadata$Full)
sum_value = sum_sample_hydro_mean + sum_sample_redox_mean

sum_sample_hydro = data.frame(sample_info, Enzyme_type, value = sum_sample_hydro_mean, error = sum_sample_hydro_sd, sum_value)

# Subset and calculate sum for oxidase enzyme type
data2_redox = subset(data2, Enzyme_type1 == "Oxidase")
sum_sample_redox = colSums(data2_redox)

sum_sample_redox_mean = aggregate(sum_sample_redox, list(metadata$Full), mean)
sum_sample_redox_mean = column_to_rownames(sum_sample_redox_mean, 'Group.1')
sum_sample_redox_mean = sum_sample_redox_mean[unique(metadata$Full), ]

sum_sample_redox_sd = aggregate(sum_sample_redox, list(metadata$Full), sd)
sum_sample_redox_sd = column_to_rownames(sum_sample_redox_sd, 'Group.1')
sum_sample_redox_sd = sum_sample_redox_sd[unique(metadata$Full), ]

Enzyme_type = rep('Oxidase', 11)
sum_sample_redox = data.frame(sample_info, Enzyme_type, value = sum_sample_redox_mean, error = sum_sample_redox_sd, sum_value)

# Combine hydrolyase and oxidase data
sum_sample_cece = rbind(sum_sample_hydro, sum_sample_redox)

# Create barplot for enzyme abundance
b = ggplot(sum_sample_cece, aes(x = sample_info, y = value, fill = Enzyme_type)) +
  geom_bar(stat = "identity", colour = "black", width = 0.6) +
  geom_errorbar(aes(ymin = sum_value - error, ymax = sum_value + error), width = 0.2) +
  scale_x_discrete(limits = unique(metadata$Full)) +
  labs(x = 'Group', y = 'Abundance (tpm)') +
  theme_prism(border = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
b

# Save barplot as PowerPoint
topptx(b, '240721EnzymeAbundanceBarplot.pptx', width = 4.8, height = 4.5)

# Define color group for hit samples
colorGroup = c('ino', "ino", "ino", 'ino', 'conincu', 'pe', 'pe', 'pp', 'pp', 'ps', 'ps')

# Create data frame for hit sample analysis
hit_sample_cece = data.frame(sample_info, mean = hit_sample_mean, sd = hit_sample_sd, colorGroup)

# Create barplot for annotated enzymes hit
c = ggplot(hit_sample_cece, aes(x = sample_info, y = mean, fill = colorGroup)) +
  geom_bar(aes(color = colorGroup), stat = "identity", width = 0.6, alpha = 0.5, lwd = 1) +
  geom_errorbar(color = 'black', aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, lwd = 0.5) +
  scale_x_discrete(limits = unique(metadata$Full)) +
  scale_fill_manual(values = c("grey", "#8896ab", "#ff6c67", "#00bd0f", "#4c9eff")) +
  scale_color_manual(values = c("grey", "#8896ab", "#ff6c67", "#00bd0f", "#4c9eff")) +
  labs(x = 'Group', y = 'Count of annotated enzymes') +
  theme_prism(border = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
c

# Save hit sample barplot as PowerPoint
topptx(c, '240721EnzymeHitBarplot.pptx', width = 4.8, height = 4.5)

# Pairwise comparison for all conditions (BH adjustment)
TestData = data.frame(sum_sample)
TestData$Type = metadata$Full
Res.ttest = pairwise.t.test(TestData[, 1], metadata$Full, p.adjust.method = "BH")

# LSD test for comparisons between groups
model = aov(TestData[, 1] ~ Type, data = TestData)
TukeyHSD(model)
out.none = LSD.test(model, "Type", p.adj = "none")
out = LSD.test(model, "Type", p.adj = "BH")

TestRes1 = out[["groups"]][unique(metadata$Full), ]
write.csv(TestRes1, "240721AllPlasticDBEachConditionBH.csv")
TestRes2 = out.none[["groups"]][unique(metadata$Full), ]
write.csv(TestRes2, "240721AllPlasticDBEachConditionNone.csv")

# Pairwise comparison for hydrolase conditions (BH adjustment)
TestData = data.frame(sum_sample_hydro)
TestData$Type = metadata$Full
Res.ttest = pairwise.t.test(TestData[, 1], metadata$Full, p.adjust.method = "BH")

# LSD test for hydrolase comparisons
model = aov(TestData[, 1] ~ Type, data = TestData)
TukeyHSD(model)
out.none = LSD.test(model, "Type", p.adj = "none")
out = LSD.test(model, "Type", p.adj = "BH")

TestRes1 = out[["groups"]][unique(metadata$Full), ]
write.csv(TestRes1, "240721HydrolasePlasticDBEachConditionBH.csv")
TestRes2 = out.none[["groups"]][unique(metadata$Full), ]
write.csv(TestRes2, "240721HydrolasePlasticDBEachConditionNone.csv")

# Pairwise comparison for redox conditions (BH adjustment)
TestData = data.frame(sum_sample_redox)
TestData$Type = metadata$Full
Res.ttest = pairwise.t.test(TestData[, 1], metadata$Full, p.adjust.method = "BH")

#####LSD#####
TestRes=vector()
model <- aov(TestData[,1] ~ Type, data = TestData)
TukeyHSD(model)
out.none <- LSD.test(model, "Type", p.adj = "none")
out<- LSD.test(model, "Type", p.adj = "BH")

TestRes1=out[["groups"]][unique(metadata$Full),]
write.csv(TestRes1,"240721RedoxPlasticDBEachConditionBH.csv")
TestRes2=out.none[["groups"]][unique(metadata$Full),]
write.csv(TestRes2,"240721RedoxPlasticDBEachConditionNone.csv")
