# Install required packages
install.packages("devtools")
install.packages("spaa")
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
devtools::install_github("GuillemSalazar/EcolUtils")

# Load required libraries
library(EcolUtils)
library(spaa)
library(devtools)
library(ggplot2)
library(vegan)
library(ggprism)
library(ggalt)
library(pairwiseAdonis)
library(readxl)
library(SQMtools)

# Load metadata and KEGG data
metadataForMetaseq <- read_excel("D:/230408 Metagenomic analysis/240311 metagenome/metaForMetagenomeUpdated.xlsx")
wormPlastic <- loadSQM("D:/230408 Metagenomic analysis/wormPlastic2/")

keggMetaseq <- wormPlastic[["functions"]][["KEGG"]][["tpm"]]
keggMetaseq <- t(keggMetaseq)

# Perform NMDS analysis
nMDSOri <- metaMDS(keggMetaseq, distance = "bray", k = 2)
nMDSOri$stress  # Check NMDS stress values
nMDSOriP <- nMDSOri$points

# Merge NMDS results with metadata
nMDSOriP <- merge(nMDSOriP, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")

# Define color palette
ColorContainGrey <- c('#8997ab', "#ff6c67", "#00bd0f", "#4c9eff")

# Create NMDS plot
nMDS_ori <- ggplot(nMDSOriP, aes(x = MDS1, y = MDS2, colour = Type, shape = Position)) +
  geom_point(size = 1.5, stroke = 2) +
  stat_ellipse(type = "t", level = 0.95, aes(linetype = Group, group = Group), color = 'grey60', lwd = 1) +
  scale_shape_manual(values = c(1, 17, 19), breaks = c("Inoculum", "Planktonic", 'Biofilm')) +
  scale_color_manual(values = ColorContainGrey) +
  theme_prism(border = TRUE) +
  theme(legend.text = element_text(size = 14)) +
  scale_fill_manual(values = ColorContainGrey) +
  geom_encircle(aes(fill = Type), alpha = 0.2, show.legend = FALSE) +
  guides(linetype = 'none') +
  ggtitle(paste0('KEGG ', round(nMDSOri$stress, 4), adonis_result[["aov.tab"]]$`Pr(>F)`[1]))

# Save NMDS plot
nMDS_ori
topptx(nMDS_ori, "240806nMDS_keeg_all.pptx", width = 5.625, height = 3.75)
ggsave(file = "nMDS_keeg_all.jpg", plot = nMDS_ori, dpi = 300, width = 5.625, height = 3.75)

# Perform PERMANOVA analysis
keggMetaseq0 <- merge(keggMetaseq, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")
keggMetaseq0 <- keggMetaseq0[, 1:10962]
row.names(keggMetaseq0) <- keggMetaseq0[, 1]
keggMetaseq0 <- keggMetaseq0[, -1]

group <- metadataForMetaseq[, c("SampleName", "Group")]
adonis_result <- adonis(keggMetaseq0 ~ Group, group, method = 'bray', permutations = 9999)
adonis_result

# Perform pairwise PERMANOVA
otu.pairwise.adonis <- pairwise.adonisKXW(x = keggMetaseq0, factors = group$Group, 
                                           sim.method = 'bray', p.adjust.m = 'BH')
otu.pairwise.adonis

# Save PERMANOVA results
prefix <- 'InoculumIncubatedKegg'
write.csv(adonis_result[["aov.tab"]], paste0(prefix, '.permanova.csv'))
write.csv(otu.pairwise.adonis, paste0(prefix, '.Pairwisepermanova.csv'))

# Perform ANOSIM analysis
anosim_result <- anosim(keggMetaseq0, group$Group, distance = 'bray', permutations = 999)
anosim_result

# Define function for pairwise.adonisKXW
pairwise.adonisKXW <- function(x, factors, sim.method, p.adjust.m) {
  library(vegan)
  co <- as.matrix(combn(unique(factors), 2))
  pairs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  for (elem in 1:ncol(co)) {
    ad <- adonis(x[factors %in% c(as.character(co[1, elem]), as.character(co[2, elem])), ] ~ 
                   factors[factors %in% c(as.character(co[1, elem]), as.character(co[2, elem]))], method = sim.method)
    pairs <- c(pairs, paste(co[1, elem], 'vs', co[2, elem]))
    F.Model <- c(F.Model, ad$aov.tab[1, 4])
    R2 <- c(R2, ad$aov.tab[1, 5])
    p.value <- c(p.value, ad$aov.tab[1, 6])
  }
  p.adjusted <- p.adjust(p.value, method = p.adjust.m)
  pairw.res <- data.frame(pairs, F.Model, R2, p.value, p.adjusted)
  return(pairw.res)
}

# Perform ANOSIM analysis
anosim_result <- anosim(keggOri, group$class, distance = 'bray', permutations = 999)
print(anosim_result)

# Perform PERMANOVA and pairwise analysis within inocula
keggMetaseq0 <- merge(keggMetaseq, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")
keggMetaseq0 <- subset(keggMetaseq0, Stage == "Inoculum")
group <- keggMetaseq0[, c("Row.names", "Ifp")]
keggMetaseq0 <- keggMetaseq0[, 1:10962]
row.names(keggMetaseq0) <- keggMetaseq0[, 1]
keggMetaseq0 <- keggMetaseq0[, -1]

adonis_result <- adonis(keggMetaseq0 ~ Ifp, group, method = 'bray', permutations = 9999)
print(adonis_result)

otu.pairwise.adonis <- pairwise.adonis(x = keggMetaseq0, factors = group$Ifp, sim.function = "vegdist", 
                                       sim.method = 'bray', p.adjust.m = 'BH', reduce = NULL, perm = 1000)
print(otu.pairwise.adonis)

prefix <- 'InoculumIfpKegg'
write.csv(adonis_result[["aov.tab"]], paste0(prefix, '.permanova.csv'))
write.csv(otu.pairwise.adonis, paste0(prefix, '.Pairwisepermanova.csv'))

# Define pairwise ANOSIM function
pairwise.anosim <- function(x, factors, sim.method, p.adjust.m) {
  library(vegan)
  co <- as.matrix(combn(unique(factors), 2))
  pairs <- c()
  R <- c()
  p.value <- c()
  
  for (elem in 1:ncol(co)) {
    ad <- anosim(x[factors %in% c(as.character(co[1, elem]), as.character(co[2, elem])), ],                
                factors[factors %in% c(as.character(co[1, elem]), as.character(co[2, elem]))], 
                permutations = 999, distance = "bray")
    pairs <- c(pairs, paste(co[1, elem], 'vs', co[2, elem]))    
    R <- c(R, ad$statistic)   
    p.value <- c(p.value, ad$signif)  
  }
  
  p.adjusted <- p.adjust(p.value, method = p.adjust.m) 
  pairw.res <- data.frame(pairs, R, p.value, p.adjusted) 
  return(pairw.res)
}

pairwise.anosim(keggMetaseq0, group$Full, sim.method = "bray", p.adjust.m = "BH")

# Perform NMDS for the "Late" stage
keggMetaseq0 <- merge(keggMetaseq, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")
keggMetaseq0 <- subset(keggMetaseq0, Stage == "Late")
keggMetaseq0 <- keggMetaseq0[, 1:10962]
row.names(keggMetaseq0) <- keggMetaseq0[, 1]
keggMetaseq0 <- keggMetaseq0[, -1]

nMDSOri <- metaMDS(keggMetaseq0, distance = "bray", k = 2)
print(nMDSOri$stress)

nMDSOriP <- nMDSOri$points
nMDSOriP <- merge(nMDSOriP, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")

ColorContainGrey <- c('#8997ab', "#ff6c67", "#00bd0f", "#4c9eff")

nMDS_ori <- ggplot(nMDSOriP, aes(x = MDS1, y = MDS2, colour = Type, shape = Position)) +
  geom_point(size = 1.5, stroke = 2) +
  scale_shape_manual(values = c(1, 17, 19), breaks = c("Source", "Planktonic", 'Biofilm')) +
  scale_color_manual(values = ColorContainGrey) +
  theme_prism(border = TRUE) +
  theme(legend.text = element_text(size = 14)) +
  scale_fill_manual(values = ColorContainGrey) +
  geom_encircle(aes(fill = Type), alpha = 0.2, show.legend = FALSE) +
  guides(linetype = 'none') +
  xlim(-0.45, 0.45) +
  ylim(-0.4, 0.4)

print(nMDS_ori)
ggsave(file = "nMDS_keeg_late.jpg", plot = nMDS_ori, dpi = 300, width = 5.625, height = 3.75)

# Test between biofilm and planktonic groups
keggMetaseq0 <- merge(keggMetaseq, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")
keggMetaseq0 <- subset(keggMetaseq0, Stage == "Late")
group <- keggMetaseq0[, c("Row.names", "Type")]
keggMetaseq0 <- keggMetaseq0[, 1:10962]
row.names(keggMetaseq0) <- keggMetaseq0[, 1]
keggMetaseq0 <- keggMetaseq0[, -1]

adonis_result <- adonis(keggMetaseq0 ~ Type, group, method = 'bray', permutations = 9999)
print(adonis_result)

otu.pairwise.adonis <- pairwise.adonis(x = keggMetaseq0, factors = group$Type, sim.function = "vegdist", 
                                       sim.method = 'bray', p.adjust.m = 'BH', reduce = NULL, perm = 1000)
print(otu.pairwise.adonis)

prefix <- 'LateTypeKegg'
write.csv(adonis_result[["aov.tab"]], paste0(prefix, '.permanova.csv'))
write.csv(otu.pairwise.adonis, paste0(prefix, '.Pairwisepermanova.csv'))

# Perform NMDS for "Inoculum" stage
keggMetaseq0 <- merge(keggMetaseq, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")
keggMetaseq0 <- subset(keggMetaseq0, Stage == "Inoculum")
keggMetaseq0 <- keggMetaseq0[, 1:10962]
row.names(keggMetaseq0) <- keggMetaseq0[, 1]
keggMetaseq0 <- keggMetaseq0[, -1]

nMDSOri <- metaMDS(keggMetaseq0, distance = "bray", k = 2)
print(nMDSOri$stress)

nMDSOriP <- nMDSOri$points
nMDSOriP <- merge(nMDSOriP, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")

nMDS_ori <- ggplot(nMDSOriP, aes(x = MDS1, y = MDS2, colour = Type, shape = Position)) +
  geom_point(size = 1.5, stroke = 2) +
  scale_shape_manual(values = c(1, 17, 19), breaks = c("Inoculum", "Planktonic", 'Biofilm')) +
  scale_color_manual(values = ColorContainGrey) +
  theme_prism(border = TRUE) +
  theme(legend.text = element_text(size = 14)) +
  scale_fill_manual(values = ColorContainGrey) +
  geom_encircle(aes(fill = Type), alpha = 0.2, show.legend = FALSE) +
  guides(linetype = 'none')

print(nMDS_ori)
ggsave(file = "nMDS_keeg_source.jpg", plot = nMDS_ori, dpi = 300, width = 5.625, height = 3.75)

# Set color schemes for comparison among samples
prefix <- "HDPE"

color <- list()
color[[1]] <- c("#000000", '#8296A0')
color[[2]] <- c("#9A1207", '#FF6C67', "#FBB0AA")
color[[3]] <- c("#287B00", '#00BD0F', "#6CFF26")
color[[4]] <- c("#193754", '#4C9EFF', "#7DDBFF")
names(color) <- c("Control", "HDPE", 'PP', "PS")

colorUsed <- color[[which(names(color) == prefix)]]

##### Taxonomy Analysis #####

# NMDS Analysis for All Samples
taxonMetaseq <- wormPlastic[["taxa"]][["genus"]][["percent"]]
taxonMetaseq <- t(taxonMetaseq)
taxonMetaseq[,601] <- taxonMetaseq[,601] + taxonMetaseq[,817]
taxonMetaseq <- taxonMetaseq[,-817]

taxonMetaseq0 <- merge(taxonMetaseq, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")
taxonMetaseq0 <- taxonMetaseq0[,1:847]
row.names(taxonMetaseq0) <- taxonMetaseq0[,1]
taxonMetaseq0 <- taxonMetaseq0[,-1]

nMDSTaxon <- metaMDS(taxonMetaseq0, distance = "bray", k = 2)
nMDSTaxon$stress
scores(nMDSTaxon)
nMDSTaxonP <- nMDSTaxon$points

group <- metadataForMetaseq[,c("SampleName", "Type")]
adonis_result <- adonis(taxonMetaseq0 ~ Type, group, method = 'bray', permutations = 9999)
adonis_result["aov.tab"]
nMDSTaxonP <- merge(nMDSTaxonP, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")

# Color configuration
ColorContainGrey <- c('#8997ab', "#ff6c67", "#00bd0f", "#4c9eff")

# Plot NMDS
title <- paste0('PlasticDB ', round(nMDSTaxon$stress, 4), adonis_result[["aov.tab"]]$`Pr(>F)`[1])
nMDS_ori <- ggplot(nMDSTaxonP, aes(x = MDS1, y = MDS2, colour = Type, shape = Position)) +
  geom_point(size = 1.5, stroke = 2) +
  stat_ellipse(type = "t", level = 0.95, aes(x = MDS1, y = MDS2, linetype = Group, group = Group), color = 'grey60', lwd = 1) +
  scale_shape_manual(values = c(1, 17, 19), breaks = c("Inoculum", "Planktonic", 'Biofilm')) +
  scale_color_manual(values = ColorContainGrey) +
  theme_prism(border = TRUE) +
  theme(legend.text = element_text(size = 14)) +
  scale_fill_manual(values = ColorContainGrey) +
  geom_encircle(aes(fill = Type), alpha = 0.2, show.legend = FALSE) +
  guides(linetype = 'none') +
  ggtitle(title)
nMDS_ori

# Save results
output_name <- '240806nMDS_genus_all'
topptx(nMDS_ori, paste0(output_name, '.pptx'), width = 5.625, height = 3.75)
ggsave(file = paste0(output_name, '.jpg'), plot = nMDS_ori, dpi = 300, width = 5.625, height = 3.75)

# PERMANOVA analysis for different groups
group <- metadataForMetaseq[, c("SampleName", "Group")]
adonis_result <- adonis(taxonMetaseq0 ~ Group, group, method = 'bray', permutations = 9999)
adonis_result["aov.tab"]

otu.pairwise.adonis <- pairwise.adonisKXW(x = taxonMetaseq0, factors = group$Group, sim.method = 'bray', p.adjust.m = 'BH')
prefix <- 'InoculumGenus'
write.csv(adonis_result[["aov.tab"]], paste0(prefix, '.permanova.csv'))
write.csv(otu.pairwise.adonis, paste0(prefix, '.Pairwisepermanova.csv'))

# NMDS Analysis Between Biofilm and Planktonic
taxonMetaseq <- wormPlastic[["taxa"]][["genus"]][["percent"]]
taxonMetaseq <- t(taxonMetaseq)
taxonMetaseq[,601] <- taxonMetaseq[,601] + taxonMetaseq[,817]
taxonMetaseq <- taxonMetaseq[,-817]

taxonMetaseq0 <- merge(taxonMetaseq, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")
taxonMetaseq0 <- subset(taxonMetaseq0, Stage == "Late")
taxonMetaseq0 <- taxonMetaseq0[,1:847]
row.names(taxonMetaseq0) <- taxonMetaseq0[,1]
taxonMetaseq0 <- taxonMetaseq0[,-1]

group <- metadataForMetaseq[, c("SampleName", "Type")]
adonis_result <- adonis(taxonMetaseq0 ~ Type, group, method = 'bray', permutations = 9999)
adonis_result["aov.tab"]

otu.pairwise.adonis <- pairwise.adonisKXW(x = taxonMetaseq0, factors = group$Group, sim.method = 'bray', p.adjust.m = 'BH')
prefix <- 'LateTypeGenus'
write.csv(adonis_result[["aov.tab"]], paste0(prefix, '.permanova.csv'))
write.csv(otu.pairwise.adonis, paste0(prefix, '.Pairwisepermanova.csv'))

# Perform NMDS analysis for taxonomic data
nMDSTaxon <- metaMDS(taxonMetaseq0, distance = "bray", k = 2) # Bray-Curtis dissimilarity
nMDSTaxon$stress
scores(nMDSTaxon)
nMDSTaxonP = nMDSTaxon$points
nMDSTaxonP = merge(nMDSTaxonP, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")

# Define colors for visualization
ColorContainGrey = c('#8997ab', "#ff6c67", "#00bd0f", "#4c9eff")

# Plot NMDS results
nMDS_ori = ggplot(nMDSTaxonP, aes(x = MDS1, y = MDS2, colour = Type, shape = Position)) +
  geom_point(size = 1.5, stroke = 2) +
  scale_shape_manual(values = c(1, 17, 19), breaks = c("Inoculum", "Planktonic", "Biofilm")) +
  scale_color_manual(values = ColorContainGrey) +
  theme_prism(border = TRUE) +
  theme(legend.text = element_text(size = 14)) +
  scale_fill_manual(values = ColorContainGrey) +
  geom_encircle(aes(fill = Type), alpha = 0.2, show.legend = FALSE) +
  guides(linetype = 'none')

nMDS_ori
ggsave(file = "nMDS_genus_late.jpg", plot = nMDS_ori, dpi = 300, width = 5.625, height = 3.75)

# Perform analysis for the inoculum stage
taxonMetaseq = wormPlastic[["taxa"]][["genus"]][["percent"]]
taxonMetaseq = t(taxonMetaseq)
taxonMetaseq[, 601] = taxonMetaseq[, 601] + taxonMetaseq[, 817]
taxonMetaseq = taxonMetaseq[, -817]

taxonMetaseq0 = merge(taxonMetaseq, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")
taxonMetaseq0 = subset(taxonMetaseq0, Stage == "Inoculum")
taxonMetaseq0 = taxonMetaseq0[, 1:847]
row.names(taxonMetaseq0) = taxonMetaseq0[, 1]
taxonMetaseq0 = taxonMetaseq0[, -1]

group = metadataForMetaseq[, c("SampleName", "Ifp")]
adonis_result = adonis(taxonMetaseq0 ~ Ifp, group, method = 'bray', permutations = 9999) # PERMANOVA
adonis_result

# Perform pairwise adonis analysis
otu.pairwise.adonis = pairwise.adonisKXW(
  x = taxonMetaseq0,
  factors = group$Group,
  sim.method = 'bray',
  p.adjust.m = 'BH'
)
otu.pairwise.adonis

# Save results
prefix = 'InoculumIfpGenus'
write.csv(adonis_result[["aov.tab"]], paste0(prefix, '.permanova.csv'))
write.csv(otu.pairwise.adonis, paste0(prefix, '.Pairwisepermanova.csv'))

# Repeat NMDS and visualization for COG data
cogMetaseq = t(wormPlastic[["functions"]][["COG"]][["tpm"]])
cogMetaseq0 = merge(cogMetaseq, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")
cogMetaseq0 = cogMetaseq0[, 1:33536]
row.names(cogMetaseq0) = cogMetaseq0[, 1]
cogMetaseq0 = cogMetaseq0[, -1]

group = metadataForMetaseq[, c("SampleName", "Group")]
adonis_result = adonis(cogMetaseq0 ~ Group, group, method = 'bray', permutations = 9999)
adonis_result

otu.pairwise.adonis = pairwise.adonisKXW(
  x = cogMetaseq0,
  factors = group$Group,
  sim.method = 'bray',
  p.adjust.m = 'BH'
)
otu.pairwise.adonis

prefix = 'AllGroupCog'
write.csv(adonis_result[["aov.tab"]], paste0(prefix, '.permanova.csv'))
write.csv(otu.pairwise.adonis, paste0(prefix, '.Pairwisepermanova.csv'))

nMDSTaxon <- metaMDS(cogMetaseq, distance = "bray", k = 2)
nMDSTaxon$stress
scores(nMDSTaxon)
nMDSTaxonP = nMDSTaxon$points
nMDSTaxonP = merge(nMDSTaxonP, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")

nMDS_ori = ggplot(nMDSTaxonP, aes(x = MDS1, y = MDS2, colour = Type, shape = Position)) +
  geom_point(size = 1.5, stroke = 2) +
  stat_ellipse(type = "t", level = 0.95, aes(x = MDS1, y = MDS2, linetype = Group, group = Group), color = 'grey60', lwd = 1) +
  scale_shape_manual(values = c(1, 17, 19), breaks = c("Inoculum", "Planktonic", "Biofilm")) +
  scale_color_manual(values = ColorContainGrey) +
  theme_prism(border = TRUE) +
  theme(legend.text = element_text(size = 14)) +
  scale_fill_manual(values = ColorContainGrey) +
  geom_encircle(aes(fill = Type), alpha = 0.2, show.legend = FALSE) +
  guides(linetype = 'none') +
  ggtitle(paste0('COG ', round(nMDSTaxon$stress, 4), adonis_result[["aov.tab"]]$`Pr(>F)`[1]))

nMDS_ori
topptx(nMDS_ori, "240806nMDS_cog_all.pptx", width = 5.625, height = 3.75)
ggsave(file = "nMDS_cog_all.jpg", plot = nMDS_ori, dpi = 300, width = 5.625, height = 3.75)

# NMDS Analysis for COG and PFAM Data

# Load COG data and preprocess for NMDS
cogMetaseq = t(wormPlastic[["functions"]][["COG"]][["tpm"]])
cogMetaseq0 = merge(cogMetaseq, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")
cogMetaseq0 = subset(cogMetaseq0, Stage == "Late")
cogMetaseq0 = cogMetaseq0[, 1:33536]
row.names(cogMetaseq0) = cogMetaseq0[, 1]
cogMetaseq0 = cogMetaseq0[, -1]

# Group metadata for analysis
group = metadataForMetaseq[, c("SampleName", "Type")]

# PERMANOVA Analysis
adonis_result = adonis(cogMetaseq0 ~ Type, group, method = 'bray', permutations = 9999)
adonis_result

# Pairwise PERMANOVA Analysis
otu.pairwise.adonis = pairwise.adonisKXW(
  x = cogMetaseq0,
  factors = group$Type,
  sim.method = 'bray',
  p.adjust.m = 'BH'
)
write.csv(adonis_result[["aov.tab"]], "LateTypeCog.permanova.csv")
write.csv(otu.pairwise.adonis, "LateTypeCog.Pairwisepermanova.csv")

# NMDS Visualization for COG Data
nMDSTaxon = metaMDS(cogMetaseq0, distance = "bray", k = 2)
colors = c('#8997ab', "#ff6c67", "#00bd0f", "#4c9eff")

data_for_plot = merge(nMDSTaxon$points, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")

nMDS_ori = ggplot(data_for_plot, aes(x = MDS1, y = MDS2, colour = Type, shape = Position)) +
  geom_point(size = 1.5, stroke = 2) +
  scale_shape_manual(values = c(1, 17, 19), breaks = c("Inoculum", "Planktonic", "Biofilm")) +
  scale_color_manual(values = colors) +
  theme_prism(border = TRUE) +
  theme(legend.text = element_text(size = 14)) +
  geom_encircle(aes(fill = Type), alpha = 0.2, show.legend = FALSE) +
  guides(linetype = 'none') +
  ylim(-0.35, 0.35)

nMDS_ori
ggsave("nMDS_cog_late.jpg", plot = nMDS_ori, dpi = 300, width = 5.625, height = 3.75)

# Repeat for Inoculum Stage
cogMetaseq0 = subset(cogMetaseq0, Stage == "Inoculum")

adonis_result = adonis(cogMetaseq0 ~ Ifp, group, method = 'bray', permutations = 9999)
adonis_result
otu.pairwise.adonis = pairwise.adonisKXW(
  x = cogMetaseq0,
  factors = group$Ifp,
  sim.method = 'bray',
  p.adjust.m = 'BH'
)
write.csv(adonis_result[["aov.tab"]], "InoculumIfpCog.permanova.csv")
write.csv(otu.pairwise.adonis, "InoculumIfpCog.Pairwisepermanova.csv")

nMDSTaxon = metaMDS(cogMetaseq0, distance = "bray", k = 2)
data_for_plot = merge(nMDSTaxon$points, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")

nMDS_ori = ggplot(data_for_plot, aes(x = MDS1, y = MDS2, colour = Type, shape = Position)) +
  geom_point(size = 1.5, stroke = 2) +
  scale_shape_manual(values = c(1, 17, 19), breaks = c("Inoculum", "Planktonic", "Biofilm")) +
  scale_color_manual(values = colors) +
  theme_prism(border = TRUE) +
  theme(legend.text = element_text(size = 14)) +
  geom_encircle(aes(fill = Type), alpha = 0.2, show.legend = FALSE) +
  xlim(-0.8, 0.8) +
  ylim(-0.5, 0.6)

nMDS_ori
ggsave("nMDS_cog_inocula.jpg", plot = nMDS_ori, dpi = 300, width = 5.625, height = 3.75)

# NMDS Analysis for PFAM
pfamMetaseq = t(wormPlastic[["functions"]][["PFAM"]][["tpm"]])
nMDSPfam = metaMDS(pfamMetaseq, distance = "bray", k = 2)

pfam_plot_data = merge(nMDSPfam$points, metadataForMetaseq, by.x = "row.names", by.y = "SampleName")

nMDS_Pfam = ggplot(pfam_plot_data, aes(x = MDS1, y = MDS2, colour = Type, shape = Position)) +
  geom_point(size = 3, stroke = 2) +
  stat_ellipse(type = "t", level = 0.68, aes(linetype = Group, group = Group), color = 'grey60', lwd = 1) +
  scale_shape_manual(values = c(1, 17, 19), breaks = c("Source", "Planktonic", "Biofilm")) +
  scale_color_manual(values = colors) +
  theme_prism(border = TRUE) +
  theme(legend.text = element_text(size = 14)) +
  guides(linetype = 'none')

nMDS_Pfam
ggsave("nMDS_pfam.jpg", plot = nMDS_Pfam, dpi = 300, width = 5.625, height = 3.75)

