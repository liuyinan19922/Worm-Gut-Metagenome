# Load required libraries
library(tidyfst)
library(eoffice)
library(ggrepel)
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(grid)
library(ggClusterNet)
library(plyr)
library(pulsar)
library(devtools)
library(EasyStat)

# Set working directory
setwd("D:/230408 Metagenomic analysis/240311 Metagenome/240513 bin network analysis")

# Define path to study samples
Samplestudy <- "D:/230408 Metagenomic analysis/240311 Metagenome/240513 bin network analysis/240529NetworkAnalysis"
list.files(Samplestudy)  # List files in the study directory

# Read OTU and taxonomy tables
otu.16s <- read.csv('CpmMediumBins.csv', header = TRUE)
tax.16s <- read.csv('TaxonTableMediumBins.csv', header = TRUE)

# Set row names for OTU and taxonomy tables
otu.16s <- column_to_rownames(otu.16s, "X")
tax.16s <- column_to_rownames(tax.16s, "X")

# Read sample metadata
sample.16s <- read.csv('metaForMetagenomeUpdated.csv', header = TRUE)
rownames(sample.16s) <- sample.16s$SampleName  # Set row names to SampleName column

# Convert OTU and taxonomy dataframes to matrices
otu.mat.16s <- as.matrix(otu.16s)
tax.mat.16s <- as.matrix(tax.16s)

# Create phyloseq objects for OTU table, taxonomy table, and sample data
asv.16s <- otu_table(otu.mat.16s, taxa_are_rows = TRUE)
classification.16s <- tax_table(tax.mat.16s)
sample.info.16s <- sample_data(sample.16s)

# Create a phyloseq object combining OTU table, taxonomy table, and sample data
phy.16s <- phyloseq(asv.16s, classification.16s, sample.info.16s)

# Sort and filter phyloseq object by sample sums
sort(sample_sums(phy.16s))  # Sort samples by total abundance
phy.16s.filt <- prune_taxa(taxa_sums(phy.16s) > 1, phy.16s)  # Filter out taxa with abundance <= 1

# Subset samples based on type and group
TotalWOC <- subset_samples(phy.16s, Type != "Control")  # Subset incubated samples
colnames(TotalWOC@sam_data)[5] <- 'OldGroup'
colnames(TotalWOC@sam_data)[3] <- 'Group'

# Network stability analysis and module comparison
tab.r <- network.pip(
  ps = TotalWOC,
  N = 500,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.8,
  p.threshold = 0.05,
  maxnode = 5,
  method = "spearman",
  label = FALSE,
  lab = "elements",
  group = "Group",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R = 10,
  ncpus = 1
)

# Save results to avoid recomputing
# saveRDS(tab.r,"network.pip.sparcc.rds")
# tab.r <- readRDS("./network.pip.sparcc.rds")

# Extract correlation matrix data
dat <- tab.r[[2]]
node <- dat$net.cor.matrix$node
edge <- dat$net.cor.matrix$edge
cortab <- dat$net.cor.matrix$cortab

# Compare modules using external data (e.g., pazy data)
pazy <- read.csv("D:/230408 Metagenomic analysis/240311 Metagenome/240715 pazy的显著性检验 含有pazy的新菌株/240715allPazy.csv")
pazy_mod <- merge(res[[2]], pazy, by.x = 'ID', by.y = 'X')
pazy_mod$sumAno <- pazy_mod$Annotated_hydrolase + pazy_mod$Annotated_oxidase
pazy_mod$sumPre <- pazy_mod$Predicted_oxidase + pazy_mod$Predicted_hydrolase
pazy_mod <- merge(pazy_mod, tax.16s[, c(6, 7)], by.x = "ID", by.y = 0)
pazy_mod <- pazy_mod[, -c(3, 4, 5, 6, 7)]

# Plot scatter for different groups
pazy_mod_ino <- subset(pazy_mod, (group %in% c("Inoculummodel_1", "Inoculummodel_2", "Inoculummodel_5")))
ggplot(pazy_mod_ino, aes(x = sumAno, y = sumPre, color = group)) + 
  geom_point(size = 4, position = position_jitter(width = 0.2, height = 0.2))

pazy_mod_bio <- subset(pazy_mod, (group %in% c("Biofilmmodel_1", "Biofilmmodel_2", "Biofilmmodel_3", 'Biofilmmodel_4', "Biofilmmodel_5", "Biofilmmodel_6")))
ggplot(pazy_mod_bio, aes(x = sumAno, y = sumPre, color = group)) + 
  geom_point(size = 4, position = position_jitter(width = 0.2, height = 0.2))

pazy_mod_pla <- subset(pazy_mod, (group %in% c("Planktonicmodel_1", "Planktonicmodel_4", "Planktonicmodel_5", 'Planktonicmodel_6')))
ggplot(pazy_mod_pla, aes(x = sumAno, y = sumPre, color = group)) + 
  geom_point(size = 4, position = position_jitter(width = 0.2, height = 0.2))

# Compare modules
res <- module.compare.m(
  corg = cortab,
  Top = 500,
  degree = TRUE,
  zipi = FALSE,
  r.threshold = 0.6,
  p.threshold = 0.05,
  method = "spearman",
  padj = FALSE,
  n = 3
)

# Visualize module similarity
p1 <- res[[1]]
p1

# Extract OTU and module info for the comparison
dat1 <- res[[2]]
head(dat1)

# Save results to PowerPoint and CSV
topptx(p1, '24061308cpmmodulesimilarityBFPK.pptx', width = 5, height = 4)
write.csv(res[[2]], '24061308cpmmoduleBFPK.csv')

# Create bar plot of similar modules
dat2 <- res[[3]]
dat2$m1 <- dat2$module1 %>% strsplit("model") %>%
  sapply(`[`, 1)
dat2$m2 <- dat2$module2 %>% strsplit("model") %>%
  sapply(`[`, 1)
dat2$cross <- paste(dat2$m1, dat2$m2, sep = "_Vs_")
dat2 <- dat2 %>% filter(module1 != "none")

p2 <- ggplot(dat2) + 
  geom_bar(aes(x = cross, fill = cross)) +
  labs(x = "", y = "Numbers of Similar Modules") + 
  theme_classic()

# Save the bar plot to PowerPoint and CSV
p2
topptx(p2, '24061308cpmNumSimilarModuleBFPK.pptx', width = 5, height = 4)
write.csv(dat2, '24061308cpmNumSimilarModuleBFPK.csv')

# -- Random Node Removal - Network Robustness ------------------------
res = Robustness.Random.removal(ps = TotalWOC,
                                Top = 500,
                                r.threshold = 0.6,
                                p.threshold = 0.05,
                                method = "spearman")
p = res[[1]]
p
topptx(p, '24061308cpmRubustnessRandomBFPK.pptx', width = 10, height = 4)
write.csv(res[[2]], '24061308cpmRubustnessRandomBFPK.csv')


# -- Key Node Removal - Network Robustness -------------------------
res2 = Robustness.Targeted.removal(ps = TotalWOC,
                                   corg = cortab,
                                   degree = TRUE,
                                   zipi = FALSE)
p3 = res2[[1]]
p3
# Extract data
dat4 = res2[[2]]
topptx(p3, '24061308cpmRubustnessKeyBFPK.pptx', width = 10, height = 4)
write.csv(res2[[2]], '24061308cpmRubustnessKeyBFPK.csv')


# -- Calculate Proportion of Negative Correlations ------------------
res4 = negative.correlation.ratio(ps = TotalWOC,
                                  corg = cortab,
                                  degree = TRUE,
                                  zipi = FALSE)

p5 = res4[[1]]
library(ggprism)
p51 = p5 + theme_prism(border = TRUE) +
  labs(x = "Network", y = 'Ratio of Negative Correlations') +
  scale_x_discrete(limit = c("Inoculum", 'Biofilm', 'Planktonic'), breaks = c("Inoculum", 'Biofilm', 'Planktonic'))
dat6 = res4[[2]]
# Negative correlation ratio data
head(dat6)
topptx(p51, '24071208cpmNegEdgeRatioBFPK.pptx', width = 5, height = 4)
write.csv(dat6, '24061308cpmNegEdgeRatioBFPK.csv')


# -- Network Significance Analysis ---------------------------------
dat = module.compare.net.pip(
  ps = NULL,
  corg = cortab,
  degree = TRUE,
  zipi = FALSE,
  r.threshold = 0.6,
  p.threshold = 0.05,
  method = "spearman",
  padj = FALSE,
  n = 2)
res = dat[[1]]
head(res)
write.csv(res, '24061308cpmSignificantTestBFPK.csv')


# -- Network Resilience --------------------------------------------
res6 = natural.con.microp(
  ps = TotalWOC,
  corg = cortab,
  norm = TRUE,
  end = 80, # Less than the number of nodes in the network
  start = 0
)
p7 = res6[[1]]
p7
dat8 = res6[[2]]
topptx(p7, '24061308cpmResistenceBFPK.pptx', width = 5, height = 4)
write.csv(dat8, '24061308cpmResistenceBFPK.csv')


# -- Network Modularity Analysis Based on Output Matrix -----------
select.mod = "no"  # Option to select which modules to display
plots = list()  # Save network modularity results
plots2 = list()  # Save modularity results
alldat = list()  # Store module data
plots3 = list()  # Store species composition of modules
plots4 = list()  # Store module diversity
alldat2 = list()  # Store species composition data
alldat3 = list()  # Store module diversity data

for (i in 1:length(names(cortab))) {
  # Specify modules for display
  select.mod = switch(i, c("model_1", "model_2", "model_5"), c("model_1", "model_2", "model_3", "model_4", "model_5", "model_6"), c("model_1", "model_4", "model_5", "model_6"))
  
  resu = module_display.2(
    pst = TotalWOC,
    method.clu = 'cluster_fast_greedy',
    corg = cortab[[names(cortab)[i]]],
    select.mod = 'no',  # Select specific modules for visualization
    num = 5,  # Do not display modules with fewer than 5 OTUs
    leg.col = 3
  )
  
  # Display all modules
  p1 = resu[[1]] + labs(title = names(cortab)[i]) +
    theme(
      plot.title = element_text(hjust = 0.5),  # Center title
      plot.margin = unit(c(3, 3, 4, 4), "lines")  # Set plot margin
    )
  p1
  
  # Extract module information: OTUs and network edge count for each module
  dat = resu[[5]]
  head(dat)
  
  # Extract network edge and node tables
  dat2 = resu[[6]]
  
  # Display network with modules colored, edges between modules in gray
  p2 = resu[[3]] + labs(title = names(cortab)[i]) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  p2
  
  alldat[[names(cortab)[i]]] = dat2
  plots[[names(cortab)[i]]] = p1
  plots2[[names(cortab)[i]]] = p2
  
  # Extract module information, including OTUs and module classification
  select.mod = switch(i, c("model_1", "model_2", "model_5"), c("model_1", "model_2", "model_3", "model_4", "model_5", "model_6"), c("model_1", "model_4", "model_5", "model_6"))
}

# Filter the modules based on selected ones
mod1 = resu$mod.groups %>%
  filter(group %in% select.mod)
head(mod1)

# Calculate species composition for the selected module
pst = TotalWOC %>%
  subset_samples.wt("Group", names(cortab)[i]) %>%
  subset_taxa.wt("OTU", colnames(cortab[[names(cortab)[i]]])) %>%
  filter_taxa(function(x) sum(x) > 0, TRUE) %>%
  scale_micro("rela")

# Generate the microbial composition at the genus level
res = module_composition1(pst = pst, mod1 = mod1, j = "genus")
p3 = res[[1]]
p3

# Define color palette for the plot
Yanse = c('#00B8DF','#FF6C65','#00B800','#9ed2ff','#feb1ad','#9dc39d','#D39100','#619BFA','#DB72F6','#c6c6c6','#dac597','#a6b7fb','#e8b2f7','#156f82','#a2524f','#177817','#765714','#4e6da1','#9c5eac')

# Customize the plot for the microbial composition
P31 = p3[[2]] +
  scale_fill_manual(values = Yanse) +
  theme_prism(border = TRUE) +
  scale_x_discrete(labels = c("Module 1", 'Module 2', 'Module 3', "Module 4", 'Module 5', 'Module 6')) +
  theme(axis.text.x = element_text(angle = 330, hjust = 0.3))

P31

# Save plot and export data
topptx(P31, '24072608cpmCompositionBIofilm.pptx', width = 6.5, height = 4)
write.csv(P31[["data"]], '240726ModuleCompositionBIofilm.csv')

plots3[[names(cortab)[i]]] = p3

# Extract and save the species composition data
ps.t = res[[3]]

otu = ps.t %>% vegan_otu() %>% t() %>% as.data.frame()
tax = ps.t %>% vegan_tax() %>% as.data.frame()
tab.res = cbind(otu, tax)
alldat2[[names(cortab)[i]]] = tab.res

# Save standardized abundance data
tab.res.3 = res[[4]]$relaabundance
write.csv(tab.res.3, '240726BFModuleComposition.csv')

# Retrieve module information for the selected module
mod1 = resu$mod.groups %>%
  dplyr::select(ID, group) %>%
  dplyr::filter(group %in% select.mod)
head(mod1)

# Calculate module alpha diversity
result = module_alpha(ps = TotalWOC, mod1 = mod1)
p5 = result[[1]]
p5

plots4[[names(cortab)[i]]] = p5

# Extract alpha diversity data
plotd = result[[4]]$alpha
alldat3[[names(cortab)[i]]] = plotd

# Save the alpha diversity significance table
sigtab = result[[4]]$sigtab

# If it's the last iteration, save the results
if (i == length(names(cortab))) {
  saveRDS(alldat, paste0('24061308cpm', names(cortab)[i], 'ModuledataBFPK.rds'))
  write.csv(alldat2, paste0('24061308cpm', names(cortab)[i], 'ModuleCompositionBFPK.csv'))
  write.csv(alldat3, paste0('24061308cpm', names(cortab)[i], 'ModuleDiversityBFPK.csv'))
}

# Arrange the plots in grid format
p1 = ggpubr::ggarrange(plotlist = plots, common.legend = FALSE, legend = "right", ncol = 3, nrow = 1)
p1

p2 = ggpubr::ggarrange(plotlist = plots2, common.legend = FALSE, legend = "right", ncol = 3, nrow = 1)
p2

p3 = ggpubr::ggarrange(plotlist = plots3, common.legend = FALSE, legend = "right", ncol = 1, nrow = 3)
p3

p4 = ggpubr::ggarrange(plotlist = plots4, common.legend = FALSE, legend = "right", ncol = 1, nrow = 3)
p4

# Save the final plots
topptx(p1, '24061308cpmModuleBFPK.pptx', width = 25, height = 7)
topptx(p2, '24061308cpmModule2BFPK.pptx', width = 25, height = 7)
topptx(p3, '24061308cpmModuleCompositionBFPK.pptx', width = 15, height = 15)
topptx(p4, '24061308cpmModuleDiversityBFPK.pptx', width = 15, height = 15)

# Module composition function
module_composition1 = function(pst = pst, mod1 = mod1, j = "Family") {
  tem = mod1$group %>% table() %>% as.data.frame() %>% dplyr::arrange(desc(Freq))
  colnames(tem) = c("Model", "OTU.num")
  otu = NULL
  map = NULL
  for (i in 1:length(tem$Model)) {
    id.s = tem$Model %>% as.character()
    id.t = mod1 %>% filter(group %in% id.s[i]) %>% .$ID
    ps.tem = subset_taxa.wt(pst, "OTU", id.t)
    otu = ps.tem %>% vegan_otu() %>% t() %>% as.data.frame()
    colnames(otu) = paste(id.s[i], colnames(otu), sep = "_")
    map = data.frame(row.names = colnames(otu), ID = colnames(otu), Group = id.s[i])
    if (i == 1) {
      otu.f = otu
      map.f = map
    } else {
      otu$ID = row.names(otu)
      otu.f$ID = row.names(otu.f)
      tem.2 = otu.f %>% full_join(otu) %>% select(ID, everything()) %>% as.data.frame()
      row.names(tem.2) = tem.2$ID
      tem.2$ID = NULL
      tem.2[is.na(tem.2)] = 0
      otu.f = tem.2
      map.f = rbind(map.f, map)
    }
  }
  pst.2 = phyloseq(otu_table(as.matrix(otu.f), taxa_are_rows = TRUE), sample_data(map.f), tax_table(pst))
  result = barMainplot(ps = pst.2, tran = F, j = j, label = FALSE, sd = FALSE, Top = 15)
  p4_1 <- result[[1]] + scale_fill_hue() + theme_classic()
  p4_2 <- result[[3]] + scale_fill_hue() + theme_classic()
  tem1 = result[[2]]
  result = barMainplot(ps = pst.2, tran = T, j = j, label = FALSE, sd = FALSE, Top = 15)
  p3_1 <- result[[1]] + scale_fill_hue() + theme_classic()
  p3_2 <- result[[3]] + scale_fill_hue() + theme_classic()
  tem2 = result[[2]]
  library(patchwork)
  p00 = p4_1 | p3_1
  p01 = p4_2 | p3_2
  plotdat = list(bundance = tem1, relaabundance = tem2)
  return(list(p00, p01, pst.2, plotdat = plotdat))
}
