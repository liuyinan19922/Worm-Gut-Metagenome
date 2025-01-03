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
library(ggprism)
library(grid)
library(eoffice)

# Set working directory
setwd("D:/230408 Metagenomic analysis/240311 Metagenome/240322ClusterprofilerKEGG/")

# Load pathway classification data
KeggList = read.csv("D:/230408 Metagenomic analysis/240311 Metagenome/pathwayClassification.csv")

# List all RData files in the directory
FileName1 = list.files(getwd(), pattern='*.RData')
FileName1

# Iterate through each RData file
for (i in 1:length(FileName1)) {

  # Load the data file
  load(FileName1[i])
  
  # Filter pathways with adjusted p-value < 0.05
  kk2_down@result = subset(kk2_down@result, p.adjust < 0.05)
  kk2_up@result = subset(kk2_up@result, p.adjust < 0.05)

  # Initialize vectors to store class and gene ratio
  ClassDown = vector()
  ClassUp = vector()
  a = vector()

  # Assign secondary categories and calculate GeneRatio for down-regulated pathways
  for (j in 1:nrow(kk2_down@result)) {
    ClassDown[j] = KeggList$Secondary.Category[which(KeggList$Name == kk2_down@result$Description[j])]
    a[j] = parse(text = kk2_down@result[["GeneRatio"]][j])
    kk2_down@result$GeneRatio[j] = eval(a[j])
  }

  # Assign secondary categories and calculate GeneRatio for up-regulated pathways
  for (j in 1:nrow(kk2_up@result)) {
    ClassUp[j] = KeggList$Secondary.Category[which(KeggList$Name == kk2_up@result$Description[j])]
    a[j] = parse(text = kk2_up@result[["GeneRatio"]][j])
    kk2_up@result$GeneRatio[j] = eval(a[j])
  }

  # Add category and GeneRatio columns to results
  kk2_down@result$ClassDown = ClassDown
  kk2_down@result$GeneRatio = as.numeric(kk2_down@result$GeneRatio)
  kk2_up@result$ClassUp = ClassUp
  kk2_up@result$GeneRatio = as.numeric(kk2_up@result$GeneRatio)

  # Filter and arrange down-regulated pathways
  kk2_down@result = subset(kk2_down@result, ClassDown != 'Global and overview maps')
  kk2_down@result = subset(kk2_down@result, ClassDown %in% c("Carbohydrate metabolism", "Glycan biosynthesis and metabolism", 'Nucleotide metabolism', 'Membrane transport', 'Metabolism of other amino acids'))
  kk2_down@result = arrange(kk2_down@result, GeneRatio)
  kk2_down@result = arrange(kk2_down@result, ClassDown)
  kk2_down@result$Description = factor(kk2_down@result$Description, levels = unique(kk2_down@result$Description))

  # Filter and arrange up-regulated pathways
  kk2_up@result = subset(kk2_up@result, ClassUp != 'Global and overview maps')
  kk2_up@result = subset(kk2_up@result, ClassUp %in% c('Xenobiotics biodegradation and metabolism', "Carbohydrate metabolism", 'Amino acid metabolism', 'Cellular community - prokaryotes', 'Lipid metabolism', 'Membrane transport'))
  kk2_up@result = arrange(kk2_up@result, GeneRatio)
  kk2_up@result = arrange(kk2_up@result, ClassUp)
  kk2_up@result$Description = factor(kk2_up@result$Description, levels = unique(kk2_up@result$Description))

  # Remove specific pathways from the results
  kk2_down@result = subset(kk2_down@result, !(Description %in% c('Propanoate metabolism', 'Butanoate metabolism')))
  kk2_up@result = subset(kk2_up@result, !(Description %in% c('Fructose and mannose metabolism')))

  # Generate plot for up-regulated pathways
  pUp = ggplot(kk2_up@result, aes(y = GeneRatio, x = Description, fill = p.adjust)) +
    ggtitle(Outputname[[i]]) +
    theme_prism(border = TRUE) +
    theme(axis.text.y = element_text(face = "plain", size = 12, lineheight = 0.6), legend.title = element_text(size = 12)) +
    geom_bar(stat = "identity", width = 0.6) + coord_flip() +
    xlab('KEGG pathway') +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 35)) +
    scale_fill_gradient(low = '#fdb7ba', high = "#920079", name = 'P.adj', limits = c(0, 0.05), breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05))

  # Save up-regulated plot
  topptx(pUp, paste0(Outputname[[i]], '_Up.pptx'), height = 6, width = 7.5)

  # Generate plot for down-regulated pathways
  pDown = ggplot(kk2_down@result, aes(y = GeneRatio, x = Description, fill = p.adjust)) +
    ggtitle(Outputname[[i]]) +
    theme_prism(border = TRUE) +
    theme(axis.text.y = element_text(face = "plain", size = 12, lineheight = 0.6), legend.title = element_text(size = 12)) +
    geom_bar(stat = "identity", width = 0.6) + coord_flip() +
    xlab('KEGG pathway') +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 35)) +
    scale_fill_gradient(low = '#92c3de', high = "#14549f", name = 'P.adj', limits = c(0, 0.05), breaks = c(0, 0.01, 0.02, 0.03, 0.04, 0.05))

  # Save down-regulated plot
  topptx(pDown, paste0(Outputname[[i]], '_Down.pptx'), height = 6, width = 7.5)
}
