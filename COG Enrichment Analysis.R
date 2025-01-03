# Install and load necessary packages
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("MicrobiomeProfiler")
library(enrichplot)
library(ggplot2)
library(MicrobiomeProfiler)
library(ggprism)
library(dplyr)
library(stringr)
library(eoffice)

# Set working directory and get file names
setwd("D:/230408 Metagenomic analysis/240311 Metagenome/240311Deseq2/cog/cog enrichment")
FileName = list.files(getwd(), pattern = '*.csv')
Outputname = lapply(strsplit(FileName, ".csv"), "[", 1)

# Loop through files and perform COG enrichment analysis
for (i in 1:26) {
  # For testing, set i = 12
  i = 12
  FileName[i]  
  
  cog_deseq <- read.csv(FileName[i], header = TRUE)
  cog_name = cog_deseq$X
  
  # Remove "Unclassified" entries from COG names
  if ("Unclassified" %in% cog_name) {
    cog_name = cog_name[-which(cog_name == "Unclassified")]
  }

  # Perform COG enrichment analysis
  cog_enrichment = enrichCOG(
    cog_name,
    dtype = "category",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2
  )
  
  # Plot results
  dotplot(cog_enrichment, showCategory = 20) + 
    ggtitle(Outputname[[i]]) + 
    theme_prism(border = TRUE)
  # Save results (optional)
  # write.csv(cog_enrichment@result, paste0("CogEnrichment_", FileName[i]))
  # save.image(paste0(Outputname[[i]], '.RData'))
}

# Combine graphs for pairs of files
for (i in 1:13) {
  FileName[2 * i - 1]; FileName[2 * i]
  
  # Process first file in the pair
  cog_deseq1 <- read.csv(FileName[2 * i - 1], header = TRUE)
  cog_name1 = cog_deseq1$X
  if ("Unclassified" %in% cog_name1) {
    cog_name1 = cog_name1[-which(cog_name1 == "Unclassified")]
  }

  # Perform COG enrichment analysis for first file
  cog_enrichment1 = enrichCOG(
    cog_name1,
    dtype = "category",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2
  )

  # Process second file in the pair
  cog_deseq2 <- read.csv(FileName[2 * i], header = TRUE)
  cog_name2 = cog_deseq2$X
  if ("Unclassified" %in% cog_name2) {
    cog_name2 = cog_name2[-which(cog_name2 == "Unclassified")]
  }

  # Perform COG enrichment analysis for second file
  cog_enrichment2 = enrichCOG(
    cog_name2,
    dtype = "category",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2
  )

  # Adjust row names for both enrichment results
  row.names(cog_enrichment2@result) = seq(1, nrow(cog_enrichment2@result), by = 1)
  row.names(cog_enrichment1@result) = seq(nrow(cog_enrichment2@result) + 1, 
                                          nrow(cog_enrichment1@result) + nrow(cog_enrichment2@result), by = 1)
  
  # Modify GeneRatio for first enrichment result
  minus = rep('-', times = nrow(cog_enrichment1@result))
  cog_enrichment1@result[["GeneRatio"]] = paste0(minus, cog_enrichment1@result[["GeneRatio"]])
  
  # Combine both enrichment results
  cog_enrichment2@result = rbind(cog_enrichment2@result, cog_enrichment1@result)

  # Calculate GeneRatio
  a = vector()
  b = vector()
  for (j in 1:nrow(cog_enrichment2@result)) {
    a[j] = parse(text = cog_enrichment2@result[["GeneRatio"]][j])
    b[j] = eval(a[j])
  }
  cog_enrichment2@result[["GeneRatio"]] = b

  # Filter and arrange results based on p-value and GeneRatio
  cog_enrichment2@result = subset(cog_enrichment2@result, p.adjust < 0.05)
  cog_enrichment2@result = arrange(cog_enrichment2@result, GeneRatio)
  cog_enrichment2@result$Description = factor(cog_enrichment2@result$Description,
                                              levels = unique(cog_enrichment2@result$Description))

  # Plot results
  p = ggplot(cog_enrichment2@result, aes(y = GeneRatio, x = Description)) +
    geom_bar(stat = "identity", width = 0.4, lwd = 1, 
             fill = ifelse(cog_enrichment2@result$GeneRatio > 0, '#ffb5b3', '#a5ceff')) +
    geom_hline(yintercept = 0, lty = 1, color = 'grey', lwd = 1) +
    ggtitle(Outputname[[2 * i]]) + 
    theme_prism(border = TRUE) +
    theme(axis.text.y = element_text(face = "plain", size = 10), 
          legend.title = element_text(size = 12)) +
    coord_flip() +
    xlab('COG functional category') +
    guides(size = guide_legend(order = 2)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
    scale_y_continuous(limits = c(-0.16, 0.16)) +
    labs(colour = 'P.adjusted', size = 'Gene count')
  
  # Save plot to PowerPoint
  topptx(p, paste0(Outputname[[2 * i]], '_bar.pptx'), height = 3.7, width = 5.5)
}

# Plot the number of significantly changed COGs
NoCog = read.table('clipboard', header = TRUE, sep = "\t")
NoCog$Group = as.factor(NoCog$Group)
levels(NoCog$Group) = c('Biofilm vs planktonic', 'Planktonic vs inoculum', 'Biofilm vs inoculum', 'Plastic-fed vs control worm gut')

p = ggplot(NoCog, aes(y = No_of_differentially_abundant_COG, x = Group)) +
  geom_bar(stat = "identity", width = 0.4, lwd = 1, 
           fill = ifelse(NoCog$No_of_differentially_abundant_COG > 0, '#ffb5b3', '#a5ceff')) +
  geom_hline(yintercept = 0, lty = 1, color = 'grey', lwd = 1) +
  theme_prism(border = TRUE) +
  theme(axis.text.y = element_text(face = "plain", size = 14), 
        legend.title = element_text(size = 12)) +
  coord_flip() +
  xlab('Group') +
  guides(size = guide_legend(order = 2)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_y_continuous(limits = c(-9500, 9500))

# Save plot to PowerPoint
topptx(p, 'Cog_bar.pptx', height = 3.7, width = 5.5)

# Old manuscript using dot plots
p = ggplot(cog_enrichment2@result, aes(y = GeneRatio, x = Description)) +
  geom_segment(aes(x = Description, xend = Description, y = 0, yend = GeneRatio), color = "grey", lwd = 0.8) +
  geom_point(aes(size = Count, color = p.adjust)) +
  ggtitle(Outputname[[2 * i]]) +
  theme_prism(border = TRUE) +
  theme(axis.text.y = element_text(face = "plain", size = 12), 
        legend.title = element_text(size = 12)) +
  geom_hline(yintercept = 0, lty = 2, color = 'grey', lwd = 0.8) +
  coord_flip() +
  xlab('COG functional category') +
  guides(size = guide_legend(order = 2)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_y_continuous(limits = c(-0.16, 0.08)) +
  labs(colour = 'P.adjusted', size = 'Gene count') +
  scale_colour_gradient(low = "red", high = "blue", labels = scales::scientific_format(), name = 'P.adjusted')

