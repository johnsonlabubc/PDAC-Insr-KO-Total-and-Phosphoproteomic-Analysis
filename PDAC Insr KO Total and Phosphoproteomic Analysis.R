library(clusterProfiler)
library(ComplexHeatmap)
library(dplyr)
library(ggh4x)
library(ggplot2)
library(ggrepel)
library(limma)
library(org.Mm.eg.db)
library(scales)
library(stringi)
library(stringr)
library(this.path)
library(tidyr)


#set working directory to location of this R script
setwd(dirname(this.path::this.path()))


### Figure 5 ###

## Fig 5A ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import files
df_total_proteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_total_protein.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_total_proteomics) <- df_total_proteomics[1,]
df_total_proteomics <- df_total_proteomics[-1,]

#selecting and reordering columns
df_total_proteomics <- df_total_proteomics[, c(3:4, 21:32, 59:78)]
df_total_proteomics <- cbind(df_total_proteomics[, c(1:4, 6:7, 9:10, 12:13)],
                             df_total_proteomics[, c(15:34)][, grepl("WT_G_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("KO_G_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("WT_G12D_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("Het_G12D_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("KO_G12D_", colnames(df_total_proteomics)[c(15:34)])])

#covert numbers to numeric
df_total_proteomics[, c(3:30)] <- sapply(df_total_proteomics[, c(3:30)], as.numeric)

#calculate fold change
df_total_proteomics$`(G12D_WT) / (G_WT)` <- rowMedians(as.matrix(df_total_proteomics[, grepl("WT_G12D_", colnames(df_total_proteomics))]))/rowMedians(as.matrix(df_total_proteomics[, grepl("WT_G_", colnames(df_total_proteomics))]))

#select needed columns
df_total_proteomics <- df_total_proteomics[, c("Gene Symbol", "(G12D_WT) / (G_WT)", "Adj. P-Value: (G12D_WT) / (G_WT)")] %>% na.omit()

#get log2 of fold change
df_total_proteomics$log2FoldChange <- log2(df_total_proteomics[, 2])

#rename columns
colnames(df_total_proteomics) <- c("Gene", "FoldChange", "padj", "log2FoldChange")

#label whether protein is differentially expressed
df_total_proteomics$DiffExpression <- "ns"
df_total_proteomics$DiffExpression[df_total_proteomics$padj < 0.05 & df_total_proteomics$log2FoldChange > 0] <- "Up in G12D"
df_total_proteomics$DiffExpression[df_total_proteomics$padj < 0.05 & df_total_proteomics$log2FoldChange < 0] <- "Down in G12D"

#denote significance and fold change with colour
df_total_proteomics$Colour[df_total_proteomics$padj < 0.05 & df_total_proteomics$log2FoldChange > 0] <- scales::hue_pal()(2)[1]
df_total_proteomics$Colour[df_total_proteomics$padj < 0.05 & df_total_proteomics$log2FoldChange < 0] <- scales::hue_pal()(2)[2]

#get -log10 of padj
df_total_proteomics$transf.log.10.padj <- -log10(df_total_proteomics$padj)

#mark proteins for labeling
df_total_proteomics$Label[df_total_proteomics$transf.log.10.padj > 10] <- "Yes"

#set expression variable factor level order
df_total_proteomics$DiffExpression <- factor(df_total_proteomics$DiffExpression, levels = sort(unique(df_total_proteomics$DiffExpression))[c(3, 1, 2)])

#create plot
ggplot2::ggplot(df_total_proteomics, ggplot2::aes(x = log2FoldChange, y = transf.log.10.padj)) +
  
  ggplot2::ggtitle(expression(italic("Ptf1a")^CreER*";"*italic("Insr")^"w/w"*";nTnG vs PK-"*italic("Insr")^"w/w")) + 
  
  ggplot2::geom_point(ggplot2::aes(colour = DiffExpression)) +
  
  ggrepel::geom_text_repel(ggplot2::aes(label = Gene),
                           seed = 123,
                           max.overlaps = 200,
                           min.segment.length = 0.2,
                           segment.alpha = 0.3,
                           segment.colour = "black",
                           colour = df_total_proteomics$Colour[!is.na(df_total_proteomics$Label)],
                           bg.color = "white",
                           bg.r = 0.05,
                           nudge_y = -0.1,
                           force = 20,
                           force_pull = 0,
                           data = df_total_proteomics %>% na.omit()) +
  
  ggplot2::coord_cartesian(xlim = c(-8, 8), 
                           ylim = c(-1, 16.5), 
                           expand = FALSE, clip = "off") +
  
  ggplot2::scale_color_manual(values = c(scales::hue_pal()(2)[1], scales::hue_pal()(2)[2], "grey"),
                              labels = c(expression("Up in PK-"*italic("Insr")^"w/w"),
                                         expression("Down in PK-"*italic("Insr")^"w/w"),
                                         "ns")) +
  
  ggplot2::labs(y = bquote(-log[10]*"(Adj. "*italic(P)*"-value)"), x = bquote(log[2]*"(Fold Change)"), colour = "Regulation") +
  
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(shape = 16))) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, vjust = -1.5,size = 15, margin = ggplot2::margin(b = 0.5, unit = "cm", )),
                 plot.margin = ggplot2::margin(t = 1, r = 0.5, b = 0.5, l = 0.5, "cm"),
                 axis.text.x = ggplot2::element_text(size = 12, colour = "black"),
                 axis.text.y = ggplot2::element_text(size = 11, colour = "black"),
                 axis.title.x = ggplot2::element_text(size = 12),
                 panel.grid.minor = ggplot2::element_blank(),
                 legend.title.align = 0,
                 legend.text.align = 0,
                 legend.direction = "vertical",
                 legend.box.just = "center",
                 legend.key.width = grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = 1/1)


## Fig 5B ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import files
df_total_proteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_total_protein.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_total_proteomics) <- df_total_proteomics[1,]
df_total_proteomics <- df_total_proteomics[-1,]

#selecting and reordering columns
df_total_proteomics <- df_total_proteomics[, c(3:4, 21:32, 59:78)]
df_total_proteomics <- cbind(df_total_proteomics[, c(1:4, 6:7, 9:10, 12:13)],
                             df_total_proteomics[, c(15:34)][, grepl("WT_G_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("KO_G_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("WT_G12D_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("Het_G12D_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("KO_G12D_", colnames(df_total_proteomics)[c(15:34)])])

#covert numbers to numeric
df_total_proteomics[, c(3:30)] <- sapply(df_total_proteomics[, c(3:30)], as.numeric)

#calculate fold change
df_total_proteomics$`(G_KO) / (G_WT)` <- rowMedians(as.matrix(df_total_proteomics[, grepl("KO_G_", colnames(df_total_proteomics))]))/rowMedians(as.matrix(df_total_proteomics[, grepl("WT_G_", colnames(df_total_proteomics))]))

#select needed columns
df_total_proteomics <- df_total_proteomics[, c("Gene Symbol", "(G_KO) / (G_WT)", "Adj. P-Value: (G_KO) / (G_WT)")] %>% na.omit()

#get log2 of fold change
df_total_proteomics$log2FoldChange <- log2(df_total_proteomics[, 2])

#rename columns
colnames(df_total_proteomics) <- c("Gene", "FoldChange", "padj", "log2FoldChange")

#label whether protein is differentially expressed
df_total_proteomics$DiffExpression <- "ns"
df_total_proteomics$DiffExpression[df_total_proteomics$padj < 0.05 & df_total_proteomics$log2FoldChange > 0] <- "Up in KO"
df_total_proteomics$DiffExpression[df_total_proteomics$padj < 0.05 & df_total_proteomics$log2FoldChange < 0] <- "Down in KO"

#denote significance and fold change with colour
df_total_proteomics$Colour[df_total_proteomics$padj < 0.05 & df_total_proteomics$log2FoldChange > 0] <- scales::hue_pal()(2)[1]
df_total_proteomics$Colour[df_total_proteomics$padj < 0.05 & df_total_proteomics$log2FoldChange < 0] <- scales::hue_pal()(2)[2]

#get -log10 of padj
df_total_proteomics$transf.log.10.padj <- -log10(df_total_proteomics$padj)

#mark proteins for labeling
df_total_proteomics$Label[df_total_proteomics$transf.log.10.padj > 10] <- "Yes"

#set expression variable factor level order
df_total_proteomics$DiffExpression <- factor(df_total_proteomics$DiffExpression, levels = sort(unique(df_total_proteomics$DiffExpression))[c(3, 1, 2)])

#create plot
ggplot2::ggplot(df_total_proteomics, ggplot2::aes(x = log2FoldChange, y = transf.log.10.padj)) +
  
  ggplot2::ggtitle(expression(italic("Ptf1a")^CreER*";"*italic("Insr")^"w/w"*";nTnG vs "*italic("Ptf1a")^CreER*";"*italic("Insr")^"f/f"*";nTnG")) + 
  
  ggplot2::geom_point(ggplot2::aes(colour = DiffExpression)) +
  
  ggrepel::geom_text_repel(ggplot2::aes(label = Gene),
                           seed = 123,
                           max.overlaps = 200,
                           min.segment.length = 0.2,
                           segment.alpha = 0.3,
                           segment.colour = "black",
                           colour = df_total_proteomics$Colour[!is.na(df_total_proteomics$Label)],
                           bg.color = "white",
                           bg.r = 0.05,
                           nudge_y = -0.1,
                           force = 20,
                           force_pull = 0,
                           data = df_total_proteomics %>% na.omit()) +
  
  ggplot2::coord_cartesian(xlim = c(-8, 8), 
                           ylim = c(-1, 16.5), 
                           expand = FALSE, clip = "off") +
  
  ggplot2::scale_color_manual(values = c(scales::hue_pal()(2)[1], scales::hue_pal()(2)[2], "grey"),
                              labels = c(expression("Up in "*italic("Ptf1a")^CreER*";"*italic("Insr")^"f/f"*";nTnG"),
                                         expression("Down in "*italic("Ptf1a")^CreER*";"*italic("Insr")^"f/f"*";nTnG"),
                                         "ns")) +
  
  ggplot2::labs(y = bquote(-log[10]*"(Adj. "*italic(P)*"-value)"), x = bquote(log[2]*"(Fold Change)"), colour = "Regulation") +
  
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(shape = 16))) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, vjust = -1.5,size = 15, margin = ggplot2::margin(b = 0.5, unit = "cm", )),
                 plot.margin = ggplot2::margin(t = 1, r = 0.5, b = 0.5, l = 0.5, "cm"),
                 axis.text.x = ggplot2::element_text(size = 12, colour = "black"),
                 axis.text.y = ggplot2::element_text(size = 11, colour = "black"),
                 axis.title.x = ggplot2::element_text(size = 12),
                 panel.grid.minor = ggplot2::element_blank(),
                 legend.title.align = 0,
                 legend.text.align = 0,
                 legend.direction = "vertical",
                 legend.box.just = "center",
                 legend.key.width = grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = 1/1)


## Fig 5C ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import files
df_total_proteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_total_protein.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_total_proteomics) <- df_total_proteomics[1,]
df_total_proteomics <- df_total_proteomics[-1,]

#selecting and reordering columns
df_total_proteomics <- df_total_proteomics[, c(3:4, 21:32, 59:78)]
df_total_proteomics <- cbind(df_total_proteomics[, c(1:4, 6:7, 9:10, 12:13)],
                             df_total_proteomics[, c(15:34)][, grepl("WT_G_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("KO_G_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("WT_G12D_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("Het_G12D_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("KO_G12D_", colnames(df_total_proteomics)[c(15:34)])])

#covert numbers to numeric
df_total_proteomics[, c(3:30)] <- sapply(df_total_proteomics[, c(3:30)], as.numeric)

#calculate fold change
df_total_proteomics$`(G12D_KO) / (G12D_WT)` <- rowMedians(as.matrix(df_total_proteomics[, grepl("KO_G12D_", colnames(df_total_proteomics))]))/rowMedians(as.matrix(df_total_proteomics[, grepl("WT_G12D_", colnames(df_total_proteomics))]))

#select needed columns
df_total_proteomics <- df_total_proteomics[, c("Gene Symbol", "(G12D_KO) / (G12D_WT)", "Adj. P-Value: (G12D_KO) / (G12D_WT)")] %>% na.omit()

#get log2 of fold change
df_total_proteomics$log2FoldChange <- log2(df_total_proteomics[, 2])

#rename columns
colnames(df_total_proteomics) <- c("Gene", "FoldChange", "padj", "log2FoldChange")

#label whether protein is differentially expressed
df_total_proteomics$DiffExpression <- "ns"
df_total_proteomics$DiffExpression[df_total_proteomics$padj < 0.05 & df_total_proteomics$log2FoldChange > 0] <- "Up in KO"
df_total_proteomics$DiffExpression[df_total_proteomics$padj < 0.05 & df_total_proteomics$log2FoldChange < 0] <- "Down in KO"

#denote significance and fold change with colour
df_total_proteomics$Colour[df_total_proteomics$padj < 0.05 & df_total_proteomics$log2FoldChange > 0] <- scales::hue_pal()(2)[1]
df_total_proteomics$Colour[df_total_proteomics$padj < 0.05 & df_total_proteomics$log2FoldChange < 0] <- scales::hue_pal()(2)[2]

#get -log10 of padj
df_total_proteomics$transf.log.10.padj <- -log10(df_total_proteomics$padj)

#mark proteins for labeling
df_total_proteomics$Label[df_total_proteomics$transf.log.10.padj > 10] <- "Yes"

#set expression variable factor level order
df_total_proteomics$DiffExpression <- factor(df_total_proteomics$DiffExpression, levels = sort(unique(df_total_proteomics$DiffExpression))[c(3, 1, 2)])

#create plot
ggplot2::ggplot(df_total_proteomics, ggplot2::aes(x = log2FoldChange, y = transf.log.10.padj)) +
  
  ggplot2::ggtitle(expression("PK-"*italic("Insr")^"w/w"*" vs "*"PK-"*italic("Insr")^"f/f")) + 
  
  ggplot2::geom_point(ggplot2::aes(colour = DiffExpression)) +
  
  ggrepel::geom_text_repel(ggplot2::aes(label = Gene),
                           seed = 123,
                           max.overlaps = 200,
                           min.segment.length = 0.2,
                           segment.alpha = 0.3,
                           segment.colour = "black",
                           colour = df_total_proteomics$Colour[!is.na(df_total_proteomics$Label)],
                           bg.color = "white",
                           bg.r = 0.05,
                           nudge_y = -0.1,
                           force = 20,
                           force_pull = 0,
                           data = df_total_proteomics %>% na.omit()) +
  
  ggplot2::coord_cartesian(xlim = c(-8, 8), 
                           ylim = c(-1, 16.5), 
                           expand = FALSE, clip = "off") +
  
  ggplot2::scale_color_manual(values = c(scales::hue_pal()(2)[1], scales::hue_pal()(2)[2], "grey"),
                              labels = c(expression("Up in PK-"*italic("Insr")^"f/f"),
                                         expression("Down in PK-"*italic("Insr")^"f/f"),
                                         "ns")) +
  
  ggplot2::labs(y = bquote(-log[10]*"(Adj. "*italic(P)*"-value)"), x = bquote(log[2]*"(Fold Change)"), colour = "Regulation") +
  
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(shape = 16))) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, vjust = -1.5, size = 15, margin = ggplot2::margin(b = 0.5, unit = "cm")),
                 plot.margin = ggplot2::margin(t = 1, r = 0.5, b = 0.5, l = 0.5, "cm"),
                 axis.text.x = ggplot2::element_text(size = 12, colour = "black"),
                 axis.text.y = ggplot2::element_text(size = 11, colour = "black"),
                 axis.title.x = ggplot2::element_text(size = 12),
                 panel.grid.minor = ggplot2::element_blank(),
                 legend.title.align = 0,
                 legend.text.align = 0,
                 legend.direction = "vertical",
                 legend.box.just = "center",
                 legend.key.width = grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = 1/1)


## Fig 5D ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import files
df_total_proteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_total_protein.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_total_proteomics) <- df_total_proteomics[1,]
df_total_proteomics <- df_total_proteomics[-1,]

#selecting and reordering columns
df_total_proteomics <- df_total_proteomics[, c(3:4, 21:32, 59:78)]

df_total_proteomics <- cbind(df_total_proteomics[, c(1:4, 6:7, 9:10, 12:13)],
                             df_total_proteomics[, c(15:34)][, grepl("WT_G_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("KO_G_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("WT_G12D_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("Het_G12D_", colnames(df_total_proteomics)[c(15:34)])],
                             df_total_proteomics[, c(15:34)][, grepl("KO_G12D_", colnames(df_total_proteomics)[c(15:34)])])

#covert numbers to numeric
df_total_proteomics[, c(3:30)] <- sapply(df_total_proteomics[, c(3:30)], as.numeric)

#filter for significant proteins
df_total_proteomics <- subset(df_total_proteomics, `Adj. P-Value: (G12D_KO) / (G12D_WT)` < 0.05,) %>% na.omit()

#get mean for each sample of each protein
df_protein_levels <- data.frame(`WT G` = rowMeans(df_total_proteomics[, c(11:30)][, grepl("WT_G_", colnames(df_total_proteomics)[c(11:30)])]),
                                `KO G` = rowMeans(df_total_proteomics[, c(11:30)][, grepl("KO_G_", colnames(df_total_proteomics)[c(11:30)])]),
                                `WT G12D` = rowMeans(df_total_proteomics[, c(11:30)][, grepl("WT_G12D_", colnames(df_total_proteomics)[c(11:30)])]),
                                `Het G12D` = rowMeans(df_total_proteomics[, c(11:30)][, grepl("Het_G12D_", colnames(df_total_proteomics)[c(11:30)])]),
                                `KO G12D` = rowMeans(df_total_proteomics[, c(11:30)][, grepl("KO_G12D_", colnames(df_total_proteomics)[c(11:30)])]),
                                check.names = FALSE)

#convert to Z-score matrix
df_protein_levels$mean_all <- rowMeans(df_protein_levels[, c(1:5)])

df_protein_levels$sd_all <- apply(df_protein_levels[, c(1:5)], 1, sd)

Z_score <- function(x) with(df_protein_levels, (x - mean_all)/sd_all)

df_protein_Z_score <- cbind(apply(df_protein_levels[, c(1:5)], 2, Z_score))

rownames(df_protein_Z_score) <- df_total_proteomics$`Gene Symbol`

matrix_Z_score_total <- as.matrix(df_protein_Z_score)

#create heatmap
heatmap_list_total <- ComplexHeatmap::Heatmap(matrix_Z_score_total,
                                              name = "matrix_Z_score_total",
                                              show_row_names = FALSE,
                                              show_column_names = FALSE,
                                              row_names_side = "right",
                                              row_names_gp = grid::gpar(fontface = "italic", fontsize = 12),
                                              column_names_side = "top",
                                              row_dend_side = "left",
                                              column_dend_side = "bottom",
                                              row_km = 8,
                                              col = circlize::colorRamp2(c(-ceiling(max(abs(matrix_Z_score_total), na.rm = T)), 0, ceiling(max(abs(matrix_Z_score_total), na.rm = T))), c(scales::hue_pal()(2)[2], "white", c(scales::hue_pal()(2)[1]))),
                                              top_annotation = ComplexHeatmap::columnAnnotation(empty = ComplexHeatmap::anno_empty(border = FALSE, height = grid::unit(12, "mm"))),
                                              column_order = 1:ncol(matrix_Z_score_total),
                                              height = grid::unit(360, "mm"),
                                              width = ncol(matrix_Z_score_total) * grid::unit(20, "mm"),
                                              border_gp = grid::gpar(col = "black"), 
                                              heatmap_legend_param = list(
                                                title = NULL,
                                                title_position = "topcenter",
                                                legend_height = grid::unit(6, "cm"),
                                                at = c(-4, -2, 0, 2, 4),
                                                direction = c("horizontal")))

#draw heatmap
heatmap_list_total <- ComplexHeatmap::draw(heatmap_list_total, heatmap_legend_side = "bottom")

#select viewport and define dimensions
grid::seekViewport("annotation_empty_1")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create top labels
grid::grid.text("Insr",
                x = (loc2$x - loc1$x) * -0.0075,
                y = 0.25,
                gp = grid::gpar(fontsize = 14,
                                fontface = "italic"),
                hjust = 1)

grid::grid.text(c(expression("w/w"),
                  expression("f/f"),
                  expression("w/w"),
                  expression("w/f"),
                  expression("f/f")),
                x = c((1/2)/ncol(matrix_Z_score_total),
                      (1 + 1/2)/ncol(matrix_Z_score_total),
                      (2 + 1/2)/ncol(matrix_Z_score_total),
                      (3 + 1/2)/ncol(matrix_Z_score_total),
                      (4 + 1/2)/ncol(matrix_Z_score_total)),
                y = 0.25,
                just = c("center", "center"),
                gp = grid::gpar(fontsize = 14,
                                col = c("#0e6e80",
                                        "#ce1e68",
                                        "#3b54a3",
                                        "#f47d23",
                                        "#bbbd32")))

grid::grid.rect(x = 0.7,
                y = 0.6,
                height = grid::unit(0.01, "inches"),
                width = 3/5*(loc2$x - loc1$x),
                just = c("center", "center"),
                gp = grid::gpar(fill = "black", 
                                col = "black", 
                                lwd = 1))

grid::grid.text(expression(italic("Kras")^"G12D"),
                x = (loc2$x - loc1$x) * 0.7,
                y = 1,
                gp = grid::gpar(fontsize = 12),
                hjust = 0.5)

#select viewport and define dimensions
grid::seekViewport("heatmap_legend")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create legend title
grid::grid.text("Protein Z-Score",
                x = 0.5,
                y = -0.5,
                gp = grid::gpar(fontsize = 11))


### Figure 6 ###

## Fig 6A ##

#clear environment
rm(list = ls())

#import file
df_total_proteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_total_protein.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_total_proteomics) <- df_total_proteomics[1,]
df_total_proteomics <- df_total_proteomics[-1,]

#select columns
df_total_proteomics <- df_total_proteomics[, c("Gene Symbol", "Accession", "(G12D_KO) / (G12D_WT)", "Adj. P-Value: (G12D_KO) / (G12D_WT)")]

#convert numbers to numeric
df_total_proteomics[, c("(G12D_KO) / (G12D_WT)", "Adj. P-Value: (G12D_KO) / (G12D_WT)")] <- as.data.frame(sapply(df_total_proteomics[, c("(G12D_KO) / (G12D_WT)", "Adj. P-Value: (G12D_KO) / (G12D_WT)")], as.numeric))

#get log2 of fold change
df_total_proteomics$`log2 (G12D_KO) / (G12D_WT)` <- log2(df_total_proteomics$`(G12D_KO) / (G12D_WT)`)

#get -log10 of adj. p-value
df_total_proteomics$`-log10 Adj. P-Value: (G12D_KO) / (G12D_WT)` <- -log10(df_total_proteomics$`Adj. P-Value: (G12D_KO) / (G12D_WT)`)

#subset for significance and select proteins
df_total_proteomics <- subset(df_total_proteomics, `Adj. P-Value: (G12D_KO) / (G12D_WT)` < 0.05 | df_total_proteomics$`Gene Symbol` %in% c("Akt1", "Akt2"))

#add network anchors
anchors <- data.frame(`Gene Symbol` = c("Akt3", "Kras", "Igf1r", "Insr", "mTor", "Raf1", "Mapk1"),
                      a = c("Q9WUA6", "P32883", "Q60751", "P15208", "Q9JLN9", "Q99N57", "P63085"), 
                      b = "",
                      c = "",
                      d = "",
                      e = "")

colnames(anchors) <- colnames(df_total_proteomics)

df_total_proteomics <- rbind(df_total_proteomics, 
                             anchors)

#export node table to total_protein_G12D_WT_vs_KO sig logFC and padj.tsv
write.table(df_total_proteomics)

#import STRING protein mapping and edge table output files
df_string_network <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_protein_G12D_WT_vs_KO sig string_interactions_uniprot.tsv", check.names = FALSE, sep = "\t")
df_string_mapping <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_protein_G12D_WT_vs_KO sig string_mapping_uniprot.tsv", check.names = FALSE, sep = "\t")

#select protein mapping columns
df_string_mapping <- df_string_mapping[, c("queryItem", "preferredName")]

#rename column names
colnames(df_string_mapping) <- c("queryItem2", "node2")

#merge mapping to network
df_string_network <- merge(df_string_network, df_string_mapping, by = c("node2"), all.x = TRUE)

#rename column names
colnames(df_string_mapping) <- c("queryItem", "#node1")

#merge mapping to network
df_string_network <- merge(df_string_network, df_string_mapping, by = c("#node1"), all.x = TRUE)

#replace mapped values with input values
df_string_network$`#node1` <- df_string_network$queryItem
df_string_network$node2 <- df_string_network$queryItem2

#get protein nodes
df_string_input <- df_total_proteomics[, c("Gene Symbol","Accession")] 

#rename column names
colnames(df_string_input) <- c("queryItem2_gene", "node2")

#merge input to network
df_string_network <- merge(df_string_network, df_string_input[, c("queryItem2_gene", "node2")], by = "node2")

#rename column names
colnames(df_string_input) <- c("queryItem_gene", "#node1")

#merge input to network
df_string_network <- merge(df_string_network, df_string_input[, c("queryItem_gene", "#node1")], by = "#node1")

#replace mapped values with input values
df_string_network$`#node1` <- df_string_network$queryItem_gene
df_string_network$node2 <- df_string_network$queryItem2_gene

#appending missing rows
df_string_network[nrow(df_string_network) + 1,] <- c("Psmd11",	"Psme2", rep("", 10), 0.937, rep("", 4))
df_string_network[nrow(df_string_network) + 1,] <- c("Psmd7",	"Psme2", rep("", 10), 0.938, rep("", 4))

#add empty columns
df_string_input[c(2:ncol(df_string_network))] <- ""

#add matching column names
colnames(df_string_input) <- colnames(df_string_network)

#remove proteins that already plotted in network
df_string_input <- subset(df_string_input, !`#node1` %in% c(df_string_network$`#node1`, df_string_network$`node2`))

#combine connected and unconnected protein
df_string_network <- rbind(df_string_network, df_string_input)

#export edge table to SCPC#46_03_Jim_Anni_pancreatic_cancer_total_protein_G12D_WT_vs_KO sig string_interactions remapped.tsv
write.table(df_string_network)


## Fig 6B ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import files
df_total_proteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_total_protein.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_total_proteomics) <- df_total_proteomics[1,]
df_total_proteomics <- df_total_proteomics[-1,]

#selecting and reordering columns
df_total_proteomics <- df_total_proteomics[, c(3:4, 59:78)]
df_total_proteomics <- cbind(df_total_proteomics[, c(1:2)],
                             df_total_proteomics[, c(3:22)][, grepl("WT_G_", colnames(df_total_proteomics)[c(3:22)])],
                             df_total_proteomics[, c(3:22)][, grepl("KO_G_", colnames(df_total_proteomics)[c(3:22)])],
                             df_total_proteomics[, c(3:22)][, grepl("WT_G12D_", colnames(df_total_proteomics)[c(3:22)])],
                             df_total_proteomics[, c(3:22)][, grepl("Het_G12D_", colnames(df_total_proteomics)[c(3:22)])],
                             df_total_proteomics[, c(3:22)][, grepl("KO_G12D_", colnames(df_total_proteomics)[c(3:22)])])

#covert numbers to numeric
df_total_proteomics[, c(3:22)] <- sapply(df_total_proteomics[, c(3:22)], as.numeric)

#select specific digestive enzymes
df_total_proteomics <- df_total_proteomics[df_total_proteomics$`Gene Symbol` %in% c("Rnase1", "Cuzd1", "Cpa2", "Klk1", "Pnliprp2", "Cela3b", "Prss2", "Cela1", "Sycn", "Amy1", "Ctrb1", "Cel", "Amy2", "Cpa1", "Pnlip", "Pnliprp1", "Cela2a", "Pla2g1b", "Clps", "Ctrc", "Rnase4", "Spink1"),]

#select abundance columns
df_protein_levels <- df_total_proteomics[, c(3:22)]

#convert to Z-score matrix
df_protein_levels$mean_all <- rowMeans(df_protein_levels[, c(1:20)])

df_protein_levels$sd_all <- apply(df_protein_levels[, c(1:20)], 1, sd) 

Z_score <- function(x) with(df_protein_levels, (x - mean_all)/sd_all)

df_protein_Z_score <- cbind(apply(df_protein_levels[, c(1:(ncol(df_protein_levels) - 2))], 2, Z_score))

rownames(df_protein_Z_score) <- df_total_proteomics$`Gene Symbol`

matrix_Z_score_total <- as.matrix(df_protein_Z_score)

#create heatmap
heatmap_list_total <- ComplexHeatmap::Heatmap(matrix_Z_score_total,
                                              name = "matrix_Z_score_total",
                                              show_column_names = FALSE,
                                              row_names_side = "right",
                                              row_names_gp = grid::gpar(fontface = "plain", fontsize = 12),
                                              column_names_side = "top",
                                              row_dend_side = "left",
                                              column_dend_side = "bottom",
                                              layer_fun = function(j, i, x, y, width, height, fill) {
                                                mat = ComplexHeatmap::restore_matrix(j, i, x, y)
                                                ind = unique(c(mat[, c(3, 5, 8, 14)]))
                                                grid.rect(x = x[ind] + grid::unit(0.5/ncol(matrix_Z_score_total), "npc"), 
                                                          y = y[ind], 
                                                          width = grid::unit(0.03, "inches"), 
                                                          height = grid::unit(1/nrow(matrix_Z_score_total), "npc"),
                                                          gp = grid::gpar(col = "white")
                                                )
                                              },
                                              col = circlize::colorRamp2(c(-ceiling(max(abs(matrix_Z_score_total), na.rm = T)), 0, ceiling(max(abs(matrix_Z_score_total), na.rm = T))), c(scales::hue_pal()(2)[2], "white", c(scales::hue_pal()(2)[1]))),
                                              top_annotation = ComplexHeatmap::columnAnnotation(empty = ComplexHeatmap::anno_empty(border = FALSE, height = grid::unit(12, "mm"))),
                                              column_order = 1:ncol(matrix_Z_score_total),
                                              height = nrow(matrix_Z_score_total) * grid::unit(4.2, "mm"),
                                              width = ncol(matrix_Z_score_total) * grid::unit(6, "mm"),
                                              border_gp = grid::gpar(col = "black"), 
                                              heatmap_legend_param = list(
                                                title = NULL,
                                                title_position = "topcenter",
                                                legend_height = grid::unit(6, "cm"),
                                                direction = c("horizontal")))

#draw heatmap
ComplexHeatmap::draw(heatmap_list_total, heatmap_legend_side = "bottom")

#select viewport and define dimensions
grid::seekViewport("annotation_empty_1")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create top labels
grid::grid.text("Insr", 
                x = (loc2$x - loc1$x) * -0.0075,
                y = 0.25,
                gp = grid::gpar(fontsize = 14,
                                fontface = "italic"), 
                hjust = 1)

grid::grid.text(c(expression("w/w"), 
                  expression("f/f"),
                  expression("w/w"),
                  expression("w/f"),
                  expression("f/f")),
                x = c((3/2)/ncol(matrix_Z_score_total),
                      (3 + 2/2)/ncol(matrix_Z_score_total),
                      (3 + 2 + 3/2)/ncol(matrix_Z_score_total),
                      (3 + 2 + 3 + 6/2)/ncol(matrix_Z_score_total),
                      (3 + 2 + 3 + 6 + 6/2)/ncol(matrix_Z_score_total)),
                y = 0.25,
                just = c("center", "center"),
                gp = grid::gpar(fontsize = 14,
                                col = c("#0e6e80",
                                        "#ce1e68",
                                        "#3b54a3",
                                        "#f47d23",
                                        "#bbbd32")))

grid::grid.rect(x = 12.5/20,
                y = 0.6,
                height = grid::unit(0.01, "inches"),
                width = 15/20*(loc2$x - loc1$x),
                just = c("center", "center"),
                gp = grid::gpar(fill = "black", 
                                col = "black", 
                                lwd = 1))

grid::grid.text(expression(italic("Kras")^"G12D"), 
                x = (loc2$x - loc1$x) * 12.5/20,
                y = 1,
                gp = grid::gpar(fontsize = 12), 
                hjust = 0.5)

#select viewport and define dimensions
grid::seekViewport("heatmap_legend")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create legend title
grid::grid.text("Protein Z-Score",
                x = 0.5,
                y = -0.5,
                gp = grid::gpar(fontsize = 11))


## Fig 6C ##

#clear environment
rm(list = ls())

#import file
df_total_proteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_total_protein.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_total_proteomics) <- df_total_proteomics[1,]
df_total_proteomics <- df_total_proteomics[-1,]

#select columns
df_total_proteomics <- df_total_proteomics[, c(3:4, 21:32, 59:78)]

#order columns
df_total_proteomics <- df_total_proteomics[, c(1:14, 17:19, 15:16, 32:34, 20:25, 26:31)]

#convert numbers to numeric
df_total_proteomics[, c(3:4, 6:7, 9:10, 12:13, 15:34)] <- sapply(df_total_proteomics[, c(3:4, 6:7, 9:10, 12:13, 15:34)], as.numeric)

#create digestive enzyme data frame
df_digestive_enzymes <- df_total_proteomics[df_total_proteomics$`Gene Symbol` %in% c("Rnase1", "Cuzd1", "Cpa2", "Klk1", "Pnliprp2", "Cela3b", "Prss2", "Cela1", "Sycn", "Amy1", "Ctrb1", "Cel", "Amy2", "Cpa1", "Pnlip", "Pnliprp1", "Cela2a", "Pla2g1b", "Clps", "Ctrc", "Rnase4", "Spink1"),]

#sum abundance of all rows for each sample
df_total_proteomics_sum <- aggregate(df_total_proteomics[, c(15:34)], by = list(rep(1, nrow(df_total_proteomics[, c(15:34)]))), FUN = sum)[, -1]
df_digestive_enzymes_sum <- aggregate(df_digestive_enzymes[, c(15:34)], by = list(rep(1, nrow(df_digestive_enzymes[, c(15:34)]))), FUN = sum)[, -1]

#calculate percent of digestive enzyme abundance in total abundance
df_percentages <- df_digestive_enzymes_sum/df_total_proteomics_sum * 100

#gather values into longer format
df_plot <- gather(df_percentages, key = "Sample", value = "Percent Digestive Enzymes out of Total")

#format condition labels
df_plot <- df_plot %>%
  dplyr::mutate(
    Condition = dplyr::case_when(
      grepl("WT_G_", df_plot$Sample) ~ "G WT",
      grepl("KO_G_", df_plot$Sample) ~ "G KO",
      grepl("WT_G12D_", df_plot$Sample) ~ "G12D WT",
      grepl("Het_G12D_", df_plot$Sample) ~ "G12D Het",
      grepl("KO_G12D_", df_plot$Sample) ~ "G12D KO",
    )
  )

#attach colours to condition
df_plot <- df_plot %>%
  dplyr::mutate(
    Colour = dplyr::case_when(
      grepl("WT_G_", df_plot$Sample) ~ "#0e6e80",
      grepl("KO_G_", df_plot$Sample) ~ "#ce1e68",
      grepl("WT_G12D_", df_plot$Sample) ~ "#3b54a3",
      grepl("Het_G12D_", df_plot$Sample) ~ "#f47d23",
      grepl("KO_G12D_", df_plot$Sample) ~ "#bbbd32",
    )
  )

#set condition factor level order
df_plot$Condition <- factor(df_plot$Condition, levels = sort(unique(df_plot$Condition))[c(2, 1, 5, 3, 4)])

#convert condition factor as numeric
df_plot$x_position <- as.numeric(factor(df_plot$Condition))

#calculate se and mean for each condition
data_matrix_se <- df_plot %>%
  dplyr::group_by(Condition) %>%
  dplyr::summarise(se = sd(`Percent Digestive Enzymes out of Total`)/sqrt(length(`Percent Digestive Enzymes out of Total`)),
                   mean = mean(`Percent Digestive Enzymes out of Total`),
                   x_position = x_position
  ) %>%
  unique()

#create plot
ggplot2::ggplot(df_plot, ggplot2::aes(colour = Condition)) +
  
  ggplot2::coord_cartesian(ylim = c(0, 20), expand = FALSE, clip = "off") +
  
  ggplot2::geom_bar(ggplot2::aes(y = `Percent Digestive Enzymes out of Total`, x = x_position),
                    stat = "summary", 
                    fun = "mean",
                    alpha =  0,
                    fill = NA,
                    colour = NA) +
  
  ggplot2::geom_bar(ggplot2::aes(y = `Percent Digestive Enzymes out of Total`, x = x_position),
                    stat = "summary", 
                    fun = "mean", 
                    width = 0.7, 
                    size = 1,
                    colour = c("#0e6e80",
                               "#ce1e68",
                               "#3b54a3",
                               "#f47d23",
                               "#bbbd32"),
                    fill = NA) +
  
  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, 
                                      ymax = mean + se, 
                                      x = x_position,
                                      colour = Condition), 
                         width = 0.4, 
                         size = 1,
                         data = data_matrix_se) +
  
  ggplot2::geom_jitter(ggplot2::aes(y = `Percent Digestive Enzymes out of Total`, x = x_position, shape = Condition, colour = Condition),
                       size = 1.5,
                       fill = NA,
                       stroke = 1.5,
                       width = 0.15) + 
  
  ggplot2::geom_text(ggplot2::aes(x = (x + xend)/2, y = y + 0.2, label = annotation),
                     colour = "black",
                     size = 6,
                     vjust = -0.4,
                     data = data.frame(x = c(1, 3, 3), xend = c(5, 4, 5),
                                       y = c(28, 24.5, 21), annotation = c("✱✱", "✱", "✱✱"))) +
  
  ggplot2::geom_segment(ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                        colour = "black",
                        size = 1,
                        data = data.frame(x = c(1, 3, 3), xend = c(5, 4, 5),
                                          y = c(28, 24.5, 21))) +
  
  ggplot2::labs(y = "Digestive enzyme abundance\n(% Total proteins)") +
  
  ggplot2::scale_colour_manual(labels = c(expression(italic("Ptf1a")^CreER*";"*italic("Insr")^"w/w"),
                                          expression(italic("Ptf1a")^CreER*";"*italic("Insr")^"f/f"),
                                          expression("PK-"*italic("Insr")^"w/w"),
                                          expression("PK-"*italic("Insr")^"w/f"),
                                          expression("PK-"*italic("Insr")^"f/f")),
                               values = c("#0e6e80",
                                          "#ce1e68",
                                          "#3b54a3",
                                          "#f47d23",
                                          "#bbbd32")) +
  
  ggplot2::scale_shape_manual(labels = c(expression(italic("Ptf1a")^CreER*";"*italic("Insr")^"w/w"),
                                         expression(italic("Ptf1a")^CreER*";"*italic("Insr")^"f/f"),
                                         expression("PK-"*italic("Insr")^"w/w"),
                                         expression("PK-"*italic("Insr")^"w/f"),
                                         expression("PK-"*italic("Insr")^"f/f")),
                              values = c(1, 0, 2, 5, 6)) +
  
  ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(c(1, 0, 2, 5, 6), 
                                                                    colour = c("#0e6e80",
                                                                               "#ce1e68",
                                                                               "#3b54a3",
                                                                               "#f47d23",
                                                                               "#bbbd32"),
                                                                    fill = "white",
                                                                    stroke = 1.5,
                                                                    size = 1.5)),
                  colour = "none",
                  alpha = "none") +
  
  ggplot2::scale_x_continuous(limits = c(0, 6)) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 0.5, unit = "cm")),
    plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0.25, "cm"),
    axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_line(colour = "black", lineend = "square", size = 1),
    axis.line.x = ggplot2::element_line(colour = "black", lineend = "square", size = 1),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 12),
    panel.border = ggplot2::element_blank(),
    axis.line.y.left = ggplot2::element_line(colour = "black", lineend = "square", size = 1),
    axis.title.x = ggplot2::element_blank(),
    axis.text.x.bottom = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.text.align = 0,
    text = ggplot2::element_text(size = 12, colour = "black"),
    aspect.ratio = 1 / 1
  )


### Figure 7 ###

## Fig 7A ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import file
df_pproteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_phosphoproteomics.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_pproteomics) <- df_pproteomics[1,]
df_pproteomics <- df_pproteomics[-1,]

#selecting and reordering columns
df_pproteomics <- df_pproteomics[, c(1:2, 16, 20, 57:75)]
df_pproteomics <- cbind(df_pproteomics[, c(1:4)],
                        df_pproteomics[, c(5:23)][, grepl("Sample, G_WT", colnames(df_pproteomics)[c(5:23)])],
                        df_pproteomics[, c(5:23)][, grepl("Sample, G_KO", colnames(df_pproteomics)[c(5:23)])],
                        df_pproteomics[, c(5:23)][, grepl("Sample, G12D_WT", colnames(df_pproteomics)[c(5:23)])],
                        df_pproteomics[, c(5:23)][, grepl("Sample, G12D_HET", colnames(df_pproteomics)[c(5:23)])],
                        df_pproteomics[, c(5:23)][, grepl("Sample, G12D_KO", colnames(df_pproteomics)[c(5:23)])])

#covert numbers to numeric
df_pproteomics[, c(4:23)] <- sapply(df_pproteomics[, c(4:23)], as.numeric)

#calculate fold change
df_pproteomics$`(G12D_KO) / (G12D_WT)` <- rowMedians(as.matrix(df_pproteomics[, grepl("Sample, G12D_KO", colnames(df_pproteomics))]))/rowMedians(as.matrix(df_pproteomics[, grepl("Sample, G12D_WT", colnames(df_pproteomics))]))

#select needed columns
df_pproteomics <- df_pproteomics[, c("Gene", "(G12D_KO) / (G12D_WT)", "Adj. P-Value: (G12D_KO) / (G12D_WT)")] %>% na.omit()

#get log2 of fold change
df_pproteomics$log2FoldChange <- log2(df_pproteomics[, 2])

#rename columns
colnames(df_pproteomics) <- c("Gene", "FoldChange", "padj", "log2FoldChange")

#label whether protein is differentially expressed
df_pproteomics$DiffExpression <- "ns"
df_pproteomics$DiffExpression[df_pproteomics$padj < 0.05 & df_pproteomics$log2FoldChange > 0] <- "Up in KO"
df_pproteomics$DiffExpression[df_pproteomics$padj < 0.05 & df_pproteomics$log2FoldChange < 0] <- "Down in KO"

#denote significance and fold change with colour
df_pproteomics$Colour[df_pproteomics$padj < 0.05 & df_pproteomics$log2FoldChange > 0] <- "red"
df_pproteomics$Colour[df_pproteomics$padj < 0.05 & df_pproteomics$log2FoldChange < 0] <- "blue"

#get -log10 of padj
df_pproteomics$transf.log.10.padj <- -log10(df_pproteomics$padj)

#mark proteins for labeling
df_pproteomics$Label[df_pproteomics$transf.log.10.padj > 10] <- "Yes"

#set expression variable factor level order
df_pproteomics$DiffExpression <- factor(df_pproteomics$DiffExpression, levels = sort(unique(df_pproteomics$DiffExpression))[c(3, 1, 2)])

#create plot
ggplot(df_pproteomics, aes(x = log2FoldChange, y = transf.log.10.padj)) +
  
  ggtitle(expression("PK-"*italic("Insr")^"w/w"*" vs "*"PK-"*italic("Insr")^"f/f")) + 
  
  geom_point(aes(colour = DiffExpression)) +
  
  geom_text_repel(aes(label = Gene),
                  seed = 123,
                  max.overlaps = 200,
                  min.segment.length = 0.2,
                  segment.alpha = 0.3,
                  segment.colour = "black",
                  colour = df_pproteomics$Colour[!is.na(df_pproteomics$Label)],
                  bg.color = "white",
                  bg.r = 0.05,
                  nudge_y = -0.1,
                  force = 20,
                  force_pull = 0,
                  data = df_pproteomics %>% na.omit()) +
  
  coord_cartesian(xlim = c(-10, 10), 
                  ylim = c(-1, 16.5), 
                  expand = FALSE, clip = "off") +
  
  scale_color_manual(values = c("red", "blue", "grey"), #c("#44CC44", "#4B0092", "grey"),
                     labels = c(expression("Up in PK-"*italic("Insr")^"f/f"),
                                expression("Down in PK-"*italic("Insr")^"f/f"),
                                "ns")) +
  
  labs(y = bquote(-log[10]*"(Adj. "*italic(P)*"-value)"), x = bquote(log[2]*"(Fold Change)"), colour = "Regulation") +
  
  guides(colour = guide_legend(override.aes = list(shape = 16))) +
  
  theme_bw() +
  
  theme(plot.title = element_text(hjust = 0.5, vjust = -1.5, size = 15, margin = margin(b = 0.5, unit = "cm", )),
        plot.margin = margin(t = 1, r = 0.5, b = 0.5, l = 0.5, "cm"),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.title.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        legend.title.align = 0,
        legend.text.align = 0,
        legend.direction = "vertical",
        legend.box.just = "center",
        legend.key.width = unit(0.75, "cm"),
        legend.margin = margin(t = 0.5, r = 0, b = 0, l = 0, "cm"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, "cm"),
        aspect.ratio = 1/1)


## Fig 7B ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import total protein file
df_total_proteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_total_protein.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_total_proteomics) <- df_total_proteomics[1,]
df_total_proteomics <- df_total_proteomics[-1,]

#select needed columns
df_total_proteomics <- df_total_proteomics[, c(3:4, 21:22)]

#covert numbers to numeric
df_total_proteomics[, c(3:4)] <- as.data.frame(sapply(df_total_proteomics[, c(3:4)], as.numeric))

#import phosphoproteomics protein file
df_pproteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_phosphoproteomics_all.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_pproteomics) <- df_pproteomics[1,]
df_pproteomics <- df_pproteomics[-1,]

#select needed columns
df_pproteomics <- df_pproteomics[, c(4:5, 12, 14:15, 17:18)]

#remove rows with no protein
df_pproteomics <- subset(df_pproteomics, !df_pproteomics$`Master Protein Accessions` == "")

#covert numbers to numeric
df_pproteomics[, c(6:7)] <- as.data.frame(sapply(df_pproteomics[, c(6:7)], as.numeric))

#separate rows with multiple proteins
df_pproteomics <- tidyr::separate_rows(df_pproteomics, `Master Protein Accessions`, sep = "; ")

#separate rows with multiple peptide position ranges
df_pproteomics <- tidyr::separate_rows(df_pproteomics, `Positions in Master Proteins`, sep = "; ")

#convert Uniprot to gene name
df_pproteomics$Gene <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = df_pproteomics$`Master Protein Accessions`, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")

#extract peptide position range
df_pproteomics$Position <- stringi::stri_join_list(regmatches(df_pproteomics$`Positions in Master Proteins`, gregexpr("(?<=\\[).+?(?=\\])", df_pproteomics$`Positions in Master Proteins`, perl = T)))

#extract peptide starting position
df_pproteomics$Position_Start <- as.numeric(stringi::stri_join_list(regmatches(df_pproteomics$`Positions in Master Proteins`, gregexpr("(?<=\\[).+?(?=\\-)", df_pproteomics$`Positions in Master Proteins`, perl = T))))

#extract phospho-site location relative to peptide
df_pproteomics$Last_Phospho_Location <- stringi::stri_locate_last_fixed(df_pproteomics$Modifications, "xPhospho")

#extract phospho modification
df_pproteomics$Last_Phospho <- substr(df_pproteomics$Modifications, df_pproteomics$Last_Phospho_Location - 1, nchar(df_pproteomics$Modifications))

#extract phospho-site
df_pproteomics$Last_Phospho_Info <- stringi::stri_join_list(regmatches(df_pproteomics$Last_Phospho, gregexpr("(?<=\\[).+?(?=\\])", df_pproteomics$Last_Phospho, perl = T)))

#remove confidence from phospho-site
df_pproteomics$Last_Phospho_Info <- gsub("\\((.*?))", "", df_pproteomics$Last_Phospho_Info)

#separate rows with multiple phospho-sites
df_pproteomics <- tidyr::separate_rows(df_pproteomics, Last_Phospho_Info, sep = "; ")

#extract phospho-site in protein
df_pproteomics$Phospho_Numeric <- ifelse(grepl("\\d", df_pproteomics$Last_Phospho_Info), paste0(stringr::str_extract(df_pproteomics$Last_Phospho_Info, "[^\\d]+"), df_pproteomics$Position_Start + as.numeric(stringr::str_extract(df_pproteomics$Last_Phospho_Info, "[\\d^]+")) - 1), df_pproteomics$Last_Phospho_Info)

#create Uniprot with phospho-site column
df_pproteomics$Sum_Index_Unique <- ifelse(grepl("\\d", df_pproteomics$Last_Phospho_Info), paste(df_pproteomics$`Master Protein Accessions`, df_pproteomics$Phospho_Numeric), paste(df_pproteomics$`Positions in Master Proteins`, df_pproteomics$Last_Phospho_Info))

#create gene name with phospho-site column
df_pproteomics$Sum_Index_Unique_Gene <- ifelse(grepl("\\d", df_pproteomics$Last_Phospho_Info), paste(df_pproteomics$Gene, df_pproteomics$Phospho_Numeric), paste(df_pproteomics$Gene, df_pproteomics$`Positions in Master Proteins`, df_pproteomics$Last_Phospho_Info))

#only include phospho-sites with defined numerical positions
df_pproteomics <- df_pproteomics[grepl("\\d", df_pproteomics$Last_Phospho_Info),]

#annotate total protein columns
colnames(df_total_proteomics)[3:4] <- paste(colnames(df_total_proteomics)[3:4], "total", sep = "_")

#annotation phospho-site columns
colnames(df_pproteomics)[6:7] <- paste(colnames(df_pproteomics)[6:7], "phospho", sep = "_") 

#rename total protein column
colnames(df_total_proteomics)[1] <- "Master Protein Accessions"

#merge phospho-site and total protein by protein
df_merge <- merge(df_pproteomics[, c(3, 6:8, 14, 16)], df_total_proteomics, by = "Master Protein Accessions")

#get log2 of total fold change
df_merge$`log2 (G12D_KO) / (G12D_WT)_total` <- log2(df_merge$`(G12D_KO) / (G12D_WT)_total`)

#get -log10 of total padj
df_merge$`-log10 Adj. P-Value: (G12D_KO) / (G12D_WT)_total` <- -log10(df_merge$`Adj. P-Value: (G12D_KO) / (G12D_WT)_total`)

#get log2 of phospho-site fold change
df_merge$`log2 (G12D_KO) / (G12D_WT)_phospho` <- log2(df_merge$`(G12D_KO) / (G12D_WT)_phospho`)

#get -log10 of phospho-site padj
df_merge$`-log10 Adj. P-Value: (G12D_KO) / (G12D_WT)_phospho` <- -log10(df_merge$`Adj. P-Value: (G12D_KO) / (G12D_WT)_phospho`)

#replace NA with 0 for -log10 of phsopho Adj p-values
df_merge$`-log10 Adj. P-Value: (G12D_KO) / (G12D_WT)_phospho`[is.na(df_merge$`-log10 Adj. P-Value: (G12D_KO) / (G12D_WT)_phospho`)] <- 0

#initialize label column
df_merge$Label <- ""

#add genes with abs(log2FC) > 5.7 to labels
df_merge$Label[abs(df_merge$`log2 (G12D_KO) / (G12D_WT)_phospho`) > 5.7] <- df_merge$`Gene Symbol`[abs(df_merge$`log2 (G12D_KO) / (G12D_WT)_phospho`) > 5.7]

#remove duplicate labels for log2FC > 5.7
df_merge$Label[df_merge$`log2 (G12D_KO) / (G12D_WT)_phospho` > 5.7][duplicated(df_merge$Label[df_merge$`log2 (G12D_KO) / (G12D_WT)_phospho` > 5.7])] <- ""

#remove duplicate labels for log2FC < 5.7
df_merge$Label[df_merge$`log2 (G12D_KO) / (G12D_WT)_phospho` < -5.7][duplicated(df_merge$Label[df_merge$`log2 (G12D_KO) / (G12D_WT)_phospho` < -5.7])] <- ""

#create plot
ggplot2::ggplot(df_merge, ggplot2::aes(x = `log2 (G12D_KO) / (G12D_WT)_total`, y = `log2 (G12D_KO) / (G12D_WT)_phospho`)) +
  
  ggplot2::geom_point(ggplot2::aes(fill = `-log10 Adj. P-Value: (G12D_KO) / (G12D_WT)_phospho`),
                      shape = 21,
                      colour = "darkgrey") +
  
  ggrepel::geom_text_repel(ggplot2::aes(label = Label),
                           seed = 123,
                           max.overlaps = Inf,
                           min.segment.length = 0.2,
                           segment.alpha = 0.3,
                           segment.colour = "black",
                           bg.color = "white",
                           bg.r = 0.05,
                           nudge_y = -0.1,
                           ylim = c(-4, -8),
                           force = 21,
                           force_pull = 0,
                           data = df_merge[df_merge$`log2 (G12D_KO) / (G12D_WT)_phospho` < 0,]) +
  
  ggrepel::geom_text_repel(ggplot2::aes(label = Label),
                           seed = 123,
                           max.overlaps = Inf,
                           min.segment.length = 0.2,
                           segment.alpha = 0.3,
                           segment.colour = "black",
                           bg.color = "white",
                           bg.r = 0.05,
                           nudge_y = -0.1,
                           ylim = c(5, NA),
                           force = 21,
                           force_pull = 0,
                           data = df_merge[df_merge$`log2 (G12D_KO) / (G12D_WT)_phospho` > 0,]) +
  
  ggplot2::coord_cartesian(xlim = c(-6, 6), 
                           ylim = c(-10, 10), 
                           expand = FALSE, clip = "off") +
  
  ggplot2::scale_fill_continuous(high = "red", low = "white", limits = c(0, 15.2)) +
  
  ggplot2::labs(y = expression(log[2]*"(Fold Change) Phosphosite"), x = expression(log[2]*"(Fold Change) Total Protein"), fill = expression(atop(-log[10]*"(Adj. "*italic(P)*"-value)", "Phosphosite"))) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, vjust = -1.5, size = 15, margin = ggplot2::margin(b = 0.5, unit = "cm", )),
                 plot.margin = ggplot2::margin(t = 1, r = 0.5, b = 0.5, l = 0.5, "cm"),
                 axis.text.x = ggplot2::element_text(size = 12, colour = "black"),
                 axis.text.y = ggplot2::element_text(size = 11, colour = "black"),
                 axis.title.x = ggplot2::element_text(size = 12),
                 panel.grid.minor = ggplot2::element_blank(),
                 legend.title.align = 0,
                 legend.text.align = 0,
                 legend.direction = "vertical",
                 legend.box.just = "center",
                 legend.key.width = grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = 1/1)


## Fig 7C ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import file
df_pproteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_phosphoproteomics.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_pproteomics) <- df_pproteomics[1,]
df_pproteomics <- df_pproteomics[-1,]

#select needed columns
df_pproteomics <- df_pproteomics[, c(1:2, 16, 18, 20, 22, 24, 26, 28, 30, 32, 57:75)]

#convert numbers to numeric
df_pproteomics[, c(4:30)] <- sapply(df_pproteomics[, c(4:30)], as.numeric)

#filter for phospho-sites on Cckar
df_pproteomics <- df_pproteomics[df_pproteomics$Gene == "Cckar",]

#calculate log2FC for each comparison
df_pproteomics$`log2 Fold Change (G12D_KO) / (G12D_WT)` <- log2(df_pproteomics$`(G12D_KO) / (G12D_WT)`)
df_pproteomics$`log2 Fold Change (G12D_HET) / (G12D_WT)` <- log2(df_pproteomics$`(G12D_HET) / (G12D_WT)`)
df_pproteomics$`log2 Fold Change (G_KO) / (G_WT)` <- log2(df_pproteomics$`(G_KO) / (G_WT)`)
df_pproteomics$`log2 Fold Change (G12D_WT) / (G_WT)` <- log2(df_pproteomics$`(G12D_WT) / (G_WT)`)

#remove rows with no significant comparisons
df_pproteomics <- df_pproteomics[rowSums(df_pproteomics[, grepl("Adj. P-Value", colnames(df_pproteomics))] < 0.05) > 0,]

#remove unneeded annotations
df_pproteomics$`Modifications in Master Proteins` <- gsub("(.*)xPhospho ", "", df_pproteomics$`Modifications in Master Proteins`)

#extract phospho-site
df_pproteomics$`Modifications in Master Proteins` <- stringi::stri_join_list(regmatches(df_pproteomics$`Modifications in Master Proteins`, gregexpr("(?<=\\[).+?(?=\\])", df_pproteomics$`Modifications in Master Proteins`, perl = T)))

#remove confidence from phospho-site
df_pproteomics$`Modifications in Master Proteins` <- gsub("\\((.*?))", "", df_pproteomics$`Modifications in Master Proteins`)

#select needed columns
df_plot <- df_pproteomics[, grepl("\\(Normalized\\): ", colnames(df_pproteomics)) | colnames(df_pproteomics) == "Modifications in Master Proteins"]

#lengthen plot values into single column
df_plot <- tidyr::pivot_longer(df_plot, names_to = "Sample", values_to = "Abundance", cols = -`Modifications in Master Proteins`)

#annotate condition
df_plot <- df_plot %>%
  dplyr::mutate(
    Condition = dplyr::case_when(
      grepl("G_WT", df_plot$Sample) ~ "G WT",
      grepl("G_KO", df_plot$Sample) ~ "G KO",
      grepl("G12D_WT", df_plot$Sample) ~ "G12D WT",
      grepl("G12D_HET", df_plot$Sample) ~ "G12D Het",
      grepl("G12D_KO", df_plot$Sample) ~ "G12D KO",
    )
  )

#annotate graph colour
df_plot <- df_plot %>%
  dplyr::mutate(
    Colour = dplyr::case_when(
      grepl("G_WT", df_plot$Sample) ~ "#0e6e80",
      grepl("G_KO", df_plot$Sample) ~ "#ce1e68",
      grepl("G12D_WT", df_plot$Sample) ~ "#3b54a3",
      grepl("G12D_HET", df_plot$Sample) ~ "#f47d23",
      grepl("G12D_KO", df_plot$Sample) ~ "#bbbd32",
    )
  )

#convert adundances to numeric
df_plot$Abundance <- as.numeric(df_plot$Abundance)

#change dot position slightly to ensure no direct overlap
df_plot$Abundance <- df_plot$Abundance + c(1:nrow(df_plot))/100

#set condition factor level order
df_plot$Condition <- factor(df_plot$Condition, levels = sort(unique(df_plot$Condition))[c(2, 1, 5, 3, 4)])

#create numerical x-axis column
df_plot$x_position <- as.numeric(factor(df_plot$Condition))

#calculate stats
data_matrix_se <- df_plot %>%
  dplyr::group_by(Condition, `Modifications in Master Proteins`) %>%
  dplyr::summarise(se = sd(Abundance)/sqrt(length(Abundance)),
                   mean = mean(Abundance),
                   x_position = x_position)

#create comparison significance annotations
anno_text_1 <- data.frame(x1 = c(1, 1, 1, 1, 3, 3), 
                          y = c(1.05 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S256; S/Y"]),
                                0.8 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S256; S265"]),
                                1.03 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S295; S296"]),
                                0.89 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S420; S425"]),
                                1.05 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S422; S425"]),
                                1.11 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S425; S430; S/Y/T"])),
                          x2 = c(3, 3, 3, 2, 4, 4), 
                          xstar = c(2, 2, 2, 1.5, 3.5, 3.5), 
                          ystar = c(1.07 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S256; S/Y"]),
                                    0.82 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S256; S265"]),
                                    1.05 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S295; S296"]),
                                    0.91 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S420; S425"]),
                                    1.07 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S422; S425"]),
                                    1.13 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S425; S430; S/Y/T"])),
                          alpha = 1,
                          lab = c("✱✱✱", "✱✱✱✱", "✱✱", "✱✱", "✱", "✱"),
                          `Modifications in Master Proteins` = sort(unique(df_plot$`Modifications in Master Proteins`)),
                          check.names = F)

anno_text_2 <- data.frame(x1 = c(1, 3, 1, 1, 3, 3), 
                          y = c(1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S256; S/Y"]),
                                1.05 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S256; S265"]),
                                1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S295; S296"]),
                                1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S420; S425"]),
                                0.8 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S422; S425"]),
                                1.04 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S425; S430; S/Y/T"])),
                          x2 = c(1, 5, 1, 1, 5, 5), 
                          xstar = c(1, 4, 1, 1, 4, 4), 
                          ystar = c(1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S256; S/Y"]),
                                    1.07 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S256; S265"]),
                                    1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S295; S296"]),
                                    1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S420; S425"]),
                                    0.82 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S422; S425"]),
                                    1.06 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S425; S430; S/Y/T"])),
                          alpha = c(0, 1, 0, 0, 1, 1),
                          lab = c("", "✱✱✱✱", "", "", "✱✱✱✱", "✱"),
                          `Modifications in Master Proteins` = sort(unique(df_plot$`Modifications in Master Proteins`)),
                          check.names = F)

anno_text_3 <- data.frame(x1 = c(1, 1, 1, 1, 1, 1), 
                          y = c(1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S256; S/Y"]),
                                0.73 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S256; S265"]),
                                1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S295; S296"]),
                                1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S420; S425"]),
                                1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S422; S425"]),
                                1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S425; S430; S/Y/T"])),
                          x2 = c(1, 2, 1, 1, 1, 1), 
                          xstar = c(1, 1.5, 1, 1, 1, 1), 
                          ystar = c(1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S256; S/Y"]),
                                    0.75 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S256; S265"]),
                                    1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S295; S296"]),
                                    1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S420; S425"]),
                                    1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S422; S425"]),
                                    1 * max(df_plot$Abundance[df_plot$`Modifications in Master Proteins` == "S425; S430; S/Y/T"])),
                          alpha = c(0, 1, 0, 0, 0, 0),
                          lab = c("", "✱", "", "", "", ""),
                          `Modifications in Master Proteins` = sort(unique(df_plot$`Modifications in Master Proteins`)),
                          check.names = F)

#create plot
ggplot2::ggplot(df_plot, ggplot2::aes(colour = Condition)) +
  
  ggplot2::geom_bar(ggplot2::aes(y = Abundance, x = x_position),
                    stat = "summary", 
                    fun = "mean",
                    alpha =  0,
                    fill = NA,
                    colour = NA) +
  
  ggplot2::geom_bar(ggplot2::aes(y = Abundance, x = x_position),
                    stat = "summary", 
                    fun = "mean", 
                    width = 0.93, 
                    colour = rep(c("#0e6e80",
                                   "#ce1e68",
                                   "#3b54a3",
                                   "#f47d23",
                                   "#bbbd32"), 6),
                    fill = "white") +
  
  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, 
                                      ymax = mean + se, 
                                      x = x_position,
                                      colour = Condition), 
                         width = 0.4, 
                         data = data_matrix_se) +
  
  ggplot2::geom_jitter(ggplot2::aes(y = Abundance, x = x_position, shape = Condition, colour = Condition),
                       size = 1,
                       fill = NA,
                       stroke = 0.75,
                       width = 0.15) + 
  
  ggplot2::geom_text(data = anno_text_1, 
                     ggplot2::aes(x = (x1 + x2)/2,  
                                  y = ystar, 
                                  label = lab),
                     colour = "black",
                     hjust = 0.5,
                     vjust = 0,
                     size = 3) +
  
  ggplot2::geom_segment(data = anno_text_1, 
                        ggplot2::aes(x = x1, 
                                     xend = x2, 
                                     y = y, 
                                     yend = y),
                        colour = "black") +
  
  ggplot2::geom_text(data = anno_text_2, 
                     ggplot2::aes(x = (x1 + x2)/2,  
                                  y = ystar, 
                                  label = lab, 
                                  alpha = alpha),
                     colour = "black",
                     hjust = 0.5,
                     vjust = 0,
                     size = 3) +
  
  ggplot2::geom_segment(data = anno_text_2, 
                        ggplot2::aes(x = x1, 
                                     xend = x2, 
                                     y = y, 
                                     yend = y, 
                                     alpha = alpha),
                        colour = "black") +
  
  ggplot2::geom_text(data = anno_text_3, 
                     ggplot2::aes(x = (x1 + x2)/2,  
                                  y = ystar, 
                                  label = lab, 
                                  alpha = alpha),
                     colour = "black",
                     hjust = 0.5,
                     vjust = 0,
                     size = 3) +
  
  ggplot2::geom_segment(data = anno_text_3, 
                        ggplot2::aes(x = x1, 
                                     xend = x2, 
                                     y = y, 
                                     yend = y, 
                                     alpha = alpha), 
                        colour = "black") +
  
  ggh4x::facet_wrap2(.~`Modifications in Master Proteins`, 
                     strip = strip_vanilla(clip = "off"),
                     strip.position = "bottom",
                     scales = "free_y",
                     ncol = 6, 
                     labeller = ggplot2::label_wrap_gen(width = 30)) +
  
  ggplot2::scale_colour_manual(labels = c(expression(italic("Ptf1a")^CreER*";"*italic("Insr")^"w/w"),
                                          expression(italic("Ptf1a")^CreER*";"*italic("Insr")^"f/f"),
                                          expression("PK-"*italic("Insr")^"w/w"),
                                          expression("PK-"*italic("Insr")^"w/f"),
                                          expression("PK-"*italic("Insr")^"f/f")),
                               values = c("#0e6e80",
                                          "#ce1e68",
                                          "#3b54a3",
                                          "#f47d23",
                                          "#bbbd32")) +
  
  ggplot2::scale_shape_manual(labels = c(expression(italic("Ptf1a")^CreER*";"*italic("Insr")^"w/w"),
                                         expression(italic("Ptf1a")^CreER*";"*italic("Insr")^"f/f"),
                                         expression("PK-"*italic("Insr")^"w/w"),
                                         expression("PK-"*italic("Insr")^"w/f"),
                                         expression("PK-"*italic("Insr")^"f/f")),
                              values = c(1, 0, 2, 5, 6)) +
  
  ggplot2::labs(y = "Phospho-peptide Intensity (AU)", shape = "Condition ") +
  
  ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(c(1, 0, 2, 5, 6), 
                                                                    colour = c("#0e6e80",
                                                                               "#ce1e68",
                                                                               "#3b54a3",
                                                                               "#f47d23",
                                                                               "#bbbd32"),
                                                                    fill = "white",
                                                                    size = 1)),
                  colour = "none",
                  alpha = "none") +
  
  ggplot2::scale_x_continuous(limits = c(0, 6)) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(
    plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0.25, "cm"),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 12),
    panel.border = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.text.x.bottom = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.text.align = 0,
    text = ggplot2::element_text(size = 12),
    aspect.ratio = 3 / 1)

## Fig 7D ##

#clear environment
rm(list = ls())

#import file
df_pproteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_phosphoproteomics.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_pproteomics) <- df_pproteomics[1,]
df_pproteomics <- df_pproteomics[-1,]

#select columns
df_pproteomics <- df_pproteomics[, c("Gene", 
                                     "Annotated Sequence", 
                                     "Modifications", 
                                     "Modifications in Master Proteins", 
                                     "Master Protein Accessions", 
                                     "Positions in Master Proteins", 
                                     "(G12D_KO) / (G12D_WT)", 
                                     "Adj. P-Value: (G12D_KO) / (G12D_WT)")]

#convert numbers to numeric
df_pproteomics[, c("(G12D_KO) / (G12D_WT)", "Adj. P-Value: (G12D_KO) / (G12D_WT)")] <- as.data.frame(sapply(df_pproteomics[, c("(G12D_KO) / (G12D_WT)", "Adj. P-Value: (G12D_KO) / (G12D_WT)")], as.numeric))

#get log2 of fold change
df_pproteomics$`log2 (G12D_KO) / (G12D_WT)` <- log2(df_pproteomics$`(G12D_KO) / (G12D_WT)`)

#get -log10 of adj. p-value
df_pproteomics$`-log10 Adj. P-Value: (G12D_KO) / (G12D_WT)` <- -log10(df_pproteomics$`Adj. P-Value: (G12D_KO) / (G12D_WT)`)

#subset for significance and select proteins
df_pproteomics <- subset(df_pproteomics, df_pproteomics$`Adj. P-Value: (G12D_KO) / (G12D_WT)` < 0.05)

#extract individual phospho-sites
df_pproteomics$`Modifications in Master Proteins` <- gsub("(.*)xPhospho ", "", df_pproteomics$`Modifications in Master Proteins`)

df_pproteomics$`Modifications in Master Proteins` <- stringi::stri_join_list(regmatches(df_pproteomics$`Modifications in Master Proteins`, gregexpr("(?<=\\[).+?(?=\\])", df_pproteomics$`Modifications in Master Proteins`, perl = T)))

df_pproteomics$`Modifications in Master Proteins` <- gsub("\\((.*?))", "", df_pproteomics$`Modifications in Master Proteins`)

#create new row per new phospho-site
df_pproteomics <- tidyr::separate_rows(df_pproteomics, `Modifications in Master Proteins`, sep = "; ")

#remove non-numbered phospho-sites
df_pproteomics <- df_pproteomics[grepl("\\d", df_pproteomics$`Modifications in Master Proteins`),]

#combine gene and phospho-site
df_pproteomics$Sum_Index_Unique_Gene <- paste(df_pproteomics$Gene, df_pproteomics$`Modifications in Master Proteins`)

#remove rows with empty gene
df_pproteomics <- subset(df_pproteomics, !Gene == "")

#order by Adj. P-value and abs(log2FC)
df_pproteomics <- df_pproteomics[with(df_pproteomics, order(df_pproteomics$Sum_Index_Unique_Gene, df_pproteomics$`Adj. P-Value: (G12D_KO) / (G12D_WT)`, -abs(as.numeric(df_pproteomics$`log2 (G12D_KO) / (G12D_WT)`)))),]

#remove duplicates
df_pproteomics <- df_pproteomics[!duplicated(df_pproteomics$Sum_Index_Unique_Gene),]

#add anchors
anchors <- df_pproteomics[0,]

anchor_genes <- c("Akt3", "Kras", "Igf1r", "Insr", "mTor", "Raf1", "Mapk1")

anchors <- anchors[1:length(anchor_genes),]

anchors[, 1:ncol(anchors)] <- ""

anchors[, 1] <- anchor_genes

df_pproteomics <- rbind(df_pproteomics, 
                        anchors)

#rename gene column name
colnames(df_pproteomics)[1] <- "#node1"

#export node table to SCPC#46_03_Jim_Anni_pancreatic_cancer_phosphoproteomics_G12D_WT_vs_KO sig logFC and padj.tsv
write.table(df_pproteomics)

#import STRING protein mapping and edge table output files
df_string_network <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_phosphoproteomics_G12D_WT_vs_KO sig string_interactions.tsv", check.names = FALSE, sep = "\t")
df_string_mapping <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_phosphoproteomics_G12D_WT_vs_KO sig string_mapping.tsv", check.names = FALSE, sep = "\t")

#select name mapping columns
df_string_mapping <- df_string_mapping[, c("queryItem", "preferredName")]

#rename column names
colnames(df_string_mapping) <- c("queryItem2", "node2")

#merge mapping to network
df_string_network <- merge(df_string_network, df_string_mapping, by = c("node2"), all.x = TRUE)

#rename column names
colnames(df_string_mapping) <- c("queryItem", "#node1")

#merge mapping to network
df_string_network <- merge(df_string_network, df_string_mapping, by = c("#node1"), all.x = TRUE)

#replace mapped values with input values
df_string_network$`#node1` <- df_string_network$queryItem
df_string_network$node2 <- df_string_network$queryItem2

#create new nodes dataframe
df_pproteomics_nodes <- df_pproteomics

#create site column
df_pproteomics_nodes$Site <- df_pproteomics_nodes$`Modifications in Master Proteins`

#create new network for phospho-sites
df_string_network_sites <- df_string_network

df_string_network_sites[c(1:15)] <- ""

df_string_network_sites <- df_string_network_sites[1:nrow(df_pproteomics_nodes),]

#create edge for every phospho-site to respective protein 
df_string_network_sites$`#node1` <- df_pproteomics_nodes$`#node1`
df_string_network_sites$node2 <- df_pproteomics_nodes$Sum_Index_Unique_Gene

#remove empty rows
df_string_network_sites <- subset(df_string_network_sites, !node2 == "")

#categorize interaction
df_string_network_sites$Interaction_Type <- "Phospho-site"
df_string_network$Interaction_Type <- "PPI"

#add phospho-site and protein edges to PPI data frame
df_string_network <- rbind(df_string_network, df_string_network_sites)

#remove duplicated rows
df_string_network <- df_string_network[!duplicated(df_string_network),]

#export edge table to SCPC#46_03_Jim_Anni_pancreatic_cancer_phosphoproteomics_G12D_WT_vs_KO sig string_interactions remapped extra site nodes.tsv
write.table(df_string_network)

#create phospho-site node data frame
df_pproteomics_nodes_sites <- df_pproteomics_nodes

#fill extra columns with blanks
df_pproteomics_nodes[c(2:ncol(df_pproteomics_nodes))] <- ""

#categorize node types
df_pproteomics_nodes$Node_Type <- "Protein"
df_pproteomics_nodes_sites$Node_Type <- "Phospho-site"

#add node label column
df_pproteomics_nodes$Node_Label <- df_pproteomics_nodes$`#node1`
df_pproteomics_nodes_sites$Node_Label <- df_pproteomics_nodes_sites$Site

#remove duplicate rows
df_pproteomics_nodes <- df_pproteomics_nodes[!duplicated(df_pproteomics_nodes),]

#replace node column with gene and site 
df_pproteomics_nodes_sites$`#node1` <- df_pproteomics_nodes_sites$Sum_Index_Unique_Gene

#remove empty rows
df_pproteomics_nodes_sites <- subset(df_pproteomics_nodes_sites, !`#node1` == "")

#add phospho-site and protein nodes to protein data frame
df_pproteomics_nodes <- rbind(df_pproteomics_nodes_sites, df_pproteomics_nodes)

#export node table to SCPC#46_03_Jim_Anni_pancreatic_cancer_phosphoproteomics_G12D_WT_vs_KO sig node table with extra site nodes.tsv
write.table(df_pproteomics_nodes)


### Figure S4 ###

## Fig S4A ##

#generate mouse input for Reactome

#import files
df_total_proteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_total_protein.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_total_proteomics) <- df_total_proteomics[1,]
df_total_proteomics <- df_total_proteomics[-1,]

#selecting and reordering columns
df_total_proteomics <- df_total_proteomics[, c(3:4, 21:32, 59:78)]
df_total_proteomics <- df_total_proteomics[, c(1:4, 6:7, 9:10, 12:13, 15:34)]

#covert numbers to numeric
df_total_proteomics[, c(3:30)] <- sapply(df_total_proteomics[, c(3:30)], as.numeric)

#calculate fold change
df_total_proteomics$`(G12D_WT) / (G_WT)` <- rowMedians(as.matrix(df_total_proteomics[, grepl("WT_G12D_", colnames(df_total_proteomics))]))/rowMedians(as.matrix(df_total_proteomics[, grepl("WT_G_", colnames(df_total_proteomics))]))

#select needed columns
df_total_proteomics <- df_total_proteomics[, c("Gene Symbol", "(G12D_WT) / (G_WT)", "Adj. P-Value: (G12D_WT) / (G_WT)")] %>% na.omit()

#get log2 of fold change
df_total_proteomics$log2FoldChange <- log2(df_total_proteomics[, 2])

#rename columns
colnames(df_total_proteomics) <- c("Gene", "FoldChange", "padj", "log2FoldChange")

#create gene list for Adj. P-Value < 0.05 used for Reactome input
df_total_proteomics <- subset(df_total_proteomics, padj < 0.05)


#generate human input for Reactome

#import files
df_normal <- read.csv(".\\Data\\proteomics_gene_level_MD_abundance_normal.cct.txt", check.names = FALSE, header = TRUE, sep = "\t")
df_tumor <- read.csv(".\\Data\\proteomics_gene_level_MD_abundance_tumor.cct.txt", check.names = FALSE, header = TRUE, sep = "\t")

#extract data columns
dat <- cbind(df_tumor[, 2:141], df_normal[, 2:76])

#set genes as row names
rownames(dat) <- df_tumor[, 1]

#define comparison groups
groups <- factor(rep(c("tumor", "normal"), c(140, 75)))

#create design matrix
design <- model.matrix(~ 0 + groups)
colnames(design) <- gsub("groups", "", colnames(design))

#perform differential expression analysis
contrast <-  limma::makeContrasts(contrasts = c("tumor-normal"), levels = design)
fit1 <- limma::lmFit(dat, design)
fit2 <- limma::contrasts.fit(fit1, contrasts = contrast)
fit3 <- limma::eBayes(fit2)

#get DE results
limma.results <- limma::topTable(fit3, number = nrow(dat))

#add genes as column
limma.results$Index <- rownames(limma.results)

#rename gene column in original data
colnames(df_tumor)[1] <- "Index"

#merge DE results with data
df_results <- merge(limma.results, cbind(df_tumor[, 1:141], df_normal[, 2:76]), by = "Index")

#add log2FC column
df_results <- cbind(log2(rowMeans(df_results[, 8:143], na.rm = T)/rowMeans(df_results[, 144:222], na.rm = T)), df_results)

#create gene list for Adj. P-Value < 1e-19 used for Reactome input
df_results <- subset(df_results, adj.P.Val < 1e-19)


#visualize Reactome results

#import Reactome result files
df_reactome_mouse <- read.csv(".\\Data\\reactome_result_mouse_protein 0.05.csv", check.names = FALSE, header = TRUE, sep = ",")
df_reactome_human <- read.csv(".\\Data\\reactome_result_human_protein 1e-19.csv", check.names = FALSE, header = TRUE, sep = ",")

#calculate gene ratio
df_reactome_mouse$`Gene Ratio` <- df_reactome_mouse$`#Entities found`/df_reactome_mouse$`#Entities total`
df_reactome_human$`Gene Ratio` <- df_reactome_human$`#Entities found`/df_reactome_human$`#Entities total`

#get top 20 mouse pathways
df_reactome_mouse_top <- head(df_reactome_mouse, 20)[, c("Pathway name", "Gene Ratio", "Entities FDR")]

#label mouse data columns
colnames(df_reactome_mouse_top)[2:3] <- paste0(colnames(df_reactome_mouse_top)[2:3], "_mouse")

#combine mouse and human Reactome result data frames
df_reactome_mouse_top <- merge(df_reactome_mouse_top,
                               df_reactome_human[df_reactome_human$`Pathway name` %in% df_reactome_mouse_top$`Pathway name`, c("Pathway name", "Gene Ratio", "Entities FDR")],
                               by = "Pathway name",
                               all = TRUE)

#label human data columns
colnames(df_reactome_mouse_top)[4:5] <- paste0(colnames(df_reactome_mouse_top)[4:5], "_human")

#arrange by min padj and max gene ratio between mouse and human
df_reactome_mouse_top <- arrange(df_reactome_mouse_top, desc(do.call(pmax, c(df_reactome_mouse_top[, c(2, 4)], list(na.rm = TRUE)))))
df_reactome_mouse_top <- arrange(df_reactome_mouse_top, desc(do.call(pmin, c(df_reactome_mouse_top[, c(3, 5)], list(na.rm = TRUE)))))

#set pathway name factor level order
df_reactome_mouse_top$`Pathway name` <- factor(df_reactome_mouse_top$`Pathway name`, levels = df_reactome_mouse_top$`Pathway name`)

#set x-axis limits
left_limit <- -0.5
right_limit <- 15

#get number of rows
num_row <- nrow(df_reactome_mouse_top)

#create plot
ggplot2::ggplot(df_reactome_mouse_top, ggplot2::aes(y = `Pathway name`)) +
  
  ggplot2::annotate(geom = "text", x = (left_limit + right_limit)/2, y = (num_row)*1.045 + 0.5, label = c("Top 20 Total Protein\nPathways in Mouse"), hjust = 0.5) +
  
  ggplot2::annotate(geom = "segment", linetype = "dashed", x = -log10(0.05), xend = -log10(0.05), y = 0.5, yend = num_row + 0.5, size = 0.5) +
  
  ggplot2::geom_point(ggplot2::aes(x = -log10(`Entities FDR_mouse`), 
                                   colour = "Mouse", 
                                   fill = "Mouse", 
                                   size = `Gene Ratio_mouse`)) +
  
  ggplot2::geom_point(ggplot2::aes(x = -log10(`Entities FDR_human`), 
                                   colour = "Human", 
                                   fill = "Human", 
                                   alpha = 0.5, 
                                   size = `Gene Ratio_human`)) +
  
  ggplot2::coord_cartesian(xlim = c(left_limit, right_limit), ylim = c(0.5, num_row + 0.5), expand = FALSE, clip = "off") +
  
  ggplot2::scale_y_discrete(labels = sapply(as.character(df_reactome_mouse_top$`Pathway name`), function(y) {ifelse(nchar(y) > 45, stringr::str_wrap(y, width = nchar(y)/1.7), y)}), drop = FALSE) +
  
  ggplot2::scale_alpha(guide = "none") +
  
  ggplot2::labs(y = NULL, 
                x = bquote(-log[10]*"(FDR)"), 
                colour = "Species", 
                shape = "Species", 
                fill = "Species", 
                size = "Gene Ratio") +
  
  ggplot2::guides(
    fill = ggplot2::guide_legend(reverse = T, 
                                 override.aes = list(size = 5), 
                                 order = 1), 
    colour = ggplot2::guide_legend(reverse = T, 
                                   order = 1), 
    shape = ggplot2::guide_legend(reverse = T, 
                                  order = 1), 
    size = ggplot2::guide_legend(reverse = T, 
                                 order = 2)) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 0.5, unit = "cm")),
                 plot.margin = ggplot2::margin(t = 2, r = 1, b = 1, l = 1, "cm"),
                 axis.text.y = ggplot2::element_text(size = 11, colour = "black"),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.grid = ggplot2::element_blank(),
                 legend.title.align = 0.5,
                 legend.direction = "vertical",
                 legend.box.just = "center",
                 legend.key.width = grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0.5, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = (nrow(df_reactome_mouse_top) + 1)/5)


## Fig S4B ##

#generate mouse input for Reactome

#import file
df_pproteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_phosphoproteomics.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_pproteomics) <- df_pproteomics[1,]
df_pproteomics <- df_pproteomics[-1,]

#selecting and reordering columns
df_pproteomics <- df_pproteomics[, c(1:2, 16, 32, 57:75)]

#covert numbers to numeric
df_pproteomics[, c(4:23)] <- sapply(df_pproteomics[, c(4:23)], as.numeric)

#calculate fold change
df_pproteomics$`(G12D_WT) / (G_WT)` <- rowMedians(as.matrix(df_pproteomics[, grepl("Sample, G12D_WT", colnames(df_pproteomics))]))/rowMedians(as.matrix(df_pproteomics[, grepl("Sample, G_WT", colnames(df_pproteomics))]))

#select needed columns
df_pproteomics <- df_pproteomics[, c("Gene", "(G12D_WT) / (G_WT)", "Adj. P-Value: (G12D_WT) / (G_WT)")] %>% na.omit()

#get log2 of fold change
df_pproteomics$log2FoldChange <- log2(df_pproteomics[, 2])

#rename columns
colnames(df_pproteomics) <- c("Gene", "FoldChange", "padj", "log2FoldChange")

#create gene list for Adj. P-Value < 0.05 used for Reactome input
df_pproteomics <- subset(df_pproteomics, padj < 0.05)


#generate human input for Reactome

#import files
df_normal <- read.csv(".\\Data\\phosphoproteomics_site_level_MD_abundance_normal.cct.txt", check.names = FALSE, header = TRUE, sep = "\t")
df_tumor <- read.csv(".\\Data\\phosphoproteomics_site_level_MD_abundance_tumor.cct.txt", check.names = FALSE, header = TRUE, sep = "\t")

#extract data columns
dat <- cbind(df_tumor[, 4:143], df_normal[, 4:78])

#set genes as row names
rownames(dat) <- df_tumor[,1]

#define comparison groups
groups <- factor(rep(c("tumor", "normal"), c(140, 75)))

#create design matrix
design <- model.matrix(~ 0 + groups)
colnames(design) <- gsub("groups", "", colnames(design))

#perform differential expression analysis
contrast <- limma::makeContrasts(contrasts = c("tumor-normal"), levels = design)
fit1 <- limma::lmFit(dat, design)
fit2 <- limma::contrasts.fit(fit1, contrasts = contrast)
fit3 <- limma::eBayes(fit2)

#get DE results
limma.results <- limma::topTable(fit3, number = nrow(dat))

#add genes as column
limma.results$Index <- rownames(limma.results)

#rename gene column in original data
colnames(df_tumor)[1] <- "Index"

#merge DE results with data
df_results <- merge(limma.results, cbind(df_tumor[, 1:143], df_normal[, 4:78]), by = "Index")

#add log2FC column
df_results <- cbind(log2(rowMeans(df_results[, 10:145], na.rm = T)/rowMeans(df_results[, 146:224], na.rm = T)), df_results)

#extract phospho-site
df_results$Site <- paste(df_results$Gene, gsub("NP(.+)_", "", df_results$Index))

#create gene list for Adj. P-Value < 1e-12 used for Reactome input
df_results <- subset(df_results, adj.P.Val < 1e-12)


#visualize Reactome results

#import Reactome result files
df_reactome_mouse <- read.csv(".\\Data\\reactome_result_mouse_pproteomics 0.05.csv", check.names = FALSE, header = TRUE, sep = ",")
df_reactome_human <- read.csv(".\\Data\\reactome_result_human_pproteomics 1e-12.csv", check.names = FALSE, header = TRUE, sep = ",")

#calculate gene ratio
df_reactome_mouse$`Gene Ratio` <- df_reactome_mouse$`#Entities found`/df_reactome_mouse$`#Entities total`
df_reactome_human$`Gene Ratio` <- df_reactome_human$`#Entities found`/df_reactome_human$`#Entities total`

#get top 20 mouse pathways
df_reactome_mouse_top <- head(df_reactome_mouse, 20)[, c("Pathway name", "Gene Ratio", "Entities FDR")]

#label mouse data columns
colnames(df_reactome_mouse_top)[2:3] <- paste0(colnames(df_reactome_mouse_top)[2:3], "_mouse")

#combine mouse and human Reactome result data frames
df_reactome_mouse_top <- merge(df_reactome_mouse_top,
                               df_reactome_human[df_reactome_human$`Pathway name` %in% df_reactome_mouse_top$`Pathway name`, c("Pathway name", "Gene Ratio", "Entities FDR")],
                               by = "Pathway name",
                               all = TRUE)

#label human data columns
colnames(df_reactome_mouse_top)[4:5] <- paste0(colnames(df_reactome_mouse_top)[4:5], "_human")

#arrange by min padj and max gene ratio between mouse and human
df_reactome_mouse_top <- arrange(df_reactome_mouse_top, desc(do.call(pmax, c(df_reactome_mouse_top[, c(2, 4)], list(na.rm = TRUE)))))
df_reactome_mouse_top <- arrange(df_reactome_mouse_top, desc(do.call(pmin, c(df_reactome_mouse_top[, c(3, 5)], list(na.rm = TRUE)))))

#set pathway name factor level order
df_reactome_mouse_top$`Pathway name` <- factor(df_reactome_mouse_top$`Pathway name`, levels = df_reactome_mouse_top$`Pathway name`)

#set x-axis limits
left_limit <- -0.5
right_limit <- 7

#get number of rows
num_row <- nrow(df_reactome_mouse_top)

#create plot
ggplot2::ggplot(df_reactome_mouse_top, ggplot2::aes(y = `Pathway name`)) +
  
  ggplot2::annotate(geom = "text", x = (left_limit + right_limit)/2, y = (num_row)*1.045 + 0.5, label = c("Top 20 Phospho-proteomics\nPathways in Mouse"), hjust = 0.5) +
  
  ggplot2::annotate(geom = "segment", linetype = "dashed", x = -log10(0.05), xend = -log10(0.05), y = 0.5, yend = num_row + 0.5, size = 0.5) +
  
  ggplot2::geom_point(ggplot2::aes(x = -log10(`Entities FDR_mouse`), 
                                   colour = "Mouse", 
                                   fill = "Mouse", 
                                   size = `Gene Ratio_mouse`)) +
  
  ggplot2::geom_point(ggplot2::aes(x = -log10(`Entities FDR_human`), 
                                   colour = "Human", 
                                   fill = "Human", 
                                   alpha = 0.5, 
                                   size = `Gene Ratio_human`)) +
  
  ggplot2::coord_cartesian(xlim = c(left_limit, right_limit), ylim = c(0.5, num_row + 0.5), expand = FALSE, clip = "off") +
  
  ggplot2::scale_y_discrete(labels = sapply(as.character(df_reactome_mouse_top$`Pathway name`), function(y) {ifelse(nchar(y) > 45, stringr::str_wrap(y, width = nchar(y)/1.7), y)}), drop = FALSE) +
  
  ggplot2::scale_alpha(guide = "none") +
  
  ggplot2::labs(y = NULL, 
                x = bquote(-log[10]*"(FDR)"), 
                colour = "Species", 
                shape = "Species", 
                fill = "Species", 
                size = "Gene Ratio") +
  
  ggplot2::guides(
    fill = ggplot2::guide_legend(reverse = T, 
                                 override.aes = list(size = 5), 
                                 order = 1), 
    colour = ggplot2::guide_legend(reverse = T, 
                                   order = 1), 
    shape = ggplot2::guide_legend(reverse = T, 
                                  order = 1), 
    size = ggplot2::guide_legend(reverse = T, 
                                 order = 2)) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 0.5, unit = "cm")),
                 plot.margin = ggplot2::margin(t = 2, r = 1, b = 1, l = 1, "cm"),
                 axis.text.y = ggplot2::element_text(size = 11, colour = "black"),
                 panel.grid.minor = ggplot2::element_blank(),
                 panel.grid = ggplot2::element_blank(),
                 legend.title.align = 0.5,
                 legend.direction = "vertical",
                 legend.box.just = "center",
                 legend.key.width = grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0.5, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = (nrow(df_reactome_mouse_top) + 1)/5)
