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


### Figure 4 ###

## Fig 4A ##

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
                 legend.key.width = grid::grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = 1/1)


## Fig 4B ##

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
df_total_proteomics <- df_total_proteomics[df_total_proteomics$`Gene Symbol` %in% c("Rnase1", "Cuzd1", "Cpa2", "Klk1", "Pnliprp2", "Cela3b", "Prss2", "Cela1", "Sycn", "Amy1", "Ctrb1", "Cel", "Amy2", "Cpa1", "Pnlip", "Pnliprp1", "Cela2a", "Pla2g1b", "Clps", "Ctrc"),]

#select abundance columns
df_protein_levels <- df_total_proteomics[, c(3:22)] %>% na.omit()

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


## Fig 4C ##

#clear environment
rm(list = ls())

#import file
df_total_proteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_total_protein.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_total_proteomics) <- df_total_proteomics[1,]
df_total_proteomics <- df_total_proteomics[-1,]

#select columns
df_total_proteomics <- df_total_proteomics[, c("Gene Symbol", "(G12D_KO) / (G12D_WT)", "Adj. P-Value: (G12D_KO) / (G12D_WT)")]

#convert numbers to numeric
df_total_proteomics[, c(2:3)] <- as.data.frame(sapply(df_total_proteomics[, c(2:3)], as.numeric))

#get log2 of fold change
df_total_proteomics$`log2 (G12D_KO) / (G12D_WT)` <- log2(df_total_proteomics$`(G12D_KO) / (G12D_WT)`)

#get -log10 of adj. p-value
df_total_proteomics$`-log10 Adj. P-Value: (G12D_KO) / (G12D_WT)` <- -log10(df_total_proteomics$`Adj. P-Value: (G12D_KO) / (G12D_WT)`)

#subset for significance and select proteins
df_total_proteomics <- subset(df_total_proteomics, `Adj. P-Value: (G12D_KO) / (G12D_WT)` < 0.05 & abs(`log2 (G12D_KO) / (G12D_WT)`) > 2.5 | df_total_proteomics$`Gene Symbol` %in% c("Akt1", "Akt2"))

#add network anchors
df_total_proteomics_anchors <- df_total_proteomics

df_total_proteomics_anchors[, c(1:ncol(df_total_proteomics_anchors))] <- ""

df_total_proteomics_anchors <- df_total_proteomics_anchors[c(1:6),]

df_total_proteomics_anchors[, 1] <- c("Kras", "Igf1r", "Insr", "mTor", "Raf1", "Mapk1")

df_total_proteomics <- rbind(df_total_proteomics, df_total_proteomics_anchors)

#rename gene column name
colnames(df_total_proteomics)[1] <- "#node1"

#export node table to SCPC#46_03_Jim_Anni_pancreatic_cancer_total_protein_G12D_WT_vs_KO sig logFC and padj.tsv
write.table(df_total_proteomics)

#import STRING protein mapping and edge table output files
df_string_network <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_protein_G12D_WT_vs_KO sig string_interactions.tsv", check.names = FALSE, sep = "\t")
df_string_mapping <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_protein_G12D_WT_vs_KO sig string_mapping.tsv", check.names = FALSE, sep = "\t")

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
df_string_input <- data.frame(`#node1` = df_total_proteomics$`#node1`, check.names = FALSE)

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


### Figure 5 ###

## Fig 5A ##

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
df_pproteomics$`(G12D_WT) / (G_WT)` <- rowMedians(as.matrix(df_pproteomics[, grepl("Sample, G12D_WT", colnames(df_pproteomics))]))/rowMedians(as.matrix(df_pproteomics[, grepl("Sample, G_WT", colnames(df_pproteomics))]))

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


## Fig 5B ##

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
  
  ggplot2::labs(y = expression(log[2]*"(Fold Change) Phosphosite"), x = expression(log[2]*"(Fold Change) Total Protein"), fill = expression(atop(-log[10]*"(Adj. "*italic(p)*"-value)", "Phosphosite"))) +
  
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


## Fig 5C ##

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
                          lab = c("***", "****", "**", "**", "*", "*"),
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
                          lab = c("", "****", "", "", "****", "*"),
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
                          lab = c("", "*", "", "", "", ""),
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
  
  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, 
                                      ymax = mean + se, 
                                      x = x_position,
                                      colour = Condition), 
                         width = 0.2, 
                         data = data_matrix_se) +
  
  ggplot2::geom_bar(ggplot2::aes(y = Abundance, x = x_position),
                    stat = "summary", 
                    fun = "mean", 
                    width = 1, 
                    colour = NA,
                    fill = rep(c("#0e6e80",
                                 "#ce1e68",
                                 "#3b54a3",
                                 "#f47d23",
                                 "#bbbd32"), 6)) +
  
  ggplot2::geom_jitter(ggplot2::aes(y = Abundance, x = x_position),
                       size = 2,
                       colour = "black",
                       fill = "white",
                       shape = 21,
                       stroke = 0.75,
                       width = 0.15) + 
  
  ggplot2::geom_text(data = anno_text_1, 
                     ggplot2::aes(x = (x1 + x2)/2,  
                                  y = ystar, 
                                  label = lab),
                     colour = "black",
                     hjust = 0.5) +
  
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
                     hjust = 0.5) +
  
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
                     hjust = 0.5) +
  
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
  
  ggplot2::labs(y = "Phospho-peptide Intensity (AU)") +
  
  ggplot2::guides(colour = guide_legend(override.aes = list(shape = 15, 
                                                            colour = c("#0e6e80",
                                                                       "#ce1e68",
                                                                       "#3b54a3",
                                                                       "#f47d23",
                                                                       "#bbbd32"),
                                                            fill = c("#0e6e80",
                                                                     "#ce1e68",
                                                                     "#3b54a3",
                                                                     "#f47d23",
                                                                     "#bbbd32"),
                                                            size = 8)),
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


## Fig 5D ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import phosphoproteomics protein file
df_pproteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_phosphoproteomics_all.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_pproteomics) <- df_pproteomics[1,]
df_pproteomics <- df_pproteomics[-1,]

#select needed columns
df_pproteomics <- df_pproteomics[, c(4:5, 12, 14:15, 29:47)]

#remove rows with no protein
df_pproteomics <- subset(df_pproteomics, !df_pproteomics$`Master Protein Accessions` == "")

#replace empty abundances with 0
df_pproteomics[df_pproteomics == ""] <- 0

#covert numbers to numeric
df_pproteomics[, c(6:24)] <- as.data.frame(sapply(df_pproteomics[, c(6:24)], as.numeric))

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

#capitalize gene names as proteins
df_pproteomics$Gene <- toupper(df_pproteomics$Gene)

#create Uniprot with phospho-site column
df_pproteomics$Sum_Index_Unique <- ifelse(grepl("\\d", df_pproteomics$Last_Phospho_Info), paste(df_pproteomics$`Master Protein Accessions`, df_pproteomics$Phospho_Numeric), paste(df_pproteomics$`Positions in Master Proteins`, df_pproteomics$Last_Phospho_Info))

#create gene name with phospho-site column
df_pproteomics$Sum_Index_Unique_Gene <- ifelse(grepl("\\d", df_pproteomics$Last_Phospho_Info), paste(df_pproteomics$Gene, df_pproteomics$Phospho_Numeric), paste(df_pproteomics$Gene, df_pproteomics$`Positions in Master Proteins`, df_pproteomics$Last_Phospho_Info))

#only include phospho-sites with defined numerical positions
df_pproteomics <- df_pproteomics[grepl("\\d", df_pproteomics$Last_Phospho_Info),]

#select columns needed
df_pproteomics <- cbind(df_pproteomics[, c(6:24)], df_pproteomics[, c("Master Protein Accessions", 
                                                                      "Gene",
                                                                      "Sum_Index_Unique",
                                                                      "Sum_Index_Unique_Gene")])

#remove duplicate rows
df_pproteomics <- df_pproteomics[!duplicated(df_pproteomics),]

#aggregate phospho-site by sum
df_pproteomics <- aggregate(df_pproteomics[, c(1:19)], by = list(df_pproteomics$`Master Protein Accessions`,
                                                                 df_pproteomics$Gene,
                                                                 df_pproteomics$Sum_Index_Unique,
                                                                 df_pproteomics$Sum_Index_Unique_Gene), FUN = sum)

#rename ID column names
colnames(df_pproteomics)[1:4] <- c("Master Protein Accessions", "Gene", "Sum_Index_Unique", "Sum_Index_Unique_Gene")

#annotation phospho-site columns
colnames(df_pproteomics)[5:23] <- paste(colnames(df_pproteomics)[5:23], rep("pproteomics", ncol(df_pproteomics[5:23])), sep = "_")

#import total protein file
df_total_proteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_total_protein.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_total_proteomics) <- df_total_proteomics[1,]
df_total_proteomics <- df_total_proteomics[-1,]

#select needed columns
df_total_proteomics <- df_total_proteomics[, c(3, 33:39, 41:52)]

#filter for only rows with all total protein sample abundances present
df_total_proteomics <- subset(df_total_proteomics, !rowSums(df_total_proteomics[, c(2:20)] == "") > 0)

#rename total protein column
colnames(df_total_proteomics)[1] <- c("Master Protein Accessions")

#annotation phospho-site columns
colnames(df_total_proteomics)[2:20] <- paste(colnames(df_total_proteomics)[2:20], rep("total", ncol(df_total_proteomics[2:20])), sep = "_")

#merge phospho-site and total protein by protein
df_pproteomics <- merge(df_pproteomics, df_total_proteomics, by = "Master Protein Accessions", all = TRUE)

#filter out non-merged rows
df_pproteomics <- subset(df_pproteomics, !rowSums(is.na(df_pproteomics[, c(3, 24)])) > 0)

#select pproteomics columns
df_pproteomics_for_norm <- df_pproteomics[, grepl("_pproteomics", colnames(df_pproteomics))]

#select total proteomics columns
df_total_proteomics_for_norm <- df_pproteomics[, grepl("_total", colnames(df_pproteomics))]

#normalize pproteomics to total
df_pproteomics_norm_to_total <- data.frame(sapply(df_pproteomics_for_norm, function(x) as.numeric(as.character(x))))/data.frame(sapply(df_total_proteomics_for_norm, function(x) as.numeric(as.character(x))))

#reattach ID columns
df_pproteomics_norm_to_total <- cbind(df_pproteomics[, c(1:4)], df_pproteomics_norm_to_total)

#create unique ID per row
df_pproteomics_norm_to_total$uniqueID <- paste(1:nrow(df_pproteomics_norm_to_total), df_pproteomics_norm_to_total$Gene)

#remove pproteomics annotation
colnames(df_pproteomics_norm_to_total) <- gsub("_pproteomics", "", colnames(df_pproteomics_norm_to_total))

#select G12D KO and G12D WT columns
dat <- df_pproteomics_norm_to_total[, c(15:23)]

#set unique ID column as row name
rownames(dat) <- df_pproteomics_norm_to_total[, "uniqueID"]

#DE analysis
groups <- factor(rep(c("G12D_KO", "G12D_WT"), c(6, 3)))
design = model.matrix(~ 0 + groups)
colnames(design) = gsub("groups", "", colnames(design))
contrast =  limma::makeContrasts(contrasts = c("G12D_KO-G12D_WT"), levels = design)
fit1 <- limma::lmFit(dat, design)
fit2 <- limma::contrasts.fit(fit1, contrasts = contrast)
fit3 <- limma::eBayes(fit2)
limma.results <- limma::topTable(fit3, number = nrow(dat))

#add unique ID as column
limma.results$uniqueID <- rownames(limma.results)

#merge DE stats with abundance data frame
df_pproteomics_norm_to_total <- merge(limma.results, df_pproteomics_norm_to_total, by = "uniqueID")

#filter for top 21 phospho-sites or Adj P-value < 0.75
df_pproteomics_norm_to_total <- head(df_pproteomics_norm_to_total[with(df_pproteomics_norm_to_total, order(adj.P.Val)),], 21)

#select and order abundance columns
df_protein_levels <- cbind(df_pproteomics_norm_to_total[, c(12:30)][, grepl("WT_G_", colnames(df_pproteomics_norm_to_total)[c(12:30)])],
                           df_pproteomics_norm_to_total[, c(12:30)][, grepl("KO_G_", colnames(df_pproteomics_norm_to_total)[c(12:30)])],
                           df_pproteomics_norm_to_total[, c(12:30)][, grepl("WT_G12D_", colnames(df_pproteomics_norm_to_total)[c(12:30)])],
                           df_pproteomics_norm_to_total[, c(12:30)][, grepl("Het_G12D_", colnames(df_pproteomics_norm_to_total)[c(12:30)])],
                           df_pproteomics_norm_to_total[, c(12:30)][, grepl("KO_G12D_", colnames(df_pproteomics_norm_to_total)[c(12:30)])])

#convert to Z-score matrix
df_protein_levels$mean_all <- rowMeans(df_protein_levels[, c(1:19)])

df_protein_levels$sd_all <- apply(df_protein_levels[, c(1:19)], 1, sd) 

Z_score <- function(x) with(df_protein_levels, (x - mean_all)/sd_all)

df_protein_Z_score <- cbind(apply(df_protein_levels[, c(1:(ncol(df_protein_levels) - 2))], 2, Z_score))

rownames(df_protein_Z_score) <- df_pproteomics_norm_to_total$Sum_Index_Unique_Gene

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
                                                ind = unique(c(mat[, c(3, 5, 8, 13)]))
                                                grid.rect(x = x[ind] + grid::unit(0.5/ncol(matrix_Z_score_total), "npc"), 
                                                          y = y[ind], 
                                                          width = grid::unit(0.03, "inches"), 
                                                          height = grid::unit(1/nrow(matrix_Z_score_total), "npc"),
                                                          gp = grid::gpar(col = "white")
                                                )
                                              },
                                              col = circlize::colorRamp2(c(-ceiling(max(abs(matrix_Z_score_total), na.rm = T)), 0, ceiling(max(abs(matrix_Z_score_total), na.rm = T))), c("blue", "white", c("red"))),
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
                      (3 + 2 + 3 + 5/2)/ncol(matrix_Z_score_total),
                      (3 + 2 + 3 + 5 + 6/2)/ncol(matrix_Z_score_total)),
                y = 0.25,
                just = c("center", "center"),
                gp = grid::gpar(fontsize = 14,
                                col = c("#0e6e80",
                                        "#ce1e68",
                                        "#3b54a3",
                                        "#f47d23",
                                        "#bbbd32")))

grid::grid.rect(x = 12/19,
                y = 0.6,
                height = grid::unit(0.01, "inches"),
                width = 14/19*(loc2$x - loc1$x),
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


## Fig 5E ##

#clear environment
rm(list = ls())

#import file
df_pproteomics <- read.csv(".\\Data\\SCPC#46_03_Jim_Anni_pancreatic_cancer_total_phospho_combined_LFQ_phosphoproteomics.csv", check.names = FALSE, header = FALSE)[-c(1:2),]

#fix column names
colnames(df_pproteomics) <- df_pproteomics[1,]
df_pproteomics <- df_pproteomics[-1,]

#select columns
df_pproteomics <- df_pproteomics[, c("Gene", 
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

#rename gene column name
colnames(df_pproteomics)[1] <- "#node1"

#order by Adj. P-value and abs(log2FC)
df_pproteomics <- df_pproteomics[with(df_pproteomics, order(df_pproteomics$`Adj. P-Value: (G12D_KO) / (G12D_WT)`, -abs(df_pproteomics$`log2 (G12D_KO) / (G12D_WT)`))),]

#remove duplicates
df_pproteomics <- df_pproteomics[!duplicated(df_pproteomics$`#node1`),]

#export node table to SCPC#46_03_Jim_Anni_pancreatic_cancer_phosphoproteomics_G12D_WT_vs_KO sig logFC and padj removed duplicates.tsv
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

#export edge table to SCPC#46_03_Jim_Anni_pancreatic_cancer_phosphoproteomics_G12D_WT_vs_KO sig string_interactions remapped.tsv
write.table(df_string_network)


### Extended Data 9 ###

## Extended Data 9A ##

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
                 legend.key.width = grid::grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = 1/1)


## Extended Data 9B ##

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
                 legend.key.width = grid::grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = 1/1)


## Extended Data 9C ##

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



