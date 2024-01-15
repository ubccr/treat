library(dplyr)
library(tidyr)
library(ggplot2)
library(EnhancedVolcano)

# Figure 6B

# Load DRBD18-RIP data
setwd("C:/Users/tylec/Read lab/DRBD18/RIPseq")
RIPdata <- read.csv("C:/Users/tylec/Read lab/DRBD18/RIPseq/diffexpr-results-shrunk-20230503.csv", sep = ",", header=TRUE)

# Deeptools Patterns
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/10_new_DeepToolsAnalysis")
deeptools <- read.csv("DeeptoolsPatternsGeneList.csv", header = T, sep = ",")
deeptools <- deeptools %>% select(X.chrom, start, end, score, strand, thickStart, thickEnd, deepTools_group, name, geneName, Description)

dtA <- deeptools %>% filter(deepTools_group=="pattern A") %>% pull(name)
dtB <- deeptools %>% filter(deepTools_group=="pattern B") %>% pull(name)
dtC <- deeptools %>% filter(deepTools_group=="pattern C") %>% pull(name)
dtD <- deeptools %>% filter(deepTools_group=="pattern D") %>% pull(name)
dtE <- deeptools %>% filter(deepTools_group=="pattern E") %>% pull(name)

# Set directory for saving files
setwd("C:/Users/tylec/Read lab/DRBD18/RIPseq/VolcanoPlots_Dec16_2023")


# Pattern A ---------------------------------------------------------------

x <- as.character(RIPdata$Gene) %in% dtA
RIPdata <- rbind(RIPdata[!x,], RIPdata[x,])
keyvals <- ifelse(RIPdata$Gene %in% dtA, 'red','gray46')
names(keyvals)[keyvals == 'red'] <- 'Pattern A'
names(keyvals)[keyvals == 'gray46'] <- 'Other'

# Full size plot
EnhancedVolcano(RIPdata,
                lab = NA,
                x = 'log2FoldChange',labSize = 3,
                y = 'padj',pCutoffCol='padj',pCutoff = .05,
                FCcutoff = 0.5849,
                #xlim = c(-4,6),
                #ylim = c(-1,15),
                xlab = "Enrichment",
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                colCustom = keyvals,
                axisLabSize = 24,
                legendLabSize = 30,
                cutoffLineType = "blank",
                vline = 0.5849) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1)) + theme(legend.position = c(0.85,0.94))
ggsave(filename="Fig6B_Deeptools_A_RIP_plot.pdf", width=9, height=6)

# Zoomed in plot
EnhancedVolcano(RIPdata,
                lab = NA,
                x = 'log2FoldChange',labSize = 3,
                y = 'padj',pCutoffCol='padj',pCutoff = .05,
                FCcutoff = 0.5849,
                #xlim = c(-4,6),
                ylim = c(0,150),
                xlab = "Enrichment",
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                colCustom = keyvals,
                axisLabSize = 24,
                legendLabSize = 24,
                cutoffLineType = "blank",
                vline = 0.5849) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1)) + theme(legend.position = c(0.88,0.98))
ggsave(filename="Deeptools_A_RIP_plot_zoom.png", width=9, height=6)


# Pattern B ---------------------------------------------------------------

x <- as.character(RIPdata$Gene) %in% dtB
RIPdata <- rbind(RIPdata[!x,], RIPdata[x,])
keyvals <- ifelse(RIPdata$Gene %in% dtB, 'red','gray46')
names(keyvals)[keyvals == 'red'] <- 'Pattern B'
names(keyvals)[keyvals == 'gray46'] <- 'Other'

# Full size plot
EnhancedVolcano(RIPdata,
                lab = NA,
                x = 'log2FoldChange',labSize = 3,
                y = 'padj',pCutoffCol='padj',pCutoff = .05,
                FCcutoff = 0.5849,
                #xlim = c(-4,6),
                #ylim = c(-1,15),
                xlab = "Enrichment",
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                colCustom = keyvals,
                axisLabSize = 24,
                legendLabSize = 30,
                cutoffLineType = "blank",
                vline = 0.5849) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1)) + theme(legend.position = c(0.85,0.94))
ggsave(filename="Fig6B_Deeptools_B_RIP_plot.pdf", width=9, height=6)

# Zoomed in plot
EnhancedVolcano(RIPdata,
                lab = NA,
                x = 'log2FoldChange',labSize = 3,
                y = 'padj',pCutoffCol='padj',pCutoff = .05,
                FCcutoff = 0.5849,
                #xlim = c(-4,6),
                ylim = c(0,150),
                xlab = "Enrichment",
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                colCustom = keyvals,
                axisLabSize = 24,
                legendLabSize = 24,
                cutoffLineType = "blank",
                vline = 0.5849) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1)) + theme(legend.position = c(0.88,0.98))
ggsave(filename="Deeptools_B_RIP_plot_zoom.png", width=9, height=6)


# Pattern C ---------------------------------------------------------------

x <- as.character(RIPdata$Gene) %in% dtC
RIPdata <- rbind(RIPdata[!x,], RIPdata[x,])
keyvals <- ifelse(RIPdata$Gene %in% dtC, 'red','gray46')
names(keyvals)[keyvals == 'red'] <- 'Pattern C'
names(keyvals)[keyvals == 'gray46'] <- 'Other'

# Full size plot
EnhancedVolcano(RIPdata,
                lab = NA,
                x = 'log2FoldChange',labSize = 3,
                y = 'padj',pCutoffCol='padj',pCutoff = .05,
                FCcutoff = 0.5849,
                #xlim = c(-4,6),
                #ylim = c(-1,15),
                xlab = "Enrichment",
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                colCustom = keyvals,
                axisLabSize = 24,
                legendLabSize = 30,
                cutoffLineType = "blank",
                vline = 0.5849) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1)) + theme(legend.position = c(0.85,0.94))
ggsave(filename="Fig6B_Deeptools_C_RIP_plot.pdf", width=9, height=6)

# Zoomed in plot
EnhancedVolcano(RIPdata,
                lab = NA,
                x = 'log2FoldChange',labSize = 3,
                y = 'padj',pCutoffCol='padj',pCutoff = .05,
                FCcutoff = 0.5849,
                #xlim = c(-4,6),
                ylim = c(0,150),
                xlab = "Enrichment",
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                colCustom = keyvals,
                axisLabSize = 24,
                legendLabSize = 24,
                cutoffLineType = "blank",
                vline = 0.5849) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1)) + theme(legend.position = c(0.88,0.98))
ggsave(filename="Deeptools_C_RIP_plot_zoom.png", width=9, height=6)


# Pattern D ---------------------------------------------------------------

x <- as.character(RIPdata$Gene) %in% dtD
RIPdata <- rbind(RIPdata[!x,], RIPdata[x,])
keyvals <- ifelse(RIPdata$Gene %in% dtD, 'red','gray46')
names(keyvals)[keyvals == 'red'] <- 'Pattern D'
names(keyvals)[keyvals == 'gray46'] <- 'Other'

# Full size plot
EnhancedVolcano(RIPdata,
                lab = NA,
                x = 'log2FoldChange',labSize = 3,
                y = 'padj',pCutoffCol='padj',pCutoff = .05,
                FCcutoff = 0.5849,
                #xlim = c(-4,6),
                #ylim = c(-1,15),
                xlab = "Enrichment",
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                colCustom = keyvals,
                axisLabSize = 24,
                legendLabSize = 30,
                cutoffLineType = "blank",
                vline = 0.5849) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1)) + theme(legend.position = c(0.85,0.94))
ggsave(filename="Fig6B_Deeptools_D_RIP_plot.pdf", width=9, height=6)

# Zoomed in plot
EnhancedVolcano(RIPdata,
                lab = NA,
                x = 'log2FoldChange',labSize = 3,
                y = 'padj',pCutoffCol='padj',pCutoff = .05,
                FCcutoff = 0.5849,
                #xlim = c(-4,6),
                ylim = c(0,150),
                xlab = "Enrichment",
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                colCustom = keyvals,
                axisLabSize = 24,
                legendLabSize = 24,
                cutoffLineType = "blank",
                vline = 0.5849) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1)) + theme(legend.position = c(0.88,0.98))
ggsave(filename="Deeptools_D_RIP_plot_zoom.png", width=9, height=6)


# Pattern E ---------------------------------------------------------------

x <- as.character(RIPdata$Gene) %in% dtE
RIPdata <- rbind(RIPdata[!x,], RIPdata[x,])
keyvals <- ifelse(RIPdata$Gene %in% dtE, 'red','gray46')
names(keyvals)[keyvals == 'red'] <- 'Pattern E'
names(keyvals)[keyvals == 'gray46'] <- 'Other'

# Full size plot
EnhancedVolcano(RIPdata,
                lab = NA,
                x = 'log2FoldChange',labSize = 3,
                y = 'padj',pCutoffCol='padj',pCutoff = .05,
                FCcutoff = 0.5849,
                #xlim = c(-4,6),
                #ylim = c(-1,15),
                xlab = "Enrichment",
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                colCustom = keyvals,
                axisLabSize = 24,
                legendLabSize = 30,
                cutoffLineType = "blank",
                vline = 0.5849) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1)) + theme(legend.position = c(0.85,0.94))
ggsave(filename="Fig6B_Deeptools_E_RIP_plot.pdf", width=9, height=6)

# Zoomed in plot
EnhancedVolcano(RIPdata,
                lab = NA,
                x = 'log2FoldChange',labSize = 3,
                y = 'padj',pCutoffCol='padj',pCutoff = .05,
                FCcutoff = 0.5849,
                #xlim = c(-4,6),
                ylim = c(0,150),
                xlab = "Enrichment",
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                colCustom = keyvals,
                axisLabSize = 24,
                legendLabSize = 24,
                cutoffLineType = "blank",
                vline = 0.5849) +
  ggplot2::scale_x_continuous(
    breaks=seq(-6,6, 1)) + theme(legend.position = c(0.88,0.98))
ggsave(filename="Deeptools_E_RIP_plot_zoom.png", width=9, height=6)
