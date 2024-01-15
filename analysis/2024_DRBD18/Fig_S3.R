library(Seurat)
library(dplyr)
library(tidyr)
library(EnhancedVolcano)

# Load data
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis")
mydata <- readRDS(file = "updates.5k.rds")
DefaultAssay(mydata) <- "RNA"
Idents(mydata) <- "integrated_snn_res.0.4"

# Load gene annotation file
gaf <- read.table("C:/Users/tylec/Read lab/TbruceiTREU927Annotations.csv", header = T, sep = ",")

# Find marker genes for cluster 0 vs cluster 7
c0v7 <- FindMarkers(mydata, ident.1="0", ident.2="7")
c0v7$gene <- rownames(c0v7)
c0v7annotated <- left_join(c0v7, gaf, by=c('gene' = 'geneID'))
write.table(c0v7annotated, "K5_C0vsC7_MarkersAnnotated.txt", sep="\t", row.names=FALSE)

# Load list of VSG genes
VSGs <- read.table("C:/Users/tylec/Read lab/DRBD18_Singlecell/11_GeneLists/VSG_Genes.csv", header = T, sep = ",")

GeneLabels <- VSGs$Gene.ID

x <- as.character(c0v7$gene) %in% GeneLabels
c0v7 <- rbind(c0v7[!x,], c0v7[x,])

keyvals <- ifelse(rownames(c0v7) %in% GeneLabels, 'red', 'gray46')

names(keyvals)[keyvals == 'gray46'] <- 'Other'
names(keyvals)[keyvals == 'red'] <- 'VSG'

EnhancedVolcano(c0v7, lab = c0v7$gene, 
                x = "avg_log2FC", y = "p_val_adj",
                title = 'Cluster 0 vs Cluster 7 Markers',
                subtitle = NULL,
                selectLab = GeneLabels,
                colCustom = keyvals,
                drawConnectors = TRUE,
                max.overlaps=30,
                pointSize = 3) +
  ggplot2::scale_x_continuous(
    breaks=seq(-3,3, 1),
    lim=c(-3.5,3.5))
ggsave(filename="C0vsC7_volcano.tiff", width=14, height=9, dpi=600, compression = "lzw")
