library(tidyr)
library(dplyr)

RIPdata <- read.csv("C:/Users/tylec/Read lab/DRBD18/RIPseq/diffexpr-results-shrunk-20230503.csv", sep = ",", header=TRUE)
RIPbound <- RIPdata %>% filter(padj < 0.05 & log2FoldChange >= 0.584963 & AddFltr==T) %>% pull(Gene)

overlaptest <- function (SampleA, SampleB) {
  q <- length(intersect(SampleA, SampleB)) - 1
  m <- length(SampleA)
  n <- 11626 - m
  k <- length(SampleB)
  o <- length(intersect(SampleA, SampleB))
  print(noquote(paste("Overlap:", o)))
  print(noquote(paste("p-value:", phyper(q,m,n,k,lower.tail=FALSE))))
}


# Fig S6A -----------------------------------------------------------------
# Life cycle lists
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/11_GeneLists")

specific_slender_table <- read.table("Roditi_2018_slender.csv", header = T, sep = ",")
specific_slender <- specific_slender_table$Gene.ID

specific_stumpy_table <- read.table("Roditi_2018_stumpy.csv", header = T, sep = ",")
specific_stumpy <- specific_stumpy_table$ID 

specific_earlyPF_table <- read.table("Roditi_2018_earlyPF.csv", header = T, sep = ",")
specific_earlyPF <- specific_earlyPF_table$Gene.ID

specific_latePF_table <- read.table("Roditi_2018_latePF.csv", header = T, sep = ",")
specific_latePF <- specific_latePF_table$Gene.ID

specific_MFup_table <- read.table("TableS6_DESeq2_UPinMF.csv", header = T, sep = ",")
specific_MFup <- specific_MFup_table$Gene.ID

specific_MFdown_table <- read.table("TableS6_DESeq2_DOWNinMF.csv", header = T, sep = ",")
specific_MFdown <- specific_MFdown_table$Gene.ID


overlaptest(RIPbound, specific_slender)
# Overlap: 1
# p-value: 0.946016032796203

overlaptest(RIPbound, specific_stumpy)
# Overlap: 32
# p-value: 1.41766759485173e-17

overlaptest(RIPbound, specific_earlyPF)
# Overlap: 0
# p-value: 1

overlaptest(RIPbound, specific_latePF)
# Overlap: 0
# p-value: 1

overlaptest(RIPbound, specific_MFup)
# Overlap: 141
# p-value: 3.39484524841009e-50

overlaptest(RIPbound, specific_MFdown)
# Overlap: 13
# p-value: 0.999926696372707


# Fig S6B -----------------------------------------------------------------
# Export impaired
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/11_GeneLists")

# Down in cyto and no change in whole cell (export impaired in DRBD18 RNAi)
export_impaired <- read.table("AMAR_export_impaired.csv", header = T, sep = ",")
expimp <- export_impaired$transcript

# TE lists
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/11_GeneLists")

TE_up <- read.table("TE_up.csv", header = F, sep = ",") %>% pull(V1)
TE_down <- read.table("TE_down.csv", header = F, sep = ",") %>% pull(V1)

# RNAi bulk RNAseq lists
RNAseqdata <- read.csv("C:/Users/tylec/Read lab/DRBD18_Singlecell/1_Read_10X-Pilot/Plus_vs_Minus-RNASeq-Description.csv", header = T, sep = ",")

RNAsequp <- RNAseqdata %>% filter(padj < 0.05 & log2FoldChange >= 0.584963) %>% pull(X)
RNAseqdown <- RNAseqdata %>% filter(padj < 0.05 & log2FoldChange <= -0.584963) %>% pull(X)


overlaptest(RIPbound, expimp)
# Overlap: 2
# p-value: 0.953944992648034

overlaptest(RIPbound, TE_up)
# Overlap: 10
# p-value: 2.39747215010508e-05

overlaptest(RIPbound, TE_down)
# Overlap: 0
# p-value: 1

overlaptest(RIPbound, RNAsequp)
# Overlap: 218
# p-value: 2.40557612115259e-134

overlaptest(RIPbound, RNAseqdown)
# Overlap: 15
# p-value: 0.807825709721741


# Fig S6C -----------------------------------------------------------------
# Cluster Marker genes
c0markers <- read.table("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis/K5_highresClusters_MarkerGenes_RNA_Annotated_C0.txt", header = T, sep = "\t")
c0markersup <- c0markers %>% filter(avg_log2FC > 0)
c0up <- c0markersup$gene

c1markers <- read.table("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis/K5_highresClusters_MarkerGenes_RNA_Annotated_C1.txt", header = T, sep = "\t")
c1markersup <- c1markers %>% filter(avg_log2FC > 0)
c1up <- c1markersup$gene

c2markers <- read.table("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis/K5_highresClusters_MarkerGenes_RNA_Annotated_C2.txt", header = T, sep = "\t")
c2markersup <- c2markers %>% filter(avg_log2FC > 0)
c2up <- c2markersup$gene

c3markers <- read.table("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis/K5_highresClusters_MarkerGenes_RNA_Annotated_C3.txt", header = T, sep = "\t")
c3markersup <- c3markers %>% filter(avg_log2FC > 0)
c3up <- c3markersup$gene

c4markers <- read.table("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis/K5_highresClusters_MarkerGenes_RNA_Annotated_C4.txt", header = T, sep = "\t")
c4markersup <- c4markers %>% filter(avg_log2FC > 0)
c4up <- c4markersup$gene

c5markers <- read.table("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis/K5_highresClusters_MarkerGenes_RNA_Annotated_C5.txt", header = T, sep = "\t")
c5markersup <- c5markers %>% filter(avg_log2FC > 0)
c5up <- c5markersup$gene

c6markers <- read.table("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis/K5_highresClusters_MarkerGenes_RNA_Annotated_C6.txt", header = T, sep = "\t")
c6markersup <- c6markers %>% filter(avg_log2FC > 0)
c6up <- c6markersup$gene

c7markers <- read.table("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis/K5_highresClusters_MarkerGenes_RNA_Annotated_C7.txt", header = T, sep = "\t")
c7markersup <- c7markers %>% filter(avg_log2FC > 0)
c7up <- c7markersup$gene


overlaptest(RIPbound, c0up)
# Overlap: 117
# p-value: 1.01538458730063e-90

overlaptest(RIPbound, c1up)
# Overlap: 0
# p-value: 1

overlaptest(RIPbound, c2up)
# Overlap: 6
# p-value: 0.240148623544208

overlaptest(RIPbound, c3up)
# Overlap: 0
# p-value: 1

overlaptest(RIPbound, c4up)
# Overlap: 8
# p-value: 0.0302843004071452

overlaptest(RIPbound, c5up)
# Overlap: 7
# p-value: 0.999839372252715

overlaptest(RIPbound, c6up)
# Overlap: 5
# p-value: 0.94063023614791

overlaptest(RIPbound, c7up)
# Overlap: 57
# p-value: 5.72769850988503e-33
