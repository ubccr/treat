# S5B ---------------------------------------------------------------------

library(tidyr)
library(dplyr)

# Marker genes
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



# Deeptools Patterns
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/10_new_DeepToolsAnalysis")
deeptools <- read.csv("DeeptoolsPatternsGeneList.csv", header = T, sep = ",")
deeptools <- deeptools %>% select(X.chrom, start, end, score, strand, thickStart, thickEnd, deepTools_group, name, geneName, Description)

length(unique(deeptools$name))
# 11439

dtA <- deeptools %>% filter(deepTools_group=="pattern A") %>% pull(name)
dtB <- deeptools %>% filter(deepTools_group=="pattern B") %>% pull(name)
dtC <- deeptools %>% filter(deepTools_group=="pattern C") %>% pull(name)
dtD <- deeptools %>% filter(deepTools_group=="pattern D") %>% pull(name)
dtE <- deeptools %>% filter(deepTools_group=="pattern E") %>% pull(name)



overlaptest <- function (SampleA, SampleB) {
  q <- length(intersect(SampleA, SampleB)) - 1
  m <- length(SampleA)
  n <- 11626 - m
  k <- length(SampleB)
  o <- length(intersect(SampleA, SampleB))
  print(noquote(paste("Overlap:", o)))
  print(noquote(paste("p-value:", phyper(q,m,n,k,lower.tail=FALSE))))
}



# Overlap tests begin here

overlaptest(c0up,dtA)
# Overlap: 2
# p-value: 1

overlaptest(c0up,dtB)
# Overlap: 7
# p-value: 1

overlaptest(c0up,dtC)
# Overlap: 25
# p-value: 0.99017

overlaptest(c0up,dtD)
# Overlap: 65
# p-value: 4.35458759521077e-17
c0_dtD_overlap <- intersect(c0up,dtD)
c0_dtD_overlap_table <- c0markersup %>% filter(gene %in% c0_dtD_overlap)
write.table(c0_dtD_overlap_table, "Cluster0_DeeptoolsD_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(c0up,dtE)
# Overlap: 133
# p-value: 1.21596949179157e-124
c0_dtE_overlap <- intersect(c0up,dtE)
c0_dtE_overlap_table <- c0markersup %>% filter(gene %in% c0_dtE_overlap)
write.table(c0_dtE_overlap_table, "Cluster0_DeeptoolsE_Overlap.csv", sep=",", row.names=FALSE)


overlaptest(c1up,dtA)
# Overlap: 2
# p-value: 0.03202

overlaptest(c1up,dtB)
# Overlap: 0
# p-value: 1

overlaptest(c1up,dtC)
# Overlap: 0
# p-value: 1

overlaptest(c1up,dtD)
# Overlap: 0
# p-value: 1

overlaptest(c1up,dtE)
# Overlap: 0
# p-value: 1


overlaptest(c2up,dtA)
# Overlap: 43
# p-value: 3.75130090428883e-12

overlaptest(c2up,dtB)
# Overlap: 21
# p-value: 0.9999

overlaptest(c2up,dtC)
# Overlap: 10
# p-value: 0.8642

overlaptest(c2up,dtD)
# Overlap: 8
# p-value: 0.4800

overlaptest(c2up,dtE)
# Overlap: 2
# p-value: 0.8733


overlaptest(c3up,dtA)
# Overlap: 1
# p-value: 0.3259

overlaptest(c3up,dtB)
# Overlap: 1
# p-value: 0.7672

overlaptest(c3up,dtC)
# Overlap: 0
# p-value: 1

overlaptest(c3up,dtD)
# Overlap: 0
# p-value: 1

overlaptest(c3up,dtE)
# Overlap: 0
# p-value: 1


overlaptest(c4up,dtA)
# Overlap: 36
# p-value: 1.32264637052072e-09

overlaptest(c4up,dtB)
# Overlap: 21
# p-value: 0.9999

overlaptest(c4up,dtC)
# Overlap: 7
# p-value: 0.9556

overlaptest(c4up,dtD)
# Overlap: 8
# p-value: 0.3407

overlaptest(c4up,dtE)
# Overlap: 2
# p-value: 0.8233


overlaptest(c5up,dtA)
# Overlap: 93
# p-value: 0.004069

overlaptest(c5up,dtB)
# Overlap: 158
# p-value: 0.9999

overlaptest(c5up,dtC)
# Overlap: 79
# p-value: 0.01575

overlaptest(c5up,dtD)
# Overlap: 55
# p-value: 0.0008999

overlaptest(c5up,dtE)
# Overlap: 17
# p-value: 0.5233


overlaptest(c6up,dtA)
# Overlap: 56
# p-value: 3.95639825008471e-06

overlaptest(c6up,dtB)
# Overlap: 68
# p-value: 0.9997

overlaptest(c6up,dtC)
# Overlap: 28
# p-value: 0.4756

overlaptest(c6up,dtD)
# Overlap: 13
# p-value: 0.7961

overlaptest(c6up,dtE)
# Overlap: 7
# p-value: 0.6059


overlaptest(c7up,dtA)
# Overlap: 5
# p-value: 0.9999

overlaptest(c7up,dtB)
# Overlap: 37
# p-value: 0.9999

overlaptest(c7up,dtC)
# Overlap: 22
# p-value: 0.8368

overlaptest(c7up,dtD)
# Overlap: 37
# p-value: 1.57939492681163e-07
c7_dtD_overlap <- intersect(c7up,dtD)
c7_dtD_overlap_table <- c7markersup %>% filter(gene %in% c7_dtD_overlap)
write.table(c7_dtD_overlap_table, "Cluster7_DeeptoolsD_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(c7up,dtE)
# Overlap: 64
# p-value: 3.13114368187263e-45
c7_dtE_overlap <- intersect(c7up,dtE)
c7_dtE_overlap_table <- c7markersup %>% filter(gene %in% c7_dtE_overlap)
write.table(c7_dtE_overlap_table, "Cluster7_DeeptoolsE_Overlap.csv", sep=",", row.names=FALSE)
