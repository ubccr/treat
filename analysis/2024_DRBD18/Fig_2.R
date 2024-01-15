library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis")
mydata <- readRDS(file = "updates.5k.rds")


# 2A ----------------------------------------------------------------------

setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis")
mydata <- readRDS(file = "updates.5k.rds")
DefaultAssay(mydata) <- "RNA"
Idents(mydata) <- "integrated_snn_res.0.4"

allmarkers <- FindAllMarkers(mydata)

c0markers <- allmarkersannotated %>% filter(cluster==0) %>% filter(avg_log2FC > 0)
c0up <- c0markers$gene

c1markers <- allmarkersannotated %>% filter(cluster==1) %>% filter(avg_log2FC > 0)
c1up <- c1markers$gene

c2markers <- allmarkersannotated %>% filter(cluster==2) %>% filter(avg_log2FC > 0)
c2up <- c2markers$gene

c3markers <- allmarkersannotated %>% filter(cluster==3) %>% filter(avg_log2FC > 0)
c3up <- c3markers$gene

c4markers <- allmarkersannotated %>% filter(cluster==4) %>% filter(avg_log2FC > 0)
c4up <- c4markers$gene

c5markers <- allmarkersannotated %>% filter(cluster==5) %>% filter(avg_log2FC > 0)
c5up <- c5markers$gene

c6markers <- allmarkersannotated %>% filter(cluster==6) %>% filter(avg_log2FC > 0)
c6up <- c6markers$gene

c7markers <- allmarkersannotated %>% filter(cluster==7) %>% filter(avg_log2FC > 0)
c7up <- c7markers$gene

# Roditi 2018 #
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/11_GeneLists")

# 57 genes #
specific_slender_table <- read.table("Roditi_2018_slender.csv", header = T, sep = ",")
specific_slender <- specific_slender_table$Gene.ID

# 103 genes #
specific_stumpy_table <- read.table("Roditi_2018_stumpy.csv", header = T, sep = ",")
specific_stumpy <- specific_stumpy_table$ID 

# 34 genes #
specific_earlyPF_table <- read.table("Roditi_2018_earlyPF.csv", header = T, sep = ",")
specific_earlyPF <- specific_earlyPF_table$Gene.ID

# 13 genes #
specific_latePF_table <- read.table("Roditi_2018_latePF.csv", header = T, sep = ",")
specific_latePF <- specific_latePF_table$Gene.ID

# Tschudi 2017 #

specific_MFup_table <- read.table("TableS6_DESeq2_UPinMF.csv", header = T, sep = ",")
specific_MFup <- specific_MFup_table$Gene.ID

specific_MFdown_table <- read.table("TableS6_DESeq2_DOWNinMF.csv", header = T, sep = ",")
specific_MFdown <- specific_MFdown_table$Gene.ID

overlaptest <- function (SampleA, SampleB) {
  q <- length(intersect(SampleA, SampleB)) - 1
  m <- length(SampleA)
  n <- 11626 - m
  k <- length(SampleB)
  o <- length(intersect(SampleA, SampleB))
  print(noquote(paste("Overlap:", o)))
  print(noquote(paste("p-value:", phyper(q,m,n,k,lower.tail=FALSE))))
}


# Overlap tests begin here #

overlaptest(c0up,specific_slender)
# Overlap: 0
# p-value: 1

overlaptest(c0up,specific_stumpy)
# Overlap: 20
# p-value: 1.46178092553897e-14
c0_stumpy_overlap <- intersect(c0up,specific_stumpy)
c0_stumpy_overlap_list <- c0markersup %>% filter(gene %in% c0_stumpy_overlap)
write.table(c0_stumpy_overlap_list, "Cluster0_Stumpy_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(c0up,specific_earlyPF)
# Overlap: 1
# p-value: 0.5025

overlaptest(c0up,specific_latePF)
# Overlap: 1
# p-value: 0.2341

overlaptest(c0up,specific_MFup)
# Overlap: 60
# p-value: 2.45614243644275e-22
c0_MFup_overlap <- intersect(c0up,specific_MFup)
c0_MFup_overlap_list <- c0markersup %>% filter(gene %in% c0_MFup_overlap)
write.table(c0_MFup_overlap_list, "Cluster0_MFup_Overlap.csv", sep=",", row.names=FALSE)

c0_MFup_overlap <- intersect(c0up,specific_MFup)
c0_MFup_overlap_table <- specific_MFup_table %>% filter (Gene.ID %in% c0_MFup_overlap)
write.table(c0_MFup_overlap_table, "Cluster0_MFup_Overlap_table.csv", sep=",", row.names=FALSE)


overlaptest(c0up,specific_MFdown)
# Overlap: 10
# p-value: 0.7977


overlaptest(c1up,specific_slender)
# Overlap: 0
# p-value: 1

overlaptest(c1up,specific_stumpy)
# Overlap: 0
# p-value: 1

overlaptest(c1up,specific_earlyPF)
# Overlap: 0
# p-value: 1

overlaptest(c1up,specific_latePF)
# Overlap: 1
# p-value: 0.002235
c1_latePF_overlap <- intersect(c1up,specific_latePF)
c1_latePF_overlap_list <- c1markersup %>% filter(gene %in% c1_latePF_overlap)
write.table(c1_latePF_overlap_list, "Cluster1_latePF_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(c1up,specific_MFup)
# Overlap: 0
# p-value: 1

overlaptest(c1up,specific_MFdown)
# Overlap: 1
# p-value: 0.1020


overlaptest(c2up,specific_slender)
# Overlap: 0
# p-value: 1

overlaptest(c2up,specific_stumpy)
# Overlap: 1
# p-value: 0.5277

overlaptest(c2up,specific_earlyPF)
# Overlap: 1
# p-value: 0.2187

overlaptest(c2up,specific_latePF)
# Overlap: 2
# p-value: 0.003820
c2_latePF_overlap <- intersect(c2up,specific_latePF)
c2_latePF_overlap_list <- c2markersup %>% filter(gene %in% c2_latePF_overlap)
write.table(c2_latePF_overlap_list, "Cluster2_latePF_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(c2up,specific_MFup)
# Overlap: 2
# p-value: 0.9664

overlaptest(c2up,specific_MFdown)
# Overlap: 26
# p-value: 5.8840205133867e-14
c2_MFdown_overlap <- intersect(c2up,specific_MFdown)
c2_MFdown_overlap_list <- c2markersup %>% filter(gene %in% c2_MFdown_overlap)
write.table(c2_MFdown_overlap_list, "Cluster2_MFdown_Overlap.csv", sep=",", row.names=FALSE)


overlaptest(c3up,specific_slender)
# Overlap: 0
# p-value: 1

overlaptest(c3up,specific_stumpy)
# Overlap: 0
# p-value: 1

overlaptest(c3up,specific_earlyPF)
# Overlap: 0
# p-value: 1

overlaptest(c3up,specific_latePF)
# Overlap: 0
# p-value: 1

overlaptest(c3up,specific_MFup)
# Overlap: 0
# p-value: 1

overlaptest(c3up,specific_MFdown)
# Overlap: 1
# p-value: 0.1020


overlaptest(c4up,specific_slender)
# Overlap: 1
# p-value: 0.3056

overlaptest(c4up,specific_stumpy)
# Overlap: 1
# p-value: 0.4834

overlaptest(c4up,specific_earlyPF)
# Overlap: 1
# p-value: 0.1954

overlaptest(c4up,specific_latePF)
# Overlap: 2
# p-value: 0.002979
c4_latePF_overlap <- intersect(c4up,specific_latePF)
c4_latePF_overlap_list <- c4markersup %>% filter(gene %in% c4_latePF_overlap)
write.table(c4_latePF_overlap_list, "Cluster4_latePF_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(c4up,specific_MFup)
# Overlap: 3
# p-value: 0.8331

overlaptest(c4up,specific_MFdown)
# Overlap: 26
# p-value: 1.8264e-15
c4_MFdown_overlap <- intersect(c4up,specific_MFdown)
c4_MFdown_overlap_list <- c4markersup %>% filter(gene %in% c4_MFdown_overlap)
write.table(c4_MFdown_overlap_list, "Cluster4_MFdown_Overlap.csv", sep=",", row.names=FALSE)


overlaptest(c5up,specific_slender)
# Overlap: 0
# p-value: 1

overlaptest(c5up,specific_stumpy)
# Overlap: 1
# p-value: 0.9737

overlaptest(c5up,specific_earlyPF)
# Overlap: 6
# p-value: 0.0009730
c5_earlyPF_overlap <- intersect(c5up,specific_earlyPF)
c5_earlyPF_overlap_list <- c5markersup %>% filter(gene %in% c5_earlyPF_overlap)
write.table(c5_earlyPF_overlap_list, "Cluster5_earlyPF_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(c5up,specific_latePF)
# Overlap: 2
# p-value: 0.07234

overlaptest(c5up,specific_MFup)
# Overlap: 9
# p-value: 0.9999

overlaptest(c5up,specific_MFdown)
# Overlap: 123
# p-value: 7.81561039358276e-63
c5_MFdown_overlap <- intersect(c5up,specific_MFdown)
c5_MFdown_overlap_list <- c5markersup %>% filter(gene %in% c5_MFdown_overlap)
write.table(c5_MFdown_overlap_list, "Cluster5_MFdown_Overlap.csv", sep=",", row.names=FALSE)


overlaptest(c6up,specific_slender)
# Overlap: 1
# p-value: 0.5796

overlaptest(c6up,specific_stumpy)
# Overlap: 5
# p-value: 0.01956
c6_stumpy_overlap <- intersect(c6up,specific_stumpy)
c6_stumpy_overlap_list <- c6markersup %>% filter(gene %in% c6_stumpy_overlap)
write.table(c6_stumpy_overlap_list, "Cluster6_Stumpy_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(c6up,specific_earlyPF)
# Overlap: 14
# p-value: 1.9431566798477e-17
c6_earlyPF_overlap <- intersect(c6up,specific_earlyPF)
c6_earlyPF_overlap_list <- c6markersup %>% filter(gene %in% c6_earlyPF_overlap)
write.table(c6_earlyPF_overlap_list, "Cluster6_earlyPF_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(c6up,specific_latePF)
# Overlap: 0
# p-value: 1

overlaptest(c6up,specific_MFup)
# Overlap: 9
# p-value: 0.7386

overlaptest(c6up,specific_MFdown)
# Overlap: 55
# p-value: 1.31845640145358e-28
c6_MFdown_overlap <- intersect(c6up,specific_MFdown)
c6_MFdown_overlap_list <- c6markersup %>% filter(gene %in% c6_MFdown_overlap)
write.table(c6_MFdown_overlap_list, "Cluster6_MFdown_Overlap.csv", sep=",", row.names=FALSE)


overlaptest(c7up,specific_slender)
# Overlap: 1
# p-value: 0.5625

overlaptest(c7up,specific_stumpy)
# Overlap: 10
# p-value: 2.13487428553096e-06
c7_stumpy_overlap <- intersect(c7up,specific_stumpy)
c7_stumpy_overlap_list <- c7markersup %>% filter(gene %in% c7_stumpy_overlap)
write.table(c7_stumpy_overlap_list, "Cluster7_Stumpy_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(c7up,specific_earlyPF)
# Overlap: 4
# p-value: 0.001360
c7_earlyPF_overlap <- intersect(c7up,specific_earlyPF)
c7_earlyPF_overlap_list <- c7markersup %>% filter(gene %in% c7_earlyPF_overlap)
write.table(c7_earlyPF_overlap_list, "Cluster7_earlyPF_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(c7up,specific_latePF)
# Overlap: 0
# p-value: 1

overlaptest(c7up,specific_MFup)
# Overlap: 52
# p-value: 6.42431013698541e-24
c7_MFup_overlap <- intersect(c7up,specific_MFup)
c7_MFup_overlap_list <- c7markersup %>% filter(gene %in% c7_MFup_overlap)
write.table(c7_MFup_overlap_list, "Cluster7_MFup_Overlap.csv", sep=",", row.names=FALSE)

c7_MFup_overlap <- intersect(c7up,specific_MFup)
c7_MFup_overlap_table <- specific_MFup_table %>% filter (Gene.ID %in% c7_MFup_overlap)
write.table(c7_MFup_overlap_table, "Cluster7_MFup_Overlap_table.csv", sep=",", row.names=FALSE)


overlaptest(c7up,specific_MFdown)
# Overlap: 9
# p-value: 0.5141



# 2B ----------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
#install.packages("forcats")
library(forcats)
#install.packages("ggdist")
library(ggdist)
#install.packages("gghalves")
library(gghalves)

setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis")
updates.5k <- readRDS(file = "updates.5k.rds")
levels(updates.5k) # check the default clusters

updates.5k$groups <- updates.5k$orig.ident
Idents(updates.5k) <- "groups"
levels(updates.5k)

# Combine replicates into Plus DOX and Minus DOX
updates.5k <- RenameIdents(updates.5k,
                           "NOTET" = "- DOX",
                           "NOTET_2" = "- DOX",
                           "TET" = "+ DOX",
                           "TET_2" = "+ DOX")
updates.5k$groups <- Idents(updates.5k)
levels(updates.5k)


# Upload life cycle genes and create gene lists
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/11_GeneLists")

specific_slender_table <- read.table("Roditi_2018_slender.csv", header = T, sep = ",")
specific_slender <- specific_slender_table$Gene.ID
slender_list <- list(specific_slender)

specific_stumpy_table <- read.table("Roditi_2018_stumpy.csv", header = T, sep = ",")
specific_stumpy <- specific_stumpy_table$ID 
stumpy_list <- list(specific_stumpy)

specific_earlyPF_table <- read.table("Roditi_2018_earlyPF.csv", header = T, sep = ",")
specific_earlyPF <- specific_earlyPF_table$Gene.ID
earlyPF_list <- list(specific_earlyPF)

specific_latePF_table <- read.table("Roditi_2018_latePF.csv", header = T, sep = ",")
specific_latePF <- specific_latePF_table$Gene.ID
latePF_list <- list(specific_latePF)

specific_MFup_table <- read.table("TableS6_DESeq2_UPinMF.csv", header = T, sep = ",")
specific_MFup <- specific_MFup_table$Gene.ID
MFup_list <- list(specific_MFup)

specific_MFdown_table <- read.table("TableS6_DESeq2_DOWNinMF.csv", header = T, sep = ",")
specific_MFdown <- specific_MFdown_table$Gene.ID
MFdown_list <- list(specific_MFdown)



# Module scoring
DefaultAssay(updates.5k) <- "RNA"

updates.5k <- AddModuleScore(
  object = updates.5k,
  features = slender_list,
  name = 'slender_genes'
)

updates.5k <- AddModuleScore(
  object = updates.5k,
  features = stumpy_list,
  name = 'stumpy_genes'
)

updates.5k <- AddModuleScore(
  object = updates.5k,
  features = earlyPF_list,
  name = 'earlyPF_genes'
)

updates.5k <- AddModuleScore(
  object = updates.5k,
  features = latePF_list,
  name = 'latePF_genes'
)

updates.5k <- AddModuleScore(
  object = updates.5k,
  features = MFup_list,
  name = 'MFup_genes'
)

updates.5k <- AddModuleScore(
  object = updates.5k,
  features = MFdown_list,
  name = 'MFdown_genes'
)



# Ridge plots
# Run one dataframe line and make the graph before running next dataframe line
counts.df <- data.frame("cells" = colnames(updates.5k),"condition" = as.factor(updates.5k$groups),"clusters" = updates.5k$integrated_snn_res.0.4,"SlenderScore" =updates.5k$slender_genes1)
counts.df <- data.frame("cells" = colnames(updates.5k),"condition" = as.factor(updates.5k$groups),"clusters" = updates.5k$integrated_snn_res.0.4,"StumpyScore" =updates.5k$stumpy_genes1)
counts.df <- data.frame("cells" = colnames(updates.5k),"condition" = as.factor(updates.5k$groups),"clusters" = updates.5k$integrated_snn_res.0.4,"EarlyPFScore" =updates.5k$earlyPF_genes1)
counts.df <- data.frame("cells" = colnames(updates.5k),"condition" = as.factor(updates.5k$groups),"clusters" = updates.5k$integrated_snn_res.0.4,"LatePFScore" =updates.5k$latePF_genes1)
counts.df <- data.frame("cells" = colnames(updates.5k),"condition" = as.factor(updates.5k$groups),"clusters" = updates.5k$integrated_snn_res.0.4,"MetacyclicUpScore" =updates.5k$MFup_genes1)
counts.df <- data.frame("cells" = colnames(updates.5k),"condition" = as.factor(updates.5k$groups),"clusters" = updates.5k$integrated_snn_res.0.4,"MetacyclicDownScore" =updates.5k$MFdown_genes1)


# Set directory for figure generation
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/ModuleScoring")
# Formatting for plots
myDPI = 600
myFORMAT = "pdf"
my_base_size = 36
my_base_family = "sans"
read_pretty_plots <- function(){
  theme_minimal(base_size = my_base_size,base_family = my_base_family)
}

##Re-make counts.df dataframe with proper module score, i.e. run the correct line of code above ^
##Complete module scoring first and be sure values are stored in object
##Note to make these plots you must change as.numeric value accordingly and y axis label
##Saving Violin Plots
ggsave(filename= "MF_down_genes.pdf",dpi=myDPI, height=12, width=12,
       counts.df %>%
         mutate(clusters = fct_rev(clusters)) %>%
         #plotting for one specific data frame
         #change as.numeric value accordingly
         ggplot(aes(x=clusters,
                    as.numeric(MetacyclicDownScore),
                    fill=condition, color=condition)) + scale_fill_manual(values = c("steelblue", "firebrick3", "green", "orange")) + scale_colour_manual(values = c("steelblue", "firebrick3", "green", "orange"))+ geom_hline(yintercept=0, colour="black") +
         ## add half-violin from {ggdist} package
         ggdist::stat_halfeye(
           adjust = .8, 
           width = .8, 
           .width = 0, 
           height = 5,
           alpha = .6,
           justification = -.5, 
           point_colour = NA
         )  +
         ## add justified jitter from the {gghalves} package
         ##Change ylab accordingly
         gghalves::geom_half_point_panel(
           shape = 21, # dot type
           stroke = .4, # dot outline
           alpha = .2, # dot transparency
           size = .8, # dot size
         ) +
         ylab("MF Down Transcript Enrichment") + xlab("Clusters") + coord_flip() + theme_minimal() + read_pretty_plots() + theme(legend.title=element_blank()) + theme(legend.position = "none")
)
