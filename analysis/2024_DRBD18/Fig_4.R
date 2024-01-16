
# 4A ----------------------------------------------------------------------

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
levels(updates.5k) # default clusters

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


# Upload deeptools patterns and create gene lists

setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/10_new_DeepToolsAnalysis")
deeptools <- read.csv("DeeptoolsPatternsGeneList.csv", header = T, sep = ",")
deeptools <- deeptools %>% select(X.chrom, start, end, score, strand, thickStart, thickEnd, deepTools_group, name, geneName, Description)

dtA_list <- list(deeptools %>% filter(deepTools_group=="pattern A") %>% pull(name))
dtB_list <- list(deeptools %>% filter(deepTools_group=="pattern B") %>% pull(name))
dtC_list <- list(deeptools %>% filter(deepTools_group=="pattern C") %>% pull(name))
dtD_list <- list(deeptools %>% filter(deepTools_group=="pattern D") %>% pull(name))
dtE_list <- list(deeptools %>% filter(deepTools_group=="pattern E") %>% pull(name))


# Module scoring #
DefaultAssay(updates.5k) <- "RNA"

updates.5k <- AddModuleScore(
  object = updates.5k,
  features = dtA_list,
  name = 'dtA_genes'
)

updates.5k <- AddModuleScore(
  object = updates.5k,
  features = dtB_list,
  name = 'dtB_genes'
)

updates.5k <- AddModuleScore(
  object = updates.5k,
  features = dtC_list,
  name = 'dtC_genes'
)

updates.5k <- AddModuleScore(
  object = updates.5k,
  features = dtD_list,
  name = 'dtD_genes'
)

updates.5k <- AddModuleScore(
  object = updates.5k,
  features = dtE_list,
  name = 'dtE_genes'
)



# Ridge plots
countsA.df <- data.frame("cells" = colnames(updates.5k),
                         "condition" = as.factor(updates.5k$groups),
                         "pattern" = "A",
                         "DeeptoolsPattern" = updates.5k$dtA_genes1)
countsB.df <- data.frame("cells" = colnames(updates.5k),
                         "condition" = as.factor(updates.5k$groups),
                         "pattern" = "B",
                         "DeeptoolsPattern" = updates.5k$dtB_genes1)
countsC.df <- data.frame("cells" = colnames(updates.5k),
                         "condition" = as.factor(updates.5k$groups),
                         "pattern" = "C",
                         "DeeptoolsPattern" = updates.5k$dtC_genes1)
countsD.df <- data.frame("cells" = colnames(updates.5k),
                         "condition" = as.factor(updates.5k$groups),
                         "pattern" = "D",
                         "DeeptoolsPattern" = updates.5k$dtD_genes1)
countsE.df <- data.frame("cells" = colnames(updates.5k),
                         "condition" = as.factor(updates.5k$groups),
                         "pattern" = "E",
                         "DeeptoolsPattern" = updates.5k$dtE_genes1)

allcounts.df <- bind_rows(countsA.df, countsB.df, countsC.df, countsD.df, countsE.df)



#Set directory for figure generation
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/ModuleScoring")
# Formatting for plots
myDPI = 600
myFORMAT = "pdf"
my_base_size = 42
my_base_family = "sans"
read_pretty_plots <- function(){
  theme_minimal(base_size = my_base_size,base_family = my_base_family)
}
ggsave(filename= "ModuleScore_allPatterns.pdf",dpi=myDPI, height=12, width=12,
       allcounts.df %>%
         mutate(pattern = fct_rev(pattern)) %>%
         #plotting for one specific data frame
         #change as.numeric value accordingly
         ggplot(aes(x=pattern,
                    as.numeric(DeeptoolsPattern),
                    fill=condition, color=condition)) + 
         scale_fill_manual(values = c("steelblue", "firebrick3", "green", "orange")) + 
         scale_colour_manual(values = c("steelblue", "firebrick3", "green", "orange")) + 
         geom_hline(yintercept=0, colour="black") +
         ## add half-violin from {ggdist} package
         ggdist::stat_halfeye(
           adjust = 5, 
           width = 1, 
           .width = 0, 
           #height = 20,
           alpha = .6,
           #justification = -.5, 
           point_colour = NA
         )  +
         ## add justified jitter from the {gghalves} package
         ##Change ylab accordingly
         #gghalves::geom_half_point_panel(
         #   shape = 21,
         #   stroke = .4,
         #   alpha = .2,
         #   size = .8,
         # ) +
         ylab("Transcript Enrichment") + xlab("Pattern") + coord_flip() + theme_minimal() + read_pretty_plots() + theme(legend.title=element_blank(),
                                                                                                                        legend.position = "none")
)



# 4C ----------------------------------------------------------------------

# Deeptools gene lists
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/10_new_DeepToolsAnalysis")
deeptools <- read.csv("DeeptoolsPatternsGeneList.csv", header = T, sep = ",")
deeptools <- deeptools %>% select(X.chrom, start, end, score, strand, thickStart, thickEnd, deepTools_group, name, geneName, Description)

dtA <- deeptools %>% filter(deepTools_group=="pattern A") %>% pull(name)
dtB <- deeptools %>% filter(deepTools_group=="pattern B") %>% pull(name)
dtC <- deeptools %>% filter(deepTools_group=="pattern C") %>% pull(name)
dtD <- deeptools %>% filter(deepTools_group=="pattern D") %>% pull(name)
dtE <- deeptools %>% filter(deepTools_group=="pattern E") %>% pull(name)

# Export and TE gene lists
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/11_GeneLists")

export_impaired <- read.table("AMAR_export_impaired.csv", header = T, sep = ",")
expimp <- export_impaired$transcript

TE_up <- read.table("TE_up.csv", header = F, sep = ",") %>% pull(V1)
TE_down <- read.table("TE_down.csv", header = F, sep = ",") %>% pull(V1)


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

overlaptest(TE_up,dtA)
# Overlap: 4
# p-value: 0.9512

overlaptest(TE_up,dtB)
# Overlap: 4
# p-value: 0.9999

overlaptest(TE_up,dtC)
# Overlap: 9
# p-value: 0.1798

overlaptest(TE_up,dtD)
# Overlap: 8
# p-value: 0.02656

overlaptest(TE_up,dtE)
# Overlap: 13
# p-value: 6.42340458437634e-09


overlaptest(TE_down,dtA)
# Overlap: 17
# p-value: 0.03593

overlaptest(TE_down,dtB)
# Overlap: 28
# p-value: 0.8521

overlaptest(TE_down,dtC)
# Overlap: 9
# p-value: 0.6275

overlaptest(TE_down,dtD)
# Overlap: 0
# p-value: 1

overlaptest(TE_down,dtE)
# Overlap: 2
# p-value: 0.7322


overlaptest(expimp,dtA)
# Overlap: 41
# p-value: 8.61479412114705e-09

overlaptest(expimp,dtB)
# Overlap: 40
# p-value: 0.9769

overlaptest(expimp,dtC)
# Overlap: 10
# p-value: 0.9412

overlaptest(expimp,dtD)
# Overlap: 4
# p-value: 0.9751

overlaptest(expimp,dtE)
# Overlap: 0
# p-value: 1