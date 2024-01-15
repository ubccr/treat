library(Seurat)
library(dplyr)
library(tidyr)

### Figure 5A
# Deeptools gene lists
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/10_new_DeepToolsAnalysis")
deeptools <- read.csv("DeeptoolsPatternsGeneList.csv", header = T, sep = ",")
deeptools <- deeptools %>% select(X.chrom, start, end, score, strand, thickStart, thickEnd, deepTools_group, name, geneName, Description)

dtA <- deeptools %>% filter(deepTools_group=="pattern A") %>% pull(name)
dtB <- deeptools %>% filter(deepTools_group=="pattern B") %>% pull(name)
dtC <- deeptools %>% filter(deepTools_group=="pattern C") %>% pull(name)
dtD <- deeptools %>% filter(deepTools_group=="pattern D") %>% pull(name)
dtE <- deeptools %>% filter(deepTools_group=="pattern E") %>% pull(name)


# Life cycle gene lists
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

overlaptest(dtA,specific_slender)
# Overlap: 15
# p-value: 0.07304

overlaptest(dtA,specific_stumpy)
# Overlap: 7
# p-value: 0.9996

overlaptest(dtA,specific_earlyPF)
# Overlap: 5
# p-value: 0.7533

overlaptest(dtA,specific_latePF)
# Overlap: 1
# p-value: 0.9231

overlaptest(dtA,specific_MFup)
# Overlap: 127
# p-value: 0.4692

overlaptest(dtA,specific_MFdown)
# Overlap: 151
# p-value: 7.27986744989804e-06
dtA_MFdown_overlap <- intersect(dtA,specific_MFdown)
dtA_MFdown_overlap_list <- deeptools %>% filter(deepTools_group=="pattern A") %>% filter(name %in% dtA_MFdown_overlap)
write.table(dtA_MFdown_overlap_list, "Pattern_A_MFdown_Overlap.csv", sep=",", row.names=FALSE)


overlaptest(dtB,specific_slender)
# Overlap: 24
# p-value: 0.9445

overlaptest(dtB,specific_stumpy)
# Overlap: 26
# p-value: 0.9999

overlaptest(dtB,specific_earlyPF)
# Overlap: 11
# p-value: 0.9930

overlaptest(dtB,specific_latePF)
# Overlap: 7
# p-value: 0.5511

overlaptest(dtB,specific_MFup)
# Overlap: 248
# p-value: 1

overlaptest(dtB,specific_MFdown)
# Overlap: 216
# p-value: 1


overlaptest(dtC,specific_slender)
# Overlap: 9
# p-value: 0.5414

overlaptest(dtC,specific_stumpy)
# Overlap: 11
# p-value: 0.9420

overlaptest(dtC,specific_earlyPF)
# Overlap: 10
# p-value: 0.03064
dtC_earlyPF_overlap <- intersect(dtC,specific_earlyPF)
dtC_earlyPF_overlap_list <- deeptools %>% filter(deepTools_group=="pattern C") %>% filter(name %in% dtC_earlyPF_overlap)
write.table(dtC_earlyPF_overlap_list, "Pattern_C_earlyPF_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(dtC,specific_latePF)
# Overlap: 1
# p-value: 0.8900

overlaptest(dtC,specific_MFup)
# Overlap: 124
# p-value: 0.07182

overlaptest(dtC,specific_MFdown)
# Overlap: 135
# p-value: 7.48683660878462e-06
dtC_MFdown_overlap <- intersect(dtC,specific_MFdown)
dtC_MFdown_overlap_list <- deeptools %>% filter(deepTools_group=="pattern C") %>% filter(name %in% dtC_MFdown_overlap)
write.table(dtC_MFdown_overlap_list, "Pattern_C_MFdown_Overlap.csv", sep=",", row.names=FALSE)


overlaptest(dtD,specific_slender)
# Overlap: 3
# p-value: 0.8948

overlaptest(dtD,specific_stumpy)
# Overlap: 39
# p-value: 8.87747302428672e-16
dtD_stumpy_overlap <- intersect(dtD,specific_stumpy)
dtD_stumpy_overlap_list <- deeptools %>% filter(deepTools_group=="pattern D") %>% filter(name %in% dtD_stumpy_overlap)
write.table(dtD_stumpy_overlap_list, "Pattern_D_stumpy_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(dtD,specific_earlyPF)
# Overlap: 4
# p-value: 0.3616

overlaptest(dtD,specific_latePF)
# Overlap: 2
# p-value: 0.3261

overlaptest(dtD,specific_MFup)
# Overlap: 102
# p-value: 4.20822216148163e-07
dtD_MFup_overlap <- intersect(dtD,specific_MFup)
dtD_MFup_overlap_list <- deeptools %>% filter(deepTools_group=="pattern D") %>% filter(name %in% dtD_MFup_overlap)
write.table(dtD_MFup_overlap_list, "Pattern_D_MFup_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(dtD,specific_MFdown)
# Overlap: 64
# p-value: 0.09473


overlaptest(dtE,specific_slender)
# Overlap: 2
# p-value: 0.6971

overlaptest(dtE,specific_stumpy)
# Overlap: 17
# p-value: 1.19741967135913e-06
dtE_stumpy_overlap <- intersect(dtE,specific_stumpy)
dtE_stumpy_overlap_list <- deeptools %>% filter(deepTools_group=="pattern E") %>% filter(name %in% dtE_stumpy_overlap)
write.table(dtE_stumpy_overlap_list, "Pattern_E_stumpy_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(dtE,specific_earlyPF)
# Overlap: 2
# p-value: 0.4208

overlaptest(dtE,specific_latePF)
# Overlap: 2
# p-value: 0.1010

overlaptest(dtE,specific_MFup)
# Overlap: 96
# p-value: 5.32226681552993e-26
dtE_MFup_overlap <- intersect(dtE,specific_MFup)
dtE_MFup_overlap_list <- deeptools %>% filter(deepTools_group=="pattern E") %>% filter(name %in% dtE_MFup_overlap)
write.table(dtE_MFup_overlap_list, "Pattern_E_MFup_Overlap.csv", sep=",", row.names=FALSE)

overlaptest(dtE,specific_MFdown)
# Overlap: 43
# p-value: 0.0005491
dtE_MFdown_overlap <- intersect(dtE,specific_MFdown)
dtE_MFdown_overlap_list <- deeptools %>% filter(deepTools_group=="pattern E") %>% filter(name %in% dtE_MFdown_overlap)
write.table(dtE_MFdown_overlap_list, "Pattern_E_MFdown_Overlap.csv", sep=",", row.names=FALSE)
