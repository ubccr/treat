library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis")
mydata <- readRDS(file = "updates.5k.rds")

mydata
dim(mydata) # 9204 rows (genes) and 41712 columns (cells)

# Create meta-data field for TET condition
mydata$tet <- ifelse((mydata$orig.ident=="NOTET")|(mydata$orig.ident=="NOTET_2"),"- DOX","+ DOX")

# Create new meta-data field called celltype.stim
# which combines the cluster number and the treatment group
mydata$celltype.stim <- paste(Idents(mydata), mydata$orig.ident, sep = "_")

# table(mydata$seurat_clusters)
# table(mydata$orig.ident)
# table(mydata$celltype.stim)


# 1B ---------------------------------------------------------------------------

Idents(mydata) <- "integrated_snn_res.0.4"
DimPlot(mydata, reduction = "umap", split.by="tet",
        label=T, label.size=6)
ggsave("Fig1B.pdf", width=8, height=5, dpi=600)

# 1C Cluster proportions -------------------------------------------------------

countTable <- prop.table(table(mydata@meta.data$integrated_snn_res.0.4, 
                               mydata@meta.data$tet),2)
countTable
countTable <- t(countTable)
countTable
countTable <- as.data.frame(countTable)
colnames(countTable) <- c("Condition", "CellType", "Fraction")
countTable
ggplot(countTable, aes(x=CellType, y=Fraction, fill=Condition)) +
  geom_bar(stat = "identity", width=0.7, position=position_dodge(width=0.8), color="black") +
  ylab("Percentage of Cells") + xlab("Cluster") +
  scale_fill_manual(name = "Sample",values = c("#e6e6e6", "#7d7d7d")) +
  theme_classic() + theme(text = element_text(size=18)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0,0.6))
ggsave("Fig1C.pdf", width=7.5, height=5, dpi=600)


# 1D Cluster Marker Genes ------------------------------------------------------

setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis")
mydata <- readRDS(file = "updates.5k.rds")
DefaultAssay(mydata) <- "RNA"
Idents(mydata) <- "integrated_snn_res.0.4"

setwd("C:/Users/tylec/Read lab/")
gaf <- read.table("C:/Users/tylec/Read lab/TbruceiTREU927Annotations.csv", header = T, sep = ",")

allmarkers <- FindAllMarkers(mydata)
allmarkersannotated <- left_join(allmarkers, gaf, by=c('gene' = 'geneID'))
write.table(allmarkersannotated, "K5_highres_MarkersAnnotated.csv", sep=",", row.names=FALSE)

c0markers <- allmarkersannotated %>% filter(cluster==0)
c1markers <- allmarkersannotated %>% filter(cluster==1)
c2markers <- allmarkersannotated %>% filter(cluster==2)
c3markers <- allmarkersannotated %>% filter(cluster==3)
c4markers <- allmarkersannotated %>% filter(cluster==4)
c5markers <- allmarkersannotated %>% filter(cluster==5)
c6markers <- allmarkersannotated %>% filter(cluster==6)
c7markers <- allmarkersannotated %>% filter(cluster==7)

write.table(highres_c0markers, "K5_highresClusters_MarkerGenes_RNA_Annotated_C0.txt", sep="\t", row.names=FALSE)
write.table(highres_c1markers, "K5_highresClusters_MarkerGenes_RNA_Annotated_C1.txt", sep="\t", row.names=FALSE)
write.table(highres_c2markers, "K5_highresClusters_MarkerGenes_RNA_Annotated_C2.txt", sep="\t", row.names=FALSE)
write.table(highres_c3markers, "K5_highresClusters_MarkerGenes_RNA_Annotated_C3.txt", sep="\t", row.names=FALSE)
write.table(highres_c4markers, "K5_highresClusters_MarkerGenes_RNA_Annotated_C4.txt", sep="\t", row.names=FALSE)
write.table(highres_c5markers, "K5_highresClusters_MarkerGenes_RNA_Annotated_C5.txt", sep="\t", row.names=FALSE)
write.table(highres_c6markers, "K5_highresClusters_MarkerGenes_RNA_Annotated_C6.txt", sep="\t", row.names=FALSE)
write.table(highres_c7markers, "K5_highresClusters_MarkerGenes_RNA_Annotated_C7.txt", sep="\t", row.names=FALSE)


# Fig1E dotplot -----------------------------------------------------------

### Dot plots of top enriched marker genes for each cluster

library(stringr)
setwd("C:/Users/tylec/Read lab/DRBD18_Singlecell/9_new_scRNA_Analysis")
mydata <- readRDS(file = "updates.5k.rds")
DefaultAssay(mydata) <- "RNA"
Idents(mydata) <- "integrated_snn_res.0.4"
levels(mydata) <- c("7","6","5","4","3","2","1","0")

DotPlot(mydata, features = c("Tb11.v5.0632",
                             "Tb927.5.2000",
                             "Tb927.10.5250",
                             "Tb927.11.1480",
                             "Tb927.9.12570",
                             "Tb8.NT.11",
                             "Tb927.1.2670",
                             "Tb927.6.3180",
                             "Tb927.11.8970",
                             "Tb927.11.9730",
                             "Tb927.10.7890",
                             "Tb927.9.7470",
                             "Tb927.8.2550",
                             "Tb927.11.13007",
                             "Tb927.10.1100",
                             "Tb927.11.7500",
                             "Tb927.6.510",
                             "Tb927.10.2010",
                             "Tb927.7.2980",
                             "Tb927.11.18430",
                             "Tb927.7.6510",
                             "Tb927.9.210"),
        cols = c('blue','red'), scale.by="size") +
        labs(x = "Gene", y = "Cluster") +
        theme(axis.text = element_text(size = 16),
              axis.title = element_text(size=20),
              axis.title.x = element_text(vjust=-0.5)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
        scale_size(range = c(1,8)) +
        scale_x_discrete(labels = str_wrap(c("Sel1 repeat", "Tb927.5.2000", "ZC3H32",
                                             "GRESAG 4", "GLK1", "Tb8.NT.11",
                                             "PF16","Tb927.6.3180","Ribose 5-phosphate isomerase",
                                             "L34", "Tb927.10.7890", "NT10", "PRI1",
                                             "L40", "L9", "DUF423", "GPEET",
                                             "HK1", "Nitroreductase family", "VSG, putative",
                                             "VSG, degenerate", "VSG, pseudogene"), width = 10))
ggsave("Fig1E.pdf", width=15, height=6, dpi=600)

