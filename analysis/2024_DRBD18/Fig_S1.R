library(EnhancedVolcano)
library(DESeq2)
setwd("C:/Users/tylec/Read Lab/DRBD18_Singlecell/1_Read_10X-Pilot")
res <- readRDS(file = "results-object.rds")

# define different cell-types that will be shaded
yellow <- c("Tb927.7.5380", "Tb927.7.3730", "Tb927.1.1500", "Tb927.6.820",
            "Tb927.10.5250", "Tb927.10.12760", "Tb927.11.8470", "Tb927.3.5250",
            "Tb927.8.2780", "Tb927.11.530", "Tb927.4.4520", "Tb927.9.4080",
            "Tb927.11.6600", "Tb927.3.1910", "Tb927.10.4430", "Tb10.v4.0033",
            "Tb927.11.14220", "Tb927.10.14550", "Tb927.2.3880", "Tb927.11.12120",
            "Tb927.10.15310", "Tb927.11.16550", "Tb927.7.2680", "Tb11.v5.0730",
            "Tb927.9.10280", "Tb927.10.5150", "Tb927.3.3960", "Tb927.6.640",
            "Tb927.7.2670", "Tb927.3.3930", "Tb927.6.4770", "Tb927.10.12780",
            "Tb927.8.6650", "Tb927.8.4540", "Tb927.10.14930", "Tb927.1.4650",
            "Tb927.4.4230", "Tb927.7.4730", "Tb927.10.12800", "Tb927.10.12660",
            "Tb927.3.2930", "Tb927.4.4580", "Tb927.9.4560")

green <- c("Tb927.10.9600", "Tb927.10.14270", "Tb927.8.7110", "Tb927.2.4200",
           "Tb927.11.14070", "Tb927.11.16790", "Tb927.2.1820", "Tb927.9.11030",
           "Tb927.7.3580", "Tb927.11.7010", "Tb927.3.1630", "Tb927.6.2030",
           "Tb927.6.1800", "Tb11.v5.0390", "Tb927.2.4280", "Tb927.8.870",
           "Tb927.5.3620", "Tb927.1.1380", "Tb927.7.6680", "Tb927.7.3650",
           "Tb927.10.5310", "Tb927.9.11100", "Tb927.10.15020", "Tb11.v5.0245",
           "Tb927.4.2500", "Tb927.10.11450", "Tb927.7.6220", "Tb927.10.8050",
           "Tb927.9.6560", "Tb927.7.4020", "Tb927.8.3550", "Tb927.7.6310",
           "Tb927.7.5770", "Tb927.7.960", "Tb927.9.14120", "Tb927.11.4610")

GeneLabels <- c(yellow, green, "Tb927.11.14090")

keyvals <- ifelse(rownames(res) %in% yellow, 'orange',
                  ifelse(rownames(res) %in% green, 'green',
                         ifelse(rownames(res)=="Tb927.11.14090", 'red', 'gray46')))

names(keyvals)[keyvals == 'gray46'] <- 'Other'
names(keyvals)[keyvals == 'orange'] <- 'RBP'
names(keyvals)[keyvals == 'green'] <- 'Kinase/Phosphatase'
names(keyvals)[keyvals == 'red'] <- 'Target'

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'DRBD18 RNAi RNA-seq',
                subtitle = NULL,
                xlim = c(-3, 5),
                pCutoff = 0.05,
                #FCcutoff = 0.5,
                selectLab = GeneLabels,
                colCustom = keyvals,
                legendPosition = 'top',
                pointSize = 3,
                titleLabSize = 18,
                axisLabSize=20,
                labSize=4,
                legendLabSize = 16) +
  ggplot2::scale_x_continuous(
    breaks=seq(-3,5, 1),
    lim=c(-3,5))
ggsave(filename="BulkRNAvolcano.tiff", width=14, height=8, dpi=600, compression = "lzw")
