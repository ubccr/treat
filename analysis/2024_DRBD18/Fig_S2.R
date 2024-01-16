library(Seurat)
setwd("C:/Users/tylec/Read lab/DRBD18_SingleCell/9_new_scRNA_Analysis")
mydata <- readRDS(file = "updates.5k.rds")

# Change TET to DOX
mydata$orig.ident <- ifelse(mydata$orig.ident=="NOTET","NODOX",
                            ifelse(mydata$orig.ident=="NOTET_2", "NODOX_2",
                                   ifelse(mydata$orig.ident=="TET", "DOX",
                                          ifelse(mydata$orig.ident=="TET_2", "DOX_2", print("end")))))

matrix <- AverageExpression(mydata, group.by="orig.ident", slot="counts")$RNA
#pairs(matrix, lower.panel = NULL) # for no R value

upper.panel<-function(x, y){
  points(x,y)
  r <- round(cor(x, y), digits=5)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.75, 0.05, txt)
}

pdf(file = "C:/Users/tylec/Read lab/DRBD18_SingleCell/Figures/FigS2.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches

pairs(matrix, lower.panel = NULL, 
      upper.panel = upper.panel)

dev.off()