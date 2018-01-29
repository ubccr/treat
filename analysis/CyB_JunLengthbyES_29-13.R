library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)

CyB.WT.all <- read.table(file.choose(), header = T, sep = ",")

cbbPalette <- c("#000000", "#D55E00", "#F0E442", "#009E73", "#E69F00", "#0072B2", "#56B4E9", "#CC79A7")
labels <- c("0","1-5","6-10",">10")

# introduce junction length bins
CyB.WT.all$bin1 <- cut(CyB.WT.all$junc_len, breaks = c(-Inf,1,6,11,Inf), right = FALSE)
levels <- levels(CyB.WT.all$bin1)

WT <- CyB.WT.all %>% filter(knock_down == "29-13") %>% group_by(edit_stop, bin1) %>%
  summarize(Avg = sum(norm_count)/5) %>% mutate(KD = "29-13")

#AvgUninduced <- RPS12allCl %>% filter(KD !="29-13" & TET == "FALSE") %>% group_by(edit_stop, bin1) %>%
#  summarise(Avg = sum(norm_count)/8) %>% mutate(KD = "AvgUn")
#AvgInduced <- RPS12allCl %>% filter(TET == "TRUE") %>% group_by(edit_stop, KD, bin1) %>%
#  summarise(Avg = sum(norm_count)/2)

Totaldata <- WT

ggplot() + 
  ylab("Percentage") + xlab("Editing Site") + ggtitle("CyB Junction Lengths by Editing Site (WT)") +
  geom_bar(data=Totaldata, aes(x = edit_stop, y = Avg, fill=bin1), stat="identity", position="fill") +
  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) +
  scale_x_reverse() + theme(legend.position = "bottom")

##### Remove junction length 0 from first ES and from last ES #####

Totaldata1 <- Totaldata %>% filter(!(edit_stop==558 & bin1=="[-Inf,1)"),
                                   !(edit_stop==582 & bin1=="[-Inf,1)"))

plot1 <- ggplot() + 
  ylab("Percentage") + xlab("Editing Site") + ggtitle("CyB Junction Lengths by Editing Site (WT)") +
  geom_bar(data=Totaldata1, aes(x = edit_stop, y = Avg, fill=bin1),stat="identity", position="fill") +
  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) + 
  scale_x_reverse() + theme(legend.position = "bottom") + theme(legend.position="none")

##### Plot abundance by ES #####


# stacked...
plot2 <- ggplot() + 
  ylab("Total Sequences") + xlab("Editing Site") + ggtitle("CyB Average Transcript Levels by Editing Site") +
  geom_bar(data=Totaldata1, aes(x = edit_stop, y = Avg, fill = bin1), stat="identity", position="stack") +
  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) + 
  scale_x_reverse() + coord_trans(y = "sqrt") +
  theme(legend.position = "bottom") + theme(legend.position="none")

################################################################################
# raw abundance in a single color...
#Totaldata3 <- Totaldata2 %>% group_by(KD, edit_stop) %>% summarize(sum=sum(Avg))

#plot2 <- ggplot() + 
#  ylab("Total Sequences") + xlab("Editing Site") + ggtitle("RPS12 Transcript Levels by Editing Site") +
#  geom_bar(data=Totaldata3, aes(x = edit_stop, y = sum), stat="identity") +
#  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) + 
#  scale_x_reverse(breaks=c(150,140,130,120,110,100,90,80,70,60,50,40,30,20,10)) +
#  facet_grid(KD~.) + coord_trans(y = "sqrt")

# to stack two plots on top of one another
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))
