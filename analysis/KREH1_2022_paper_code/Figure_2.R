## Copyright 2015 TREAT Authors. All rights reserved.
##
## This file is part of TREAT.
##
## TREAT is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## TREAT is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of 
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with TREAT.  If not, see <http://www.gnu.org/licenses/>.
##
## --------------------------------------------------------------------------
##
## This document contains all R code necessary to recapitulate analysis in this study.
## Additional code can be found in file "StatisticalProcessing.m" for determination of EPSs.


##### Figure 2A #####

# This code will create a table showing pre-edited, fully edited, and partially
# edited mRNA levels for average WT and for KREH1 KO samples.

library(dplyr)
library(tidyr)
library(ggplot2)

# Load the files you're going to need
KREH1KO <- read.csv(file.choose(), sep = ",", header=TRUE)
A63WT <- read.csv(file.choose(), sep = ",", header=TRUE)

# Change WT knock_down column to character class
# (R interprets the WT 29-13 cell line name as an integer)
A63WT$knock_down <- as.character(A63WT$knock_down)

# Check to make sure the correct files were input
unique(KREH1KO$sample) # should see 2 samples
unique(A63WT$sample) # should see 5 samples

# Combine tables into one
A6 <- bind_rows(KREH1KO, A63WT)

# Labeling Pre, Fully, and Partially edited reads
A6 <- A6 %>%
  mutate(status=ifelse((edit_stop==24 & junc_len==0), 'Pre',
                       ifelse((edit_stop==115 & junc_len==0), 'Full',
                              'Partial')))

# Separate and average WT samples
A6.WT <- A6 %>% filter(knock_down=="2913")

# Sum partially edited before calculating mean,sd
A6.pefe.WT <- A6.WT %>% filter(status=='Pre' | status=='Full') %>%
  select(sample, knock_down, replicate, norm_count, status)
A6.partial.WT <- A6.WT %>% filter(status=='Partial') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, norm_count, status) %>%
  distinct(., .keep_all = TRUE)

# Group back together and calculate mean, sd
A6.summary.WT <- bind_rows(A6.pefe.WT, A6.partial.WT) %>%
  group_by(status) %>%
  summarise(mean=mean(norm_count),stdev=sd(norm_count), knock_down="WT")

# Now separate and average knockout Samples
A6.KO <- A6 %>% filter(knock_down=="REH1KO")

# Sum partially edited before calculating mean, sd
A6.pefe.KO <- A6.KO %>% filter(status=='Pre' | status=='Full') %>%
  select(sample, knock_down, replicate, norm_count, status)
A6.partial.KO <- A6.KO %>% filter(status=='Partial') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, norm_count, status) %>%
  distinct(., .keep_all = TRUE)

# Group back together and calculate mean, sd
A6.summary.KO <- bind_rows(A6.pefe.KO, A6.partial.KO) %>%
  group_by(status) %>%
  summarise(mean=mean(norm_count),stdev=sd(norm_count), knock_down="KREH1 KO")


A6forfigure <- bind_rows(A6.summary.WT, A6.summary.KO)
A6forfigure$status <- factor(A6forfigure$status, levels=c('Pre','Partial','Full'))
A6forfigure$knock_down <- factor(A6forfigure$knock_down, levels=c('WT', 'KREH1 KO'))
Palette <- c("#bdbdbd", "#FF0000")

ggplot(A6forfigure, aes(x=status, y=mean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity") +
  ggtitle("A6") +
  ylab("Total sequences (Average)") + xlab("Transcript type") +
  theme_bw() + scale_y_continuous(limits = c(0, 90000)) +
  scale_fill_manual(values=Palette, name="Sample") +
  theme(legend.position="bottom", axis.title.x=element_blank(),
        plot.title = element_text(size=32, hjust = 0.5),
        axis.text = element_text(size=24),
        axis.title = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=24)) +
  geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev),
                width=.15, size=1,
                position=position_dodge(.9))

# Tests of significance
A6wt <- bind_rows(A6.pefe.WT, A6.partial.WT) %>% mutate(knock_down="WT")
A6ko <- bind_rows(A6.pefe.KO, A6.partial.KO)
A6data <- bind_rows(A6wt,A6ko)

A6test <- A6data %>% 
  group_by(knock_down, status) %>% 
  summarise(norm_count = list(norm_count)) %>% 
  spread(knock_down, norm_count) %>% 
  group_by(status) %>% 
  mutate(p_value = t.test(unlist(WT), unlist(REH1KO), var.equal = TRUE)$p.value)
A6test <- A6test[order(A6test$p_value),]
A6test$adj <- p.adjust(A6test$p_value, method = "BY")


##### Figure 2B #####

library(dplyr)
library(tidyr)
library(ggplot2)

# Our EPS calculation is performed by Runpu Chen (see author list in
# Smith et al. 2020, NAR; or McAdams et al. 2019, RNA).
# Generation of files sent to Runpu and analysis of results file are shown here.

# Load the files you're going to need
KREH1KO <- read.csv(file.choose(), sep = ",", header=TRUE)
A63WT <- read.csv(file.choose(), sep = ",", header=TRUE)

# Change WT knock_down column to character class
# (R interprets the WT 29-13 cell line name as an integer)
A63WT$knock_down <- as.character(A63WT$knock_down)

# Check to make sure the correct files were input
unique(KREH1KO$sample) # should see 2 samples
unique(A63WT$sample) # should see 5 samples

# Change KO samples to tetracycline = true
KREH1KO$tetracycline <- "true"

# Combine tables into one
A6 <- bind_rows(KREH1KO, A63WT)

# Set editing sites corresponding to Pre-edited, Fully edited reads
pre <- 24
fe <- 115

# Remove pre-edited and fully edited transcripts and set tetracycline column as logical
A6nopre <- A6 %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))

# Determine the total norm count of sequences at each editing stop site
ESSct <- A6nopre %>% group_by(sample, edit_stop) %>% 
  mutate(nm = sum(norm_count))
ESSctbyES <- ESSct %>% select(sample, edit_stop, tetracycline, replicate, knock_down, nm) %>% 
  distinct(.keep_all = TRUE)

# Save the file you wish to send to Runpu as follows:
write.table(ESSctbyES, "A6KREH1KOEPSforRunpu.csv", sep=",", row.names=FALSE)
###

###
# After results file is received from Runpu, convert to .csv and import:
epsall <- read.table(file.choose(), header=T, sep = ",")

# Take Runpu's data frame and generate the EPS table
avg1 <- epsall %>% group_by(edit_stop) %>% filter(tetracycline==FALSE) %>% 
  mutate(avg = mean(nm)) %>% select(edit_stop, avg) # averaging uninduced samples
epsandavg <- inner_join(avg1, epsall) # rejoin to dataframe
siteseps <- epsandavg %>% group_by(sample, edit_stop) %>% filter(tetracycline==TRUE) %>%
  mutate(epsbyrep = ifelse(((nm>avg)&(pVal<0.05)&(qVal<0.05)),TRUE, FALSE)) # checking whether induced has significantly changed
sites1 <- siteseps %>% group_by(knock_down, edit_stop) %>% 
  mutate(eps = ifelse(epsbyrep==TRUE, TRUE, FALSE)) 
sites2 <- sites1 %>% select(knock_down, edit_stop, replicate, epsbyrep) %>% 
  distinct(knock_down, edit_stop, replicate, .keep_all=TRUE) %>% 
  mutate(replicate = sub("^", "R", replicate)) %>% 
  spread(replicate, epsbyrep,convert=TRUE) %>% rowwise() %>%
  mutate(trueeps = ifelse(R1==TRUE&R2==TRUE, TRUE, FALSE)) # can add any additional replicates here with &
EPSdeterm <- sites2 %>% mutate(EPS = ifelse(trueeps==TRUE,as.character(edit_stop),'false')) %>%
  select(-trueeps) %>% select(knock_down, edit_stop, EPS)
DatafromRunpu <- sites1 %>% select(edit_stop, sample, knock_down, tetracycline, replicate, nm, pVal, qVal, avg)
EPSfromRunpu <- inner_join(EPSdeterm, DatafromRunpu) %>% 
  select(edit_stop, sample, knock_down, tetracycline, replicate, nm, pVal, qVal, avg, EPS) %>%
  distinct(edit_stop, sample, knock_down, tetracycline, replicate, nm, pVal, qVal, avg, EPS, .keep_all=TRUE)

# Save EPS table
write.table(EPSfromRunpu, "A6_EPS_KREH1KO.csv", sep=",", row.names=FALSE)
###

###
# Now create bar charts to quantify norm count changes for each EPS

EPStable <- EPSfromRunpu  # If you have the object already defined in R
# OR #
# Load EPS table
EPStable <- read.table(file.choose(), header = T, sep = ",")

# Need the raw data as well
A6data <- A6nopre  # If you have the object already defined in R
# OR # 
# Load export data
KREH1KO <- read.csv(file.choose(), sep = ",", header=TRUE)
A63WT <- read.csv(file.choose(), sep = ",", header=TRUE)
A63WT$knock_down <- as.character(A63WT$knock_down)
A6 <- bind_rows(KREH1KO, A63WT)
A6nopre <- A6 %>% 
  filter(!(junc_len==0&edit_stop==24),!(junc_len==0&edit_stop==115)) %>%
  mutate(tetracycline = as.logical(tetracycline))
A6data <- A6nopre

# Create table of average counts for each EPS
EPStrue <- EPStable %>% filter(EPS!="FALSE") %>% group_by(knock_down, edit_stop) %>%
  mutate(EPSavg=mean(nm)) %>% select(edit_stop, knock_down, avg, EPSavg) %>%
  distinct(.keep_all = TRUE) %>% gather(sample, normcount, c(avg,EPSavg))
EPStrue[EPStrue$sample=="avg", 2] <- "WT"
EPStrue <- EPStrue %>% select(-sample) %>% distinct(.keep_all = TRUE)
EPSlist <- unique(EPStrue$edit_stop)

# Calculate standard deviation of counts for each EPS from raw data
A6data2 <- A6data %>% filter(edit_stop %in% EPSlist)
A6data2[A6data2$knock_down=="2913",4] <- "WT"
A6data3 <- A6data2 %>% group_by(sample, edit_stop) %>%
  mutate(normcount=sum(norm_count)) %>%
  select(sample, knock_down, edit_stop, normcount) %>%
  distinct(.keep_all = TRUE)
A6data4 <- A6data3 %>% group_by(edit_stop, knock_down) %>%
  summarise(stdev=sd(normcount))

# Combine and prepare data for graph
EPStotal <- left_join(EPStrue,A6data4, by=c('edit_stop', 'knock_down'))
EPStotal[EPStotal$knock_down=="REH1KO",2] <- "KO"
EPStotal$knock_down <- factor(EPStotal$knock_down, 
                              levels=c("WT","KO"))
EPStotal$edit_stop <- factor(EPStotal$edit_stop, levels=c("46", "31", "29", "25", "24"))
Palette <- c("#bdbdbd", "#FF0000")

ggplot() + 
  ylab("Average Norm Count") + ggtitle("A6 Exacerbated Pause Sites: KREH1 KO") +
  geom_bar(data=EPStotal, aes(x = knock_down, y = normcount, fill=knock_down), stat="identity") + theme_bw() + 
  scale_fill_manual(values=Palette) + 
  theme(plot.title = element_text(size=32, hjust=0.5),
        axis.text = element_text(size=24),
        axis.title = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=24),
        axis.title.x = element_blank(),
        legend.position = "none") +
  geom_errorbar(data=EPStotal, aes(x = knock_down, ymin=normcount-stdev, ymax=normcount+stdev),
                width=.15, size=1,
                position=position_dodge(.9)) +
  facet_wrap(~edit_stop, scales="free", nrow=1)


##### Figure 2C and S1 #####

# Block 1 = ESS 41  and JL0
# Block 2 = ESS 48  and JL0
# Block 3 = ESS 59  and JL0
# Block 4 = ESS 73  and JL0  covers gRNAs 4 and 5
# Block 5 = ESS 88  and JL0  covers gRNA 6
# Block 6 = ESS 97  and JL0  covers gRNA 7
# Block 7 = ESS 105 and JL0  covers gRNA 8
# Block 8 = ESS 112 and JL0  covers gRNAs 9 and 10
# Block 9 = ESS 119 and JL0  covers gRNA 11

# Load the files you're going to need
setwd("/Users/tylec/Read lab/REH1 Paper/TREAT raw export files")
A63WT <- read.csv("A6 PF WT export.csv", sep = ",", header=TRUE)
KREH1KO <- read.csv("A6 KREH1 KO export.csv", sep = ",", header=TRUE)

# Fig 2C Blocks 1-3 #
WT1 <- A63WT %>% filter(edit_stop == 41 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="WT", block=1)
WT2 <- A63WT %>% filter(edit_stop == 48 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="WT", block=2)
WT3 <- A63WT %>% filter(edit_stop == 59 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="WT", block=3)

KO1 <- KREH1KO %>% filter(edit_stop == 41 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="REH1KO", block=1)
KO2 <- KREH1KO %>% filter(edit_stop == 48 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="REH1KO", block=2)
KO3 <- KREH1KO %>% filter(edit_stop == 59 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="REH1KO", block=3)

makegraph <- bind_rows(WT1,KO1,WT2,KO2,WT3,KO3) %>% group_by(knock_down, block) %>%
  summarise(mean=mean(total),sd=sd(total))
makegraph[makegraph$knock_down=="REH1KO",1] <- "KREH1 KO"
makegraph$knock_down <- factor(makegraph$knock_down, levels=c('WT', 'KREH1 KO'))

# make ggplot here #
KOPalette <- c("#bdbdbd", "#FF0000")
ggplot(makegraph, aes(x=block, y=mean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity") +
  ggtitle("A6") +
  ylab("Total sequences (Average)") + xlab("Block") +
  theme_bw() + scale_y_continuous(limits = c(0, 3000)) +
  scale_fill_manual(values=KOPalette, name="Sample") +
  theme(legend.position="bottom",
        plot.title = element_text(size=32, hjust = 0.5),
        axis.text = element_text(size=24),
        axis.title = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=24)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                width=.15, size=1,
                position=position_dodge(.9))

combinedtest <- bind_rows(WT1,KO1,WT2,KO2,WT3,KO3)
A6test <- combinedtest %>% 
  group_by(knock_down, block) %>% 
  summarise(total = list(total)) %>% 
  spread(knock_down, total) %>% 
  group_by(block) %>% 
  mutate(p_value = t.test(unlist(WT), unlist(REH1KO), var.equal = TRUE)$p.value)
A6test <- A6test[order(A6test$p_value),]
A6test$adj <- p.adjust(A6test$p_value, method = "BY")


# Fig S1 Blocks 1-9 #
WT1 <- A63WT %>% filter(edit_stop == 41 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="WT", block=1)
WT2 <- A63WT %>% filter(edit_stop == 48 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="WT", block=2)
WT3 <- A63WT %>% filter(edit_stop == 59 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="WT", block=3)
WT4 <- A63WT %>% filter(edit_stop == 73 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="WT", block=4)
WT5 <- A63WT %>% filter(edit_stop == 88 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="WT", block=5)
WT6 <- A63WT %>% filter(edit_stop == 97 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="WT", block=6)
WT7 <- A63WT %>% filter(edit_stop == 105 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="WT", block=7)
WT8 <- A63WT %>% filter(edit_stop == 112 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="WT", block=8)
WT9 <- A63WT %>% filter(edit_stop == 119 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="WT", block=9)

KO1 <- KREH1KO %>% filter(edit_stop == 41 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="REH1KO", block=1)
KO2 <- KREH1KO %>% filter(edit_stop == 48 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="REH1KO", block=2)
KO3 <- KREH1KO %>% filter(edit_stop == 59 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="REH1KO", block=3)
KO4 <- KREH1KO %>% filter(edit_stop == 73 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="REH1KO", block=4)
KO5 <- KREH1KO %>% filter(edit_stop == 88 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="REH1KO", block=5)
KO6 <- KREH1KO %>% filter(edit_stop == 97 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="REH1KO", block=6)
KO7 <- KREH1KO %>% filter(edit_stop == 105 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="REH1KO", block=7)
KO8 <- KREH1KO %>% filter(edit_stop == 112 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="REH1KO", block=8)
KO9 <- KREH1KO %>% filter(edit_stop == 119 & junc_len == 0) %>% group_by(sample) %>%
  summarise(total=sum(norm_count)) %>%
  mutate(knock_down="REH1KO", block=9)

makegraph <- bind_rows(WT1,KO1,WT2,KO2,WT3,KO3,WT4,KO4,WT5,KO5,
                       WT6,KO6,WT7,KO7,WT8,KO8,WT9,KO9) %>%
  group_by(knock_down, block) %>%
  summarise(mean=mean(total),sd=sd(total))
makegraph[makegraph$knock_down=="REH1KO",1] <- "KREH1 KO"
makegraph$knock_down <- factor(makegraph$knock_down, levels=c('WT', 'KREH1 KO'))

# make ggplot here #
KOPalette <- c("#bdbdbd", "#FF0000")
ggplot(makegraph, aes(x=block, y=mean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity") +
  ggtitle("A6") +
  ylab("Total sequences (Average)") + xlab("Block") +
  theme_bw() + scale_y_continuous(limits = c(0, 3000)) +
  scale_fill_manual(values=KOPalette, name="Sample") +
  theme(legend.position="bottom",
        plot.title = element_text(size=32, hjust = 0.5),
        axis.text = element_text(size=24),
        axis.title = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=24)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                width=.15, size=1,
                position=position_dodge(.9))

combinedtest <- bind_rows(WT1,KO1,WT2,KO2,WT3,KO3,WT4,KO4,WT5,KO5,
                          WT6,KO6,WT7,KO7,WT8,KO8,WT9,KO9)
A6test <- combinedtest %>% 
  group_by(knock_down, block) %>% 
  summarise(total = list(total)) %>% 
  spread(knock_down, total) %>% 
  group_by(block) %>% 
  mutate(p_value = t.test(unlist(WT), unlist(REH1KO), var.equal = TRUE)$p.value)
A6test <- A6test[order(A6test$p_value),]
A6test$adj <- p.adjust(A6test$p_value, method = "BY")
