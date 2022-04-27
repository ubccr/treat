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

### This code contains the analysis for Figure 5C and 5D (KREH1 OE in RPS12)

##### Fig 5C #####
library(dplyr)
library(tidyr)
library(ggplot2)

# Load the KREH1 OE RPS12 TREAT export file
KREH1OERPS12 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(KREH1OERPS12$sample) #should see 12 samples
unique(KREH1OERPS12$tetracycline)

# Label each sequence as Pre-edited, fully edited, or partially edited
KREH1OERPS12 <- KREH1OERPS12 %>%
  mutate(status=ifelse((edit_stop==9 & junc_len==0), 'Pre',
                       ifelse((edit_stop==147 & junc_len==0), 'Full',
                              'Partial')))

# Separate and average the uninduced samples
RPS12.un <- KREH1OERPS12 %>% filter(tetracycline=="false")
RPS12.un[,4] <- "Uninduced"
RPS12.pefe.un <- RPS12.un %>% filter(status=='Pre' | status=='Full') %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status)
RPS12.partial.un <- RPS12.un %>% filter(status=='Partial') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
RPS12.totals.un <- bind_rows(RPS12.pefe.un, RPS12.partial.un)
RPS12.avg.un <- RPS12.totals.un %>% group_by(status) %>%
  summarise(mean=sum(norm_count)/6, stdv=sd(norm_count)) %>%
  mutate(knock_down="Uninduced")

# Separate and average the induced samples
RPS12.in <- KREH1OERPS12 %>% filter(tetracycline=="true")
RPS12.in[RPS12.in$knock_down=="REH1_WT_OE",4] <- "WT-OE"
RPS12.in[RPS12.in$knock_down=="REH1_K168A",4] <- "KA-OE"
RPS12.in[RPS12.in$knock_down=="REH1_E268Q",4] <- "EQ-OE"
RPS12.pefe.in <- RPS12.in %>% filter(status=='Pre' | status=='Full') %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status)
RPS12.partial.in <- RPS12.in %>% filter(status=='Partial') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
RPS12.totals.in <- bind_rows(RPS12.pefe.in, RPS12.partial.in)
RPS12.avg.in <- RPS12.totals.in %>% group_by(knock_down, status) %>%
  summarise(mean=sum(norm_count)/2, stdv=sd(norm_count))

# Combine tables and create graph
RPS12forfigure <- bind_rows(RPS12.avg.un, RPS12.avg.in)
RPS12forfigure$status <- factor(RPS12forfigure$status, levels=c('Pre','Partial','Full'))
RPS12forfigure$knock_down <- factor(RPS12forfigure$knock_down, 
                                 levels=c('Uninduced', 'WT-OE', 'KA-OE', 'EQ-OE'))
OEPalette <- c("#bdbdbd", "#118811", "#90BFF9", "#C000C0")

ggplot(RPS12forfigure, aes(x=status, y=mean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity") +
  ggtitle("RPS12") +
  ylab("Total sequences (Average)") + xlab("Transcript type") +
  theme_bw() + scale_y_continuous(limits = c(0, 90000)) +
  scale_fill_manual(values=OEPalette, name="Sample") +
  theme(legend.position="bottom", axis.title.x=element_blank(),
        plot.title = element_text(size=32, hjust = 0.5),
        axis.text = element_text(size=24),
        axis.title = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=24)) +
  geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv),
                width=.15, size=1,
                position=position_dodge(.9))

# Significance tests #
RPS12.totals <- bind_rows(RPS12.totals.un, RPS12.totals.in)
RPS12.totals[RPS12.totals$knock_down=="WT-OE",2] <- "WTOE"
RPS12.totals[RPS12.totals$knock_down=="KA-OE",2] <- "KAOE"
RPS12.totals[RPS12.totals$knock_down=="EQ-OE",2] <- "EQOE"

WTtest1 <- RPS12.totals %>% filter(knock_down=="Uninduced" | knock_down=="WTOE")
WTtest2 <- WTtest1 %>% 
  group_by(knock_down, status) %>% 
  summarise(norm_count = list(norm_count)) %>% 
  spread(knock_down, norm_count) %>% ungroup() %>%
  group_by(status) %>% 
  mutate(p_value = tryCatch(t.test(unlist(Uninduced), unlist(WTOE), var.equal = FALSE)$p.value,
                            error = function(err) return (NA)))
WTtest2 <- WTtest2[order(WTtest2$p_value),]
WTtest2$adj <- p.adjust(WTtest2$p_value, method = "BY")


KAtest1 <- RPS12.totals %>% filter(knock_down=="Uninduced" | knock_down=="KAOE")
KAtest2 <- KAtest1 %>% 
  group_by(knock_down, status) %>% 
  summarise(norm_count = list(norm_count)) %>% 
  spread(knock_down, norm_count) %>% ungroup() %>%
  group_by(status) %>% 
  mutate(p_value = tryCatch(t.test(unlist(Uninduced), unlist(KAOE), var.equal = FALSE)$p.value,
                            error = function(err) return (NA)))
KAtest2 <- KAtest2[order(KAtest2$p_value),]
KAtest2$adj <- p.adjust(KAtest2$p_value, method = "BY")


EQtest1 <- RPS12.totals %>% filter(knock_down=="Uninduced" | knock_down=="EQOE")
EQtest2 <- EQtest1 %>% 
  group_by(knock_down, status) %>% 
  summarise(norm_count = list(norm_count)) %>% 
  spread(knock_down, norm_count) %>% ungroup() %>%
  group_by(status) %>% 
  mutate(p_value = tryCatch(t.test(unlist(Uninduced), unlist(EQOE), var.equal = FALSE)$p.value,
                            error = function(err) return (NA)))
EQtest2 <- EQtest2[order(EQtest2$p_value),]
EQtest2$adj <- p.adjust(EQtest2$p_value, method = "BY")


##### Fig 5D #####

library(dplyr)
library(tidyr)
library(ggplot2)

# Our EPS calculation is performed by Runpu Chen (see author list in
# Smith et al. 2020, NAR; or McAdams et al. 2019, RNA).
# Generation of files sent to Runpu and analysis of results file are shown here.

# Load the files you're going to need
KREH1OERPS12 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(KREH1OERPS12$sample) # Should have 12 samples

# Set editing sites corresponding to Pre-edited, Fully edited reads
pre <- 9
fe <- 147

# Remove pre-edited and fully edited transcripts and set tetracycline column as logical
RPS12nopre <- KREH1OERPS12 %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))

# Determine the total norm count of sequences at each editing stop site
ESSct <- RPS12nopre %>% group_by(sample, edit_stop) %>% 
  mutate(nm = sum(norm_count))
ESSctbyES <- ESSct %>% select(sample, edit_stop, tetracycline, replicate, knock_down, nm) %>% 
  distinct(.keep_all = TRUE)

# Save the file you wish to send to Runpu as follows:
write.table(ESSctbyES, "RPS12KREH1OEEPSforRunpu.csv", sep=",", row.names=FALSE)
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
write.table(EPSfromRunpu, "RPS12_EPS_KREH1KO.csv", sep=",", row.names=FALSE)
###


###
# Now create bar charts to quantify norm count changes for each EPS

EPStable <- EPSfromRunpu  # If you have the object already defined in R
# OR #
# Load EPS table
EPStable <- read.table(file.choose(), header = T, sep = ",")

# Need the raw data as well
RPS12data <- RPS12nopre  # If you have the object already defined in R
# OR # 
# Load export data
KREH1OERPS12 <- read.csv(file.choose(), sep = ",", header=TRUE)
RPS12nopre <- KREH1OERPS12 %>% 
  filter(!(junc_len==0&edit_stop==9),!(junc_len==0&edit_stop==147)) %>%
  mutate(tetracycline = as.logical(tetracycline))
RPS12data <- RPS12nopre

# Create table of average counts for each EPS
EPStrue <- EPStable %>% filter(EPS!="FALSE") %>% group_by(knock_down, edit_stop) %>%
  mutate(EPSavg=mean(nm)) %>% select(edit_stop, knock_down, avg, EPSavg) %>%
  distinct(.keep_all = TRUE) %>% gather(sample, normcount, c(avg,EPSavg))
EPStrue[EPStrue$sample=="avg", 2] <- "AvgUn"
EPStrue <- EPStrue %>% select(-sample) %>% distinct(.keep_all = TRUE)
EPSlist <- unique(EPStrue$edit_stop)

# Calculate standard deviation of counts for each EPS from raw data
RPS12data2 <- RPS12data %>% filter(edit_stop %in% EPSlist)
RPS12data2[RPS12data2$tetracycline=="FALSE",4] <- "AvgUn"
RPS12data3 <- RPS12data2 %>% group_by(sample, edit_stop) %>%
  mutate(normcount=sum(norm_count)) %>%
  select(sample, knock_down, edit_stop, normcount) %>%
  distinct(.keep_all = TRUE)
RPS12data4 <- RPS12data3 %>% group_by(edit_stop, knock_down) %>%
  summarise(stdev=sd(normcount))

# Combine and perpare data for graph
EPStotal <- left_join(EPStrue,RPS12data4, by=c('edit_stop', 'knock_down'))
EPStotal[EPStotal$knock_down=="REH1_WT_OE", 2] <- "KREH1-WT"
EPStotal[EPStotal$knock_down=="REH1_K168A", 2] <- "KREH1-KA"
EPStotal[EPStotal$knock_down=="REH1_E268Q", 2] <- "KREH1-EQ"
EPStotal$knock_down <- factor(EPStotal$knock_down, 
                              levels=c("AvgUn","KREH1-WT","KREH1-KA","KREH1-EQ"))
ColorPalette <- c("#bdbdbd", "#118811", "#90BFF9", "#C000C0")

ggplot() + 
  ylab("Average Norm Count") + ggtitle("RPS12 Exacerbated Pause Sites: KREH1 OE") +
  geom_bar(data=EPStotal, aes(x = knock_down, y = normcount, fill=knock_down), stat="identity") + theme_bw() + 
  scale_fill_manual(values=ColorPalette) + 
  theme(plot.title = element_text(size=24, hjust=0.5),
        axis.text = element_text(size=15),
        axis.title = element_text(size=22),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=18),
        axis.title.x = element_blank(),
        legend.position = "none") +
  #labs(fill="Sample") +
  geom_errorbar(data=EPStotal, aes(x = knock_down, ymin=normcount-stdev, ymax=normcount+stdev),
                width=.15, size=1,
                position=position_dodge(.9)) +
  facet_wrap(~edit_stop, scales="free")
