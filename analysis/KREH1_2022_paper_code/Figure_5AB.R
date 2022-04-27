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

### This code contains the analysis for Figure 5A and 5B (KREH1 OE in A6)

##### Fig 5A #####
library(dplyr)
library(tidyr)
library(ggplot2)

# Load the KREH1 OE A6 TREAT export file
KREH1OEA6 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(KREH1OEA6$sample) #should see 12 samples
unique(KREH1OEA6$tetracycline)

# Label each sequence as Pre-edited, fully edited, or partially edited
KREH1OEA6 <- KREH1OEA6 %>%
  mutate(status=ifelse((edit_stop==24 & junc_len==0), 'Pre',
                       ifelse((edit_stop==115 & junc_len==0), 'Full',
                              'Partial')))

# Separate and average the uninduced samples
A6.un <- KREH1OEA6 %>% filter(tetracycline=="false")
A6.un[,4] <- "Uninduced"
A6.pefe.un <- A6.un %>% filter(status=='Pre' | status=='Full') %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status)
A6.partial.un <- A6.un %>% filter(status=='Partial') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
A6.totals.un <- bind_rows(A6.pefe.un, A6.partial.un)
A6.avg.un <- A6.totals.un %>% group_by(status) %>%
  summarise(mean=sum(norm_count)/6, stdv=sd(norm_count)) %>%
  mutate(knock_down="Uninduced")
  
# Separate and average the induced samples
A6.in <- KREH1OEA6 %>% filter(tetracycline=="true")
A6.in[A6.in$knock_down=="REH1_WT_OE",4] <- "WT-OE"
A6.in[A6.in$knock_down=="REH1_K168A",4] <- "KA-OE"
A6.in[A6.in$knock_down=="REH1_E268Q",4] <- "EQ-OE"
A6.pefe.in <- A6.in %>% filter(status=='Pre' | status=='Full') %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status)
A6.partial.in <- A6.in %>% filter(status=='Partial') %>%
  group_by(sample) %>% mutate(norm_count=sum(norm_count)) %>%
  select(sample, knock_down, replicate, tetracycline, norm_count, status) %>%
  distinct(., .keep_all = TRUE)
A6.totals.in <- bind_rows(A6.pefe.in, A6.partial.in)
A6.avg.in <- A6.totals.in %>% group_by(knock_down, status) %>%
  summarise(mean=sum(norm_count)/2, stdv=sd(norm_count))

# Combine tables and create graph
A6forfigure <- bind_rows(A6.avg.un, A6.avg.in)
A6forfigure$status <- factor(A6forfigure$status, levels=c('Pre','Partial','Full'))
A6forfigure$knock_down <- factor(A6forfigure$knock_down, 
                                 levels=c('Uninduced', 'WT-OE', 'KA-OE', 'EQ-OE'))
OEPalette <- c("#bdbdbd", "#118811", "#90BFF9", "#C000C0")

ggplot(A6forfigure, aes(x=status, y=mean, fill=knock_down)) +
  geom_bar(aes(fill=knock_down), position="dodge", stat="identity") +
  ggtitle("A6") +
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
A6.totals <- bind_rows(A6.totals.un, A6.totals.in)
A6.totals[A6.totals$knock_down=="WT-OE",2] <- "WTOE"
A6.totals[A6.totals$knock_down=="KA-OE",2] <- "KAOE"
A6.totals[A6.totals$knock_down=="EQ-OE",2] <- "EQOE"

WTtest1 <- A6.totals %>% filter(knock_down=="Uninduced" | knock_down=="WTOE")
WTtest2 <- WTtest1 %>% 
  group_by(knock_down, status) %>% 
  summarise(norm_count = list(norm_count)) %>% 
  spread(knock_down, norm_count) %>% ungroup() %>%
  group_by(status) %>% 
  mutate(p_value = tryCatch(t.test(unlist(Uninduced), unlist(WTOE), var.equal = FALSE)$p.value,
                            error = function(err) return (NA)))
WTtest2 <- WTtest2[order(WTtest2$p_value),]
WTtest2$adj <- p.adjust(WTtest2$p_value, method = "BY")


KAtest1 <- A6.totals %>% filter(knock_down=="Uninduced" | knock_down=="KAOE")
KAtest2 <- KAtest1 %>% 
  group_by(knock_down, status) %>% 
  summarise(norm_count = list(norm_count)) %>% 
  spread(knock_down, norm_count) %>% ungroup() %>%
  group_by(status) %>% 
  mutate(p_value = tryCatch(t.test(unlist(Uninduced), unlist(KAOE), var.equal = FALSE)$p.value,
                            error = function(err) return (NA)))
KAtest2 <- KAtest2[order(KAtest2$p_value),]
KAtest2$adj <- p.adjust(KAtest2$p_value, method = "BY")


EQtest1 <- A6.totals %>% filter(knock_down=="Uninduced" | knock_down=="EQOE")
EQtest2 <- EQtest1 %>% 
  group_by(knock_down, status) %>% 
  summarise(norm_count = list(norm_count)) %>% 
  spread(knock_down, norm_count) %>% ungroup() %>%
  group_by(status) %>% 
  mutate(p_value = tryCatch(t.test(unlist(Uninduced), unlist(EQOE), var.equal = FALSE)$p.value,
                            error = function(err) return (NA)))
EQtest2 <- EQtest2[order(EQtest2$p_value),]
EQtest2$adj <- p.adjust(EQtest2$p_value, method = "BY")



##### Fig 5B #####

library(dplyr)
library(tidyr)
library(ggplot2)

# Our EPS calculation is performed by Runpu Chen (see author list in
# Smith et al. 2020, NAR; or McAdams et al. 2019, RNA).
# Generation of files sent to Runpu and analysis of results file are shown here.

# Load the files you're going to need
KREH1OEA6 <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(KREH1OEA6$sample) # Should have 12 samples

# Set editing sites corresponding to Pre-edited, Fully edited reads
pre <- 24
fe <- 115

# Remove pre-edited and fully edited transcripts and set tetracycline column as logical
A6nopre <- KREH1OEA6 %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))

# Determine the total norm count of sequences at each editing stop site
ESSct <- A6nopre %>% group_by(sample, edit_stop) %>% 
  mutate(nm = sum(norm_count))
ESSctbyES <- ESSct %>% select(sample, edit_stop, tetracycline, replicate, knock_down, nm) %>% 
  distinct(.keep_all = TRUE)

# Save the file you wish to send to Runpu as follows:
write.table(ESSctbyES, "A6KREH1OEEPSforRunpu.csv", sep=",", row.names=FALSE)
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
KREH1OEA6 <- read.csv(file.choose(), sep = ",", header=TRUE)
A6nopre <- KREH1OEA6 %>% 
  filter(!(junc_len==0&edit_stop==24),!(junc_len==0&edit_stop==115)) %>%
  mutate(tetracycline = as.logical(tetracycline))
A6data <- A6nopre

# Create table of average counts for each EPS
EPStrue <- EPStable %>% filter(EPS!="FALSE") %>% group_by(knock_down, edit_stop) %>%
  mutate(EPSavg=mean(nm)) %>% select(edit_stop, knock_down, avg, EPSavg) %>%
  distinct(.keep_all = TRUE) %>% gather(sample, normcount, c(avg,EPSavg))
EPStrue[EPStrue$sample=="avg", 2] <- "AvgUn"
EPStrue <- EPStrue %>% select(-sample) %>% distinct(.keep_all = TRUE)
EPSlist <- unique(EPStrue$edit_stop)

# Calculate standard deviation of counts for each EPS from raw data
A6data2 <- A6data %>% filter(edit_stop %in% EPSlist)
A6data2[A6data2$tetracycline=="FALSE",4] <- "AvgUn"
A6data3 <- A6data2 %>% group_by(sample, edit_stop) %>%
  mutate(normcount=sum(norm_count)) %>%
  select(sample, knock_down, edit_stop, normcount) %>%
  distinct(.keep_all = TRUE)
A6data4 <- A6data3 %>% group_by(edit_stop, knock_down) %>%
  summarise(stdev=sd(normcount))

# Combine and perpare data for graph
EPStotal <- left_join(EPStrue,A6data4, by=c('edit_stop', 'knock_down'))
EPStotal[EPStotal$knock_down=="REH1_WT_OE", 2] <- "KREH1-WT"
EPStotal[EPStotal$knock_down=="REH1_K168A", 2] <- "KREH1-KA"
EPStotal[EPStotal$knock_down=="REH1_E268Q", 2] <- "KREH1-EQ"
EPStotal$knock_down <- factor(EPStotal$knock_down, 
                              levels=c("AvgUn","KREH1-WT","KREH1-KA","KREH1-EQ"))
ColorPalette <- c("#bdbdbd", "#118811", "#90BFF9", "#C000C0")

ggplot() + 
  ylab("Average Norm Count") + ggtitle("A6 Exacerbated Pause Sites: KREH1 OE") +
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


