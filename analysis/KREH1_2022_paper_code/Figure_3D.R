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

### This code recreates the "gRNA skipping" analysis shown in Figure 3D.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
library(Biostrings)
library(dplyr)
library(tidyr)
library(ggplot2)

# Upload A6 WT and KREH1 KO export files
WT <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(WT$replicate)
WT$knock_down <- as.character(WT$knock_down)
KO <- read.csv(file.choose(), sep = ",", header=TRUE)
unique(KO$replicate)

# Remove pre-edited reads
WTnopre <- WT %>% filter(!(junc_len==0&edit_stop==24))
KOnopre <- KO %>% filter(!(junc_len==0&edit_stop==24))

# Subset reads with ESS 24
WT24 <- WTnopre %>% filter(edit_stop==24)
KO24 <- KOnopre %>% filter(edit_stop==24)
unique(WT24$edit_stop)
unique(KO24$edit_stop)

# Perform the junction analysis each sample individually
WT24_1 <- WT24 %>% filter(replicate==1)
WT24_2 <- WT24 %>% filter(replicate==2)
WT24_3 <- WT24 %>% filter(replicate==3)
WT24_4 <- WT24 %>% filter(replicate==4)
WT24_5 <- WT24 %>% filter(replicate==5)

KO24_1 <- KO24 %>% filter(replicate==1)
KO24_2 <- KO24 %>% filter(replicate==2)

# Use this code to assign one sample to the A6ESS24 variable at a time
# Then run the analysis below
A6ESS24 <- as.data.frame(WT24_1)
A6ESS24 <- as.data.frame(WT24_2)
A6ESS24 <- as.data.frame(WT24_3)
A6ESS24 <- as.data.frame(WT24_4)
A6ESS24 <- as.data.frame(WT24_5)

A6ESS24 <- as.data.frame(KO24_1)
A6ESS24 <- as.data.frame(KO24_2)

# Code for the gRNA skipping analysis starts here #

# Create a table to store the analysis results
results <- data.frame(matrix(NA, nrow = nrow(A6ESS24),
                             ncol = length(seq(1:6))))
pre <- "AAAACACCCAUUUUUAGGAGGAUAAGAGGGGAGAAAAGGGGAAAUGGAAUUGGGAAUUGCCUUUGCCAAACUUUUAGAAGAAAGAGCAGGAAAGGUUAGGGGGAGGAGAGAAGAAAGGGAAAGUUGUGAUUUU"
pres <- strsplit(reverse(pre),"A|G|C") %>% unlist()
nchar(pres)
preES <- nchar(pres)

for (rowIdx in 1:nrow(results)) {
  results[rowIdx, 1] <- A6ESS24[rowIdx,"junc_seq"]
  results[rowIdx, 2] <- A6ESS24[rowIdx,"norm_count"]
  results[rowIdx, 3] <- A6ESS24[rowIdx,"junc_len"]
  results[rowIdx, 4] <- A6ESS24[rowIdx,"junc_end"]
  
  # Split junction into editing sites
  s<-strsplit(reverse(A6ESS24[rowIdx,"junc_seq"]),"A|G|C") %>% unlist()
  
  # Count number of T's at each editing site
  scount <- nchar(s)
  
  # Subset pre-edited sequence to match number of editing sites in junction
  subpre <- preES[c(1:length(scount))]
  
  # Test whether each editing site in junction matches pre-edited number of T's
  test <- scount==subpre
  
  # Record the number of editing sites in junction not matching pre-edited sequence
  results[rowIdx, 5] <- length(which(test=="FALSE"))
  
  # Count how many of the first 15 editing sites match the pre-edited sequence
  totalpre <- sum(test[1:15])
  if(is.na(totalpre)) {totalpre <- "NA"}
  
  # Finally: check if all subset editing sites are pre-edited (equal to 15)
  if(totalpre==15) {
    results[rowIdx, 6] <- "TRUE"
  } else {results[rowIdx, 6] <- "FALSE"}
  
  if(rowIdx==1) {print(noquote("Progress:"))}
  if(rowIdx %% 1000 == 0){
    prog <- round((rowIdx/nrow(results)*100), digits=2)
    print(noquote(c(prog,"%")))}
}

# Name the columns of the results table and save as csv file #
resultstest <- setNames(results, c("seq", "normcount", "junc_len", "junc_end", "ESdiff", "PregRNA1"))
write.table(resultstest, "A6REH1KOESS24_2_junctionanalysis.csv", sep=",", row.names=FALSE)

# Save results table as a unique object for each sample
WT24_1_JA <- resultstest %>% mutate(kd="WT", sample="WT1")
WT24_2_JA <- resultstest %>% mutate(kd="WT", sample="WT2")
WT24_3_JA <- resultstest %>% mutate(kd="WT", sample="WT3")
WT24_4_JA <- resultstest %>% mutate(kd="WT", sample="WT4")
WT24_5_JA <- resultstest %>% mutate(kd="WT", sample="WT5")

KO24_1_JA <- resultstest %>% mutate(kd="KO", sample="KO1")
KO24_2_JA <- resultstest %>% mutate(kd="KO", sample="KO2")

# Combine tables
combined <- bind_rows(WT24_1_JA, WT24_2_JA, WT24_3_JA, WT24_4_JA, WT24_5_JA,
                      KO24_1_JA, KO24_2_JA)

# next... recreate graph in figure 3D
# save each output as a separate object
# join objects and calculate mean, sd for KO and WT
# join KO and WT, make ggplot

combinedtotalES24 <- combined %>% group_by(sample, kd) %>% 
             summarise(total=sum(normcount))

combinedtrue <- combined %>% filter(PregRNA1==TRUE) %>% group_by(sample, kd) %>%
                summarise(pretrue=sum(normcount))

combinedgraph <- full_join(combinedtotalES24,combinedtrue) %>%
                 mutate(ratio=pretrue/total*100) %>% group_by(kd) %>%
                 summarise(mean=mean(ratio),stdev=sd(ratio))
combinedgraph$kd <- factor(combinedgraph$kd, levels=c('WT', 'KO'))
Palette <- c("#bdbdbd", "#FF0000")

ggplot(combinedgraph, aes(x=kd, y=mean, fill=kd)) +
  geom_bar(aes(fill=kd), position="dodge", stat="identity") +
  ggtitle("A6 ESS 24 gRNA1 Pre-edited") +
  ylab("Percent of ESS24 Reads") + xlab("Sample") +
  theme_bw() +
  scale_fill_manual(values=Palette, name="Sample") +
  scale_y_continuous(labels=function(x) paste0(x,"%")) +
  theme(legend.position="none", axis.title.x=element_blank(),
        plot.title = element_text(size=32, hjust = 0.5),
        axis.text = element_text(size=24),
        axis.title = element_text(size=24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=24),
        strip.text.x = element_text(size=24)) +
  geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev),
                width=.15, size=1,
                position=position_dodge(.9))

combinedtest <- full_join(combinedtotalES24,combinedtrue) %>%
  mutate(ratio=pretrue/total*100) %>% group_by(kd) %>%
  summarise(Perc = list(ratio)) %>% spread(kd, Perc) %>%
  mutate(p_value = t.test(unlist(WT), unlist(KO), var.equal = TRUE)$p.value)



#### This code creates a fold change table from average counts in gRNA skipping analysis

WTJAtable <- bind_rows(WT24_1_JA, WT24_2_JA, WT24_3_JA, WT24_4_JA, WT24_5_JA) %>%
  group_by(junc_len, junc_end, seq) %>% mutate(avgnorm=sum(normcount)/5) %>%
  select(seq, avgnorm, junc_len, junc_end, ESdiff, PregRNA1, kd) %>%
  distinct(., .keep_all = TRUE) # distinct keeps only unique rows of a table

KOJAtable <- bind_rows(KO24_1_JA, KO24_2_JA) %>%
  group_by(junc_len, junc_end, seq) %>% mutate(avgnorm=sum(normcount)/2) %>%
  select(seq, avgnorm, junc_len, junc_end, ESdiff, PregRNA1, kd) %>%
  distinct(., .keep_all = TRUE) # distinct keeps only unique rows of a table

A6JAfc <- left_join(KOJAtable, WTJAtable, by=c('junc_len', 'seq')) %>%
  mutate(FCKO=avgnorm.x/avgnorm.y, FCWT=avgnorm.y/avgnorm.x) %>%
  select(seq, junc_len, junc_end.x, ESdiff.x, PregRNA1.x, avgnorm.x, avgnorm.y, FCKO, FCWT)
colnames(A6JAfc) <- c("seq", "junc_len", "junc_end", "ESdiff", "PregRNA1", 
                      "KOavg", "WTavg", "FCKO", "FCWT")
write.table(A6JAfc, "A6REH1KOESS24_JAfc.csv", sep=",", row.names=FALSE)
