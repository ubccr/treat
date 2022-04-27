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

## Code for the gRNA skipping analysis shown in Figure 6B (KREH1 OE for RPS12)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
library(Biostrings)
library(dplyr)
library(tidyr)
library(ggplot2)

# Upload TREAT export file (raw data)
KREH1OERPS12 <- read.csv(file.choose(), sep = ",", header=TRUE)

# Remove pre-edited reads
KREH1OERPS12nopre <- KREH1OERPS12 %>% filter(!(junc_len==0&edit_stop==9))

# Subset reads with ESS 24
KREH1OERPS12_9 <- KREH1OERPS12nopre %>% filter(edit_stop==9)

# Perform the junction analysis each sample individually
OEUN_WT1 <- KREH1OERPS12_9 %>% filter(sample=="REH1_MinTet_Rep1_RPS12_OE")
OEUN_WT2 <- KREH1OERPS12_9 %>% filter(sample=="REH1_MinTet_Rep2_RPS12_OE")
OEUN_KA1 <- KREH1OERPS12_9 %>% filter(sample=="REH1_MinTet_Rep1_RPS12_M1")
OEUN_KA2 <- KREH1OERPS12_9 %>% filter(sample=="REH1_MinTet_Rep2_RPS12_M1")
OEUN_EQ1 <- KREH1OERPS12_9 %>% filter(sample=="REH1_MinTet_Rep1_RPS12_E1")
OEUN_EQ2 <- KREH1OERPS12_9 %>% filter(sample=="REH1_MinTet_Rep2_RPS12_E1")

OEIN_WT1 <- KREH1OERPS12_9 %>% filter(sample=="REH1_PlusTet_Rep1_RPS12_OE")
OEIN_WT2 <- KREH1OERPS12_9 %>% filter(sample=="REH1_PlusTet_Rep2_RPS12_OE")
OEIN_KA1 <- KREH1OERPS12_9 %>% filter(sample=="REH1_PlusTet_Rep1_RPS12_M1")
OEIN_KA2 <- KREH1OERPS12_9 %>% filter(sample=="REH1_PlusTet_Rep2_RPS12_M1")
OEIN_EQ1 <- KREH1OERPS12_9 %>% filter(sample=="REH1_PlusTet_Rep1_RPS12_E1")
OEIN_EQ2 <- KREH1OERPS12_9 %>% filter(sample=="REH1_PlusTet_Rep2_RPS12_E1")

# Use this code to assign one sample to the RPS12ESS24 variable at a time
# Then run the analysis below
RPS12ESS9 <- OEUN_WT1
RPS12ESS9 <- OEUN_WT2
RPS12ESS9 <- OEUN_KA1
RPS12ESS9 <- OEUN_KA2
RPS12ESS9 <- OEUN_EQ1
RPS12ESS9 <- OEUN_EQ2

RPS12ESS9 <- OEIN_WT1
RPS12ESS9 <- OEIN_WT2
RPS12ESS9 <- OEIN_KA1
RPS12ESS9 <- OEIN_KA2
RPS12ESS9 <- OEIN_EQ1
RPS12ESS9 <- OEIN_EQ2


# Code for the gRNA skipping analysis starts here #

# Create a table to store the analysis results
results <- data.frame(matrix(NA, nrow = nrow(RPS12ESS9),
                             ncol = length(seq(1:6))))
pre <- "UAAAGUAAAAAGGCGAGGAUUUUUUGAGUGGGACUGGAGAGAAAGAGCCGUUCGAGCCCAGCCGGAACCGACGGAGAGCUUCUUUUGAAUAAAAGGGAGGCGGGGAGGAGAGUUUCAAAAAGAUUUGGGUGGGGGGAACCCUUUGUUUUGGUUAAAGAAACAUCGUUUAGAAGAGAUUUU"
pres <- strsplit(reverse(pre),"A|G|C") %>% unlist()
nchar(pres)
preES <- nchar(pres)

for (rowIdx in 1:nrow(results)) {
  results[rowIdx, 1] <- RPS12ESS9[rowIdx,"junc_seq"]
  results[rowIdx, 2] <- RPS12ESS9[rowIdx,"norm_count"]
  results[rowIdx, 3] <- RPS12ESS9[rowIdx,"junc_len"]
  results[rowIdx, 4] <- RPS12ESS9[rowIdx,"junc_end"]
  
  # Split junction into editing sites
  s<-strsplit(reverse(RPS12ESS9[rowIdx,"junc_seq"]),"A|G|C") %>% unlist()
  
  # Count number of T's at each editing site
  scount <- nchar(s)
  
  # Subset pre-edited sequence to match number of editing sites in junction
  subpre <- preES[c(1:length(scount))]
  
  # Test whether each editing site in junction matches pre-edited number of T's
  test <- scount==subpre
  
  # Record the number of editing sites in junction not matching pre-edited sequence
  results[rowIdx, 5] <- length(which(test=="FALSE"))
  
  # Count how many of the first 15 editing sites match the pre-edited sequence
  totalpre <- sum(test[1:13])
  if(is.na(totalpre)) {totalpre <- "NA"}
  
  # Finally: check if all subset editing sites are pre-edited (equal to 13)
  if(totalpre==13) {
    results[rowIdx, 6] <- "TRUE"
  } else {results[rowIdx, 6] <- "FALSE"}
  
  if(rowIdx==1) {print(noquote("Progress:"))}
  if(rowIdx %% 1000 == 0){
    prog <- round((rowIdx/nrow(results)*100), digits=2)
    print(noquote(c(prog,"%")))}
}

# Name the columns of the results table and save as csv file #
resultstest <- setNames(results, c("seq", "normcount", "junc_len", "junc_end", "ESdiff", "PregRNA1"))
write.table(resultstest, "RPS12REH1_WTOE_minTet_Rep2_ESS9_junctionanalysis.csv", sep=",", row.names=FALSE)

# Save results table as a unique object for each sample
OEUN_WT1_JA <- resultstest %>% mutate(kd="AvgUn", sample="UnWT1")
OEUN_WT2_JA <- resultstest %>% mutate(kd="AvgUn", sample="UnWT2")
OEUN_KA1_JA <- resultstest %>% mutate(kd="AvgUn", sample="UnKA1")
OEUN_KA2_JA <- resultstest %>% mutate(kd="AvgUn", sample="UnKA2")
OEUN_EQ1_JA <- resultstest %>% mutate(kd="AvgUn", sample="UnEQ1")
OEUN_EQ2_JA <- resultstest %>% mutate(kd="AvgUn", sample="UnEQ2")

OEIN_WT1_JA <- resultstest %>% mutate(kd="WTOE", sample="InWT1")
OEIN_WT2_JA <- resultstest %>% mutate(kd="WTOE", sample="InWT2")
OEIN_KA1_JA <- resultstest %>% mutate(kd="KAOE", sample="InKA1")
OEIN_KA2_JA <- resultstest %>% mutate(kd="KAOE", sample="InKA2")
OEIN_EQ1_JA <- resultstest %>% mutate(kd="EQOE", sample="InEQ1")
OEIN_EQ2_JA <- resultstest %>% mutate(kd="EQOE", sample="InEQ2")

# Combine tables
combined <- bind_rows(OEUN_WT1_JA, OEUN_WT2_JA, OEUN_KA1_JA, OEUN_KA2_JA, OEUN_EQ1_JA, OEUN_EQ2s_JA,
                      OEIN_WT1_JA, OEIN_WT2_JA, OEIN_KA1_JA, OEIN_KA2_JA, OEIN_EQ1_JA, OEIN_EQ2_JA)

# next... recreate graph in figure 6A
# save each output as a separate object
# join objects and calculate mean, sd for KO and WT
# join KO and WT, make ggplot

combinedtotalES9 <- combined %>% group_by(sample, kd) %>% 
  summarise(total=sum(normcount))

combinedtrue <- combined %>% filter(PregRNA1==TRUE) %>% group_by(sample, kd) %>%
  summarise(pretrue=sum(normcount))

combinedgraph <- full_join(combinedtotalES9,combinedtrue) %>%
  mutate(ratio=pretrue/total*100) %>% group_by(kd) %>%
  summarise(mean=mean(ratio),stdev=sd(ratio))

combinedgraph$kd <- factor(combinedgraph$kd, levels=c('AvgUn', 'WTOE', 'KAOE', 'EQOE'))

OEPalette <- c("#bdbdbd", "#118811", "#90BFF9", "#C000C0")

ggplot(combinedgraph, aes(x=kd, y=mean, fill=kd)) +
  geom_bar(aes(fill=kd), position="dodge", stat="identity") +
  ggtitle("RPS12 gRNA1 Region Pre-edited") +
  ylab("Percent of ESS9 Reads") + xlab("Sample") +
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

# Significance tests
WTtest1 <- combinedtotalES9 %>% filter(kd=="AvgUn" | kd=="WTOE")
WTtest2 <- WTtest1 %>% 
  group_by(kd, sample) %>% 
  summarise(total = list(total)) %>% 
  spread(kd, total) %>% ungroup() %>%
  group_by(kd) %>% 
  mutate(p_value = tryCatch(t.test(unlist(AvgUn), unlist(WTOE), var.equal = FALSE)$p.value,
                            error = function(err) return (NA)))
WTtest2 <- WTtest2[order(WTtest2$p_value),]
WTtest2$adj <- p.adjust(WTtest2$p_value, method = "BY")

KAtest1 <- combinedtotalES9 %>% filter(kd=="AvgUn" | kd=="KAOE")
KAtest2 <- KAtest1 %>% 
  group_by(kd, sample) %>% 
  summarise(total = list(total)) %>% 
  spread(kd, total) %>% ungroup() %>%
  group_by(kd) %>% 
  mutate(p_value = tryCatch(t.test(unlist(AvgUn), unlist(KAOE), var.equal = FALSE)$p.value,
                            error = function(err) return (NA)))
KAtest2 <- KAtest2[order(KAtest2$p_value),]
KAtest2$adj <- p.adjust(KAtest2$p_value, method = "BY")

EQtest1 <- combinedtotalES9 %>% filter(kd=="AvgUn" | kd=="EQOE")
EQtest2 <- EQtest1 %>% 
  group_by(kd, sample) %>% 
  summarise(total = list(total)) %>% 
  spread(kd, total) %>% ungroup() %>%
  group_by(kd) %>% 
  mutate(p_value = tryCatch(t.test(unlist(AvgUn), unlist(EQOE), var.equal = FALSE)$p.value,
                            error = function(err) return (NA)))
EQtest2 <- EQtest2[order(EQtest2$p_value),]
EQtest2$adj <- p.adjust(EQtest2$p_value, method = "BY")
