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
##---------------------------------------------------------------------------
## 
## This document contains all R code necessary to recapitulate analysis in this study
## additional code can be found in file "StatisticalProcessing.m" for determination of EPSs
## The majority of this code is the same regardless of what trascript one examines
## 
##
library(dplyr)
#install.packages("tidyr")
library(tidyr)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("reshape2")
library(reshape2)

####Graphing the EPSs#############################
## First combine all data sets to analyze
#######For RPS12###########
#First create a file that combines both the current data set (MRB10130 KD) and the previous uninduced samples for other REMC proteins
#RS Rep1 Uninduced samples
RPS12rep1All <- read.table(file.choose(), header = T, sep = ",")
#RS Rep2 Uninduced samples
RPS12rep2All <- read.table(file.choose(), header = T, sep = ",")
#NM MRB10130 all 4 samples
RPS12rep3All <- read.table(file.choose(), header = T, sep = ",")
RPS12cleanNMRS <-bind_rows(RPS12rep2All, RPS12rep1All, RPS12rep3All)
####Warning messages: this is ok!!!!
#1: In bind_rows_(x, .id) : Unequal factor levels: coercing to character
#2: In bind_rows_(x, .id) : Unequal factor levels: coercing to character
#3: In bind_rows_(x, .id) : Unequal factor levels: coercing to character

write.csv(RPS12cleanNMRS, file = "RPS12_10UNIN_10130_NM.csv")



#open the clean table "RPS12_10UNIN_10130_NM.csv" organized by KD and induction (this will include 2 MRB800 and induced & 10 uninduced samples)
RPS12allCl <- read.csv(file.choose(), header = TRUE, sep = ",")

#Need to average all ten uninduced samples
AvgUninducedRPS12 <- RPS12allCl %>% filter(tetracycline == "false") %>% group_by(edit_stop, gene) %>%
  summarize(Avg = sum(norm_count)/10) %>% mutate(KD = "AvgUn")

#If want to separate out KD by KD  
AvgInducedRPS12 <- RPS12allCl %>% filter(knock_down == "MRB10130"& tetracycline == "true") %>% group_by(edit_stop, gene) %>%
  summarize(Avg = sum(norm_count)/2) %>% mutate(KD = "MRB10130")

#Combine data 
TotaldataRPS12 <- bind_rows(AvgUninducedRPS12,AvgInducedRPS12)

#define which sites you want to look at
EPS10130 <- c(12,21,28,29,30,31,45,46,48,50,51,53,60,61)

TotaldataRPS12b <- TotaldataRPS12 %>% filter(edit_stop%in%EPS10130)

allbar <- ggplot(TotaldataRPS12b, aes(x=KD, y=Avg)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~edit_stop, scales = "free") + scale_fill_manual(values=c("black"), name="Junction Length", labels=c("0","1-10","11-50","50+")) +
  theme(panel.margin = unit(1, "lines"))

#Final graph
allbar

#Other variations that could be used 
allbar8 <- ggplot(TotaldataRPS12b, aes(x=KD, y=Avg)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~gene + edit_stop, nrow=2, scales = "free") + scale_fill_manual(values=c("black"), name="Junction Length", labels=c("0","1-10","11-50","50+")) +
  theme(panel.margin = unit(1, "lines"))

allbar5 <- ggplot(TotaldataRPS12b, aes(x=KD, y=Avg)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~gene + edit_stop, nrow=5, scales = "free") + scale_fill_manual(values=c("black"), name="Junction Length", labels=c("0","1-10","11-50","50+")) +
  theme(panel.margin = unit(1, "lines"))

allbar5 <- ggplot(TotaldataRPS12b, aes(x=KD, y=Avg)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~edit_stop, nrow=5, scales = "free") + scale_fill_manual(values=c("black"), name="Junction Length", labels=c("0","1-10","11-50","50+")) +
  theme(panel.margin = unit(1, "lines"))


############Top Junction Sequences gRNA 1-2 ####RPS12##########################
library(dplyr)
#install.packages("tidyr")
library(tidyr)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("reshape2")
library(reshape2)

# open the clean table organized by KD and induction (this will include 2 MRB10130 induced & 10 uninduced samples)
RPS12allCl <- read.csv(file.choose(), header = TRUE, sep = ",")

#to look at the first and second gRNAs - edit_stop<40,edit_stop>10) and remove pre-edited
RPS12allg2 <- RPS12allCl %>% filter(edit_stop<40,edit_stop>10,!(edit_stop==9&junc_len==0))

#Average all the 10 uninduced samples

RPS12allg1un <- RPS12allg2 %>% filter(tetracycline == "false") %>% group_by(edit_stop, junc_len, junc_seq) %>%
  summarize(Avg = sum(norm_count)/10) %>% mutate(KD = "AvgUn")

RPS12allg1un <- RPS12allg2 %>% filter(tetracycline=="false") %>% group_by(edit_stop, junc_len, junc_seq) %>%
  mutate(AvgNm = sum(norm_count)/10) %>% select(edit_stop, junc_seq, junc_len, tetracycline, AvgNm) %>%
  mutate(KD="AvgUn") %>% distinct(., .keep_all = TRUE)
glimpse(RPS12allg1un)
head(RPS12allg1un)

# make an equivalent table with the top junction sequences
# found in the KDs, taking the average of the norm count across replicates

RPS12allg1in <- RPS12allg2 %>% filter(tetracycline=="true") %>% group_by(edit_stop, junc_len, junc_seq) %>%
  mutate(AvgNm = sum(norm_count)/2) %>% select(edit_stop, junc_seq, junc_len, tetracycline, AvgNm) %>%
  distinct(., .keep_all = TRUE)

write.table(RPS12allg1un, "RPS12gRNA1_2Unin_itopseq.csv", sep= ",")
write.table(RPS12allg1in, "RPS12gRNA1_2Induc_itopseq.csv", sep= ",")


# combine the tables together to make graphing easier

#Contains all average junction sequences for Unind and KD
RPS12allg1avg <- bind_rows(RPS12allg1un, RPS12allg1in) %>% group_by(KD) %>% 
  mutate(total = sum(AvgNm)) %>% rowwise() %>% mutate(perc = 100*(AvgNm/total))
View(RPS12allg1avg)

# put the top sequences by KD in a chart compare side by side (>100 avg copies)
RPS12g1top100 <- RPS12allg1avg %>% filter(AvgNm > 50) %>% arrange(KD, desc(AvgNm))
glimpse(RPS12g1top100)
View(RPS12g1top100)
summary(RPS12g1top100)

write.table(RPS12g1top100, "RPS12gRNA1_2topSeq_MRB10130.csv", sep= ",")

