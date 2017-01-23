## Copyright 2015 TREAT Authors. All rights reserved.
##
## This file is part of TREAT.
##
## TREAT is free software: you can redistribute it and#or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## TREAT is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with TREAT.  If not, see <http:##www.gnu.org#licenses#>.
##
##-------------------------------------------------------------------------------------------
##
## this document contains all R code necessary to recapitulate analyses in this study
## additional code can be found in file "StatisticalProcessing.m" for determination of EPSs
##  The majority of this code is the same regardless of what transcript one examines
##  when this is the case, RPS12 only is shown to avoid unnecessary duplication of code
##  The same code will work with ND7-5' data.

##  These packages are used in the following code
## if you haven't done it already, install packages and load them for data wrangling
#install.packages("stringr")
library("stringr")

#install.packages("tidyr")
library("tidyr")

#install.packages("dplyr")
library("dplyr")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("gdata")
library(gdata)

#installed.packages("grid")
library(grid)

#install.packages("stringr")
library(stringr)


##  Table "RPS12allCl" is the full database as taken from the "Search" Page on TREAT 
##  using the "export" function

## preparation of table for determination of significant changes (EPS, EJE, junction 0, etc.)

RPS12all <- RPS12allCl %>% filter(!(junc_length==0&edit_stop=9)) 
#removal of pre-edited
# other parameters can be set here (include all junc_length==0)
# can be used for other transcipts.                                                              # 
glimpse(RPS12all)

RPS12allR1un <- RPS12all %>% filter(replicate==TRUE, tet==FALSE) %>%
  select(sample, norm_count, edit_stop) %>% spread(sample, norm_count)
glimpse(RPS12allR1un)

RPS12allR1in <- RPS12all %>% filter(replicate==TRUE, tet==TRUE) %>%
  select(sample, norm_count, edit_stop) %>% spread(sample, norm_count)
RPS12allR2un <- RPS12all %>% filter(replicate==FALSE, tet==FALSE) %>% 
  select(sample, norm_count, edit_stop) %>% spread(sample, norm_count)
RPS12allR2in <- RPS12all %>% filter(replicate==FALSE, tet==TRUE) %>%
  select(sample, norm_count, edit_stop) %>% spread(sample, norm_count)

RPS12R1 <- left_join(RPS12allR1un, RPS12allR1in)
glimpse(RPS12R1)

RPS12R2 <- left_join(RPS12allR2un, RPS12allR2in)
glimpse(RPS12R2)

## these files can be saved and used in additional code developed by RC to determine p and q values
## comparison of which p and q values were significant across multiple replicates was performed manually
#########################################################################################################
#########################################################################################################

## correlation of EPSs with specific editing characteristics
## this code can be adjusted for other Editing Sites of interest (EJE, JSS etc)

##determining the overlap of sets with editing characteristics
##use to get the information for Fisher's Exact Test
##need a table with each variable to test correlation with as a column
##containing y or n for variable's presence.

##Input:  Characteristics table and column of ES subset (IPS, MJEs etc)

##output:  need to have the count of TRUE, FALSE categories for overall (ES 
##         character column) and for overlap with subset in question 
##         (e.g. IPSs)

# import editing site characteristic table (
# this table was generated manually and can be made for other transcripts
# see additional file "RPS12ESCharacteristics.csv" and "ND75ESCharacteristics.csv"
RPS12esch <- read.table(file.choose(), header=T, sep = ",")

##import data from RC code that has TRUE and FALSE for EPSs present in 
##all replicates -- see additional files "RPS12epsR1andR2.csv" and "ND75epsR1andR2.csv"

ESSall <- read.table(file.choose(), header=T, sep = ",")

##add the columns of the relevant characteristic chart 
epsA <- inner_join(ESSall, RPS12esch)

##gather data from RC
epsB <- epsA %>% gather(name, eps, 2:5)
ESset <- epsB
glimpse(ESset)

##getting all the pvalues

T2 <- ESset %>% rowwise() %>% mutate(BothAdd = isTRUE(Add&eps)) ##if TRUE in set and characteristic
T3 <- T2 %>% rowwise() %>% mutate(BothDel = isTRUE(Del&eps)) %>% mutate(BothNoAxn = isTRUE(NoAxn&eps)) %>%
  mutate(Both3A = isTRUE(X3A&eps)) %>% mutate(Both5A = isTRUE(X5A&eps)) %>% 
  mutate(Both3G = isTRUE(X3G&eps)) %>% mutate(Both5G = isTRUE(X5G&eps)) %>% mutate(Both3C = isTRUE(X3C&eps)) %>%
  mutate(Both5C = isTRUE(X5C&eps))
F1t <- T3 %>% group_by(name) %>% summarise_each(funs(sum(. ==TRUE)))
F1bf <- F1t %>% rowwise() %>% mutate(bAddf = eps-BothAdd) %>% mutate(bDelf = eps-BothDel) %>% mutate(bNoAxnf = eps-BothNoAxn) %>%
  mutate(b3Af = eps-Both3A) %>% mutate(b5Af = eps-Both5A) %>% mutate(b3Gf = eps-Both3G) %>% mutate(b5Gf = eps-Both5G) %>% 
  mutate(b3Cf = eps-Both3C) %>% mutate(b5Cf = eps-Both5C)
#count TRUE and FALSE of characteristic
F1f <- T3 %>% group_by(name) %>% summarise_each(funs(sum(. ==FALSE))) ##do something to subtract the eps count from true both count... 
F2ct <- F1t %>% select(name, 3:11) %>% gather(char, chTRUEct, -name) 
F2bt <- F1t %>% select(name, 13:21) %>% gather(both, bTRUEct, -name)
F2bf <- F1bf %>% select(name, 22:30) %>% gather(both, bFALSEct, -name)
F2cf <- F1f %>% select(name, 3:11) %>% gather(char, chFALSEct, -name) 
F3ct <- F2ct %>% select(chTRUEct)
F3bt <- F2bt %>% select(name, both, bTRUEct)
F3bf <- F2bf %>% select(bFALSEct)
F3cf <- F2cf %>% select(chFALSEct)
F3 <- bind_cols(F3bt, F3ct, F3bf, F3cf) 
F3u <- as.matrix(F3[,3:6])
FishTest <- apply(F3u,1, function(x) fisher.test(matrix(x, nrow = 2))$p.value)
Fsum <- bind_cols(F3, as.data.frame(FishTest))                                              
glimpse(Fsum)

# can write output table for future reference 
write.table(Fsum, file = "FishTestEPSr1r2.csv", sep = ",")

##################################################################################################
###################################################################################################

## Determining the action at JSS that correlate with TbRGG2 and MRB8180 RNAi caused EPSs
## limit to the EPSs requiring deletions at the next (5') site and remove pre-edited
EPSpreDel <- c(9, 19, 30, 31, 53, 78, 106)
RPreDel <-  RPS12allCl %>% filter(edit_stop %in% EPSpreDel, !(edit_stop==9 & junc_len==0)) 
summary(RPreDel)

## We will now check if the correct deletion has been performed 
## to do this, we will grab the right most "n" characters 
## using the stringr package formula str_sub(col, start = , end = )
## We'll make a column of the substring
## then we can use this column to check the prop'n with no del, add, part
## for this, I will substring the length of the U's deleted plus one nucl
## given size and locations of deletions in RPS12, -3, -4 and -5 should be good for all
#group sites I can test in a meaningful way to simplify the analysis
ESS9noDel <- c("", "TTTT", "ATTTT")
ESS9add <- c("TTTTT")
T9sites <- RPreDel %>% filter(edit_stop==9) %>% 
  rowwise() %>% mutate(delstr = str_sub(junc_seq, -5)) %>%
  mutate(nodel = ifelse(delstr %in% ESS9noDel, "ND", ifelse(delstr=="TTTTT", "ADD","PD")))
glimpse(T9sites)

# check to see if the ESS9 junctions have over deletion (would show up as PD here)
chk9PD <- T9sites %>% filter(nodel == "PD")
glimpse(chk9PD)
View(chk9PD)

ESS30_78noDel <- c("", "TTTTG", "GTTTTG", "CTTTTG")
ESS30_78add <- c("TTTTTG")
#make last five with this subset
TTTTsites <- RPreDel %>% filter(edit_stop %in% c(30, 78)) %>% 
  rowwise() %>% mutate(delstr = str_sub(junc_seq, -6)) %>%
  mutate(nodel = ifelse(delstr %in% ESS30_78noDel, "ND", ifelse(delstr=="TTTTTG", "ADD","PD")))
glimpse(TTTTsites)
chk <- TTTTsites %>% filter(nodel == "ND" & edit_stop==78)
summary(chk)

ESS31_53noDel <- c("","TTTG", "TTTC", "CTTTG", "GTTTC")
ESS31_53add <- c("TTTTG", "TTTTC")
T3sites <- RPreDel %>% filter(edit_stop==31|edit_stop==53) %>% 
  rowwise() %>% mutate(delstr = str_sub(junc_seq, -5)) %>%
  mutate(nodel = ifelse(delstr %in% ESS31_53noDel, "ND", ifelse(delstr%in%ESS31_53add, "ADD","PD")))
glimpse(T3sites)

#ESS19 and 106 require small ins.-- cut down to -2 (not -4)
# because "T" could be a partial deletion of 106 and "TT" an add for 19
# we will do these separately.
ESS106noDel <- c("", "TTC","GTTC")
ESS106add <- c("TTTC")
ESS19noDel <- c("","ATC", "TC") 
ESS19add <- c("TTC")
T2sites <- RPreDel %>% filter(edit_stop==106) %>% 
  rowwise() %>% mutate(delstr = str_sub(junc_seq, -4)) %>%
  mutate(nodel = ifelse(delstr %in% ESS106noDel, "ND", ifelse(delstr=="TTTC", "ADD","PD")))
glimpse(T2sites)

T1sites <- RPreDel %>% filter(edit_stop==19) %>% 
  rowwise() %>% mutate(delstr = str_sub(junc_seq, -3)) %>%
  mutate(nodel = ifelse(delstr %in% ESS19noDel, "ND", ifelse(delstr=="TTC", "ADD","PD")))
glimpse(T1sites)
T1chk <- T1sites %>% filter(nodel == "PD")
glimpse(T1chk)

# calculate the proportion of "no del" junctions at each ESS
# first combine the tables together
# because we have all same column names, we can use bind_rows%%%%
NoDelAll <- bind_rows(TTTTsites, T9sites, T3sites, T2sites, T1sites) %>% filter(KD!="29-13") %>%
  arrange(nodel)
glimpse(NoDelAll)

## start by just vizualizing the change in proportion across the different groups...
DelSum <- NoDelAll %>% group_by(KD, TET, edit_stop, nodel) %>% 
  summarize(sum = sum(norm_count))
View(DelSum)

DelQuickViz <- ggplot(DelSum, aes(x=TET, y=sum, fill = nodel)) + 
  geom_bar(stat = "identity", position = "fill") + facet_grid(edit_stop ~ KD, scales = "free") +
  labs(title = "Rep1 and Rep2") + scale_fill_brewer(palette="Spectral")
DelQuickViz

# check that trend is same across both replicates
R1R2Sum <- NoDelAll %>% group_by(REP1, KD, TET, edit_stop, nodel) %>% 
  summarize(sum = sum(norm_count))

R1QuickViz <- R1R2Sum %>% filter(REP1==TRUE) %>% ggplot(., aes(x=TET, y=sum, fill = nodel)) + 
  geom_bar(stat = "identity", position = "fill") + facet_grid(edit_stop ~ KD, scales = "free") +
  labs(title = "Rep1") + scale_fill_brewer(palette="Spectral")
R1QuickViz

R2QuickViz <- R1R2Sum %>% filter(REP1==FALSE) %>% ggplot(., aes(x=TET, y=sum, fill = nodel)) + 
  geom_bar(stat = "identity", position = "fill") + facet_grid(edit_stop ~ KD, scales = "free") +
  labs(title = "Rep2") + scale_fill_brewer(palette="Spectral")
R2QuickViz

###########################################################################################
###########################################################################################

## for analyses of exact junction sequences, the average of the norm_count of 
## each unique sequence is taken across all relevant replicates
##  that is, across all 8 uninduced samples and between both induced samples for each RNAi
## this calculation is incorporated into subsequent code

## Determination of number of sequences of each junction length 
## across all editing stop sites (Fig. 5a etc)

RPS12allJLun <- RPS12allCl %>% filter(TET==FALSE, !(edit_stop==9&junc_len==0), KD!="29-13") %>% 
  group_by(junc_len) %>% 
  summarise(JLtot = sum(norm_count)/8) %>% mutate(KD="AvgUn")
glimpse(RPS12allJLun)

RPS12allJLin <- RPS12allCl %>% filter(TET==TRUE, !(edit_stop==9&junc_len==0)) %>% group_by(KD, junc_len) %>% 
  summarise(JLtot = sum(norm_count)/2)

RPS12allJL <- bind_rows(RPS12allJLun, RPS12allJLin)
target = c("AvgUn", "GAP1", "81704160", "81802c", "TbRGG2")
RPS12allJL$KD <- reorder.factor(RPS12allJL$KD, new.order = target)

glimpse(RPS12allJL)

RPS12allshort <- RPS12allJL %>% filter(junc_len<6)
Tb2Palette <- c("#009E73", "#CC79A7", "#0072B2", "#56B4E9", "#E69F00")
JLplot <- ggplot(RPS12allshort, aes(x=junc_len, y=JLtot, fill=KD)) + geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values =  Tb2Palette)
JLplot

########################################################
## Pre-edited graphs were made in a similar fashion
RPS12allJLunPE <- RPS12allCl %>% filter(TET==FALSE, (edit_stop==9&junc_len==0), KD!="29-13") %>% 
  group_by(junc_len) %>% 
  summarise(JLtot = sum(norm_count)/8) %>% mutate(KD="AvgUn")
glimpse(RPS12allJLunPE)

RPS12allJLinPE <- RPS12allCl %>% filter(TET==TRUE, (edit_stop==9&junc_len==0)) %>% group_by(KD, junc_len) %>% 
  summarise(JLtot = sum(norm_count)/2)

RPS12allPE <- bind_rows(RPS12allJLunPE, RPS12allJLinPE)
target = c("AvgUn", "GAP1", "81704160", "81802c", "TbRGG2")
RPS12allPE$KD <- reorder.factor(RPS12allPE$KD, new.order = target)

glimpse(RPS12allPE)

Tb2Palette <- c("#009E73", "#CC79A7", "#0072B2", "#56B4E9", "#E69F00")
PEplot <- ggplot(RPS12allPE, aes(x=KD, y=JLtot, fill=KD)) + geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values =  Tb2Palette)
PEplot

###################################################################################
###################################################################################
## plot proportion of junction lengths at each ESS
##  This section of code was written by Brianna L Tylec (BLT)
# introduce junction length bins
RPS12allCl$bin1 <- cut(RPS12allCl$junc_len, breaks = c(-Inf,1,11,51,Inf), right = FALSE)
levels <- levels(RPS12allCl$bin1)

WT <- RPS12allCl %>% filter(KD == "29-13") %>% group_by(edit_stop, bin1) %>%
  summarize(Avg = sum(norm_count)/2) %>% mutate(KD = "29-13")

AvgUninduced <- RPS12allCl %>% filter(KD !="29-13" & TET == "FALSE") %>% group_by(edit_stop, bin1) %>%
  summarise(Avg = sum(norm_count)/8) %>% mutate(KD = "AvgUn")

AvgInduced <- RPS12allCl %>% filter(TET == "TRUE") %>% group_by(edit_stop, KD, bin1) %>%
  summarise(Avg = sum(norm_count)/2)

Totaldata <- bind_rows(WT,AvgUninduced,AvgInduced)
Totaldata$KD <- reorder.factor(Totaldata$KD, new.order = c("29-13","AvgUn","GAP1","81704160","81802c","TbRGG2"))
Totaldata1 <- Totaldata %>% arrange(KD)

GAP1EPS <- data.frame(edit_stop=c(22,23,24,40,78,81,106,119), Avg=1.05,
                      KD="GAP1")
MRB8180EPS <- data.frame(edit_stop=c(14,15,16,19,25,26,27,29,30,31,35,36,37,53,78,80),
                         Avg=1.05, KD="81802c")
TbRGG2EPS <- data.frame(edit_stop=c(12,15,16,19,30,35,36,37,53,78,83,100),
                        Avg=1.05, KD="TbRGG2")
MRB8170EPS <- data.frame(edit_stop=c(9,12,135), Avg=1.05, KD="81704160")

ggplot() + 
  ylab("Percentage") + xlab("Editing Site") + ggtitle("RPS12 Junction Lengths by Editing Site") +
  geom_bar(data=Totaldata1, aes(x = edit_stop, y = Avg, fill=bin1),stat="identity", position="fill") +
  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) + 
  scale_x_reverse(breaks=c(150,140,130,120,110,100,90,80,70,60,50,40,30,20,10)) +
  facet_grid(KD~.) +
  geom_point(data=GAP1EPS, aes(x = edit_stop, y = Avg), pch=21) +  
  geom_point(data=MRB8180EPS, aes(x = edit_stop, y = Avg), pch=21) +
  geom_point(data=TbRGG2EPS, aes(x = edit_stop, y = Avg), pch=21) +
  geom_point(data=MRB8170EPS, aes(x = edit_stop, y = Avg), pch=21)

## Remove junction length 0 from first ES

Totaldata2 <- Totaldata1 %>% filter(!(edit_stop==9 & bin1=="[-Inf,1)"))

ggplot() + 
  ylab("Percentage") + xlab("Editing Site") + ggtitle("RPS12 Junction Lengths by Editing Site") +
  geom_bar(data=Totaldata2, aes(x = edit_stop, y = Avg, fill=bin1),stat="identity", position="fill") +
  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) + 
  scale_x_reverse(breaks=c(150,140,130,120,110,100,90,80,70,60,50,40,30,20,10)) +
  facet_grid(KD~.) +
  geom_point(data=GAP1EPS, aes(x = edit_stop, y = Avg), pch=21) +  
  geom_point(data=MRB8180EPS, aes(x = edit_stop, y = Avg), pch=21) +
  geom_point(data=TbRGG2EPS, aes(x = edit_stop, y = Avg), pch=21) +
  geom_point(data=MRB8170EPS, aes(x = edit_stop, y = Avg), pch=21)

## Plot just the region with first 2 guide RNAs

Totaldata3 <- Totaldata2 %>% filter(edit_stop <= 45)

GAP1EPS2 <- data.frame(edit_stop=c(22,23,24,40), Avg=1.05, KD="GAP1")
MRB8180EPS2 <- data.frame(edit_stop=c(14,15,16,19,25,26,27,29,30,31,35,36,37), Avg=1.05, KD="81802c")
TbRGG2EPS2 <- data.frame(edit_stop=c(12,15,16,19,30,35,36,37), Avg=1.05, KD="TbRGG2")
MRB8170EPS2 <- data.frame(edit_stop=c(9,12), Avg=1.05, KD="81704160")

ggplot() + 
  ylab("Percentage") + xlab("Editing Site") + ggtitle("RPS12 Junction Lengths by Editing Site") +
  geom_bar(data=Totaldata3, aes(x = edit_stop, y = Avg, fill=bin1),stat="identity", position="fill") +
  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) + 
  scale_x_reverse(breaks=c(45,40,35,30,25,20,15,10,5)) +
  facet_grid(KD~.) +
  geom_point(data=GAP1EPS2, aes(x = edit_stop, y = Avg), pch=21) +  
  geom_point(data=MRB8180EPS2, aes(x = edit_stop, y = Avg), pch=21) +
  geom_point(data=TbRGG2EPS2, aes(x = edit_stop, y = Avg), pch=21) +
  geom_point(data=MRB8170EPS2, aes(x = edit_stop, y = Avg), pch=21)

# ND7 ---------------------------------------------------------------------

ND7rep1All <- read.table(file.choose(), header = T, sep = ",")
ND7rep2All <- read.table(file.choose(), header = T, sep = ",")
ND7all <-bind_rows(ND7rep2All, ND7rep1All)

KD<-c("29-13","81704160","81802c","GAP1","TbRGG2")
plusTET <-c("+")
dateRep1 <-c("6-3")
table<-ND7all
samples <-ND7all$sample
parameter <- ND7all$junc_len

ND7allCl <- CleanTREAT(ND7all, ND7all$sample, ND7all$junc_len, KD, plusTET, dateRep1)

ND7allCl$bin1 <- cut(ND7allCl$junc_len, breaks = c(-Inf,1,11,51,Inf), right = FALSE)
levels <- levels(ND7allCl$bin1)

WT <- ND7allCl %>% filter(KD == "29-13") %>% group_by(edit_stop, bin1) %>%
  summarize(Avg = sum(norm_count)/2) %>% mutate(KD = "29-13")

AvgUninduced <- ND7allCl %>% filter(KD !="29-13" & TET == "FALSE") %>% group_by(edit_stop, bin1) %>%
  summarise(Avg = sum(norm_count)/8) %>% mutate(KD = "AvgUn")

AvgInduced <- ND7allCl %>% filter(TET == "TRUE") %>% group_by(edit_stop, KD, bin1) %>%
  summarise(Avg = sum(norm_count)/2)

ND7Totaldata <- bind_rows(WT,AvgUninduced,AvgInduced)
ND7Totaldata$KD <- factor(ND7Totaldata$KD, levels = c("29-13","AvgUn","GAP1","81704160","81802c","TbRGG2"))
ND7Totaldata1 <- ND7Totaldata %>% arrange(KD)

ND7GAP1EPS <- data.frame(edit_stop=c(33,34), Avg=1.05, KD="GAP1")
ND7MRB8180EPS <- data.frame(edit_stop=c(30,33), Avg=1.05, KD="81802c")
ND7TbRGG2EPS <- data.frame(edit_stop=c(22,23,24,30,33), Avg=1.05, KD="TbRGG2")
ND7MRB8170EPS <- data.frame(edit_stop=c(17,36,37,38,46,47,49,50,52,55,60,61,64,67), 
                            Avg=1.05, KD="81704160")

ggplot() + 
  ylab("Percentage") + xlab("Editing Site") + ggtitle("ND7 Junction Lengths by Editing Site") +
  geom_bar(data=ND7Totaldata1, aes(x = edit_stop, y = Avg, fill=bin1),stat="identity", position="fill") +
  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) + 
  scale_x_reverse(breaks=c(100,90,80,70,60,50,40,30,20,10)) +
  facet_grid(KD~.) +
  geom_point(data=ND7GAP1EPS, aes(x = edit_stop, y = Avg), pch=21) +  
  geom_point(data=ND7MRB8180EPS, aes(x = edit_stop, y = Avg), pch=21) +
  geom_point(data=ND7TbRGG2EPS, aes(x = edit_stop, y = Avg), pch=21) +
  geom_point(data=ND7MRB8170EPS, aes(x = edit_stop, y = Avg), pch=21)

## Remove 0 length junctions from ES21...

ND7Totaldata2 <- ND7Totaldata1 %>% filter(!(edit_stop==21 & bin1=="[-Inf,1)"))

ggplot() + 
  ylab("Percentage") + xlab("Editing Site") + ggtitle("ND7 Junction Lengths by Editing Site") +
  geom_bar(data=ND7Totaldata2, aes(x = edit_stop, y = Avg, fill=bin1),stat="identity", position="fill") +
  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) + 
  scale_x_reverse(breaks=c(100,90,80,70,60,50,40,30,20,10)) +
  facet_grid(KD~.) +
  geom_point(data=ND7GAP1EPS, aes(x = edit_stop, y = Avg), pch=21) +  
  geom_point(data=ND7MRB8180EPS, aes(x = edit_stop, y = Avg), pch=21) +
  geom_point(data=ND7TbRGG2EPS, aes(x = edit_stop, y = Avg), pch=21) +
  geom_point(data=ND7MRB8170EPS, aes(x = edit_stop, y = Avg), pch=21)

# create graph for first two gRNAs

ND7Totaldata3 <- ND7Totaldata2 %>% filter(edit_stop <= 45)

ND7GAP1EPS2 <- data.frame(edit_stop=c(33,34), Avg=1.05, KD="GAP1")
ND7MRB8180EPS2 <- data.frame(edit_stop=c(30,33), Avg=1.05, KD="81802c")
ND7TbRGG2EPS2 <- data.frame(edit_stop=c(22,23,24,30,33), Avg=1.05, KD="TbRGG2")
ND7MRB8170EPS2 <- data.frame(edit_stop=c(17,36,37,38), Avg=1.05, KD="81704160")

ggplot() + 
  ylab("Percentage") + xlab("Editing Site") + ggtitle("ND7 Junction Lengths by Editing Site") +
  geom_bar(data=ND7Totaldata3, aes(x = edit_stop, y = Avg, fill=bin1),stat="identity", position="fill") +
  scale_fill_manual(values=cbbPalette, name="Junction Length", labels=labels) + 
  scale_x_reverse(breaks=c(45,40,35,30,25,20,15,10,5)) +
  facet_grid(KD~.) +
  geom_point(data=ND7GAP1EPS2, aes(x = edit_stop, y = Avg), pch=21) +  
  geom_point(data=ND7MRB8180EPS2, aes(x = edit_stop, y = Avg), pch=21) +
  geom_point(data=ND7TbRGG2EPS2, aes(x = edit_stop, y = Avg), pch=21) +
  geom_point(data=ND7MRB8170EPS2, aes(x = edit_stop, y = Avg), pch=21)

################################################################################################
################################################################################################

## To compare both abundance and proportion, do a stacked bar graph showing both
##  Fig. S2
target <- c("AvgUn","GAP1","81704160","81802c","TbRGG2")

  RPS12allCl$Bin1 <- cut(RPS12allCl$junc_len, breaks = c(-Inf,1,11,51,Inf), right = FALSE)
  JLsBinAvgUn <- RPS12allCl %>% filter(TET == "FALSE") %>% group_by(edit_stop, Bin1, TET) %>% summarise(Avg = sum(norm_count)/8) %>%
    mutate(KD = "AvgUn") %>% select(edit_stop, KD, Bin1, TET, Avg)
  JLsBinAvgIn <- RPS12allCl %>% filter(TET == "TRUE") %>% group_by(edit_stop, KD, Bin1, TET) %>% summarise(Avg = sum(norm_count)/2)
  JLsBinAvg <- bind_rows(JLsBinAvgIn, JLsBinAvgUn)
  JLsBinAvg$KD <- reorder.factor(JLsBinAvg$KD, new.order = target)
  JLsBinAvg1 <- JLsBinAvg %>% arrange(KD)

# limit only to all EPS's
RPS12allEPS <- c(22,23,24,40,78,81,106,119,14,15,16,19,25,26,27,29,30,31,35,36,37,53,78,80,83,100,9, 135)
# for ND7-5' the following EPS string is used
# ND75allEPS <- c(17,36,37,38,46,47,49,50,52,55,60,61,64,67, 30, 33,22,23,24, 34)

JLsBinAvg1EPS <- JLsBinAvg1 %>% filter(edit_stop %in% RPS12allEPS)
  
glimpse(JLsBinAvg1)

###make the bar graphs 
##color blind friendly palette
cbbPalette <- c("#000000", "#D55E00", "#F0E442", "#009E73", "#E69F00", "#0072B2", "#56B4E9", "#CC79A7")


##this is the favorite form
bars <- ggplot(JLsBinAvg1EPS, aes(x=KD, y=Avg, fill = Bin1)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~edit_stop, scales = "free") + scale_fill_manual(values=cbbPalette) 
bars

#################################################################################################
#################################################################################################
##  Determination of major junction sequences in gRNA-1

# limit to sequences whose edit_stop is within the first gRNA,
# remove 29-13 (if applicable), and remove pre-edited
head(RPS12allCl)
RPS12allg1 <- RPS12allCl %>% filter(edit_stop<23,KD!="29-13",!(edit_stop==9&junc_len==0))
summary(RPS12allg1)
check <- RPS12allg1 %>% filter(junc_len==0)
head(check)

# this is the data set of junction sequences which start in gRNA-1
# next, merge replicates 1 and 2 by creating a new avgnorm for all seq
# take uninduced together and induced together

RPS12allg1un <- RPS12allg1 %>% filter(TET==FALSE) %>% group_by(edit_stop, junc_len, junc_seq) %>%
  mutate(AvgNm = sum(norm_count)/8) %>% select(edit_stop, junc_seq, junc_len, TET, AvgNm) %>%
  mutate(KD="AvgUn") %>% distinct(., .keep_all = TRUE)
glimpse(RPS12allg1un)
head(RPS12allg1un)
chkcounts <- RPS12allg1un %>% group_by(KD) %>% summarise(sum = sum(AvgNm))
comparelen <- RPS12allg1 %>% group_by(KD, REP1, TET) %>% summarise(len = length(id))
comparelen
chkcounts

# make an equivalent table with the junction sequences
# found in the KDs, taking the average of the norm count across replicates

RPS12allg1in <- RPS12allg1 %>% filter(TET==TRUE) %>% group_by(KD, edit_stop, junc_len, junc_seq) %>%
  mutate(AvgNm = sum(norm_count)/2) %>% select(edit_stop, junc_seq, junc_len, TET, AvgNm, KD) %>%
  distinct(., .keep_all = TRUE)

head(RPS12allg1in)

# combine the tables together and add percentage

RPS12allg1avg <- bind_rows(RPS12allg1un, RPS12allg1in) %>% group_by(KD) %>% 
  mutate(total = sum(AvgNm)) %>% rowwise() %>% mutate(perc = 100*(AvgNm/total))
glimpse(RPS12allg1avg)

# put the top sequences by KD in a chart compare side by side (>100 avg copies)
RPS12g1top100 <- RPS12allg1avg %>% filter(AvgNm > 100) %>% arrange(KD, desc(AvgNm))
glimpse(RPS12g1top100)
View(RPS12g1top100)
summary(RPS12g1top100)

# make a hierarchical tree of order of action past ESS15
# first determine the number of sequences with specific actions

# each single action of editing after ESS15 can be calculated 
# start with junction length of 1, this represents action at ESS16 after 15
# junction length 2 can be first (if no action at 16) or second (if action at 16)
# junc_len 3 can be 1st (no other action), 2nd or third action
# etc.

g1ess15to16 <- RPS12allg1avg %>% filter((edit_stop==16&junc_len==0) || (edit_stop==15&junc_len==1)) %>%
  group_by(KD) %>% summarise(axn16 = sum(AvgNm), unq16 = length(junc_seq))
g1ess15to16

# for 17 alone, need to check if it is the only sequence with action, look for seq with intact pre at ES16
g1ess15to17 <- RPS12allg1avg %>% filter(edit_stop==15&junc_len==2) %>% mutate(ES16 = str_sub(junc_seq, -2)) %>%
  filter(ES16=="GA") %>% group_by(KD) %>% summarise(axn17 = sum(AvgNm), unq17 = length(junc_seq))
g1ess15to17

# use similar logic for 18 alone, the seq if ES16 and ES17 are skipped would be "AGA"
g1ess15to18 <- RPS12allg1avg %>% filter(edit_stop==15, junc_len==3) %>% mutate(ES16_17 = str_sub(junc_seq, -3)) %>%
  filter(ES16_17=="AGA") %>% group_by(KD) %>% summarise(axn18 = sum(AvgNm), unq18 = length(junc_seq))
View(g1ess15to18)

# analagous logic for 19 and up
g1ess15to19 <- RPS12allg1avg %>% filter(edit_stop==15, junc_len==4) %>% mutate(ES16_18 = str_sub(junc_seq, -7)) %>%
  filter(ES16_18=="GTTTAGA") %>% group_by(KD) %>% summarise(axn19 = sum(AvgNm), unq19 = length(junc_seq))
g1ess15to19

g1ess15to20 <- RPS12allg1avg %>% filter(edit_stop==15, junc_len==5) %>% mutate(ES16_19 = str_sub(junc_seq, -9)) %>%
  filter(ES16_19=="CGTTTAGA") %>% group_by(KD) %>% summarise(axn20 = sum(AvgNm), unq20 = length(junc_seq))
g1ess15to20

g1ess15to21 <- RPS12allg1avg %>% filter(edit_stop==15, junc_len==6) %>% mutate(ES16_20 = str_sub(junc_seq, -11)) %>%
  filter(ES16_20=="ATCGTTTAGA") %>% group_by(KD) %>% summarise(axn21 = sum(AvgNm), unq21 = length(junc_seq))
g1ess15to21

g1ess15to22 <- RPS12allg1avg %>% filter(edit_stop==15, junc_len==7) %>% mutate(ES16_21 = str_sub(junc_seq, -12)) %>%
  filter(ES16_21=="CATCGTTTAGA") %>% group_by(KD) %>% summarise(axn22 = sum(AvgNm), unq22 = length(junc_seq))
g1ess15to22

g1ess15to23 <- RPS12allg1avg %>% filter(edit_stop==15, junc_len==8) %>% 
  filter(grepl('ACATCGTTTAGA', junc_seq)) %>% group_by(KD) %>% 
  summarise(axn23 = sum(AvgNm), unq23 = length(junc_seq))
g1ess15to23

g1ess15to24 <- RPS12allg1avg %>% filter(edit_stop==15, junc_len==9) %>% 
  filter(grepl('AACATCGTTTAGA', junc_seq)) %>% group_by(KD) %>% summarise(axn24 = sum(AvgNm), unq24 = length(junc_seq))
g1ess15to24

g1ess15to25 <- RPS12allg1avg %>% filter(edit_stop==15, junc_len==10) %>% 
  filter(grepl('AAACATCGTTTAGA', junc_seq)) %>% group_by(KD) %>% summarise(axn25 = sum(AvgNm), unq25 = length(junc_seq))
g1ess15to25

# to determine changes in relative proportion of single action sites upon KDs
# put these into a single table and calculate relative ratio of each value

g1ess15singleA <- full_join(g1ess15to16, g1ess15to17, by = "KD")
g1ess15singleB <- full_join(g1ess15singleA, g1ess15to18, by ="KD")
g1ess15singleC <- full_join(g1ess15singleB, g1ess15to19, by ="KD")
g1ess15singleD <- full_join(g1ess15singleC, g1ess15to20, by ="KD") 
g1ess15singleE <- full_join(g1ess15singleD, g1ess15to21, by ="KD") 
g1ess15singleF <- full_join(g1ess15singleE, g1ess15to22, by ="KD") 
g1ess15singleG <- full_join(g1ess15singleF, g1ess15to23, by ="KD") 
g1ess15singleH <- full_join(g1ess15singleG, g1ess15to24, by ="KD") 
g1ess15single <- full_join(g1ess15singleH, g1ess15to25, by ="KD") 
g1ess15single

View(RPS12allg1avg)
g1ess15singleTalla <- g1ess15single %>% select(KD, contains("axn")) %>% 
  gather("Axn1", "AvgNm", 2:11) 
g1ess15singleTallb <- g1ess15single %>% select(KD, contains("unq")) %>%
  gather("Unq1", "UnqCt", 2:11)
g1ess15singleTall <- bind_cols(g1ess15singleTalla, g1ess15singleTallb[,2:3])
g1ess15singleTall

g1ess15singlePerc <- g1ess15singleTall %>%  group_by(KD) %>% 
  mutate(Tot = sum(AvgNm, na.rm=TRUE)) %>%
  rowwise() %>% mutate(Perc = 100*(AvgNm/Tot)) %>% select(-Unq1) %>% arrange(KD)
glimpse(g1ess15singlePerc)
View(g1ess15singlePerc)
## calculating the number of sequences with two or more modifications in this region involves a little more
# direct examination of the junction sequences...

g1ess16then17 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==2)||(edit_stop==16&junc_len==1)||(edit_stop==17&junc_len==0)) %>%
  mutate(ES16 = str_sub(junc_seq, -2)) %>% filter(ES16 != "GA") %>% group_by(KD) %>%
  summarise(axn16_17 = sum(AvgNm), unq16_17 = length(junc_seq))
g1ess16then17

# here I use grepl in the filter to remove sequences with ESS15 and modifications at ES17
g1ess16then18 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==3)||(edit_stop==16&junc_len==2)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -3)) %>% filter(ESS16_17 != "AGA") %>% 
  filter(!(edit_stop==16&ESS16_17=="ATG")) %>% filter((grepl('AG', junc_seq)&edit_stop==15)||edit_stop==16) %>%
  group_by(KD) %>% summarise(axn16_18 = sum(AvgNm), unq16_18 = length(junc_seq))
View(g1ess16then18)

g1ess16then19 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==4)||(edit_stop==16&junc_len==3)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -3)) %>% filter(ESS16_17 != "AGA") %>% 
  filter(!(edit_stop==16&ESS16_17=="ATG")) %>% filter((grepl('AG', junc_seq)&edit_stop==15)||edit_stop==16) %>%
  filter(grepl('GTTTA', junc_seq)) %>%
  group_by(KD) %>% summarise(axn16_19 = sum(AvgNm), unq16_19 = length(junc_seq))
g1ess16then19

g1ess16then20 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==5)||(edit_stop==16&junc_len==4)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -3)) %>% filter(ESS16_17 != "AGA") %>% 
  filter(!(edit_stop==16&ESS16_17=="ATG")) %>% filter((grepl('AG', junc_seq)&edit_stop==15)||edit_stop==16) %>%
  filter(grepl('CGTTTA', junc_seq)) %>%
  group_by(KD) %>% summarise(axn16_20 = sum(AvgNm), unq16_20 = length(junc_seq))
g1ess16then20

g1ess16then21 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==6)||(edit_stop==16&junc_len==5)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -3)) %>% filter(ESS16_17 != "AGA") %>% 
  filter(!(edit_stop==16&ESS16_17=="ATG")) %>% filter((grepl('AG', junc_seq)&edit_stop==15)||edit_stop==16) %>%
  filter(grepl('ATCGTTTA', junc_seq)) %>%
  group_by(KD) %>% summarise(axn16_21 = sum(AvgNm), unq16_21 = length(junc_seq))
g1ess16then21

g1ess16then22 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)||(edit_stop==16&junc_len==6)) %>%
  filter(grepl('CATCGTTTAG', junc_seq)) %>% filter(!(grepl('GA', junc_seq))) %>%
  group_by(KD) %>% summarise(axn16_22 = sum(AvgNm), unq16_22 = length(junc_seq))
g1ess16then22

g1ess16then23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)||(edit_stop==16&junc_len==7)) %>%
  filter(grepl('ACATCGTTTAG', junc_seq)) %>% filter(!(grepl('GA', junc_seq))) %>%
  group_by(KD) %>% summarise(axn16_23 = sum(AvgNm), unq16_23 = length(junc_seq))
g1ess16then23

# some pairs we won't be able to distinguish order in which they are generated
# given this, I'm going to assume that largest to largest is order
# since so few sequences had 17 and 19 as first site, I'll check the ESS18 to next ones...

g1ess18then17 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==3)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -3)) %>% filter(ESS16_17 != "AGA") %>% 
  filter((grepl('GA', ESS16_17))) %>%
  group_by(KD) %>% summarise(axn18_17 = sum(AvgNm), unq18_17 = length(junc_seq))
g1ess18then17

g1ess18then19  <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==4)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -3)) %>% filter(ESS16_17 == "AGA") %>% 
  filter(!(grepl('GTTTA', junc_seq))) %>%
  group_by(KD) %>% summarise(axn18_19 = sum(AvgNm), unq18_19 = length(junc_seq))
g1ess18then19

g1ess18then20  <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==5)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -3)) %>% filter(ESS16_17 == "AGA") %>% 
  filter(!(grepl('GTTTA', junc_seq))) %>% filter((grepl('CG', junc_seq))) %>%
  group_by(KD) %>% summarise(axn18_20 = sum(AvgNm), unq18_20 = length(junc_seq))
View(g1ess18then20)

g1ess18then21 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==6)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -3)) %>% filter(ESS16_17 == "AGA") %>% 
  filter(!(grepl('GTTTA', junc_seq))) %>% filter((grepl('ATCG', junc_seq))) %>%
  group_by(KD) %>% summarise(axn18_21 = sum(AvgNm), unq18_21 = length(junc_seq))
g1ess18then21

g1ess18then22 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -3)) %>% filter(ESS16_17 == "AGA") %>% 
  filter(!(grepl('GTTTA', junc_seq))) %>% filter((grepl('CATCG', junc_seq))) %>%
  group_by(KD) %>% summarise(axn18_22 = sum(AvgNm), unq18_22 = length(junc_seq))
View(g1ess18then22)

g1ess18then23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -3)) %>% filter(ESS16_17 == "AGA") %>% 
  filter(!(grepl('GTTTA', junc_seq))) %>% filter((grepl('ACATCG', junc_seq))) %>%
  group_by(KD) %>% summarise(axn18_23 = sum(AvgNm), unq18_23 = length(junc_seq))
View(g1ess18then23)


# fill in the missing pairwise pairs (unlikely to exist given single piece data but do for completeness)

g1ess17then19 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==4)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% 
  filter((grepl('GTTTA', junc_seq))) %>% filter(!(grepl('AGA', junc_seq))) %>%
  group_by(KD) %>% summarise(axn17_19 = sum(AvgNm), unq17_19 = length(junc_seq))
g1ess17then19

g1ess17then20 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==5)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% 
  filter((grepl('CGTTTA', junc_seq))) %>% filter(!(grepl('AG', junc_seq)))%>%
  group_by(KD) %>% summarise(axn17_20 = sum(AvgNm), unq17_20 = length(junc_seq))
g1ess17then20

g1ess17then21 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==6)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% 
  filter((grepl('ATCGTTTA', junc_seq))) %>% filter(!(grepl('AG', junc_seq))) %>%
  group_by(KD) %>% summarise(axn17_21 = sum(AvgNm), unq17_21 = length(junc_seq))
g1ess17then21

g1ess17then22 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% 
  filter((grepl('CATCGTTTA', junc_seq))) %>% filter(!(grepl('AG', junc_seq))) %>%
  group_by(KD) %>% summarise(axn17_22 = sum(AvgNm), unq17_22 = length(junc_seq))
g1ess17then22

g1ess17then23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% 
  filter((grepl('CATCGTTTA', junc_seq))) %>% filter(!(grepl('AG', junc_seq))) %>%
  filter((grepl('ACA', junc_seq))) %>%
  group_by(KD) %>% summarise(axn17_23 = sum(AvgNm), unq17_23 = length(junc_seq))
g1ess17then23

g1ess19then20 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==5)) %>%
  filter((grepl('GTTTAGA', junc_seq))) %>% filter(!(grepl('CG', junc_seq))) %>%
  group_by(KD) %>% summarise(axn19_20 = sum(AvgNm), unq19_20 = length(junc_seq))
g1ess19then20

g1ess19then21 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==6)) %>%
  filter((grepl('GTTTAGA', junc_seq))) %>% filter(!(grepl('CG', junc_seq))) %>%
  group_by(KD) %>% summarise(axn19_21 = sum(AvgNm), unq19_21 = length(junc_seq))
g1ess19then21

g1ess19then22 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  filter((grepl('GTTTAGA', junc_seq))) %>% filter(!(grepl('CG', junc_seq))) %>%
  group_by(KD) %>% summarise(axn19_22 = sum(AvgNm), unq19_22 = length(junc_seq))
g1ess19then22

g1ess19then23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  filter((grepl('GTTTAGA', junc_seq))) %>% filter(!(grepl('CG', junc_seq))) %>%
  filter((grepl('ACA', junc_seq))) %>%
  group_by(KD) %>% summarise(axn19_23 = sum(AvgNm), unq19_23 = length(junc_seq))
g1ess19then23

g1ess20then21 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==6)) %>%
  filter((grepl('CGTTTAGA', junc_seq))) %>% filter(!(grepl('CATC', junc_seq))) %>%
  group_by(KD) %>% summarise(axn20_21 = sum(AvgNm), unq20_21 = length(junc_seq))
g1ess20then21

g1ess20then22 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  filter((grepl('CGTTTAGA', junc_seq))) %>% filter(!(grepl('ATC', junc_seq))) %>%
  filter((grepl('CA', junc_seq))) %>%
  group_by(KD) %>% summarise(axn20_22 = sum(AvgNm), unq20_22 = length(junc_seq))
g1ess20then22

g1ess20then23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  filter((grepl('CGTTTAGA', junc_seq))) %>% filter(!(grepl('ATC', junc_seq))) %>%
  filter((grepl('ACA', junc_seq))) %>%
  group_by(KD) %>% summarise(axn20_23 = sum(AvgNm), unq20_23 = length(junc_seq))
g1ess20then23

g1ess21then22 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  filter((grepl('ATCGTTTAGA', junc_seq))) %>% filter(!(grepl('ACA', junc_seq))) %>% 
  filter(!(grepl('CA', junc_seq))) %>% 
  group_by(KD) %>% summarise(axn21_22 = sum(AvgNm), unq21_22 = length(junc_seq))
g1ess21then22

g1ess21then23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  filter((grepl('ATCGTTTAGA', junc_seq))) %>% filter((grepl('AC', junc_seq))) %>%
  filter(!(grepl('CA', junc_seq))) %>% 
  group_by(KD) %>% summarise(axn21_23 = sum(AvgNm), unq21_23 = length(junc_seq))
g1ess21then23

g1ess22then23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  filter((grepl('CATCGTTTAGA', junc_seq))) %>% filter(!(grepl('AC', junc_seq))) %>%
  group_by(KD) %>% summarise(axn22_23 = sum(AvgNm), unq22_23 = length(junc_seq))
g1ess22then23

## make a single table of the major second action sequences

g1ess15twoA <- full_join(g1ess16then17, g1ess16then18, by = "KD")
g1ess15twoB <- full_join(g1ess15twoA, g1ess16then19, by = "KD")
g1ess15twoC <- full_join(g1ess15twoB, g1ess16then20, by = "KD")
g1ess15twoD <- full_join(g1ess15twoC, g1ess16then21, by = "KD")
g1ess15twoE <- full_join(g1ess15twoD, g1ess16then22, by = "KD")
g1ess15twoF <- full_join(g1ess15twoE, g1ess18then17, by = "KD")
g1ess15twoG <- full_join(g1ess15twoF, g1ess18then19, by = "KD")
g1ess15twoH <- full_join(g1ess15twoG, g1ess18then20, by = "KD")
g1ess15twoI <- full_join(g1ess15twoH, g1ess18then21, by = "KD")
g1ess15twoJ <- full_join(g1ess15twoI, g1ess18then22, by = "KD")
g1ess15twoK <- full_join(g1ess15twoJ, g1ess17then19, by = "KD")
g1ess15twoL <- full_join(g1ess15twoK, g1ess17then20, by = "KD")
g1ess15twoM <- full_join(g1ess15twoL, g1ess17then21, by = "KD")
g1ess15twoN <- full_join(g1ess15twoM, g1ess19then20, by = "KD")
g1ess15twoO <- full_join(g1ess15twoN, g1ess19then21, by = "KD")
g1ess15twoP <- full_join(g1ess15twoO, g1ess19then22, by = "KD")
g1ess15twoQ <- full_join(g1ess15twoP, g1ess20then21, by = "KD")
g1ess15twoR <- full_join(g1ess15twoQ, g1ess20then22, by = "KD")
g1ess15twoS <- full_join(g1ess15twoR, g1ess16then23, by = "KD")
g1ess15twoT <- full_join(g1ess15twoS, g1ess17then23, by = "KD")
g1ess15twoU <- full_join(g1ess15twoT, g1ess18then23, by = "KD")
g1ess15twoV <- full_join(g1ess15twoU, g1ess19then23, by = "KD")
g1ess15twoW <- full_join(g1ess15twoV, g1ess20then23, by = "KD")
g1ess15twoX <- full_join(g1ess15twoW, g1ess21then23, by = "KD")
g1ess15twoY <- full_join(g1ess15twoX, g1ess22then23, by = "KD")
g1ess15twoAxn <- full_join(g1ess15twoY, g1ess21then22, by = "KD")

g1ess15twoATalla <- g1ess15twoAxn %>% select(KD, contains("axn")) %>%
  gather("Axn1", "AvgNm", contains("axn"))
g1ess15twoATallb <-  g1ess15twoAxn %>% select(KD, contains("unq")) %>%
  gather("Unq1", "UnqCt", contains("unq"))
g1ess15twoATall <- bind_cols(g1ess15twoATalla, g1ess15twoATallb[,2:3])
g1ess15twoATall

g1ess15twoAPerc <- g1ess15twoATall %>% group_by(KD) %>% 
  mutate(Tot = sum(AvgNm,na.rm=TRUE)) %>%
  rowwise() %>% mutate(Perc = 100*(AvgNm/Tot)) %>% select(-Unq1) %>% arrange(KD, Axn1)
head(g1ess15twoAPerc)
View(g1ess15twoAPerc)

write.table(g1ess15twoAPerc, "RPS12gRNA1twoAxn.csv", sep= ",")

## third step...
# this portion i taken piece-wise considering first the major steps1 and 2

g1ess15_16_17_18 <-  RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==3)||(edit_stop==16&junc_len==2)||
           (edit_stop==17&junc_len==1)|| (edit_stop==19&junc_len==0)) %>% #all junc len options
  filter(!(grepl("AGA", junc_seq))) %>%  # removes those that skip 16 and 17
  mutate(ES17 = str_sub(junc_seq, -2)) %>% filter(ES17 != "AG",ES17!="GA") %>% 
  # removes non-modified 16 and 17 when they're JSS, ES18 modified by default as JES
  group_by(KD) %>%  summarise(axn16_17_18 = sum(AvgNm), unq16_17_18 = length(junc_seq))
g1ess15_16_17_18

g1ess15_16_17_19 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==4)||(edit_stop==16&junc_len==3)||
           (edit_stop==17&junc_len==2)) %>%  
  filter(!(grepl("AGA", junc_seq))) %>% filter((grepl("GTTTA", junc_seq))) %>% 
  mutate(ES17 = str_sub(junc_seq, -2)) %>% filter(ES17 != "AG") %>% 
  filter(!(grepl("CG", junc_seq))) %>% filter(!(grepl('GA', junc_seq))) %>%
  group_by(KD) %>%  summarise(axn16_17_19 = sum(AvgNm), unq16_17_19 = length(junc_seq))
g1ess15_16_17_19

g1ess15_16_17_20 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==5)||(edit_stop==16&junc_len==4)||
           (edit_stop==17&junc_len==3)) %>%
  filter(!(grepl("AGA", junc_seq))) %>% filter((grepl("CGTTTA", junc_seq))) %>% 
  mutate(ES17 = str_sub(junc_seq, -2)) %>% filter(ES17 != "AG") %>%
  group_by(KD) %>%  summarise(axn16_17_20 = sum(AvgNm), unq16_17_20 = length(junc_seq))
g1ess15_16_17_20

g1ess15_16_17_21 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==6)||(edit_stop==16&junc_len==5)||
           (edit_stop==17&junc_len==4)) %>% filter((grepl('ATC', junc_seq))) %>%
  filter(!(grepl("AGA", junc_seq))) %>% filter((grepl("ACGTTTA", junc_seq))) %>% 
  mutate(ES17 = str_sub(junc_seq, -2)) %>% filter(ES17 != "AG") %>%
  group_by(KD) %>%  summarise(axn16_17_21 = sum(AvgNm), unq16_17_21 = length(junc_seq))
g1ess15_16_17_21

g1ess15_16_17_22 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==7)||(edit_stop==16&junc_len==6)||
           (edit_stop==17&junc_len==5)) %>%
  filter(!(grepl("AGA", junc_seq))) %>% filter((grepl("CATCGTTTA", junc_seq))) %>% 
  mutate(ES17 = str_sub(junc_seq, -2)) %>% filter(ES17 != "AG") %>%
  filter(ES17 != "GA") %>%
  group_by(KD) %>%  summarise(axn16_17_22 = sum(AvgNm), unq16_17_22 = length(junc_seq))
g1ess15_16_17_22

g1ess15_16_17_23 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==8)||(edit_stop==16&junc_len==7)||
           (edit_stop==17&junc_len==6)) %>%
  filter(!(grepl("AGA", junc_seq))) %>% filter((grepl("ACATCGTTTA", junc_seq))) %>% 
  mutate(ES17 = str_sub(junc_seq, -2)) %>% filter(ES17 != "AG") %>%
  group_by(KD) %>%  summarise(axn16_17_23 = sum(AvgNm), unq16_17_23 = length(junc_seq))
g1ess15_16_17_23

g1ess15_16_18_19 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==4)||(edit_stop==16&junc_len==3)) %>%
  filter((grepl("AG", junc_seq))) %>% filter(!(grepl("GTTTA", junc_seq))) %>% 
  filter(!(grepl("ATG", junc_seq))) %>%
  mutate(ES16 = str_sub(junc_seq, -2)) %>% filter(ES16 != "GA") %>%
  group_by(KD) %>%  summarise(axn16_18_19 = sum(AvgNm), unq16_18_19 = length(junc_seq))
g1ess15_16_18_19

g1ess15_16_18_20 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==5)||(edit_stop==16&junc_len==4)) %>%
  filter(!(grepl("AGA", junc_seq))) %>% filter(!(grepl("CGTTTA", junc_seq))) %>% 
  filter((grepl("AG", junc_seq))) %>% filter((grepl('CG', junc_seq))) %>%
  mutate(ES16 = str_sub(junc_seq, -2)) %>% filter(ES16 != "GA") %>%
  group_by(KD) %>%  summarise(axn16_18_20 = sum(AvgNm), unq16_18_20 = length(junc_seq))
g1ess15_16_18_20

g1ess15_16_18_21 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==6)||(edit_stop==16&junc_len==5)) %>% # limit to proper jxn lengths
  mutate(ES16 = str_sub(junc_seq, -2)) %>% filter(ES16 != "GA") %>% # ensures ES16 is modified
  filter((grepl("AG", junc_seq))) %>% filter((grepl('CG', junc_seq))) %>% #first ensures ES17 unmod, second 19unmod
  filter(!(grepl("CGTTTA", junc_seq))) %>% filter((grepl('ATCG', junc_seq))) %>%
  #first ES18mod (can incl CG here to diff from ES16 mod because of limits on previous line)
  # second ES20 is unmod
  group_by(KD) %>%  summarise(axn16_18_21 = sum(AvgNm), unq16_18_21 = length(junc_seq))
g1ess15_16_18_21

g1ess15_16_18_22 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==7)||(edit_stop==16&junc_len==6)) %>%
  mutate(ES16 = str_sub(junc_seq, -2)) %>% filter(ES16 != "GA") %>% # ensures ES16 is modified
  filter((grepl("AG", junc_seq))) %>% filter((grepl('CG', junc_seq))) %>% #first ensures ES17 unmod, second 19unmod
  filter(!(grepl("CGTTTA", junc_seq))) %>% filter((grepl('CATCG', junc_seq))) %>% #ES18mod, ES20 and 21 unmod
  group_by(KD) %>%  summarise(axn16_18_22 = sum(AvgNm), unq16_18_22 = length(junc_seq))
g1ess15_16_18_22

g1ess15_16_18_23 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==8)||(edit_stop==16&junc_len==7)) %>%
  mutate(ES16 = str_sub(junc_seq, -2)) %>% filter(ES16 != "GA") %>% # ensures ES16 is modified
  filter((grepl("AG", junc_seq))) %>% filter((grepl('CG', junc_seq))) %>% #first ensures ES17 unmod, second 19unmod
  filter(!(grepl("CGTTTA", junc_seq))) %>% filter((grepl('ACATCG', junc_seq))) %>% #ES18mod, ES20, 21 and 22 unmod
  group_by(KD) %>%  summarise(axn16_18_23 = sum(AvgNm), unq16_18_23 = length(junc_seq))
g1ess15_16_18_23

g1ess15_16_19_21 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==6)||(edit_stop==16&junc_len==5)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 != "GA") %>% #ES16 mod
  filter((grepl("AG", junc_seq))) %>% #ES17unmod
  filter(grepl('GTTTAG', junc_seq)) %>% #ES18mod (G diffs from ES16)
  filter(!(grepl('CG', junc_seq))) %>% #ES19 mod, 
  filter((grepl('ATG', junc_seq))) %>%  #ES20 unmod, ES21 mod by default jes
  group_by(KD) %>% summarise(axn16_19_21 = sum(AvgNm), unq16_19_21 = length(junc_seq))
g1ess15_16_19_21

g1ess15_16_20_19 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==5)||(edit_stop==16&junc_len==4)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 != "GA") %>% #ES16 mod
  filter((grepl("AG", junc_seq))) %>% #ES17unmod
  filter(!(grepl('ATG', junc_seq))) %>%  #ES20 mod
  filter(!(grepl('CG', junc_seq))) %>% #ES19 mod, 
  filter(grepl('GTTTAG', junc_seq)) %>% #ES18mod (G diffs from ES16)
  group_by(KD) %>% summarise(axn16_20_19 = sum(AvgNm), unq16_20_19 = length(junc_seq))
g1ess15_16_20_19

g1ess15_16_19_22 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==7)||(edit_stop==16&junc_len==6)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 != "GA") %>% #ES16 mod
  filter((grepl("AG", junc_seq))) %>% #ES17unmod
  filter(grepl('GTTTAG', junc_seq)) %>% #ES18mod (G diffs from ES16)
  filter(!(grepl('CG', junc_seq))) %>% #ES19 mod, 
  filter((grepl('CATC', junc_seq))) %>%  #ES20 unmod, ES21 unmod, ES22 mod by default jes
  group_by(KD) %>% summarise(axn16_19_22 = sum(AvgNm), unq16_19_22 = length(junc_seq))
g1ess15_16_19_22

g1ess15_16_19_23 <- RPS12allg1avg %>% 
  filter((edit_stop==15&junc_len==8)||(edit_stop==16&junc_len==7)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 != "GA") %>% #ES16 mod
  filter((grepl("AG", junc_seq))) %>% #ES17unmod
  filter(grepl('GTTTAG', junc_seq)) %>% #ES18mod (G diffs from ES16)
  filter(!(grepl('CG', junc_seq))) %>% #ES19 mod, 
  filter((grepl('ACATC', junc_seq))) %>%  #ES20-22 unmod, ES23 mod by default jes
  group_by(KD) %>% summarise(axn16_19_23 = sum(AvgNm), unq16_19_23 = length(junc_seq))
g1ess15_16_19_23

g1ess15_16_20_21 <- RPS12allg1avg %>%
  filter((edit_stop==15&junc_len==6)||(edit_stop==16&junc_len==5)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 != "GA") %>% # ES16 mod
  filter(grepl('CGTTTAG', junc_seq)) %>% # ES17, 18, 19 unmodified
  filter(!(grepl('ATC', junc_seq))) %>%  #ES20 mod, ES21 mod implied by junc end
  group_by(KD) %>% summarise(axn16_20_21 = sum(AvgNm), unq16_20_21 = length(junc_seq))
g1ess15_16_20_21

g1ess15_16_20_22 <- RPS12allg1avg %>%
  filter((edit_stop==15&junc_len==7)||(edit_stop==16&junc_len==6)) %>% # 
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>%  filter(ESS16_17 != "GA") %>% #ES16 modified
  filter((grepl('AG', junc_seq))) %>% #ES 17 unmod
  filter(!(grepl('ATC', junc_seq))) %>%   #ES 20 mod
  filter(grepl('CGTTTA', junc_seq)) %>% # ES 18 and 19 unmod
  filter(grepl('CA', junc_seq)) %>% #ES 21 un mod, ES22 mod by default as JES
  group_by(KD) %>% summarise(axn16_20_22 = sum(AvgNm), unq16_20_22 = length(junc_seq))
g1ess15_16_20_22

g1ess15_16_20_23 <- RPS12allg1avg %>%
  filter((edit_stop==15&junc_len==8)||(edit_stop==16&junc_len==7)) %>% # 
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>%  filter(ESS16_17 != "GA") %>% #ES16 modified
  filter(grepl('CGTTTAG', junc_seq)) %>% # ES17-19 unmod
  filter(!(grepl('ATC', junc_seq))) %>%   #ES 20 mod
  filter(grepl('ACA', junc_seq)) %>% #ES 21-22 un mod, ES23 mod by default as JES
  group_by(KD) %>% summarise(axn16_20_23 = sum(AvgNm), unq16_20_23 = length(junc_seq))
g1ess15_16_20_23

g1ess15_16_21_22 <- RPS12allg1avg %>%
  filter((edit_stop==15&junc_len==7)||(edit_stop==16&junc_len==6)) %>% # 
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>%  filter(ESS16_17 != "GA") %>% #ES16 modified
  filter(grepl('ATCGTTTAG', junc_seq)) %>% # ES17-20 unmod
  filter(!(grepl('AC', junc_seq))) %>%   #ES 21 mod
  group_by(KD) %>% summarise(axn16_21_22 = sum(AvgNm), unq16_21_22 = length(junc_seq))
g1ess15_16_21_22

g1ess15_16_21_23 <- RPS12allg1avg %>%
  filter((edit_stop==15&junc_len==8)||(edit_stop==16&junc_len==7)) %>% # 
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>%  filter(ESS16_17 != "GA") %>% #ES16 modified
  filter(grepl('ATCGTTTAG', junc_seq)) %>% # ES17-20 unmod
  filter(!(grepl('CA', junc_seq))) %>%   #ES 21 mod
  filter((grepl('AC', junc_seq))) %>% #ES22 unmod
  group_by(KD) %>% summarise(axn16_21_23 = sum(AvgNm), unq16_21_23 = length(junc_seq))
g1ess15_16_21_23

g1ess15_16_22_23 <- RPS12allg1avg %>%
  filter((edit_stop==15&junc_len==8)||(edit_stop==16&junc_len==7)) %>% # 
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>%  filter(ESS16_17 != "GA") %>% #ES16 modified
  filter(grepl('CATCGTTTAG', junc_seq)) %>% # ES17-21 unmod
  filter(!(grepl('AC', junc_seq))) %>% #ES22 mod
  group_by(KD) %>% summarise(axn16_22_23 = sum(AvgNm), unq16_22_23 = length(junc_seq))
g1ess15_16_22_23

g1ess15_17_18_19 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==4)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% # ES16 unmod
  filter(!(grepl('GTTTA', junc_seq))) %>% filter(!(grepl('AG', junc_seq))) %>% # ES18 mod, ES17mod
  group_by(KD) %>% summarise(axn17_18_19 = sum(AvgNm), unq17_18_19 = length(junc_seq))
g1ess15_17_18_19

g1ess15_18_20_17 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==5)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% # ES16 unmod
  filter(!(grepl('GTTTA', junc_seq))) %>% filter((grepl('CG', junc_seq))) %>% # ES18 mod, ES19 un mod
  filter(!(grepl('AG', junc_seq))) %>% filter(!(grepl('ATC', junc_seq))) %>% # ES17 and ES20 mod
  group_by(KD) %>% summarise(axn18_20_17 = sum(AvgNm), unq18_20_17 = length(junc_seq))
g1ess15_18_20_17

g1ess15_18_17_21 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==6)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% # ES16 unmod
  filter(!(grepl('GTTTA', junc_seq))) %>% filter((grepl('CG', junc_seq))) %>% # ES18 mod, ES19 un mod
  filter(!(grepl('AG', junc_seq))) %>% filter((grepl('ATC', junc_seq))) %>% # ES17 and ES20 unmod
  group_by(KD) %>% summarise(axn18_17_21 = sum(AvgNm), unq18_17_21 = length(junc_seq))
g1ess15_18_17_21

g1ess15_18_17_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% # ES16 unmod
  filter(!(grepl('AG', junc_seq))) %>% #ES17 mod
  filter(!(grepl('GTTTA', junc_seq))) %>% # ES18 mod
  filter((grepl('ACATCG', junc_seq))) %>% #ES19-22 unmod
  group_by(KD) %>% summarise(axn18_17_23 = sum(AvgNm), unq18_17_23 = length(junc_seq))
g1ess15_18_17_23

g1ess15_18_20_19 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==5)) %>% 
  filter(grepl('AGA', junc_seq)) %>% # ES16 an 17 unmod
  filter(!(grepl('GTTTA', junc_seq))) %>% filter(!(grepl('CG', junc_seq))) %>% #ES18 and 19 mod
  filter(!(grepl('ATC', junc_seq))) %>% # ES20 mod
  group_by(KD) %>% summarise(axn18_20_19 = sum(AvgNm), unq18_20_19 = length(junc_seq))
g1ess15_18_20_19

g1ess15_18_20_21 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==6)) %>%
  filter(grepl('AGA', junc_seq)) %>% #ES16 and ES17 unmod
  filter(!(grepl('GTTTA', junc_seq))) %>% filter((grepl('CG', junc_seq))) %>% #ES18mod, ES19 unmod
  filter(!(grepl('ATC', junc_seq))) %>% # ES20 mod, ES21 mod by default as JES
  group_by(KD) %>% summarise(axn18_20_21 = sum(AvgNm), unq18_20_21 = length(junc_seq))
g1ess15_18_20_21

g1ess15_18_20_22 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  filter(grepl('AGA', junc_seq)) %>% #ES16 and ES17 unmod
  filter(!(grepl('GTTTA', junc_seq))) %>% filter((grepl('CG', junc_seq))) %>% #ES18mod, ES19 unmod
  filter(!(grepl('ATC', junc_seq))) %>% filter((grepl('CA', junc_seq))) %>% 
  # ES20 mod, ES21 unmod, ES22 mod by default as JES
  group_by(KD) %>% summarise(axn18_20_22 = sum(AvgNm), unq18_20_22 = length(junc_seq))
g1ess15_18_20_22

g1ess15_18_20_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  filter(grepl('AGA', junc_seq)) %>% #ES16 and ES17 unmod
  filter(!(grepl('GTTTA', junc_seq))) %>% filter((grepl('CG', junc_seq))) %>% #ES18mod, ES19 unmod
  filter(!(grepl('ATC', junc_seq))) %>% filter((grepl('ACA', junc_seq))) %>% 
  # ES20 mod, ES21, 22 unmod, ES23 mod by default as JES
  group_by(KD) %>% summarise(axn18_20_23 = sum(AvgNm), unq18_20_23 = length(junc_seq))
View(g1ess15_18_20_23)

g1ess15_18_22_17 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% #ES16 unmod
  filter(!(grepl('AG',junc_seq))) %>%  #ES17 mod
  filter(!(grepl('GTTTA', junc_seq))) %>%  #ES18 mod
  filter((grepl('CATCG', junc_seq))) %>% # ES19-21 unmod, ES22 mod as JES
  group_by(KD) %>% summarise(axn18_22_17 = sum(AvgNm), unq18_22_17 = length(junc_seq))
g1ess15_18_22_17

g1ess15_18_22_19 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  filter(!(grepl("AGA", junc_seq))) %>% # ES16 and 17 unmod
  filter(!(grepl('GTTTA', junc_seq))) %>% filter((grepl('CATC', junc_seq))) %>% # ES18 mod, ES20, 21 unmod
  filter(!(grepl('CG', junc_seq))) %>%  #ES19mod
  group_by(KD) %>% summarise(axn18_22_19 = sum(AvgNm), unq18_22_19 = length(junc_seq))
g1ess15_18_22_19


g1ess15_18_22_21 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  filter((grepl("AGA",junc_seq))) %>% # ES16 and ES17 unmod
  filter(!(grepl('GTTTA', junc_seq))) %>% filter((grepl('ATCG', junc_seq))) %>% #ES18mod, ES19, 20 unmod
  filter(!(grepl('CA', junc_seq))) %>% #ES21 mod
  group_by(KD) %>% summarise(axn18_22_21 = sum(AvgNm), unq18_22_21 = length(junc_seq))
g1ess15_18_22_21

g1ess15_18_22_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  filter((grepl("AGA",junc_seq))) %>% # ES16 and ES17 unmod
  filter(!(grepl('GTTTA', junc_seq))) %>% #ES18 mod
  filter((grepl('CATCG', junc_seq))) %>% filter(!(grepl('AC', junc_seq))) %>% #ES19, 20, 21 unmod; ES22mod
  group_by(KD) %>% summarise(axn18_22_23 = sum(AvgNm), unq18_22_23 = length(junc_seq))
g1ess15_18_22_23

g1ess15_17_19_20 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==5)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% #ES16 unmod
  filter(!(grepl('AG', junc_seq))) %>% #ES17mod
  filter((grepl('GTTTA', junc_seq))) %>% #ES18 unmod
  filter(!(grepl('CG', junc_seq))) %>% #ES19mod, ES20 mod by jes
  group_by(KD) %>% summarise(axn17_19_20 = sum(AvgNm), unq17_19_20 = length(junc_seq))
g1ess15_17_19_20

g1ess15_17_19_21 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==6)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% #ES16 unmod
  filter(!(grepl('AG', junc_seq))) %>% #ES17mod
  filter((grepl('GTTTA', junc_seq))) %>% #ES18 unmod
  filter(!(grepl('CG', junc_seq))) %>% filter(grepl('ATC', junc_seq)) %>% #ES19mod, ES20 unmod, ES21 by jes
  group_by(KD) %>% summarise(axn17_19_21 = sum(AvgNm), unq17_19_21 = length(junc_seq))
g1ess15_17_19_21

g1ess15_17_19_22 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% #ES16 unmod
  filter(!(grepl('AG', junc_seq))) %>% #ES17mod
  filter((grepl('GTTTA', junc_seq))) %>% #ES18 unmod
  filter(!(grepl('CG', junc_seq))) %>% filter(grepl('CATC', junc_seq)) %>% #ES19mod, ES20-21 unmod, ES22 by jes
  group_by(KD) %>% summarise(axn17_19_22 = sum(AvgNm), unq17_19_22 = length(junc_seq))
g1ess15_17_19_22

g1ess15_17_19_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% #ES16 unmod
  filter(!(grepl('AG', junc_seq))) %>% #ES17mod
  filter((grepl('GTTTA', junc_seq))) %>% #ES18 unmod
  filter(!(grepl('CG', junc_seq))) %>% filter(grepl('ACATC', junc_seq)) %>% #ES19mod, ES20-22 unmod, ES23 by jes
  group_by(KD) %>% summarise(axn17_19_23 = sum(AvgNm), unq17_19_23 = length(junc_seq))
g1ess15_17_19_23

g1ess15_17_20_21 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==6)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% #ES16 unmod
  filter(!(grepl('AG', junc_seq))) %>% #ES17mod
  filter((grepl('GTTTA', junc_seq))) %>% #ES18 unmod
  filter((grepl('CG', junc_seq))) %>% filter(!(grepl('ATC', junc_seq))) %>% #ES19unmod, ES20unmod, ES21mod by JES
  group_by(KD) %>% summarise(axn17_20_21 = sum(AvgNm), unq17_20_21 = length(junc_seq))
g1ess15_17_20_21

g1ess15_20_22_17 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% #ES16 unmod
  filter(!(grepl('AG', junc_seq))) %>% #ES17 mod
  filter((grepl('CGTTTA', junc_seq))) %>% filter(!(grepl('ATC', junc_seq))) %>% #ES18,19unmod, ES20mod
  filter((grepl('CA', junc_seq))) %>%  #ES21unmod, ES22mod by default as JES
  group_by(KD) %>% summarise(axn20_22_17 = sum(AvgNm), unq20_22_17 = length(junc_seq)) #
g1ess15_20_22_17

g1ess15_17_20_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% #ES16 unmod
  filter(!(grepl('AG', junc_seq))) %>% #ES17mod
  filter((grepl('GTTTA', junc_seq))) %>% #ES18 unmod
  filter((grepl('CG', junc_seq))) %>% filter(!(grepl('ATC', junc_seq))) %>% #ES19unmod, ES20unmod, ES21mod by JES
  filter((grepl('ACA', junc_seq))) %>% #ES21,22 unmod
  group_by(KD) %>% summarise(axn17_20_23 = sum(AvgNm), unq17_20_23 = length(junc_seq))
g1ess15_17_20_23

g1ess15_17_21_22 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% #ES16 unmod
  filter(!(grepl('AG', junc_seq))) %>% #ES17mod
  filter((grepl('ATCGTTTA', junc_seq))) %>% #ES18-20 unmod
  filter(!(grepl('CA', junc_seq))) %>% #ES21 mod, ES22 mod by JES
  group_by(KD) %>% summarise(axn17_21_22 = sum(AvgNm), unq17_21_22 = length(junc_seq))
g1ess15_17_21_22

g1ess15_17_21_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% #ES16 unmod
  filter(!(grepl('AG', junc_seq))) %>% #ES17mod
  filter((grepl('ATCGTTTA', junc_seq))) %>% #ES18-20 unmod
  filter(!(grepl('CA', junc_seq))) %>% filter((grepl('AC', junc_seq))) %>% #ES21 mod, ES22 mod by JES
  group_by(KD) %>% summarise(axn17_21_23 = sum(AvgNm), unq17_21_23 = length(junc_seq))
g1ess15_17_21_23

g1ess15_17_22_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -2)) %>% filter(ESS16_17 == "GA") %>% #ES16 unmod
  filter(!(grepl('AG', junc_seq))) %>% #ES17mod
  filter((grepl('CATCGTTTA', junc_seq))) %>% #ES18-21 unmod
  filter(!(grepl('AC', junc_seq))) %>% #ES22 mod, ES23 mod by JES
  group_by(KD) %>% summarise(axn17_22_23 = sum(AvgNm), unq17_22_23 = length(junc_seq))
g1ess15_17_22_23

g1ess15_18_19_21 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==6)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -3)) %>% filter(ESS16_17 == "AGA") %>% #ES16,17 unmod
  filter(!(grepl('GTTTA', junc_seq))) %>% #ES18 mod
  filter(!(grepl('CG', junc_seq))) %>% filter((grepl('ATC',junc_seq))) %>% #ES19 mod, ES20 unmod, ES21mod JES
  group_by(KD) %>% summarise(axn18_19_21  = sum(AvgNm), unq18_19_21  = length(junc_seq))
g1ess15_18_19_21 

g1ess15_18_19_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -3)) %>% filter(ESS16_17 == "AGA") %>% #ES16,17 unmod
  filter(!(grepl('GTTTA', junc_seq))) %>% #ES18 mod
  filter(!(grepl('CG', junc_seq))) %>% filter((grepl('ACATC',junc_seq))) %>% #ES19 mod, ES20-22 unmod, ES23mod JES
  group_by(KD) %>% summarise(axn18_19_23  = sum(AvgNm), unq18_19_23  = length(junc_seq))
g1ess15_18_19_23 

g1ess15_18_21_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  mutate(ESS16_17 = str_sub(junc_seq, -3)) %>% filter(ESS16_17 == "AGA") %>% #ES16,17 unmod
  filter(!(grepl('GTTTA', junc_seq))) %>% #ES18 mod
  filter(!(grepl('CG', junc_seq))) %>% filter((grepl('ATC',junc_seq))) %>% #ES19 mod, ES20 unmod, 
  filter((grepl('AC',junc_seq))) %>% filter(!(grepl('CA', junc_seq))) %>% #ES22unmod ES21 mod, ES23mod by JES
  group_by(KD) %>% summarise(axn18_21_23 = sum(AvgNm), unq18_21_23  = length(junc_seq))
g1ess15_18_21_23

g1ess15_19_20_21 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==6)) %>%
  filter((grepl('GTTTAGA', junc_seq))) %>% #ES16-18 unmod
  filter(!(grepl('CG', junc_seq))) %>% filter(!(grepl('ATC',junc_seq))) %>% #ES19 mod, ES20 mod, ES21 by JES
  group_by(KD) %>% summarise(axn19_20_21 = sum(AvgNm), unq19_20_21  = length(junc_seq))
g1ess15_19_20_21

g1ess15_20_22_19 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  filter((grepl('GTTTAGA', junc_seq))) %>% filter(!(grepl('CG', junc_seq))) %>%  #ES16-18 unmod, ES19 mod
  filter(!(grepl('ATC', junc_seq))) %>% filter((grepl('CA', junc_seq))) %>% #ES20 mod, ES21unmod, 22 by JES 
  group_by(KD) %>% summarise(axn20_22_19 = sum(AvgNm), unq20_22_19 = length(junc_seq))
g1ess15_20_22_19

g1ess15_19_21_22 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  filter((grepl('GTTTAGA', junc_seq))) %>% filter(!(grepl('CG', junc_seq))) %>%  #ES16-18 unmod, ES19 mod
  filter((grepl('ATC', junc_seq))) %>% filter(!(grepl('CA', junc_seq))) %>% #ES20 unmod, ES21unmod, 22 by JES 
  group_by(KD) %>% summarise(axn19_21_22 = sum(AvgNm), unq19_21_22 = length(junc_seq))
g1ess15_19_21_22

g1ess15_19_21_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  filter((grepl('GTTTAGA', junc_seq))) %>% filter(!(grepl('CG', junc_seq))) %>%  #ES16-18 unmod, ES19 mod
  filter((grepl('ATC', junc_seq))) %>% filter(!(grepl('CA', junc_seq))) %>% #ES20 unmod, ES21unmod, 22 by JES 
  filter((grepl('AC', junc_seq))) %>% #ES22 unmod
  group_by(KD) %>% summarise(axn19_21_23 = sum(AvgNm), unq19_21_23 = length(junc_seq))
g1ess15_19_21_23

g1ess15_19_22_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  filter((grepl('GTTTAGA', junc_seq))) %>% filter(!(grepl('CG', junc_seq))) %>%  #ES16-18 unmod, ES19 mod
  filter((grepl('CATC', junc_seq))) %>% filter(!(grepl('AC', junc_seq))) %>% #ES20-21 unmod, ES22mod, 23 by JES 
  group_by(KD) %>% summarise(axn19_22_23 = sum(AvgNm), unq19_22_23 = length(junc_seq))
g1ess15_19_22_23

g1ess15_19_20_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  filter((grepl('GTTTAGA', junc_seq))) %>% #ES16-18 unmod
  filter(!(grepl('CG', junc_seq))) %>% filter(!(grepl('ATC',junc_seq))) %>% #ES19 mod, ES20 mod, ES21 by JES
  filter((grepl('ACA', junc_seq))) %>% # ES21-22 unmod
  group_by(KD) %>% summarise(axn19_20_23 = sum(AvgNm), unq19_20_23  = length(junc_seq))
g1ess15_19_20_23

g1ess15_20_22_21 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==7)) %>%
  filter((grepl('CGTTTAGA', junc_seq))) %>% #ES16-19 unmod
  filter(!(grepl('ATC', junc_seq))) %>% filter(!(grepl('CA', junc_seq))) %>% #ES20mod, ES21mod
  group_by(KD) %>% summarise(axn20_22_21 = sum(AvgNm), unq20_22_21 = length(junc_seq))
g1ess15_20_22_21

g1ess15_20_21_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  filter((grepl('CGTTTAGA', junc_seq))) %>% #ES16-19 unmod
  filter(!(grepl('ATC', junc_seq))) %>% filter(!(grepl('CA', junc_seq))) %>% #ES20mod, ES21mod
  filter((grepl('AC', junc_seq))) %>% #ES22 unmod, ES23 by JES
  group_by(KD) %>% summarise(axn20_21_23 = sum(AvgNm), unq20_21_23 = length(junc_seq))
g1ess15_20_21_23

g1ess15_20_22_23 <- RPS12allg1avg %>% filter((edit_stop==15&junc_len==8)) %>%
  filter((grepl('CGTTTAGA', junc_seq))) %>% #ES16-19 unmod
  filter(!(grepl('ATC', junc_seq))) %>% filter((grepl('CA', junc_seq))) %>% #ES20mod, ES21unmod
  filter(!(grepl('ACA', junc_seq))) %>% #ES22 mod
  group_by(KD) %>% summarise(axn20_22_23 = sum(AvgNm), unq20_22_23 = length(junc_seq))
g1ess15_20_22_23

## make a single table of the major three action sequences

g1ess15tripA <- full_join(g1ess15_16_17_18, g1ess15_16_17_19, by = "KD")
g1ess15tripB <- full_join(g1ess15tripA, g1ess15_16_17_20, by = "KD")
g1ess15tripC <- full_join(g1ess15tripB, g1ess15_16_17_21, by = "KD")
g1ess15tripD <- full_join(g1ess15tripC, g1ess15_16_17_22, by = "KD")
g1ess15tripE <- full_join(g1ess15tripD, g1ess15_16_17_23, by = "KD")
g1ess15tripF <- full_join(g1ess15tripE, g1ess15_16_18_19, by = "KD")
g1ess15tripG <- full_join(g1ess15tripF, g1ess15_16_18_20, by = "KD")
g1ess15tripH <- full_join(g1ess15tripG, g1ess15_16_18_21, by = "KD")
g1ess15tripI <- full_join(g1ess15tripH, g1ess15_16_18_22, by = "KD")
g1ess15tripJ <- full_join(g1ess15tripI, g1ess15_16_18_23, by = "KD")
g1ess15tripK <- full_join(g1ess15tripJ, g1ess15_16_20_19, by = "KD")
g1ess15tripL <- full_join(g1ess15tripK, g1ess15_16_19_21, by = "KD")
g1ess15tripM <- full_join(g1ess15tripL, g1ess15_16_19_22, by = "KD")
g1ess15tripN <- full_join(g1ess15tripM, g1ess15_16_19_23, by = "KD")
g1ess15tripO <- full_join(g1ess15tripN, g1ess15_16_20_21, by = "KD")
g1ess15tripP <- full_join(g1ess15tripO, g1ess15_16_20_22, by = "KD")
g1ess15tripQ <- full_join(g1ess15tripP, g1ess15_16_20_23, by = "KD")
g1ess15tripR <- full_join(g1ess15tripQ, g1ess15_16_21_22, by = "KD")
g1ess15tripS <- full_join(g1ess15tripR, g1ess15_16_21_23, by = "KD")
g1ess15tripT <- full_join(g1ess15tripS, g1ess15_16_22_23, by = "KD")
g1ess15tripU <- full_join(g1ess15tripT, g1ess15_17_18_19, by = "KD")
g1ess15tripV <- full_join(g1ess15tripU, g1ess15_18_20_17, by = "KD")
g1ess15tripW <- full_join(g1ess15tripV, g1ess15_18_17_21, by = "KD")
g1ess15tripX <- full_join(g1ess15tripW, g1ess15_18_22_17, by = "KD")
g1ess15tripY <- full_join(g1ess15tripX, g1ess15_18_17_23, by = "KD")
g1ess15tripAa <- full_join(g1ess15tripY, g1ess15_17_19_20, by = "KD")
g1ess15tripAb <- full_join(g1ess15tripAa, g1ess15_17_19_21, by = "KD")
g1ess15tripAc <- full_join(g1ess15tripAb, g1ess15_17_19_22, by = "KD")
g1ess15tripAd <- full_join(g1ess15tripAc, g1ess15_17_19_23, by = "KD")
g1ess15tripAe <- full_join(g1ess15tripAd, g1ess15_17_20_21, by = "KD")
g1ess15tripAf <- full_join(g1ess15tripAe, g1ess15_20_22_17, by = "KD")
g1ess15tripAg <- full_join(g1ess15tripAf, g1ess15_17_20_23, by = "KD")
g1ess15tripAh <- full_join(g1ess15tripAg, g1ess15_17_21_22, by = "KD")
g1ess15tripAi <- full_join(g1ess15tripAh, g1ess15_17_21_23, by = "KD")
g1ess15tripAj <- full_join(g1ess15tripAi, g1ess15_17_22_23, by = "KD")
g1ess15tripAk <- full_join(g1ess15tripAj, g1ess15_18_19_20, by = "KD")
g1ess15tripAl <- full_join(g1ess15tripAk, g1ess15_18_19_21, by = "KD")
g1ess15tripAm <- full_join(g1ess15tripAi, g1ess15_18_22_19, by = "KD")
g1ess15tripAn <- full_join(g1ess15tripAm, g1ess15_18_19_23, by = "KD")
g1ess15tripAo <- full_join(g1ess15tripAn, g1ess15_18_20_21, by = "KD")
g1ess15tripAp <- full_join(g1ess15tripAo, g1ess15_18_20_22, by = "KD")
g1ess15tripAq <- full_join(g1ess15tripAp, g1ess15_18_20_23, by = "KD")
g1ess15tripAr <- full_join(g1ess15tripAq, g1ess15_18_22_21, by = "KD")
g1ess15tripAs <- full_join(g1ess15tripAr, g1ess15_18_21_23, by = "KD")
g1ess15tripAt <- full_join(g1ess15tripAs, g1ess15_18_22_23, by = "KD")
g1ess15tripAu <- full_join(g1ess15tripAt, g1ess15_19_20_21, by = "KD")
g1ess15tripAv <- full_join(g1ess15tripAu, g1ess15_20_22_19, by = "KD")
g1ess15tripAw <- full_join(g1ess15tripAv, g1ess15_19_20_23, by = "KD")
g1ess15tripAx <- full_join(g1ess15tripAw, g1ess15_19_21_22, by = "KD")
g1ess15tripAy <- full_join(g1ess15tripAx, g1ess15_19_21_23, by = "KD")
g1ess15tripAz <- full_join(g1ess15tripAy, g1ess15_19_22_23, by = "KD")
g1ess15tripBa <- full_join(g1ess15tripAz, g1ess15_20_22_21, by = "KD")
g1ess15tripBb <- full_join(g1ess15tripBa, g1ess15_20_21_23, by = "KD")
g1ess15tripAxn <- full_join(g1ess15tripBb, g1ess15_20_22_23, by = "KD")
glimpse(g1ess15tripAxn)

g1ess15tripATalla <- g1ess15tripAxn %>% select(KD, contains("axn")) %>%
  gather("Axn1", "AvgNm", contains("axn"))
g1ess15tripATallb <-  g1ess15tripAxn %>% select(KD, contains("unq")) %>%
  gather("Unq1", "UnqCt", contains("unq"))
g1ess15tripATall <- bind_cols(g1ess15tripATalla, g1ess15tripATallb[,2:3])
glimpse(g1ess15tripATall)

g1ess15tripAPerc <- g1ess15tripATall %>% group_by(KD) %>% 
  mutate(Tot = sum(AvgNm,na.rm=TRUE)) %>%
  rowwise() %>% mutate(Perc = 100*(AvgNm/Tot)) %>% select(-Unq1) %>% arrange(KD, Axn1)
head(g1ess15tripAPerc)
View(g1ess15tripAPerc)
#order KDs as desired
target <- c("AvgUn", "GAP1", "81704160", "81802c", "TbRGG2")
g1ess15tripAPerc$KD <- reorder.factor(g1ess15tripAPerc$KD, new.order = target)
g1ess15tripAPercO <- g1ess15tripAPerc %>% arrange(KD)
View(g1ess15tripAPercO)

#make the two option one paper friendly
target <- c("AvgUn", "GAP1", "81704160", "81802c", "TbRGG2")
g1ess15twoAPerc$KD <- reorder.factor(g1ess15twoAPerc$KD, new.order = target)
g1ess15twoAPercO <- g1ess15twoAPerc %>% arrange(KD)
View(g1ess15twoAPercO)

View(g1ess15singlePerc)

# write all tables to files for reference
write.table(g1ess15singlePerc, "RPS12gRNA1singleAxn.csv", sep= ",")
write.table(g1ess15twoAPercO, "RPS12gRNA1twoAxn.csv", sep= ",")
write.table(g1ess15tripAPercO, "RPS12gRNA1tripleAxn.csv", sep= ",")

#########################################################################
#########################################################################

##  Examination of long junction sequences (Fig. 9)
##  mostly pre-edited versus mostly fully edited 

LJsites <- RPS12allCl %>% filter(edit_stop %in% c(9, 12, 15, 16, 19), junc_len>50, KD!=29-13)
glimpse(LJsites) # the c(a, b, c, ....) are ESS that are of interest (EPS w long junctions)

prerunES20toES48 <- c("AGATTTGGGTGGGGGGAACCCTTTGTTTTGGTTAAAGAAACA")
fullrunES20toES48 <- c("AGAGGGTGGTGGTTTTGTTGATTTACCCGGTGTAAAGTATTATACA")

LJsitePRE <- LJsites %>% filter((grepl(prerunES20toES48, junc_seq)))
LJsiteFE <- LJsites %>% filter((grepl(fullrunES20toES48, junc_seq)))

SumPREin <- LJsitePRE %>% filter(TET==TRUE) %>% group_by(edit_stop, KD) %>% 
  summarise(AVG = sum(norm_count)/2) %>%  mutate(BP = "PRE")

SumFEin <- LJsiteFE %>% filter(TET==TRUE) %>% group_by(edit_stop, KD) %>% 
  summarise(AVG = sum(norm_count)/2) %>%  mutate(BP = "FULL")

SumPREun <- LJsitePRE %>% filter(TET==FALSE) %>% group_by(edit_stop) %>% 
  summarise(AVG = sum(norm_count)/8) %>%  mutate(BP = "PRE") %>% mutate(KD = "AvgUn")

SumFEun <- LJsiteFE %>% filter(TET==FALSE) %>% group_by(edit_stop) %>% 
  summarise(AVG = sum(norm_count)/8) %>%  mutate(BP = "FULL") %>% mutate(KD = "AvgUn")

SumLJin <- bind_rows(SumPREin, SumFEin)
glimpse(SumLJin)

SumLJun <- bind_rows(SumPREun, SumFEun)
glimpse(SumLJun)

SumLJboth <- bind_rows(SumLJin, SumLJun)
target <- c("AvgUn", "GAP1", "81704160", "81802c", "TbRGG2")
SumLJboth$KD <- reorder.factor(SumLJboth$KD, new.order = target)
SumLJbothO <- SumLJboth %>% arrange(KD)

#make a graph of this...

barLJ <- ggplot(SumLJbothO, aes(x=KD, y=AVG, fill=BP)) + geom_bar(stat = "identity", position = "dodge") +
  facet_grid(edit_stop~., scales = "free")
barLJ

propLJ <- ggplot(SumLJbothO, aes(x=KD, y=AVG, fill=BP)) + geom_bar(stat = "identity", position = "fill") +
  facet_grid(edit_stop~., scales = "free")
propLJ

SumLJbothOsm <- SumLJbothO %>% filter(edit_stop %in% c(9,12))
barJLsmall <- ggplot(SumLJbothOsm, aes(x=KD, y=AVG, fill=BP)) + geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~edit_stop, scales = "free")
barJLsmall


## just to check that this trend is consistent in both replicates, graph them separately here :)

SumPREinR <- LJsitePRE %>% filter(TET==TRUE) %>% group_by(REP1, edit_stop, KD) %>% 
  summarise(AVG = sum(norm_count)) %>%  mutate(BP = "PRE")

SumFEinR <- LJsiteFE %>% filter(TET==TRUE) %>% group_by(REP1, edit_stop, KD) %>% 
  summarise(AVG = sum(norm_count)) %>%  mutate(BP = "FULL")

SumPREunR <- LJsitePRE %>% filter(TET==FALSE) %>% group_by(REP1, edit_stop) %>% 
  summarise(AVG = sum(norm_count)/4) %>%  mutate(BP = "PRE") %>% mutate(KD = "AvgUn")

SumFEunR <- LJsiteFE %>% filter(TET==FALSE) %>% group_by(REP1, edit_stop) %>% 
  summarise(AVG = sum(norm_count)/8) %>%  mutate(BP = "FULL") %>% mutate(KD = "AvgUn")

SumLJinR <- bind_rows(SumPREinR, SumFEinR)
glimpse(SumLJinR)

SumLJunR <- bind_rows(SumPREunR, SumFEunR)
glimpse(SumLJunR)

SumLJbothR <- bind_rows(SumLJinR, SumLJunR)

SumLJbothR$KD <- reorder.factor(SumLJbothR$KD, new.order = target)
SumLJbothOR <- SumLJbothR %>% arrange(KD)

SumLJbothORsm <- SumLJbothOR %>% filter(edit_stop %in% c(9,12))
barJLsmall <- ggplot(SumLJbothORsm, aes(x=KD, y=AVG, fill=BP)) + geom_bar(stat = "identity", position = "dodge") +
  facet_grid(REP1~edit_stop, scales = "free")
barJLsmall


