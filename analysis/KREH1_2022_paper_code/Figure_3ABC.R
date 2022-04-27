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

library(dplyr)
library(tidyr)
library(ggplot2)

##### Figure 3A #####

# Load the files you're going to need
KREH1KO <- read.csv(file.choose(), sep = ",", header=TRUE)
A63WT <- read.csv(file.choose(), sep = ",", header=TRUE)

# Change WT knock_down column to character class
# (R interprets the WT 29-13 cell line name as an integer)
A63WT$knock_down <- as.character(A63WT$knock_down)
# Then set KO as Tet=TRUE
KREH1KO$tetracycline <- "true"

# Check to make sure the correct files were input
unique(KREH1KO$sample) # should see 2 samples
unique(A63WT$sample) # should see 5 samples

# Combine tables into one
A6 <- bind_rows(KREH1KO, A63WT)
# Set editing sites corresponding to Pre-edited, Fully edited reads
pre <- 24
fe <- 115

# Remove pre-edited and fully edited transcripts and set tetracycline column as logical
A6nopre <- A6 %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe))

# Cut and total JL counts in each sample
A6.cut <- A6nopre %>% 
  mutate(bin = cut(junc_len,
                   breaks = c(-Inf,1,11,21,31,41,51,61,71,81,91,Inf),
                   right=FALSE)) %>%
  group_by(sample, bin) %>% mutate(total=sum(norm_count)) %>%
  select(sample, knock_down, replicate, tetracycline, total, bin) %>%
  distinct(.keep_all = TRUE)

# Calculate mean and sd for -Tet and +Tet JLs
A6.stats <- A6.cut %>% group_by(tetracycline,bin) %>%
  summarise(mean=mean(total),stdev=sd(total))

# Making the graph
A6.stats[A6.stats$tetracycline=="false",1] <- "WT"
A6.stats[A6.stats$tetracycline=="true",1] <- "KREH1 KO"
A6.stats$tetracycline <- factor(A6.stats$tetracycline, levels=c('WT', 'KREH1 KO'))
Palette <- c("#bdbdbd", "#FF0000")

ggplot(A6.stats, aes(x=bin, y=mean, fill=tetracycline)) +
  geom_bar(aes(fill=tetracycline), position="dodge", stat="identity") +
  ggtitle("A6") +
  ylab("Total sequences (Average)") + xlab("Junction Length") +
  theme_bw() +
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

A6plot2 <- A6.stats %>% filter(bin=="[51,61)" | bin=="[61,71)" |
                               bin=="[71,81)" | bin=="[91, Inf)" | bin=="[81,91)")
ggplot(A6plot2, aes(x=bin, y=mean, fill=tetracycline)) +
  geom_bar(aes(fill=tetracycline), position="dodge", stat="identity") +
  ggtitle("A6") +
  ylab("Total sequences (Average)") + xlab("Junction Length") +
  theme_bw() +
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

A6test <- A6.cut %>% 
  group_by(knock_down, bin) %>% 
  summarise(total = list(total)) %>% 
  spread(knock_down, total) %>% 
  group_by(bin) %>% 
  mutate(p_value = t.test(unlist(2913), unlist(REH1KO), var.equal = TRUE)$p.value)
A6test <- A6test[order(A6test$p_value),]
A6test$adj <- p.adjust(A6test$p_value, method = "BY")


##### Figure 3B #####

# This code will identify the list of reads with editing stop site 24 as well as
# their average counts and fold changes between WT and KREH1 KO samples.

library(dplyr)
library(tidyr)
library(ggplot2)

# Load the files you're going to need
A63WT <- read.csv(file.choose(), sep = ",", header=TRUE)
KREH1KO <- read.csv(file.choose(), sep = ",", header=TRUE)

# Check to make sure the correct files were input
unique(A63WT$sample) # should see 5 samples
unique(KREH1KO$sample) # should see 2 samples

# Remove pre-edited sequences and filter for editing stop site 24
A63WT24 <- A63WT %>% filter(!(edit_stop==24&junc_len==0)) %>%
  filter(edit_stop==24)
KREH1KO24 <- KREH1KO %>% filter(!(edit_stop==24&junc_len==0)) %>%
  filter(edit_stop==24)

# Calculate the average counts for each unique sequence
A63WTavg <- A63WT24 %>% mutate(knock_down="WT") %>%
  group_by(tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/5) %>% 
  select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE)
KREH1KOavg <- KREH1KO24 %>% mutate(knock_down="KREH1 KO") %>%
  group_by(tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/2) %>% 
  select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE)

# Calculate fold changes and sort by normcount in KO
A624fc <- left_join(KREH1KOavg, A63WTavg, by=c('edit_stop','junc_seq')) %>%
  mutate(FC=AvgNm.x/AvgNm.y) %>%
  select(edit_stop, junc_seq, junc_len.x, junc_end.x, AvgNm.x, AvgNm.y, FC)
colnames(A624fc) <- c("edit_stop", "junc_seq", "junc_len", "junc_end", "KOavg", "WTavg", "FC")
A624fc <- arrange(A624fc, desc(KOavg))


##### Figure 3C #####
library(dplyr)
library(tidyr)
library(ggplot2)

# Load the files you're going to need
A63WT <- read.csv(file.choose(), sep = ",", header=TRUE)
KREH1KO <- read.csv(file.choose(), sep = ",", header=TRUE)

# Check to make sure the correct files were input
unique(A63WT$sample) # should see 5 samples
unique(KREH1KO$sample) # should see 2 samples

# Remove pre-edited sequences and
# filter for editing stop site 24 and junction end site > 39
A63WT24 <- A63WT %>% filter(!(edit_stop==24&junc_len==0)) %>%
  filter(edit_stop==24&junc_end>39)
KREH1KO24 <- KREH1KO %>% filter(!(edit_stop==24&junc_len==0)) %>%
  filter(edit_stop==24&junc_end>39)

# Calculate the average counts for each unique sequence
A63WTavg <- A63WT24 %>% mutate(knock_down="WT") %>%
  group_by(tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/5) %>% 
  select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE)
KREH1KOavg <- KREH1KO24 %>% mutate(knock_down="KREH1 KO") %>%
  group_by(tetracycline, edit_stop, junc_seq, junc_end) %>%
  mutate(AvgNm = sum(norm_count)/2) %>% 
  select(knock_down, tetracycline, edit_stop, junc_len, junc_end, junc_seq, AvgNm) %>%
  distinct(., .keep_all = TRUE)

# Calculate fold changes and sort by normcount in KO and in WT
A624fc <- left_join(KREH1KOavg, A63WTavg, by=c('edit_stop','junc_seq')) %>%
  mutate(FCKO=AvgNm.x/AvgNm.y, FCWT=AvgNm.y/AvgNm.x) %>%
  select(edit_stop, junc_seq, junc_len.x, junc_end.x, AvgNm.x, AvgNm.y, FCKO, FCWT)
colnames(A624fc) <- c("edit_stop", "junc_seq", "junc_len", "junc_end", "KOavg", "WTavg", "FCKO", "FCWT")
A624fc <- arrange(A624fc, desc(KOavg))
write.table(A624fc, "A6REH1KOESS24JL39fc.csv", sep=",", row.names=FALSE)
