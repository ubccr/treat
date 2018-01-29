library("stringr")
library("tidyr")
library("dplyr")
library("ggplot2")
library(gdata)
library(grid)

CyB.WT.all <-read.csv(file.choose(), sep = ",", header=TRUE)

CyB.prefull  <- CyB.WT.all %>%
    filter((edit_stop==558 & junc_len==0)|(edit_stop==582 & junc_len==0)) %>%
    select(norm_count, replicate, edit_stop) %>%
    mutate(status=ifelse(edit_stop==558, 'Pre-edited', 'Fully edited'))

CyB.prefull$status = factor(CyB.prefull$status, levels=c('Pre-edited','Fully edited'))
CyB.prefull$replicate = factor(CyB.prefull$replicate)

ggplot(CyB.prefull, aes(x=replicate, y=norm_count, fill=status)) +
    geom_bar(stat="identity") +
    xlab("Replicate") + ylab("Norm Count") +
    ggtitle("CyB") + facet_grid(.~status) + theme(legend.position="none")

CyB.pfavg <- CyB.prefull %>% group_by(status) %>%
    summarise(mean=mean(norm_count),stdv=sd(norm_count))

ggplot(CyB.pfavg, aes(x=status, y=mean, fill=status)) +
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv), size=.75, width=.1) +
    xlab("edit stop") + ylab("normcount") +
    ggtitle("CyB prefully levels") + theme(legend.position="none")

ggplot(CyB.pfavg, aes(x=status, y=mean, fill=status)) + 
    geom_bar(stat="identity") +
    geom_point(data=CyB.prefull, aes(x=status, y=norm_count), size=2) +
    xlab("Transcript type") + ylab("Normalized count") +
    ggtitle("CYb pre-edited and fully edited transcript levels")  +
    theme_light(base_size = 11, base_family = "") + theme(legend.position="none")
