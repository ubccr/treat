# IPS -----
#install.packages("stringr")
library('stringr')

#install.packages("tidyr")
library('tidyr')

#install.packages("dplyr")
library('dplyr')

#install.packages("ggplot2")
library('ggplot2')

#install.packages("gdata")
library('gdata')

#install.packages("grid")
library('grid')

# data from TREAT: WT CyB 5 reps
treatdata <-read.table(file.choose(), header = T, sep = ",")

#load in the number of the editing stop sites associated with pre and fully edited
pre <- 558
fe <- 582

# Remove pre-edited and fully edited transcripts and setting tetracycline column as logical
# change the two edit_stop== inside of the filter() function to equal the editing stop site 
# that corresponds with the transcript you are examining
treatnopre <- treatdata %>% 
  filter(!(junc_len==0&edit_stop==pre),!(junc_len==0&edit_stop==fe)) %>%
  mutate(tetracycline = as.logical(tetracycline))

## to determine the total norm count of seqences at each editing stop site
##  this code is needed for both IPS and EPS determination
ESSct <-treatnopre %>% group_by(sample, edit_stop) %>% 
  mutate(nm = sum(norm_count))
ESSctbyES <- ESSct %>% select(sample, edit_stop, tetracycline, replicate, knock_down, nm) %>% 
  distinct(.keep_all = TRUE)

## Then write the file you wish to send to Runpu as follows:
#  IF NO ADDED CONTROLS: 
#write.table(ESSctbyES, "CyBESSforRunpu100917.csv", sep=",", row.names=FALSE)

## To determine the IPSs
## load the outlier threshold function
oThresh<-function(x){
  + IQR(x)*1.5 + quantile(x,prob=0.75)}

oThreshR1 <- ESSctbyES %>% group_by(sample) %>% mutate(oThr = oThresh(nm)) 
# samples each get their own threshold
IPStable<-oThreshR1 %>% group_by(sample) %>%
  mutate(IPS = ifelse(nm>oThr,edit_stop, "FALSE"))
IPStf <- oThreshR1 %>% group_by(sample, edit_stop) %>% 
  mutate(IPS = ifelse(nm>oThr,TRUE, FALSE)) %>%
  select(sample, edit_stop, tetracycline, replicate, knock_down, IPS)  %>% 
  distinct(.keep_all = TRUE)
write.table(IPStable, "CyBIPStableWT.csv", sep=",", row.names=FALSE)

# IPS graphs -----
IPSpeaks <- IPStable %>% filter(IPS != "FALSE") %>%
    mutate(difference = nm-oThr)

ggplot(IPSpeaks, aes(x=replicate, y=difference, fill=IPS)) +
    geom_bar(stat="identity") +
    xlab("WT Sample") + ylab("Difference above threshold") +
    ggtitle("CyB IPS peak levels") +
    facet_grid(.~IPS)

CyBIPSsum <- IPSpeaks %>% group_by(IPS) %>%
    summarise(mean=mean(difference), stdv=sd(difference))

ggplot(CyBIPSsum, aes(x=IPS, y = mean)) +
    geom_bar(stat="identity", fill="#2ca25f") +
    geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv),
                  width=.2, size=.75,
                  position=position_dodge(.9)) +
    xlab("IPS") + ylab("Difference above threshold") +
    ggtitle("CyB IPS peak levels")

## line graph of each ES with mean, standard deviation of 5 reps
CyBWTESline <- ESSctbyES %>% group_by(edit_stop) %>%
    summarise(mean=mean(nm),stdv=sd(nm))

ggplot(CyBWTESline, aes(x=edit_stop, y=mean)) +
    geom_line(stat="identity") +
    geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv), size=.75, width=.1) +
    xlab("edit stop") + ylab("normcount") +
    ggtitle("normcount levels") + scale_x_reverse()

# MJES ------
## The following code is analagous to the IPS/EPS code but for junction ends
## It will determine major junction ends and exacerbated junction ends
## to determine the total norm count of seqences at each editing stop site
##  this code is needed for both MJES and EJE determination

JESct <-treatnopre %>% group_by(sample, junc_end) %>% mutate(nm = sum(norm_count))
JESctbyES <- JESct %>% select(sample, junc_end, tetracycline, replicate, knock_down, nm) %>% 
  distinct(.keep_all = TRUE)

## Then write the file you wish to send to Runpu as follows:
#  IF NO ADDED CONTROLS: 
#write.table(JESctbyES, "CyBJESforRunpu100917.csv", sep=",", row.names=FALSE)

## To determine the MJESs
oThreshR1 <- JESctbyES %>% group_by(sample) %>% mutate(oThr = oThresh(nm))
MJEStable<-oThreshR1 %>% group_by(sample) %>% 
  mutate(MJES = ifelse(nm>oThr,junc_end, "FALSE"))
MJEStf <- oThreshR1 %>% group_by(sample) %>% 
  mutate(MJES = ifelse(nm>oThr,TRUE, FALSE)) %>%
  mutate(edit_stop = junc_end) %>% 
  select(sample, edit_stop, tetracycline, replicate, knock_down, MJES) %>% 
  distinct(.keep_all = TRUE)
write.table(MJEStable, "CyB-WT-MJEStable.csv", sep=",", row.names=FALSE)

ipsmje <- inner_join(IPStf, MJEStf)
write.table(ipsmje, "CyB-WT-IPSMJEStable.csv", sep=",", row.names=FALSE)

# MJLs -------------------------------
# This section of code calculates major junction lengths and exacerbated junction lengths
JLct <-treatnopre %>% group_by(sample, junc_len) %>% mutate(nm = sum(norm_count))
JLctbyES <- JLct %>% select(sample, junc_len, tetracycline, replicate, knock_down, nm) %>% 
  distinct(.keep_all = TRUE)

## Then write the file you wish to send to Runpu as follows:
#  IF NO ADDED CONTROLS: 
#write.table(JLctbyES, "CyBJLforRunpu100917.csv", sep=",", row.names=FALSE)

## To determine the MJLs
oThreshR1 <- JLctbyES %>% group_by(sample) %>% mutate(oThr = oThresh(nm))
MJLtable<-oThreshR1 %>% group_by(sample) %>% 
  mutate(MJL = ifelse(nm>oThr,junc_len, "FALSE")) %>% 
  distinct(.keep_all = TRUE)

write.table(MJLtable, "CyB-WT-MJLtable.csv", sep=",", row.names=FALSE)



# k

