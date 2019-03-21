#Look at the spring blitz versus the fall blitz


library(tidyverse)
library(lubridate)
library(reshape2)
library(vegan)
library(pwr)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(visreg)
library(readxl)
library(simr)


source("blitzclean.R")


#subset teh samples from both blitzes in 2018
bugsblitzSF <- filter(inverts2, year(Date) >= 2018 & 
                        ((month(Date) >= 3 & month(Date)<=5) | 
                           (month(Date) >= 10 & month(Date)<=11)))

#check to see if we haven't assigned targets for any samples
#bugsx = summarize(group_by(bugsblitzSF, SampleID, Station, Region, Volume, Sampletype), total = sum(TotalCount, na.rm = T))
#test = merge(bugsx, targets, all.x = T)
#test2 = test[which(is.na(test$Target)),]


#add some more info
targets <- read.csv("blitz2/targets.csv")

bugsblitzSF = merge(bugsblitzSF, targets, by = "SampleID")

#We didn't have neuston or benthic cores in the fall, and only a subset of the sites.
bugsblitzSF = filter(bugsblitzSF, Target != "benthic" & Target != "neuston")

#import site types
sitetypes = read_excel("blitz2/sites.xlsx")
bugsblitzSF = merge(bugsblitzSF, sitetypes[,c(2,6,7)])

#filter out sampled sites
bugsblitzSF = filter(bugsblitzSF, site == "Prospect" |site == "Ryer" |site == "Winter"| site == "Browns")

bugsxSF = summarize(group_by(bugsblitzSF, SampleID, Station, Region, Volume, Sampletype, Target), total = sum(TotalCount, na.rm = T))



bugsxSF$effort = rep(NA, nrow(bugsxSF))

#effort for trawls ad sweep nets is volume
bugsxSF$effort[which(bugsxSF$Sampletype=="Mysid net"|bugsxSF$Sampletype=="sweepnet")]= 
  bugsxSF$Volume[which(bugsxSF$Sampletype=="Mysid net"|bugsxSF$Sampletype=="sweepnet")]



#attach the efforts to the origional dataset
bugsblitzSF= merge(bugsblitzSF, bugsxSF)
bugsblitzSF$TotalCount[which(is.na(bugsblitzSF$TotalCount))] =0
#Assime 100% subsampled unless otherwise specified
bugsblitzSF$subsampled[which(is.na(bugsblitzSF$subsampled))] =100


# Adjust total count for subsampling
bugsblitzSF$atotal = bugsblitzSF$TotalCount*(1/(bugsblitzSF$subsampled/100))

#Calculate CPUE
bugsblitzSF$CPUE = bugsblitzSF$atotal/bugsblitzSF$effort

#add some analysis groups
zoocodes <- read_excel("blitz2/zoocodes.xlsx")
zoocodes <- zoocodes[,c(3,12,13, 14)]
bugsblitzSF = merge(bugsblitzSF, zoocodes)
bugsblitzSF = bugsblitzSF[order(bugsblitzSF$SampleID),]



#I'm going to go ahead and remove the mesozooplanktoton from all the blitz bugs samples. I may want to make some extra notes 
#on the catch of cladocera in sweep nets, but just to make stuff clearer for now...
bugsblitzSF = filter(bugsblitzSF, Analy2 != "Copepoda" & Analy2 != "Cladocera" & Analy2 != "Ostracoda")


sites = group_by(bugsblitzSF, Station, Region, Region2) %>% summarize(count = length(Station))

#import some info on the sites
sitetypes = read_excel("blitz2/sites.xlsx")
bugsblitzSF = merge(bugsblitzSF, sitetypes[,c(2,6,7)])

#put the sties in order along the estuarine/fresh gradient
bugsblitzSF$site = factor(bugsblitzSF$site, levels = c("Ryer",  "Browns",
                                                   "Winter", "Prospect"))

#Create a new colum where "targets" (habitat type) combines all the vegetation samples into one type.
bugsblitzSF$targets2 = as.character(bugsblitzSF$Target)
bugsblitzSF$targets2[which(bugsblitzSF$Target == "EAV" | bugsblitzSF$Target == "SAV"| bugsblitzSF$Target == "FAV" )] = "sweep net"

#Summarize by sample and calculate the total CPUE of all the macroinvertebrates in the sample
bugsblitzSF.1 = summarize(group_by(bugsblitzSF, SampleID, Station, Region2, Target, targets2, 
                                 Region, Sampletype, site, sitetype, Date),
                        tcount = sum(atotal, na.rm = T), tCPUE = sum(CPUE, na.rm = T), richness = length(unique(Analy)))

bugsblitzSF.1$site = as.factor(bugsblitzSF.1$site)
bugsblitzSF.1$sitetype = as.factor(bugsblitzSF.1$sitetype)
bugsblitzSF.1$targets2 = as.factor(bugsblitzSF.1$targets2)

ggplot(bugsblitzSF, aes(x=site, y = CPUE, fill = Analy2)) + 
  geom_bar(stat = "identity", position = "fill")+
  facet_grid(targets2~., scales = "free", space = "free_x") +
  scale_fill_manual(values = mypal, name = NULL) + 
  xlab("Site")+ ylab("Relative percent composition") + mytheme
