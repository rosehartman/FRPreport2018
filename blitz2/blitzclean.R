#A start on the macroinvertebrate "bliz" analyses from 2018 

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

#Upload the data
#First load the function I wrote to query the FRP database
source("querydatabase.R")

#Now specify the path to the FRP database
path = "U:/FRPA/MONITORING/Labs/Databases/FRPdata28DEC2018.accdb"

#Query the invertebrate data
inverts2 <- GetFRPdata(path, type = "inverts")

#Query the station information
stations = GetFRPdata(path, type = "stations")

#attatch station information
inverts2 = merge(inverts2, stations[,c(2,7)], all.x = T)

#assume that if "subsampled" is left blankn, it's 100%
inverts2$subsampled[which(is.na(inverts2$subsampled))]= 100

#convert the date from factor to date
inverts2$Date <- ymd(inverts2$Date)

#clean up the 'sample types' column
inverts2$Sampletype = as.character(inverts2$Sampletype)
inverts2$Sampletype[which(inverts2$Sampletype == "Mysid net (no sled)" | 
                            inverts2$Sampletype == "larval sled oblique"|
                            inverts2$Sampletype == "larval sled benthic")]=
  "Mysid net"

#subset teh samples from 2017 and 2018
bugsblitz <- filter(inverts2, year(Date) >= 2017 & (month(Date) >= 3 & month(Date)<=5))

#check to see if we haven't assigned targets for any samples
#bugsx = summarize(group_by(bugsblitz, SampleID, Station, Region, Volume, Sampletype), total = sum(TotalCount, na.rm = T))
#test = merge(bugsx, targets, all.x = T)
#test2 = test[which(is.na(test$Target)),]


#add some more info
targets <- read.csv("blitz2/targets.csv")

bugsblitz = merge(bugsblitz, targets, by = "SampleID")

bugsx = summarize(group_by(bugsblitz, SampleID, Station, Region, Volume, Sampletype, Target), total = sum(TotalCount, na.rm = T))


bugsx$effort = rep(NA, nrow(bugsx))

#effort for trawls and sweepnets is volume
bugsx$effort[which(bugsx$Sampletype=="Mysid net")]= bugsx$Volume[which(bugsx$Sampletype=="Mysid net")]
bugsx$effort[which(bugsx$Sampletype=="neuston trawl")]= bugsx$Volume[which(bugsx$Sampletype=="neuston trawl")]
bugsx$effort[which(bugsx$Sampletype=="sweepnet")]= bugsx$Volume[which(bugsx$Sampletype=="sweepnet")]


#effort for benthic cores and grabs is per unit area

bugsx$effort[which(bugsx$Sampletype=="Ponar grab")]=
  rep(0.0231, nrow(bugsx[which(bugsx$Sampletype=="Ponar grab"),]))

bugsx$effort[which(bugsx$Sampletype=="PCV core")]=
  rep(0.00811, nrow(bugsx[which(bugsx$Sampletype=="PCV core"),]))

#everything else has an effort of "1"

bugsx$effort[which(is.na(bugsx$effort))] = 1

#attach the efforts to the origional dataset
bugsblitz= merge(bugsblitz, bugsx)
bugsblitz$TotalCount[which(is.na(bugsblitz$TotalCount))] =0
bugsblitz$subsampled[which(is.na(bugsblitz$subsampled))] =100


# Adjust total count for subsampling
bugsblitz$atotal = bugsblitz$TotalCount*(1/(bugsblitz$subsampled/100))

#Calculate CPUE
bugsblitz$CPUE = bugsblitz$atotal/bugsblitz$effort

#add some analysis groups
zoocodes <- read_excel("blitz2/zoocodes.xlsx")
zoocodes <- zoocodes[,c(3,12,13, 14)]
bugsblitz = merge(bugsblitz, zoocodes)
bugsblitz = bugsblitz[order(bugsblitz$SampleID),]



#I'm going to go ahead and remove the mesozooplanktoton from all the blitz bugs samples. I may want to make some extra notes 
#on the catch of cladocera in sweep nets, but just to make stuff clearer for now...
bugsblitz = filter(bugsblitz, Analy2 != "Copepoda" & Analy2 != "Cladocera" & Analy2 != "Ostracoda")


sites = group_by(bugsblitz, Station, Region, Region2) %>% summarize(count = length(Station))

#import some info on the sites
sitetypes = read_excel("blitz2/sites.xlsx")
bugsblitz = merge(bugsblitz, sitetypes[,c(2,6,7)])

#put the sties in order along the estuarine/fresh gradient
bugsblitz$site = factor(bugsblitz$site, levels = c("Ryer", "Grizzly", "Tule Red", "Blacklock", "Bradmoor", "LHB", "Browns",
                                                     "Winter", "Broad", "Horseshoe", "Decker", "Stacys", "Miner", "Prospect",
                                                     "Liberty", "Lindsey", "Flyway"))

#Create a new colum where "targets" (habitat type) combines all the vegetation samples into one type.
bugsblitz$targets2 = as.character(bugsblitz$Target)
bugsblitz$targets2[which(bugsblitz$Target == "EAV" | bugsblitz$Target == "SAV"| bugsblitz$Target == "FAV" )] = "sweep net"


#remove samples from Dow Wetlands, since some of it's wierd and it isn't an FRP site anymore,
#plus some that weren't sampled during the Blitz
bugsblitz = filter(bugsblitz, site != "Dow" & site != "Honker" & site != "Sherman")


#Summarize by sample and calculate the total CPUE of all the macroinvertebrates in the sample
bugsblitz.1 = summarize(group_by(bugsblitz, SampleID, Station, Region2, Target, targets2, 
                                Region, Sampletype, site, sitetype, Date),
                       tcount = sum(atotal, na.rm = T), tCPUE = sum(CPUE, na.rm = T), richness = length(unique(Analy)))

bugsblitz.1$site = as.factor(bugsblitz.1$site)
bugsblitz.1$sitetype = as.factor(bugsblitz.1$sitetype)
bugsblitz.1$targets2 = as.factor(bugsblitz.1$targets2)

#summarize by Region and habitat type so I can make a quick plot

#add up all the critters that are in the same analysis group
bugssum = summarize(group_by(bugsblitz, site, targets2, 
                             Region2, Station, sitetype, Sampletype, SampleID, Analy), 
                    tCPUE = sum(CPUE))

#I need to add the zeros in for taxa that were not present in a sample so I can take an appropriate mean
#I usually do this by transfering it from "long" to "wide" and back again.

#Make it wide
ComMat.1 = dcast(bugsblitz, formula = SampleID~Analy2, value.var="CPUE", 
                 fun.aggregate = sum, fill = 0)
#do it by smaller grouping
ComMat.2 = dcast(bugsblitz, formula = SampleID~Analy, value.var="CPUE", 
                 fun.aggregate = sum, fill = 0)

#do it by common name
ComMat.CN = dcast(bugsblitz, formula = SampleID~CommonName, value.var="CPUE", 
                 fun.aggregate = sum, fill = 0)

#do it by cleaned up common name
ComMat.CN2 = dcast(bugsblitz, formula = SampleID~CN, value.var="CPUE", 
                  fun.aggregate = sum, fill = 0)

#Go back to long
bugsblitz.2x = melt(ComMat.1, id.vars = "SampleID", variable.name = "Analy2", value.name = "CPUE")
bugsblitz.2x = merge(bugsblitz.1, bugsblitz.2x)

#Means by analysis group for each location
bugssum.1 = summarize(group_by(bugsblitz.2x, site, Region2, sitetype, targets2, 
                               Region, Sampletype, Analy2), mCPUE = mean(CPUE),
                      sdCPUE = sd(CPUE), seCPUE = sd(CPUE)/length(CPUE))

#Now calculate total CPUE
bugstot = summarize(group_by(bugsblitz.2x, SampleID, Date, site, sitetype, Target, targets2, 
                             Region, Region2, Sampletype), tCPUE = sum(CPUE, na.rm = T))


bugstot$Year = year(bugstot$Date)
bugstot$Year2 = as.factor(bugstot$Year)
bugstot$logtot = log(bugstot$tCPUE +1)
bugstotNoB = filter(bugstot, targets2!= "benthic")

#Mean CPUE, sample size, and standard error for each location and habitat so I can make a pretty graph.
bugstotave = summarize(group_by(bugstot, site, sitetype, targets2, 
                                 Region2), mCPUE = mean(tCPUE, na.rm = T), 
                       sdCPUE = sd(tCPUE, na.rm = T), seCPUE = sd(tCPUE)/length(tCPUE), N = length(SampleID))

#Some sites didn't have samples for a particular habitat type, so I added rows to make it easier to graph
sitetarget = expand.grid(site=unique(bugstotave$site), targets2 = unique(bugstotave$targets2))
foo = merge(sitetypes[,c(4,6,7)], sitetarget)
foo = unique(foo)

bugstotave2 = merge(foo, bugstotave, all.x = T)
bugstotave2$mCPUE[which(is.na(bugstotave2$mCPUE))] = 0
bugstotave2$N[which(is.na(bugstotave2$N))] = 0
bugstotave2$site = factor(bugstotave2$site, levels = c("Ryer", "Grizzly", "Tule Red", "Blacklock", "Bradmoor", "LHB", "Browns",
                                                   "Winter", "Broad", "Horseshoe", "Decker", "Stacys", "Miner", "Prospect",
                                                   "Liberty", "Lindsey", "Flyway"))


###############################################################################################
#subset by habitat type

ben = filter(bugsblitz, Target == "benthic")
#put in all the zeros for clams
benMat.2 = dcast(ben, formula = SampleID~CN, value.var="CPUE", 
                 fun.aggregate = sum, fill = 0)
#take out the non-clams
benMat.2 = benMat.2[,c("SampleID", "Corbicula", "Potamocorbula", "Clam Other")]
ben.2x = melt(benMat.2, id.vars = "SampleID", variable.name = "CN", value.name = "CPUE")
ben.2x = merge(ben[,c(1,3,4,6,8,26, 34,35,36)], ben.2x, by = "SampleID")
ben.2x = ben.2x[(!duplicated(ben.2x)),]


ben.1 = summarize(group_by(ben.2x, SampleID, 
                                 Region2, site, sitetype, Date),
 tCPUE = sum(CPUE, na.rm = T))

ben.1$site = as.factor(ben.1$site)
ben.1$sitetype = as.factor(ben.1$sitetype)


bensum = summarize(group_by(ben.2x, site, targets2, Region2,
                             sitetype, SampleID, CN, Date), 
                    tCPUE = sum(CPUE))




#######################################################################################

mys = filter(bugsblitz, Target == "mysid")
mys.1 = filter(bugsblitz.1, Target == "mysid")
mys.2x = filter(bugsblitz.2x, Target == "mysid")
myssum.1 = filter(bugssum.1, targets2 == "mysid")

#######################################################################################

neu = filter(bugsblitz, Target == "neuston")
neu.1 = filter(bugsblitz.1, Target == "neuston")
neu.2x = filter(bugsblitz.2x, Target == "neuston")
neusum.1 = filter(bugssum.1, targets2 == "neuston")

#######################################################################################

sweeps = filter(bugsblitz, targets2 == "vegetation")
sweeps.1 = filter(bugsblitz.1, targets2 == "vegetation")
sweeps.2x = filter(bugsblitz.2x, targets2 == "vegetation")
sweepssum.1 = filter(bugssum.1, targets2 == "vegetation")


############################################################################################
#set up some labels and stuff to make plotting easier

#set up a color pallete
mypal = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"), "black")

#set up a theme
mytheme = theme(strip.text.x = element_text(size=14),strip.text.y = element_text(size=18),
                axis.text.x  = element_text(angle=90, hjust=1, vjust = 0.5, size=12),
                axis.title.y  = element_text(size=16),
                axis.title.x  = element_text(size = 16), axis.text.y = element_text(size = 14),
                legend.text = element_text(size = 14),
                legend.title=element_blank(), legend.background =element_blank()) 

LocLab = c("Broad Sl.", "Browns Is.", "Decker Is.", "Dow", "Horseshoe", 
           "Lindsey Sl.", "Miner Sl.", "Prospect Is.", "Ryer Is.", "Tule Red", "Stacys Is.", "Winter Is.")








