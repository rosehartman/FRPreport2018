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


source("blitz2/blitzclean.R")


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
sitetypes = read_excel("blitz2/Stations2.xlsx")
bugsblitzSF = merge(bugsblitzSF, sitetypes[,c(2,9,10,11)])

#filter out sampled sites
bugsblitzSF = filter(bugsblitzSF, site == "Prospect" |site == "Ryer" |site == "Winter"| site == "Browns")

bugsxSF = summarize(group_by(bugsblitzSF, SampleID, Station, Region, Volume, Sampletype, Target), total = sum(TotalCount, na.rm = T))



bugsxSF$effort = rep(NA, nrow(bugsxSF))

#effort for trawls ad sweep nets is volume
bugsxSF$effort[which(bugsxSF$Sampletype=="Mysid net"|bugsxSF$Sampletype=="sweep net")]= 
  bugsxSF$Volume[which(bugsxSF$Sampletype=="Mysid net"|bugsxSF$Sampletype=="sweep net")]



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


#put the sties in order along the estuarine/fresh gradient
bugsblitzSF$site = factor(bugsblitzSF$site, levels = c("Ryer",  "Browns",
                                                   "Winter", "Prospect"))

#Create a new colum where "targets" (habitat type) combines all the vegetation samples into one type.
bugsblitzSF$targets2 = as.character(bugsblitzSF$Target)
bugsblitzSF$targets2[which(bugsblitzSF$Target == "EAV" | bugsblitzSF$Target == "SAV"| bugsblitzSF$Target == "FAV" )] = "sweep net"


#add a new spring versus fall factor
bugsblitzSF$season = rep(NA, nrow(bugsblitzSF))
bugsblitzSF$season[which(month(bugsblitzSF$Date)> 9)] = "fall"
bugsblitzSF$season[which(month(bugsblitzSF$Date)< 9)] = "spring"
bugsblitzSF$season= factor(bugsblitzSF$season, levels = c("spring", "fall"))

#Summarize by sample and calculate the total CPUE of all the macroinvertebrates in the sample
bugsblitzSF.1 = summarize(group_by(bugsblitzSF, SampleID, Station, Region2, Target, targets2, 
                                 Region, Sampletype, site, sitetype, Date, season),
                        tcount = sum(atotal, na.rm = T), tCPUE = sum(CPUE, na.rm = T), richness = length(unique(Analy)))


bugsblitzSF.1$site = as.factor(bugsblitzSF.1$site)
bugsblitzSF.1$sitetype = as.factor(bugsblitzSF.1$sitetype)
bugsblitzSF.1$targets2 = as.factor(bugsblitzSF.1$targets2)

#stacked bar plot of community composition
ggplot(bugsblitzSF, aes(x=site, y = CPUE, fill = Analy2)) + 
  geom_bar(stat = "identity", position = "fill")+
  facet_grid(targets2~season, scales = "free", space = "free_x") +
  scale_fill_manual(values = mypal, name = NULL) + 
  xlab("Site")+ ylab("Relative percent composition") + mytheme

#mean and standard deviation by site and gear type
bugsblitzSF.ave = summarize(group_by(bugsblitzSF.1, Region2, targets2, 
                                     site, sitetype, season),
                            mCPUE = mean(tCPUE, na.rm = T),  sdCPUE = sd(tCPUE, na.rm = T), seCPUE = sd(tCPUE, na.rm = T)/length(tCPUE))


#graph of meanCPUE
ggplot(bugsblitzSF.ave, aes(x=site, y = mCPUE, fill = sitetype)) + geom_bar(stat = "identity") +
  facet_grid(targets2~season, scales = "free") + geom_errorbar(aes(ymin = mCPUE-seCPUE, ymax = mCPUE+seCPUE))

#graph of meanCPUE
ggplot(bugsblitzSF.ave, aes(x=site, y = log(mCPUE), fill = sitetype)) + geom_bar(stat = "identity") +
  facet_grid(targets2~season, scales = "free") 

#hmmmmm.... there was really high catch in mac6-22mar2018. that might make everything hard to interpret.
#lets take it out and see what happens
#mean and standard deviation by site and gear type
foo = summarize(group_by(filter( bugsblitzSF.1, SampleID!="MAC6-22MAR2018"), 
                                      Region2, targets2, site, sitetype, season),
                            mCPUE = mean(tCPUE, na.rm = T), mlogCPUE = mean(log(tCPUE), na.rm = T),  sdCPUE = sd(tCPUE, na.rm = T), 
                            seCPUE = sd(tCPUE, na.rm = T)/length(tCPUE), selCPUE = sd(log(tCPUE), na.rm = T)/length(tCPUE))


#graph of meanCPUE
ggplot(foo, aes(x=site, y = mCPUE, fill = sitetype)) + geom_bar(stat = "identity") +
  facet_grid(targets2~season, scales = "free") + 
  geom_errorbar(aes(ymin = mCPUE-seCPUE, ymax = mCPUE+seCPUE))

ggplot(foo, aes(x=site, y = mlogCPUE, fill = sitetype)) + geom_bar(stat = "identity") +
  facet_grid(targets2~season, scales = "free") + 
  geom_errorbar(aes(ymin = mlogCPUE-selCPUE, ymax = mlogCPUE+selCPUE)) +
  scale_fill_manual(values = c("red", "limegreen"), name = "Site type") +
  ylab("mean log CPUE") 

#try a linear model
lm1 = lm(log(tCPUE)~site +targets2+ season, data =bugsblitzSF.1[which(bugsblitzSF.1$SampleID!="MAC6-22MAR2018"),])
summary(lm1)

lm2 = lm(log(tCPUE)~site +targets2+ season, data =bugsblitzSF.1)
summary(lm2)
visreg(lm2)
#meh

lm3 = aov(log(tCPUE)~site +targets2+ season, data =bugsblitzSF.1)
summary(lm3)
TukeyHSD(lm3)

########################################################################################################

#multivariates

#set up the community matrix
sf = dcast(bugsblitzSF, formula = SampleID~Analy2, value.var="CPUE", 
                 fun.aggregate = sum, fill = 0)


row.names(sf) = sf$SampleID
sf.12 = sf[,2:ncol(sf)]
sf2 = bugsblitzSF.1[which(rowSums(sf.12 )!=0),]
sf.12 = sf.12[which(rowSums(sf.12)!=0),]
#make one based on proportion of total catch for each sample
sf.12p = sf.12/rowSums(sf.12)


#PERMANOVA with relative abundance
pn1 = adonis(sf.12p~site +  targets2 +season , data = sf2)
print(pn1)
#NMDS

NMDSsf = metaMDS(sf.12p, try = 50, trymax = 500)
print( NMDSsf)

#plot it
source("plotNMDS.R")
PlotNMDS(NMDSsf, data = sf2, group = "site")
#test significance
envfit(NMDSsf~site, data = sf2)


PlotNMDS(NMDSsf, data = sf2, group = "season")
#test significance
envfit(NMDSsf~season, data = sf2)
