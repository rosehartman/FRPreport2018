#decker time series analysis - 2018


library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(lubridate)
library(readxl)


mypal = c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"))
mytheme = theme(text = element_text(size=14), legend.text.align = 0)

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


#import some info on the sites
sitetypes = read_excel("blitz2/Stations2.xlsx")
inverts3 = merge(inverts2, sitetypes[,c(2,9,10,11)])
inverts3$year = year(inverts3$Date)


#subset teh samples from 2017 and 2018
Decker <- filter(inverts3, site == "Decker" & (year ==2017 | year == 2018))

#add some more info
targets <- read.csv("blitz2/targets.csv")

#check to see if we haven't assigned targets for any samples
#bugsx = summarize(group_by(Decker, SampleID, Station, site, Region, Volume, Sampletype), total = sum(TotalCount, na.rm = T))
#test = merge(bugsx, targets, all.x = T)
#test2 = test[which(is.na(test$Target)),]

Decker = merge(Decker, targets, by = "SampleID", all.x = T)

bugsx = summarize(group_by(Decker, SampleID, Station, Volume, Sampletype, Target), total = sum(TotalCount, na.rm = T))


bugsx$effort = rep(NA, nrow(bugsx))

#effort for trawls and sweepnets is volume
bugsx$effort[which(bugsx$Sampletype=="Mysid net")]= bugsx$Volume[which(bugsx$Sampletype=="Mysid net")]
bugsx$effort[which(bugsx$Sampletype=="neuston trawl")]= bugsx$Volume[which(bugsx$Sampletype=="neuston trawl")]
bugsx$effort[which(bugsx$Sampletype=="sweepnet")]= bugsx$Volume[which(bugsx$Sampletype=="sweepnet")]


#effort for benthic cores and grabs is per unit area

bugsx$effort[which(bugsx$Sampletype=="Ponar grab")]=
  rep(0.05226, nrow(bugsx[which(bugsx$Sampletype=="Ponar grab"),]))

bugsx$effort[which(bugsx$Sampletype=="Petite Ponar")]=
  rep(0.0231, nrow(bugsx[which(bugsx$Sampletype=="Petite Ponar"),]))


bugsx$effort[which(bugsx$Sampletype=="PCV core")]=
  rep(0.00811, nrow(bugsx[which(bugsx$Sampletype=="PCV core"),]))

#everything else has an effort of "1"

bugsx$effort[which(is.na(bugsx$effort))] = 1

#attach the efforts to the origional dataset
Decker= merge(Decker, bugsx)
Decker$TotalCount[which(is.na(Decker$TotalCount))] =0
Decker$subsampled[which(is.na(Decker$subsampled))] =100


# Adjust total count for subsampling
Decker$atotal = Decker$TotalCount*(1/(Decker$subsampled/100))

#Calculate CPUE
Decker$CPUE = Decker$atotal/Decker$effort

#add some analysis groups
zoocodes <- read_excel("blitz2/zoocodes.xlsx")
zoocodes <- zoocodes[,c(3,12,13, 14)]
Decker = merge(Decker, zoocodes)
Decker = Decker[order(Decker$SampleID),]



#I'm going to go ahead and remove the mesozooplanktoton from all the blitz bugs samples. I may want to make some extra notes 
#on the catch of cladocera in sweep nets, but just to make stuff clearer for now...
Decker = filter(Decker, Analy2 != "Copepoda" & Analy2 != "Cladocera" & Analy2 != "Ostracoda")


sites = group_by(inverts2, Station) %>% summarize(count = length(Station))

#import some info on the sites
sitetypes = read_excel("blitz2/Stations2.xlsx")
Decker = merge(Decker, sitetypes[,c(2,9,10,11)])

#put the sties in order along the estuarine/fresh gradient
Decker$site = factor(Decker$site, levels = c("Ryer", "Grizzly", "Tule Red", "Wings", "Blacklock", "Bradmoor",
                                                   "L Honker", "Browns","Winter", "Broad", "Dow",
                                                   "Sherman", "Horseshoe", "Decker", "Stacys",
                                                   "Lindsey", "Liberty", "Prospect", "Miner", "Flyway"))

#Create a new colum where "targets" (habitat type) combines all the vegetation samples into one type.
Decker$targets2 = as.character(Decker$Target)
Decker$targets2[which(Decker$Target == "EAV" | Decker$Target == "SAV"| Decker$Target == "FAV" )] = "sweep net"

#Summarize by sample and calculate the total CPUE of all the macroinvertebrates in the sample
Decker.1 = summarize(group_by(Decker, SampleID, Station, Region2, Target, targets2, 
                                 Region, Sampletype, site, sitetype, Date),
                        tcount = sum(atotal, na.rm = T), tCPUE = sum(CPUE, na.rm = T), richness = length(unique(Analy)))

Decker.1$site = as.factor(Decker.1$site)
Decker.1$sitetype = as.factor(Decker.1$sitetype)
Decker.1$targets2 = as.factor(Decker.1$targets2)

#summarize by Region and habitat type so I can make a quick plot

#add up all the critters that are in the same analysis group
bugssum = summarize(group_by(Decker, site, targets2, 
                             Region2, Station, sitetype, Sampletype, SampleID, Analy), 
                    tCPUE = sum(CPUE))

#I need to add the zeros in for taxa that were not present in a sample so I can take an appropriate mean
#I usually do this by transfering it from "long" to "wide" and back again.

#Make it wide
ComMat.1 = dcast(Decker, formula = SampleID~Analy2, value.var="CPUE", 
                 fun.aggregate = sum, fill = 0)
#do it by smaller grouping
ComMat.2 = dcast(Decker, formula = SampleID~Analy, value.var="CPUE", 
                 fun.aggregate = sum, fill = 0)

#do it by common name
ComMat.CN = dcast(Decker, formula = SampleID~CommonName, value.var="CPUE", 
                  fun.aggregate = sum, fill = 0)

#do it by cleaned up common name
ComMat.CN2 = dcast(Decker, formula = SampleID~CN, value.var="CPUE", 
                   fun.aggregate = sum, fill = 0)

#Go back to long
Decker.2x = melt(ComMat.1, id.vars = "SampleID", variable.name = "Analy2", value.name = "CPUE")
Decker.2x = merge(Decker.1, Decker.2x)

#Means by analysis group for each location
bugssum.1 = summarize(group_by(Decker.2x, site, Region2, sitetype, targets2, 
                               Region, Sampletype, Analy2), mCPUE = mean(CPUE),
                      sdCPUE = sd(CPUE), seCPUE = sd(CPUE)/length(CPUE))

#Now calculate total CPUE
bugstot = summarize(group_by(Decker.2x, SampleID, Date, site, sitetype, Target, targets2, 
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



#we just want to look at spring of 2017
Deckerbugs = filter(Deckerbugs, Date > "2017-1-1" & Date < "2017-7-1")


#take out the zooplankotn
Deckerbugs = filter(Deckerbugs,  analysisB != "Copepoda" & analysisB != "Cladocera")


#Summarize by sample and get rid of the zooplankotn
Deckerbugs.1 = summarize(group_by(Deckerbugs, miSampleID, Station, Target, 
                                    sampletype, Date, effort),
                           tcount = sum(atotal), tCPUE = sum(CPUE), richness = length(unique(analysisA)))

######################################################################
#now lets make a few quick plots
#put zeros in so we can make area plots and stuff
Deckerbugs2 = dcast(Deckerbugs, formula = miSampleID~analysisB, value.var="CPUE", 
                    fun.aggregate = sum, fill = 0) %>%
  melt(id.vars = "miSampleID", variable.name = "analysisB", value.name = "CPUE")
Deckerbugs2.1 = merge(Deckerbugs2, Deckerbugs.1)
Deckerbugs2.1$Month = month(Deckerbugs2.1$Date)

#summerize by month
Deckerbugs3 = group_by(Deckerbugs2.1, Month, sampletype, analysisB) %>% summarise(mCPUE = mean(CPUE, na.rm = T))

#plot catch over time
c1 = ggplot(Deckerbugs3[which(Deckerbugs3$sampletype != "benthic"),], 
              aes(x=Month, y = mCPUE))
c1 +  geom_area(aes(fill =analysisB)) + 
  facet_wrap(~sampletype, scales = "free_y") +
  ylab("CPUE") + 
  scale_fill_manual(values = mypal, name = "Taxon")

################################################################################
#try a GLM of total catch by date
Deckerbugs.1$yday = yday(Deckerbugs.1$Date)
m1 = glm(log(tCPUE) ~ sampletype + yday, data =Deckerbugs.1)
summary(m1)
m1l =lm(log(tCPUE) ~ sampletype + yday, data =Deckerbugs.1)
summary(m1l)

# if there are any intermediate peaks, a quadratic should fit better than a linear function

m2 = glm(log(tCPUE) ~ sampletype + I(yday^2), data =Deckerbugs.1)
summary(m2)
m2l = lm(log(tCPUE) ~ sampletype + I(yday^2), data =Deckerbugs.1)
summary(m2l)

#Not really. Linear it is!

####################################################################################
#Quick graph of other survey's smelt catch over the past 15 years

decksmelt <- read_excel("deckersmelt.xlsx",     sheet = "2002-2017")

#calculate total catch by month and mean catch by month
dssum = group_by(decksmelt, month, Survey) %>% 
  summarize(tcatch = sum(Catch), mcatch = tcatch/16)

dssum1 = dcast(dssum, month ~ Survey, value.var = "tcatch")
dssum2 = melt(dssum1, id.vars = "month", variable.name = "Survey", value.name = "tcatch")
#do total catch instead, and add zeros, just to make graphing easier
dssum2$tcatch[which(is.na(dssum2$tcatch))] = 0
dssum2$mcatch = dssum2$tcatch/16


ds3 = ggplot(dssum2, aes(x=month, y = mcatch, fill = Survey))
ds3 + geom_area(stat= "identity") + xlab("Month") + ylab("Mean Catch, 2002-2017") +
  scale_x_continuous(breaks = c(1:12))

###########################################################################################
#Chinook salmon

Chipps <- read_excel("Chipps Island Trawls CHN & POD Species 2012-2018.xlsx")
Chipps$Date = as.Date(Chipps$Date)
Chipps$month = month(Chipps$Date)

stagelabs = data.frame(Stage = factor(c("1","2","3","4","5","6", "N/P")), 
                       labels = factor(x=c("yolk-sac","fry", "par", "silver par",
                                           "smolt", "adult", "Not Provided" ), 
                                       levels = c("yolk-sac","fry", "par", "silver par",
                                                  "smolt", "adult", "Not Provided" ),
                                       ordered = T))

Chipps1 = filter(Chipps, Species == "CHN" & Date < "2018-1-1")
Chipps1 = merge(Chipps, stagelabs)
#calculate total catch by month and mean catch by month
chipsum = group_by(Chipps1, Stage, labels, month) %>% 
  summarize(tcatch = sum(Catch), mcatch = tcatch/6)

chipsum1 = dcast(chipsum, month ~ labels, value.var = "tcatch")
chipsum2 = melt(chipsum1, id.vars = "month", variable.name = "Stage", value.name = "tcatch")
chipsum2$tcatch[which(is.na(chipsum2$tcatch))] = 0
chipsum2$mcatch = chipsum2$tcatch/16


chip2 = ggplot(filter(chipsum2, Stage != "adult" & Stage != "fry"), aes(x=month, y = tcatch, fill = Stage))
chip2 + geom_area(stat= "identity") + xlab("Month") + ylab("Total Catch, 2012-2017") +
  scale_x_continuous(breaks = c(1:12))


########################################################################
#I want to plot all the different groups on the same scale so I can look at overlap.
#To do this, I'll calculate the % of total catch that falls in each month.

#First the zooplankotn
#dzooTot$month = month(dzooTot$Date)
#dzooTot2 = group_by(dzooTot, month) %>% summarise(mCPUE = mean(tCPUE))
#tot = sum(dzooTot2$mCPUE)
#dzooTot2$prop = dzooTot2$mCPUE/tot
#dzooTot2$sampletype = "zoop"
#Actually, I'm not that interested in zoops, because i"m targeting macroinvertebrate sample timing


#now macroinvertebrates
Deckerbugs.1$month = month(Deckerbugs.1$Date)
dbug = group_by(Deckerbugs.1, month, sampletype) %>% summarise(mCPUE = mean(tCPUE))
bugtot = group_by(dbug, sampletype) %>% summarise(tot = sum(mCPUE))
dbug = merge(dbug, bugtot)
dbug$prop = dbug$mCPUE/dbug$tot
dbug$tot = NULL
dbug = filter(dbug, sampletype != "benthic")

#now salmon
chipMonth = group_by(filter(Chipps1, labels == "smolt"), month) %>% 
  summarize(mCPUE = sum(Catch)/6)
tot = sum(chipMonth$mCPUE)
chipMonth$prop = chipMonth$mCPUE/tot
chipMonth$sampletype = "Chinook smolts"

#now smelt
dstot = group_by(filter(dssum2, Survey == "SKT"), month) %>% summarise(mCPUE = mean(mcatch))
tot = sum(dstot$mCPUE)
dstot$prop = dstot$mCPUE/tot
dstot$sampletype = "Delta Smelt adults"

#combine the datasets
decker.critters = rbind(dstot, chipMonth, dbug)
decker.fish = rbind(dstot, chipMonth)

#plot it
critters = ggplot(decker.critters, aes(x=month, y=prop))
critters + geom_line(aes(color = sampletype)) + geom_point(aes(color = sampletype)) +
  coord_cartesian(xlim = c(1, 7)) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), 
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul"))+
  scale_color_manual(values = mypal)

#I want the combined hight to be the greatest.
critters + geom_bar(aes(fill = sampletype), stat = "identity") + 
  coord_cartesian(xlim = c(1, 7)) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7), 
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul"))+
  scale_color_manual(values = mypal)


fish = ggplot(decker.fish, aes(x=month, y=prop))
fish + geom_bar(aes(fill = sampletype), stat = "identity") + 
 coord_cartesian(xlim = c(1, 6), ylim = c(0,0.6)) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6), 
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun"))+
  scale_fill_manual(values = mypal, name = NULL) +
  geom_smooth(aes(x=month, y = prop, color = "Macroinvertebrates"), data = dbug, method = "lm")+
  ylab("Percentage of total catch per month") +
  scale_color_manual(name = NULL, values = "blue")

