#decker time series analysis - 2018


library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(lubridate)
library(readxl)
library(MuMIn)
library(visreg)


mypal = c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"))
mytheme = theme(text = element_text(size=14), legend.text.align = 0)

#source("querydatabase.R")

#Now specify the path to the FRP database
#path = "U:/FRPA/MONITORING/Labs/Databases/FRPdata28DEC2018.accdb"

#Query the invertebrate data
#inverts2 <- GetFRPdata(path, type = "inverts")

#Query the station information
#stations = GetFRPdata(path, type = "stations")

#Load the invertebrate data
inverts2 <- read_excel("invertsqry.xlsx")

#load the station information
stations = read_excel("Stations2.xlsx")


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
inverts3$month = month(inverts3$Date)
inverts3$month2 = as.factor(inverts3$month)
inverts3$Year2 = as.factor(inverts3$year)


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

#Take out the benthic samples
Decker = filter(Decker, Target != "benthic")

#We didn't have FAV or SAV samples in all months, which I think is throwing us off
sites = group_by(inverts2, Station) %>% summarize(count = length(Station))
Decker = filter(Decker, Target != "SAV" & Target != "FAV")

#One of the january MAC samples got totally filled with vegetation
Decker = filter(Decker, SampleID != "MAC1-23JAN2018")

#import some info on the sites
sitetypes = read_excel("blitz2/Stations2.xlsx")
Decker = merge(Decker, sitetypes[,c(2,9,10,11)])

#Create a new colum where "targets" (habitat type) combines all the vegetation samples into one type.
Decker$targets2 = as.character(Decker$Target)
Decker$targets2[which(Decker$Target == "EAV" | Decker$Target == "SAV"| Decker$Target == "FAV" )] = "sweep net"

#Summarize by sample and calculate the total CPUE of all the macroinvertebrates in the sample
Decker.1 = summarize(group_by(Decker, SampleID,targets2, 
                                 Sampletype, month, month2, year, Year2, Date, Volume),
                        tcount = sum(atotal, na.rm = T), tCPUE = sum(CPUE, na.rm = T), richness = length(unique(Analy)))


Decker.1$targets2 = as.factor(Decker.1$targets2)

#summarize by Region and habitat type so I can make a quick plot

#add up all the critters that are in the same analysis group
decksum = summarize(group_by(Decker, targets2, 
                              Sampletype, SampleID, Analy, Analy2, month, month2, year, Year2), 
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

#Means by analysis group for each month
decksum.1 = summarize(group_by(Decker.2x, targets2, 
                               Analy2, month, month2, year, Year2), mCPUE = mean(CPUE),
                      sdCPUE = sd(CPUE), seCPUE = sd(CPUE)/length(CPUE))

#Now calculate total CPUE
decktot = summarize(group_by(Decker.2x, SampleID, Date, month, month2, year, Year2, targets2, 
                              Sampletype), tCPUE = sum(CPUE, na.rm = T))



decktot$logtot = log(decktot$tCPUE +1)
decktotNoB = filter(decktot, targets2!= "benthic")

#Mean CPUE, sample size, and standard error for each location and habitat so I can make a pretty graph.
decktotave = summarize(group_by(decktot, targets2, month, year, month2, Year2), 
                       mCPUE = mean(tCPUE, na.rm = T), mlogCPUE = mean(logtot), 
                       sdCPUE = sd(tCPUE, na.rm = T), seCPUE = sd(tCPUE)/length(tCPUE), N = length(SampleID))


######################################################################
#now lets make a few quick plots

#plot catch over time
c1 = ggplot(decksum.1, 
              aes(x=month, y = mCPUE))
c1 +  geom_area(aes(fill =Analy2)) + 
  facet_wrap(Year2~targets2, scales = "free_y") +
  ylab("CPUE") + 
  scale_fill_manual(values = mypal, name = "Taxon") +
  coord_cartesian(xlim = c(1,6))

################################################################################
#try a GLM of total catch by date
Decker.1$yday = yday(Decker.1$Date)
m1 = glm(log(tCPUE + 1) ~ targets2 + yday*Year2, data =Decker.1)
summary(m1)
visreg(m1)
visreg(m1, xvar = "yday", by = "Year2")
m1l =lm(log(tCPUE+1) ~ targets2 + yday, data =Decker.1)
summary(m1l)

# if there are any intermediate peaks, a quadratic should fit better than a linear function

m2 = glm(log(tCPUE+1) ~ targets2 + I(yday^2), data =Decker.1)
summary(m2)
m2l = lm(log(tCPUE+1) ~ targets2 + I(yday^2), data =Decker.1)
summary(m2l)


####################################################################################
#Quick graph of other survey's smelt catch over the past 15 years

decksmelt <- read_excel("time series/deckersmelt.xlsx",     sheet = "2002-2017")

#calculate total catch by month and mean catch by month
dssum = group_by(decksmelt, month, Survey) %>% 
  summarize(tcatch = sum(Catch), mcatch = tcatch/17)

dssum1 = dcast(dssum, month ~ Survey, value.var = "tcatch")
dssum2 = melt(dssum1, id.vars = "month", variable.name = "Survey", value.name = "tcatch")
#do total catch instead, and add zeros, just to make graphing easier
dssum2$tcatch[which(is.na(dssum2$tcatch))] = 0
dssum2$mcatch = dssum2$tcatch/18


ds3 = ggplot(dssum2, aes(x=month, y = mcatch, fill = Survey))
ds3 + geom_area(stat= "identity") + xlab("Month") + ylab("Mean Catch, 2002-2018") +
  scale_x_continuous(breaks = c(1:12))

###########################################################################################
#Chinook salmon
chipps = read.csv("time series/2002-2018_DJFMP_trawl_fish_and_water_quality_data.csv")
chipps$SampleDate = as.Date(chipps$SampleDate) 
chipps = filter(chipps, Location == "Chipps Island" & OrganismCode == "CHN" & year(SampleDate) >2002)
chipps$month = month(chipps$SampleDate)

stagelabs = data.frame(StageCode = factor(c("1","2","3","4","5","6", "n/p")), 
                       labels = factor(x=c("yolk-sac","fry", "par", "silver par",
                                           "smolt", "adult", "Not Provided" ), 
                                       levels = c("yolk-sac","fry", "par", "silver par",
                                                  "smolt", "adult", "Not Provided" ),
                                       ordered = T))

Chipps1 = merge(chipps, stagelabs)
#calculate total catch by month and mean catch by month
chipsum = group_by(Chipps1, StageCode, labels, month) %>% 
  summarize(tcatch = sum(Count), mcatch = tcatch/6)

chipsum1 = dcast(chipsum, month ~ labels, value.var = "tcatch")
chipsum2 = melt(chipsum1, id.vars = "month", variable.name = "Stage", value.name = "tcatch")
chipsum2$tcatch[which(is.na(chipsum2$tcatch))] = 0
chipsum2$mcatch = chipsum2$tcatch/16


chip2 = ggplot(filter(chipsum2, Stage != "adult"), aes(x=month, y = tcatch, fill = Stage))
chip2 + geom_area(stat= "identity") + xlab("Month") + ylab("Total Catch, 2002-2017") +
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
bugtot = group_by(decktot,month, targets2) %>% summarise(tot = sum(tCPUE))
bugtot2 = group_by(bugtot, targets2) %>% summarise(totall = sum(tot))
dbug = merge(bugtot, bugtot2)
dbug$prop = dbug$tot/dbug$totall
dbug$tot = dbug$totall
dbug$totall = NULL
dbug$sampletype = dbug$targets2
dbug$targets2 = NULL


#now salmon
chipMonth = group_by(filter(Chipps1, labels == "smolt"), month) %>% 
  summarize(tot = sum(Count)/6)
totall = sum(chipMonth$tot)
chipMonth$prop = chipMonth$tot/totall
chipMonth$sampletype = "Chinook smolts"

#now smelt
dstot = group_by(filter(dssum2, Survey == "SKT"), month) %>% summarise(tot = mean(mcatch))
totall = sum(dstot$tot)
dstot$prop = dstot$tot/totall
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


#See how EMP's mysid data compares
EMPMysid <- read_excel("time series/1972-2018MysidMatrix.xlsx", 
                       sheet = "Mysid CPUE Matrix 1972-2018 ", 
                       col_types = c("numeric", "numeric", "numeric", 
                                     "numeric", "date", "text", "text", 
                                     "text", "numeric", "text", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric", "numeric", "numeric", 
                                     "numeric"))

#filter out just the station near Decker, just the past few years
EMPDeck = filter(EMPMysid, Station == "NZ064", Year >2001)

#total mysids
EMPDeck$tCPUE = rowSums(EMPDeck[,c(17:24)])
EMPDeck$month = month(EMPDeck$SampleDate)
EMPmonth = group_by(EMPDeck, month) %>% summarize(tot = mean(tCPUE))
totall = sum(EMPmonth$tot)
EMPmonth$prop = EMPmonth$tot/totall

#try graphing it again
fish + geom_bar(aes(fill = sampletype, color = sampletype), stat = "identity") + 
  coord_cartesian(xlim = c(1, 6), ylim = c(0,0.6)) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6), 
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun"))+
  scale_fill_manual(values = mypal, name = NULL) +
  geom_smooth(aes(x=month, y = prop, color = "FRP Macroinverts", fill = "FRP Macroinverts"), data = dbug, method = "lm")+
  geom_smooth(aes(x=month, y = prop, color = "EMP Mysids", fill = "EMP Mysids"), data = EMPmonth, method = "lm")+
  ylab("Percentage of total catch per month") +
  scale_color_manual(name = NULL, values = c("red", "blue", "limegreen", "purple")) +
  scale_fill_manual(name = NULL, values = c("red", "blue", "limegreen", "purple"))


###########################################################################################
#Dan suggested trying to use flow rather than day of the year to see what the trend looks like
dayflow2017 <- read_csv("time series/dayflowCalculations2017.csv")
dayflow2018 <- read_csv("time series/dayflowCalculations2018.csv")
dayflow = rbind(dayflow2017, dayflow2018)
dayflow$Date = as.Date(dayflow$Date)

#Average sacramento river flow by month
monthflow = group_by(dayflow, Mo, Year) %>% summarize(sac = mean(SAC))
monthflow$month = monthflow$Mo
monthflow$year = monthflow$Year

deckday = merge(Decker.1, monthflow)

#try a GLM of total catch by flow
m1 = glm(log(tCPUE +1) ~ targets2 + sac + month + year + sac:year, data =deckday, na.action = na.fail)
dredge(m1)
#best model had just river flow and gear type
mbest = glm(log(tCPUE +1) ~ targets2 + sac*year, data =deckday)
summary(mbest)
visreg(mbest)
visreg(mbest,xvar = "sac", by = "year")

mfl = glm(log(tCPUE +1) ~ targets2 + sac, data =deckday)
summary(mfl)
visreg(mfl)

#just a plot of catch versus flow
ggplot(deckday, aes(x=sac, y = tCPUE)) + geom_point(aes(color = Year2)) + geom_smooth(method = "lm", color = "black")

ggplot(deckday, aes(x=sac, y = tCPUE)) + geom_point(aes(color = Year2)) +
  scale_y_log10()+ geom_smooth(method = "lm", color = "black") + xlab("Sacramento River Flow (CFS)") +
  ylab("total CPUE of macroinvertebrates per sample") +
  scale_color_manual(values = c("red", "blue"), name = NULL)

#add the fish in there
chipp1718 <- read_excel("time series/Chipps Island Trawls CHN & POD Species 2012-2019.xlsx")
chipps1718 = droplevels(filter(chipp1718, year(Date) >2016 & Species == "CHN" & (Stage == "5" | Stage == "4")))
str(chipps1718)
chipps1718$month = month(chipps1718$Date)
chipps1718$year = year(chipps1718$Date)
chipps1718 = merge(chipps1718, monthflow)
chipps1718m = group_by(chipps1718, month, year, sac) %>% summarize(tcatch = sum(Catch)) 

ggplot(chipps1718m, aes(x=sac, y = tcatch)) + geom_point(aes(color = as.factor(year)))

ggplot(chipps1718m, aes(x=sac, y = tcatch)) + geom_point(aes(color = as.factor(year))) +
  scale_y_log10()
