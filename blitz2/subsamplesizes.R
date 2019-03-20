#Let's see how often our mysid and zoop catch is way below waht we want

library(tidyverse)
library(readxl)

#read in a query of all the critters counted
#Note: I need to go back and remove cladocera and stuff.
bugs = read_xlsx("totalcatch.xlsx")


#add some analysis groups
zoocodes <- read_excel("zoocodes.xlsx")
zoocodes <- zoocodes[,c(3,12,13)]
inverts3 = merge(inverts2, zoocodes)
inverts3 = inverts3[order(inverts3$SampleID),]

#I'm going to go ahead and remove the mesozooplanktoton from all the blitz bugs samples since we
#stopped counting them in 2018
inverts3 = filter(inverts3, Analy2 != "Copepoda" & Analy2 != "Cladocera"& Analy2 != "Ostracoda")


#now calculate tthe total critters counted
bugs = group_by(inverts3, Date, SampleID, Region, Sampletype, Volume) %>% summarize(totalcount = sum(TotalCount, na.rm = T))


#look at total number of organisms counted from all macroinvertebrate samples
ggplot(bugs, aes(x = totalcount)) + geom_histogram()

#zoom in on the small number
ggplot(bugs, aes(x = totalcount)) + geom_histogram(binwidth = 5)+
  coord_cartesian(xlim = c(0,100))

#facet by sample type
ggplot(bugs, aes(x = totalcount)) + geom_histogram(binwidth = 5)+
  coord_cartesian(xlim = c(0,100)) + facet_wrap(~Sampletype)

#Just look at the mysid trawls
bugs2 = filter(bugs, Sampletype == "Mysid net")

ggplot(bugs2, aes(x = totalcount)) + geom_histogram(binwidth = 5) +
  ylab("number of trawls") +
  coord_cartesian(xlim = c(0,100))

bugssmall = filter(bugs2, totalcount < 100)
#we have 190 samples out of 264 that had under 100 critter total

bugsVsmall = filter(bugs2, totalcount < 50)
#143 samples had under 50 critters.

ggplot(bugssmall, aes(x = totalcount)) + geom_histogram(binwidth = 5) +
  facet_wrap(~Region)


#Now let's look at the zooplankton tows
zoops = read_xlsx("totalcatch.xlsx", sheet = "zoops")
ggplot(zoops, aes(x = SumOfCount)) + geom_histogram() + ylab("Number of trawls") + xlab("total critters counted")

ggplot(zoops, aes(x = SumOfCount)) + geom_histogram(binwidth = 5)+
  coord_cartesian(xlim = c(0,400))+ ylab("Number of trawls") + xlab("total critters counted")

zoopssmall = filter(zoops, SumOfCount < 400)
zoopsVsmall = filter(zoops, SumOfCount < 100)

###############################################################################
#Ok, now we don't want to increase mysid trawls too much if it means too much fish catch

#What was teh total volume of samples taken with the mysid net?
vol = sum(bugs2$Volume)

#import listed fishes
fish =read_excel("listedfishmysids.xlsx")
#so, tehre were two delta smelt and 25 longfin smelt collected in mysid nets in March of 2018. No other listed
#fish in any of the mysid trawls. 

fishx = merge(fish[,1:3], bugs2, all.y = T)
fishwide = spread(fishx, Species, Count, fill = 0)
fishlong = gather(fishwide[,1:8], key = "Species", value = "Count", DELSME, LONSME)
fishlong = fishlong[order(fishlong$Date),]

#cumulative sum of volumes
fishlong$cvol = cumsum(fishlong$Volume)

#cumulative sum of fish catch
fishlfn = filter(fishlong, Species == "LONSME")
fishlfn$csum = cumsum(fishlfn$Count)
fishdelt = filter(fishlong, Species == "DELSME")
fishdelt$csum = cumsum(fishdelt$Count)
fishlong = rbind(fishlfn, fishdelt)
#there was probably an easier way to do that

#cumulative volume versus number of smelt caught
ggplot(fishlong, aes(x=cvol, y=csum, col = Species)) + geom_line() +
  ylab("cumulative sum of fish catch") + xlab("cumulative volume sampled")

ggplot(fishlong, aes(x=Date, y=csum, col = Species)) + geom_line()+
  ylab("cumulative sum of fish catch") + xlab("date")

ggplot(fishlong, aes(x=Date, y=cvol)) + geom_line()+
  ylab("cumulative sum of volume") + xlab("date")
