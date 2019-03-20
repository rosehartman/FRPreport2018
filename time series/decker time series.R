#decker time series analysis - clean


library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(lubridate)
library(readxl)


mypal = c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"))
mytheme = theme(text = element_text(size=14), legend.text.align = 0)



#Check out macroinvertebrates
Deckerbugs <- read_excel("Deckerbugs.xlsx", 
                         sheet = "deckerbugs", col_types = c("text", 
                                                             "numeric", "numeric", "numeric", 
                                                             "text", "text", "date", "text", "text", 
                                                             "text", "text", "text", "text", "text", 
                                                             "text", "logical", "text", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric", "numeric", "numeric", 
                                                             "numeric"))

#we just want to look at spring of 2017
Deckerbugs = filter(Deckerbugs, Date > "2017-1-1" & Date < "2017-7-1")

#add some more info
targets <- read.csv("targets.csv")
Deckerbugs = merge(Deckerbugs, targets, all.x = T)

Deckerbugs.e = summarize(group_by(Deckerbugs, miSampleID, Station, Target,VegWeight,
                                  sampletype, Date, NetMeterSTart, NetMeterEnd, LatStart, LongStart,
                                  LatEnd, LongEnd),
                         mean.en= mean(NetMeterSTart))




# Calculate effort for the trawls
Deckerbugs.e$distance = rep(NA, nrow(Deckerbugs.e))

Deckerbugs.e$distance[which(Deckerbugs.e$sampletype=="neuston")]= apply(Deckerbugs.e[which(Deckerbugs.e$sampletype=="neuston"),c(9:12)], 1, function(x) {
  R <- 6371
  rads <- x * pi/180
  d <- sin(rads[1])*sin(rads[3])+cos(rads[1])*cos(rads[3])*cos(abs(rads[2]-rads[4]))
  d <- acos(d)
  R*d*1000
} )

Deckerbugs.e$effort = rep(NA, nrow(Deckerbugs.e))

# Effort for neuston trawl is surface area sampled, distance in m * width of net mouth (45cm)

# Some of the neuston trawls didn't have GPS coordinates, so I'll just average the distance
# for the rest of the 3-minute trawls.
#meandis = mean(Deckerbugs.e$distance[which(Deckerbugs.e$sampletype =="neuston")], na.rm=T)
#Deckerbugs.e$distance[which(is.na(Deckerbugs.e$distance) & Deckerbugs.e$sampletype=="neuston")] = meandis
#Deckerbugs.e$effort[which(Deckerbugs.e$sampletype=="neuston")] = Deckerbugs.e$distance[which(Deckerbugs.e$sampletype=="neuston")]*.45

# Effort for benthic and oblique trawls is based on flowmeter readings and net mouth area

Deckerbugs.e$effort[which(Deckerbugs.e$sampletype=="trawl")]= 
  (Deckerbugs.e$NetMeterEnd[which(Deckerbugs.e$sampletype=="trawl") ]-
     Deckerbugs.e$NetMeterSTart[which(Deckerbugs.e$sampletype=="trawl")])*0.026873027*0.16 #I've used the factory calibration here

Deckerbugs.e$effort[which(Deckerbugs.e$effort<0)] = ((1000000 - Deckerbugs.e$NetMeterEnd[which(Deckerbugs.e$effort<0)])+
                                                       Deckerbugs.e$NetMeterSTart[which(Deckerbugs.e$effort<0)])*0.026873027*0.16 #I've used the factory calibration here

meantrawl = mean(Deckerbugs.e$effort[which(Deckerbugs.e$sampletype=="trawl" &
                                             Deckerbugs.e$miSampleID != "MAC1-8MAY2017" &
                                             Deckerbugs.e$miSampleID != "MAC1D-7JUN2017" )], na.rm = T)

#Something's wrong with a couple of the mysid trawls
Deckerbugs.e$effort[which(Deckerbugs.e$miSampleID == "MAC1-8MAY2017" | 
                            Deckerbugs.e$miSampleID == "MAC1D-7JUN2017"|
                            Deckerbugs.e$miSampleID == "MYSID1-14feb2017"|
                            Deckerbugs.e$miSampleID == "MYSID1-25jan2017")] = meantrawl



Deckerbugs.e$effort[which(Deckerbugs.e$sampletype=="sweepnet" & Deckerbugs.e$VegWeight>0)] =
  Deckerbugs.e$VegWeight[which(Deckerbugs.e$sampletype=="sweepnet" & Deckerbugs.e$VegWeight>0)]

#everything else has an effort of "1"

Deckerbugs.e$effort[which(is.na(Deckerbugs.e$effort))] = 1

#attach the efforts to the origional dataset
Deckerbugs= merge(Deckerbugs, Deckerbugs.e)


# Adjust total count for subsampling
Deckerbugs$atotal = Deckerbugs$TotalCount*(1/(Deckerbugs$subsampled/100))

#Calculate CPUE
Deckerbugs$CPUE = Deckerbugs$atotal/Deckerbugs$effort

#add some analysis groups
zoocodes <- read.csv("zoocodes.csv", as.is=F)
zoocodes <- zoocodes[,c(1,11,12, 13)]
zoocodes$InvertCode = zoocodes$Zoo.Code
Deckerbugs = merge(Deckerbugs, zoocodes)
Deckerbugs = Deckerbugs[order(Deckerbugs$miSampleID),]

#Summarize by order, family, taxon
tax = summarize(group_by(Deckerbugs,  analysisB, Phylum, Class, Order), total = sum(TotalCount))

#take out the blank ones
Deckerbugs = Deckerbugs[-which(is.na(Deckerbugs$analysisA)),]

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

