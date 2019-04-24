#zooplankton - 2018 update
library(tidyverse)
library(lubridate)
library(ggplot2)
library(readxl)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(simr)
library(lme4)
library(lmerTest)
################################################################################################################################
#First load the function I wrote to query the FRP database
source("querydatabase.R")

#Now specify the path to the FRP database
path = "U:/FRPA/MONITORING/Labs/Databases/FRPdata28DEC2018.accdb"
zoopsfrp <- GetFRPdata(path, "zoop")
View(zoopsfrp)

#add more station information
stations = read_excel("blitz2/Stations2.xlsx")
zoopsfrp = merge(zoopsfrp, stations)

################################################################################################################################
#Calculate CPUE

#combine the slides
totalcells = group_by(zoopsfrp, SampleID) %>% 
  summarise(totcells = max(CellNumber), totCountall = sum(Count))

zoo = group_by(zoopsfrp, Station, site, sitetype, sitetype2, Region2, Date, SampleID, Dilution, CommonName,
                Volume) %>% 
  summarise(totCount = sum(Count))

zoo = merge(zoo, totalcells)

#adjust for subsampling
zoo$atotal = (zoo$totCount/zoo$totcells)*zoo$Dilution

#Calculate CPUE
zoo$CPUE = zoo$atotal/zoo$Volume


#figuer out how many critters total
tot = sum(zoo$CPUE, na.rm = T)

#Now how many critters per taxon
zooptax = group_by(zoo, CommonName) %>% summarise(totaltaxa = sum(CPUE, na.rm = T), 
                                                  prop = totaltaxa/sum(zoo$CPUE, na.rm = T), 
                                                  tax = CommonName[1])

frptaxa <- read_excel("blitz2/zoocodes.xlsx")

zoop = merge(zoo, frptaxa)

#some critters have barely any individuals, lets combine groups
zoop$Analy[which(zoop$Analy == "Platyhelminthes" | zoop$Analy == "Bivalvia" | 
                   zoop$Analy == "Cumacea" | zoop$Analy == "Arachnida" | zoop$Analy == "Cnidaria" )] = "other"
zoop$Analy[which( zoop$Analy=="terrestrial" | zoop$Analy == "Trichoptera" |
                   zoop$Analy =="Hymenoptera" |zoop$Analy == "Coleoptera" | zoop$Analy == "Hemiptera" | 
                    zoop$Analy =="Ephemeroptera" | zoop$Analy =="Diptera" |
                    zoop$Analy =="Odonata" |  zoop$Analy =="Thysanoptera")] = "Insecta"

zoop$Analy[which(zoop$Analy == "corophium" | zoop$Analy == "Gammaridae")] = "Amphipoda"

zoop$AnalyLS = paste(zoop$Analy, zoop$LifeStage)
zoop$AnalyLS[which(zoop$AnalyLS == "other adult" | zoop$AnalyLS == "other larvae" | 
                     zoop$AnalyLS == "other juvenile" | zoop$AnalyLS == "NA NA" |zoop$AnalyLS ==  "Insecta NA")] = "other"

zoop$AnalyLS[which(zoop$AnalyLS == "fish adult" | zoop$AnalyLS == "fish larvae" | zoop$AnalyLS == "fish NA")] = "fish"
zoop$AnalyLS[which(zoop$AnalyLS == "Amphipoda adult" | zoop$AnalyLS == "Amphipoda larvae" )] = "Amphipoda"
zoop$AnalyLS[which(zoop$AnalyLS == "Insecta pupae")] = "Insecta adult"



#how many of each do I got?
zooptax2 = group_by(zoop, AnalyLS) %>% summarise(totaltaxa = sum(CPUE, na.rm = T), 
                                                prop = totaltaxa/sum(zoo$CPUE, na.rm = T))



#set up some labels for graphing
zooLab = c("Amphipoda", "Annelida", "Calanoida", "Cal juv", "Cladocera", "Collembola", 
            "Cyclopoda", "Cyclo juv", "Decapoda", "fish","Gastropoda", "Harpacticoida", "Insect adult", "Insect juv",
           "Isopoda", "Mysidea", "Nematoda","Ostracoda", "other", "Rotifera", "Tanaid")
zoop$AnalyLS = factor(zoop$AnalyLS, labels=zooLab)


#Now let's look at what we got!
#############################################################################################################

#filter for just samples from the blitz
zoop$month = month(zoop$Date)
zooblitz = filter(zoop, month >2 & month <5)


#set up a color pallete
mypal = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"), "black", "white", "grey", "pink")

#put zeros back in to get average catch by species
zoowide = dcast(zooblitz, SampleID~AnalyLS, fun.aggregate = sum, value.var = CPUE)

zooave = group_by(zooblitz, Analy, CommonName, site, AnalyLS) %>% summarize(CPUE = mean(CPUE, na.rm = T))
zooave =droplevels(zooave)

#quick plot of CPUE by location
z1 = ggplot(zooave, aes(x=site, y= CPUE))
z1 + geom_bar(stat = "identity", aes(fill = AnalyLS)) + scale_fill_manual(values = mypal)+
  coord_cartesian(ylim = c(0, 4000))

#total CPUE
zootot = group_by(zoo2017, ZSampleID, date, Station, 
                  Location, Region, SiteType2) %>% summarize(CPUE = sum(CPUE, na.rm = T))

z2 = ggplot(zootot, aes(x=Location, y= CPUE))
z2+geom_boxplot()
#One major outlier at Prospect. Let's figure out what it is

pros = filter(zoo2017, Location =="Prospect" & CPUE >1000)
#The major outlier is ZOOP3-25APR2017
#I'm going to take it out to make the analysis easier. 
zoo2017 = filter(zoo2017, ZSampleID != "ZOOP3-25APR2017")

#recalculate:
#quick plot of CPUE by location
zoo2017ave = group_by(zoo2017, Analy, CommonName, Location, AnalyLS, label) %>% summarize(CPUE = mean(CPUE, na.rm = T))
zoo2017ave =droplevels(zoo2017ave)
z1 = ggplot(zoo2017ave, aes(x=Location, y= CPUE))
z1 + geom_bar(stat = "identity", aes(fill = Analy)) + scale_fill_manual(values = mypal)

#total CPUE
zootot = group_by(zoo2017, ZSampleID, date, Station, 
                  Location, Region, SiteType2, label) %>% summarize(CPUE = sum(CPUE, na.rm = T))

z2 = ggplot(zootot, aes(x=Location, y= CPUE))
z2+geom_boxplot()
#better

zoototave = group_by(zootot,
                  Location, Region, SiteType2, label) %>% 
  summarize(mCPUE = mean(CPUE, na.rm = T), seCPUE = sd(CPUE)/length(CPUE), N = length(CPUE))

z3 = ggplot(zoototave, aes(x=label, y=mCPUE))
z3 + geom_bar(stat = "identity", aes(fill = Location)) +
  facet_grid(.~Region, scales = "free", space = "free_x") + 
  geom_errorbar(aes(ymin = mCPUE - seCPUE, ymax = mCPUE + seCPUE)) +
  scale_fill_manual(values = mypal, guide = "none")+
  geom_label(aes(label = paste("n = ", N), y = mCPUE+ 300), label.padding = unit(0.1, "lines"), size = 4) +
  mytheme + ylab("CPUE") + xlab(label = NULL)

#GLMm of total CPUE
zootot$SiteType2 = as.factor(zootot$SiteType2)
zblitz = lmer(log(CPUE) ~ Region + SiteType2 + (1|Location), data = zootot)
summary(zblitz)
visreg(zblitz)
#Differences by region, but not by habitat type. Opposite of macroinvertebrates!!!!
zblitza = glm(log(CPUE) ~ Region + SiteType2, data = zootot)
summary(zblitza)
visreg(zblitza)

zblitzb = aov(log(CPUE) ~ Region + SiteType2 + Location, data = zootot)
summary(zblitzb)
TukeyHSD(zblitzb)
#################################################################################

#calculate coefficient of variation in CPUE
zCVs = summarize(group_by(zootot, Location), mCPUE = mean(CPUE), 
                sdCPUE = sd(CPUE), seCPUE = sd(CPUE)/length(CPUE), N = length(ZSampleID), CV = sdCPUE/mCPUE)

#Calculate average withing-site CV and the CV of within-site means by sampling type.
zmCVs = summarize(zCVs, mCV = mean(CV, na.rm = T), mmCPUE = mean(mCPUE), 
                 sdmCPUE= sd(mCPUE), CV2 = sdmCPUE/mmCPUE)


############################################################################################

#plot of community composition 
z4 = ggplot(zoo2017, aes(x=label, y=CPUE))
z4 + geom_bar(stat = "identity", aes(fill = AnalyLS), position = "fill") +
  facet_grid(.~Region, scales = "free", space = "free_x") + 
  scale_fill_manual(values = c(mypal, "grey"), labels = zooLab, name = NULL) +
  mytheme + ylab("Relative Abudance") +xlab(label = NULL)

#plot of community composition by region
z4 = ggplot(zoo2017, aes(x=Region, y=CPUE))
z4 + geom_bar(stat = "identity", aes(fill = AnalyLS), position = "fill") +
  
  scale_fill_manual(values = c(mypal, "grey"), labels = zooLab, name = NULL) +
  mytheme + ylab("Relative Abudance") +xlab(label = NULL)


#plot of community composition by site type
z4 = ggplot(zoo2017, aes(x=SiteType2, y=CPUE))
z4 + geom_bar(stat = "identity", aes(fill = AnalyLS), position = "fill") +
  
  scale_fill_manual(values = c(mypal, "grey"), labels = zooLab, name = NULL) +
  mytheme + ylab("Relative Abudance") +xlab(label = NULL)

#Community matrix by CPUE
zMat = dcast(zoo2017, formula = ZSampleID~AnalyLS, value.var="CPUE", 
             fun.aggregate = sum, fill = 0)
zMat = zMat[,2:19]
#now by proprotions
zMatp = zMat/rowSums(zMat)

#PerMANOVA
za1 = adonis(zMatp~ Region + SiteType2 + Location, data = zootot)
za1
#region and site type are both significant, but site type much more so than region.

zNMDS = metaMDS(zMatp, try = 50)
zNMDS

zootot$Region = as.factor(zootot$Region)
PlotNMDS2(zNMDS, data = zootot, group = "Region", xlimit = c(-.8, 0.8), ylimit = c(-.8, .8))

zootot$SiteType2 = as.factor(zootot$SiteType2)
PlotNMDS2(zNMDS, data = zootot, group = "SiteType2", xlimit = c(-.8, 3), ylimit = c(-.8, .8))

#So, Tule red was really wierd. LEt's try re-running it again without tule red

zMatp1 = zMatp[which(zootot$Location != "Tule Red"),]

#run the Permanova again
za2 = adonis(zMatp1~ SiteType2 + Region, data = filter(zootot, Location != "Tule Red"))
za2

#Run the NMDS again
zNMDS2 = metaMDS(zMatp1, try = 50)
zNMDS2

#plot it
zootot$Region = as.factor(zootot$Region)
zootot$SiteType2 = as.factor(zootot$SiteType2)
PlotNMDS(zNMDS2, data = filter(zootot, Location != "Tule Red"), group = "Region")
envfit(zNMDS2~Region, data = filter(zootot, Location != "Tule Red"))

PlotNMDS(zNMDS2, data = filter(zootot, Location != "Tule Red"), group = "SiteType2")
envfit(zNMDS2~SiteType2, data = filter(zootot, Location != "Tule Red"))

########################################################################################################
#Look at catch at Tule Red to figure out what in the "other" category was so dominant

TR = droplevels(filter(zoop, Location == "Tule Red"))
TR1 = ggplot(TR, aes(x=CommonName, y = CPUE))
TR1 + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))

TR2 = ggplot(TR, aes(x=Analy, y = CPUE, fill = CommonName))
TR2 + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))

others = filter(TR, Analy == "other")
TR3 = ggplot(others, aes(x = CommonName, y = CPUE))
TR3 + geom_bar(stat = "identity")
