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
zoopsfrpx = merge(zoopsfrp, stations)

################################################################################################################################
#Calculate CPUE

#combine the slides
totalcells = group_by(zoopsfrpx, SampleID) %>% 
  summarise(totcells = max(CellNumber), totCountall = sum(Count))

zoo = group_by(zoopsfrpx, Station, site, sitetype, sitetype2, Region2, Date, SampleID, Dilution, CommonName,
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
zooLab = c("Amphipoda", "Annelida", "Barnacle Nauplii", "Calanoida", "Cal juv", "Cladocera", "Collembola", 
            "Cyclopoda", "Cyclo juv", "Decapoda", "fish","Gastropoda", "Harpacticoida", "Insect adult", "Insect juv",
           "Isopoda", "Mysidea", "Nematoda","Ostracoda", "other", "Rotifera", "Tanaid")
zoop$foo = zoop$AnalyLS
zoop$AnalyLS = factor(zoop$AnalyLS, labels = zooLab)

test = group_by(zoop, AnalyLS, foo) %>% summarise(catch = sum(CPUE))
#Now let's look at what we got!
#############################################################################################################

#filter for just samples from the blitz
zoop$month = month(zoop$Date)
zoop$year = year(zoop$Date)
zooblitz = filter(zoop, month >2 & month <5 & year >2016)
zooblitz = filter(zooblitz, site != "Lindsey" & 
                    site != "Sherman" & site != "Dow" & site != "Wings")


#set up a color pallete
mypal = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"), "black", "white", "grey", "pink")

#put zeros back in to get average catch by species
zoowide = dcast(zooblitz, SampleID+Region2+site+Date+month+year+sitetype2~AnalyLS, fun.aggregate = sum, na.rm = T, value.var = "CPUE")
zooblitz0 = melt(zoowide, id.vars = c("SampleID", "site", "Region2","sitetype2", "year", "Date", "month"), variable.name = "AnalyLS", value.name = "CPUE" )

zooave = group_by(zooblitz0, site, sitetype2, Region2, AnalyLS, year) %>% summarize(CPUE = mean(CPUE, na.rm = T))
zooave =droplevels(zooave)

#quick plot of CPUE by location
z1 = ggplot(zooave, aes(x=site, y= CPUE))
z1 + geom_bar(stat = "identity", aes(fill = AnalyLS)) + scale_fill_manual(values = mypal) +
  coord_cartesian(ylim = c(0, 7000)) + facet_grid(year~Region2, scales = "free")

#relative abundance
z1 + geom_bar(stat = "identity", aes(fill = AnalyLS), position = "fill") + 
  scale_fill_manual(values = mypal) +
  facet_grid(year~Region2, scales = "free_x")


#total CPUE
zootot = group_by(zooblitz0, SampleID, Date, site, year, month, 
                 Region2, sitetype2) %>% summarize(CPUE = sum(CPUE, na.rm = T))

z2 = ggplot(zootot, aes(x=site, y= CPUE))
z2+geom_boxplot()
#Bradmoor had crzy high CPUE.
z2+geom_boxplot() + coord_cartesian(ylim = c(0,5000))


#average total cpue for bar plot

zoototave = group_by(zootot, site, Region2, sitetype2, year) %>% 
  summarize(mCPUE = mean(CPUE, na.rm = T), seCPUE = sd(CPUE)/length(CPUE), N = length(CPUE))

z3 = ggplot(zoototave, aes(x=site, y=mCPUE))
z3 + geom_bar(stat = "identity", aes(fill = site)) +
  facet_grid(year~Region2, scales = "free", space = "free_x") + 
  geom_errorbar(aes(ymin = mCPUE - seCPUE, ymax = mCPUE + seCPUE)) +
  scale_fill_manual(values = mypal, guide = "none")+
  geom_label(aes(label = paste("n = ", N), y = mCPUE+ 300), label.padding = unit(0.1, "lines"), size = 4) +
   ylab("CPUE") + xlab(label = NULL)+
coord_cartesian(ylim = c(0,10000))

#GLMm of total CPUE
zootot$sitetype2 = as.factor(zootot$sitetype2)
zblitz = lmer(log(CPUE) ~ Region2 + sitetype2 + (1|site), data = zootot)
summary(zblitz)
visreg(zblitz)
#Differences by site type, not region
zblitza = glm(log(CPUE) ~ Region2 + sitetype2, data = zootot)
summary(zblitza)
visreg(zblitza)

zblitzb = aov(log(CPUE) ~ Region2 + sitetype2 + site,  data = zootot)
summary(zblitzb)
TukeyHSD(zblitzb)
#################################################################################

#calculate coefficient of variation in CPUE
zCVs = summarize(group_by(zootot, site), mCPUE = mean(CPUE), 
                sdCPUE = sd(CPUE), seCPUE = sd(CPUE)/length(CPUE), N = length(SampleID), CV = sdCPUE/mCPUE)

#Calculate average withing-site CV and the CV of within-site means by sampling type.
zmCVs = summarize(zCVs, mCV = mean(CV, na.rm = T), mmCPUE = mean(mCPUE), 
                 sdmCPUE= sd(mCPUE), CV2 = sdmCPUE/mmCPUE)


############################################################################################

#plot of community composition 
z4 = ggplot(zooblitz, aes(x=site, y=CPUE))
z4 + geom_bar(stat = "identity", aes(fill = AnalyLS), position = "fill") +
  facet_grid(year~Region2, scales = "free", space = "free_x") + 
  scale_fill_manual(values = c(mypal, "grey"), name = NULL) +
   labs(y="Relative Abudance", x=NULL)

#plot of community composition by region
z4 = ggplot(zooblitz, aes(x=Region2, y=CPUE))
z4 + geom_bar(stat = "identity", aes(fill = AnalyLS), position = "fill") +
  
  scale_fill_manual(values = c(mypal, "grey"), labels = zooLab, name = NULL) +
  ylab("Relative Abudance") +xlab(label = NULL)


#plot of community composition by site type
z4 = ggplot(zooblitz, aes(x=sitetype2, y=CPUE))
z4 + geom_bar(stat = "identity", aes(fill = AnalyLS), position = "fill") +
  
  scale_fill_manual(values = c(mypal, "grey"), labels = zooLab, name = NULL) +
   ylab("Relative Abudance") +xlab(label = NULL)

#Community matrix by CPUE
zMat = zoowide[,8:28]
#now by proprotions
zMatp = zMat/rowSums(zMat)

#PerMANOVA
za1 = adonis(zMatp~ Region2 + sitetype2 + year, data = zootot)
za1
#region and site type are both significant, and year

zNMDS = metaMDS(zMatp, try = 1000)
zNMDS
#no convergent solutions

#but let's graph it anyway!
zootot$Region2 = as.factor(zootot$Region2)
source("PlotNMDS.R")
PlotNMDS(zNMDS, data = zootot, group = "Region2")

zootot$sitetype2 = as.factor(zootot$sitetype2)
PlotNMDS2(zNMDS, data = zootot, group = "sitetype2")

#Hm. Last year, Tule red through everyting off, but now it's all a mess anyway, so that's less of an issue.
#I'm just going to run the 2018 data

#######################################################################################################


#Community matrix by CPUE
zMat2018 = filter(zoowide, year != 2017)[,8:27]
#now by proprotions
zMat2018p = zMat2018/rowSums(zMat2018)

zootot2018 = filter(zootot, year != 2017)

#PerMANOVA
za118 = adonis(zMat2018p~ Region2 + sitetype2, data = zootot2018)
za118
#region and site type are both significant, 

zNMDS = metaMDS(zMat2018p, try = 1000)
zNMDS

zootot2018$Region2 = droplevels(as.factor(zootot2018$Region2))
zootot2018$sitetype2 = droplevels(as.factor(zootot2018$sitetype2))
source("PlotNMDS.R")
PlotNMDS(zNMDS, data = zootot2018, group = "sitetype2")
PlotNMDS(zNMDS, data = zootot2018, group = "Region2")



########################################################################################################
#Look at catch at Tule Red to figure out what in the "other" category was so dominant

TR = droplevels(filter(zoop, site == "Tule Red"))
TR1 = ggplot(TR, aes(x=CommonName, y = CPUE))
TR1 + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))

TR2 = ggplot(TR, aes(x=Analy, y = CPUE, fill = CommonName))
TR2 + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))

others = filter(TR, Analy == "other")
TR3 = ggplot(others, aes(x = CommonName, y = CPUE))
TR3 + geom_bar(stat = "identity")

#check Grizzly bay too
griz = droplevels(filter(zoop, site == "Grizzly"))
gr = ggplot(griz, aes(x=CommonName, y = CPUE))
gr + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))

others = filter(griz, Analy == "other")
gr3 = ggplot(others, aes(x = CommonName, y = CPUE))
gr3 + geom_bar(stat = "identity")
