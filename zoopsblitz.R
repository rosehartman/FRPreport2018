#zooplankton
library(dplyr)
library(ggplot2)
library(readxl)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(simr)
library(lme4)
library(lmerTest)
################################################################################################################################

zoopsfrp <- read_excel("zoopsfrp.xlsx", col_types = c("numeric", 
                                                      "text", "text", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", "numeric", 
                                                      "text", "text", "text", "text", "numeric", 
                                                      "numeric", "numeric", "date", "text", 
                                                      "text", "text", "text", "text", "text", 
                                                      "text", "numeric"))
View(zoopsfrp)

################################################################################################################################
#Calculate CPUE

# Effort for benthic and oblique trawls is based on flowmeter readings and net mouth area

#combine the slides
totalcells = group_by(zoopsfrp, ZSampleID) %>% 
  summarise(totcells = max(CellNumber), totCountall = sum(Count))

zoo = group_by(zoopsfrp, Station, Location, Date, ZSampleID, Dilution, ZooCode, CommonName,
                Phylum, Class, Order, Family,Genus, Species, LifeStage, MeterStart2, MeterEnd2) %>% 
  summarise(totCount = sum(Count))

zoo = merge(zoo, totalcells)

#adjust for subsampling
zoo$atotal = (zoo$totCount/zoo$totcells)*zoo$Dilution

#Volume sampled
zoo$volume = (zoo$MeterEnd2-zoo$MeterStart2)*0.026873027*0.0167


#volume for samples where the flowmeter turned over to zero.
zoo$volume[which(zoo$volume<0)] = ((1000000 - zoo$MeterStart2[which(zoo$volume<0)])+
                                     zoo$MeterEnd2[which(zoo$volume<0)])*0.026873027*0.0167 #I've used the factory calibration here

#some of the samples iddn't have flowmeters, so we'll use the average volume for those
zoo$volume[which(is.na(zoo$volume))] = mean(zoo$volume, na.rm = T)

#Calculate CPUE
zoo$CPUE = zoo$atotal/zoo$volume


#figuer out how many critters total
tot = sum(zoo$CPUE)

#Now how many critters per taxon
zooptax = group_by(zoo, CommonName) %>% summarise(totaltaxa = sum(CPUE, na.rm = T), 
                                                  prop = totaltaxa/sum(zoo$CPUE, na.rm = T), 
                                                  tax = CommonName[1])
frptaxa <- read_excel("FRP_EMPcrosswalk.xlsx", 
                      sheet = "allfrptaxa")

zoop = merge(zoo, frptaxa)

#some critters have barely any individuals, lets combine groups
zoop$Analy[which(zoop$Analy == "platyhelminthes" | zoop$Analy == "mollusca" | 
                   zoop$Analy=="terrestrial" | zoop$Analy == "cumaceans")] = "other"
zoop$Analy[which(zoop$Analy == "corophium" | zoop$Analy == "gammarid" )] = "Amphipod"

zoop$AnalyLS = paste(zoop$Analy, zoop$lifestage)
zoop$AnalyLS[which(zoop$AnalyLS == "other adult" | zoop$AnalyLS == "other juv" | zoop$AnalyLS == "insect adult")] = "other"

#set up some labels for graphing
zooLab = c("Amphipoda", "Annelida", "Calanoida", "Cal juv", "Cal nauplii","Cladocera", "Collembola", 
           "Cyclo nauplii", "Cyclopoda", "Cyclo juv", "fish", "Harpacticoida", "Insect juv", "Isopoda", "Mysidea", 
           "Ostracoda", "other", "Rotifera", "Tanaid")
zoop$AnalyLS = factor(zoop$AnalyLS, labels=c("Amphipoda", "Annelida", "Calanoida", "Cal juv", "Cal nauplii","Cladocera", "Collembola", 
                                                "Cyclo nauplii", "Cyclopoda", "Cyclo juv", "fish", "Harpacticoida", "Insect juv", "Isopoda", "Mysidea", 
                                                "Ostracoda", "other", "Rotifera", "Tanaid"))


#Now let's look at what we got!
#############################################################################################################

#filter for just samples from the blitz
samplesBlitz <- read_excel("samplesBlitz.xlsx",  sheet = "blitzzoops", 
                                      col_types = c("text", "text", "text", "text", "text", "date"))
zoo2017 = merge(select(zoop,-c(Station, Location, Date)), samplesBlitz)
zoo2017 = merge(zoo2017, sitetypes)


#set up a color pallete
mypal = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"), "black", "white", "grey", "pink")



#quick plot of CPUE by location
zoo2017ave = group_by(zoo2017, Analy, CommonName, Location, AnalyLS) %>% summarize(CPUE = mean(CPUE, na.rm = T))
zoo2017ave =droplevels(zoo2017ave)
z1 = ggplot(zoo2017ave, aes(x=Location, y= CPUE))
z1 + geom_bar(stat = "identity", aes(fill = Analy)) + scale_fill_manual(values = mypal)

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

zblitz = lmer(log(CPUE) ~ Region + SiteType2 + (1|Location), data = zootot)
summary(zblitz)
visreg(zblitz)
#Differences by region, but not by habitat type. Opposite of macroinvertebrates!!!!


###########################################################################################
#Power analysis
library(pwr)

#look at everything
#power with given sample size
snf = summary(zblitz1)$r.squared[1]/(1-summary(zblitz1)$r.squared[1])
pwr.f2.test(summary(zblitz1)$df[1], v = zblitz1$df.residual, f2 = snf, sig.level = 0.05)
#sample size at 80% power
pwr.f2.test(summary(zblitz1)$df[1],  f2 = snf, sig.level = 0.05, power = 0.8)
#I think this is telling us we need 31 samples to get some significant differences
#but we probably need to look factor by factor

zblitz3 = lm(log(CPUE) ~ Region, data = zootot)
zblitz5 = lm(log(CPUE) ~ SiteType2, data = zootot)

#differences between regions
#power with given sample size
Rf = summary(zblitz3)$r.squared[1]/(1-summary(zblitz3)$r.squared[1])
pwr.f2.test(summary(zblitz3)$df[1], v = zblitz3$df.residual, f2 = Rf, sig.level = 0.05)
#sample size at 80% power
pwr.f2.test(summary(zblitz3)$df[1],  f2 = Rf, sig.level = 0.05, power = 0.8)
#We need 36 samples for enough power to differentiate between regions

#differences between site types
#power with given sample size
stf = summary(blitz5)$r.squared[1]/(1-summary(blitz5)$r.squared[1])
pwr.f2.test(summary(blitz5)$df[1], v = blitz5$df.residual, f2 = stf, sig.level = 0.05)
#sample size at 80% power
pwr.f2.test(summary(blitz5)$df[1],  f2 = stf, sig.level = 0.05, power = 0.8)
#We need 80 samples for enough power to differentiate between site types (when there are 4 site types)

#####################################################################################################################
#try a differen tpower analysis
zpc1 = powerCurve(zblitz)
zpp1 = powerSim(zblitz, nsim = 299)

#OK, it looks like we shouldn't use the observed effect size, because results can be misleading. Instead 
# we should decide what effect we are interested in. Let's say we want to be able to tell the
#difference between 0.5 differrence in log CPUE
fixef(zblitz)
fixef(zblitz)= c(7, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5)
zpc2 = powerCurve(zblitz, nsim = 299)
plot(zpc2)
zpp2 = powerSim(zblitz, nsim = 299)
#So we are pretty good on power.
zblitz2 = extend(zblitz, within = "Region + SiteType2", n = 10)
zpc2.1 = powerCurve(zblitz2, nsim = 299, within = "Region + SiteType2", breaks = c(1,2,3,4,5,6,7,8,9,10))
zpc2.1
#nice!
zpp = powerSim(zblitz, test = fixed("SiteType2"), nsim = 200)
zpp = powerSim(zblitz, test = fixed("Region"), nsim = 200)
zpp

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
PlotNMDS(zNMDS2, data = filter(zootot, Location != "Tule Red"), group = "Region")
envfit(zNMDS2~Region, data = filter(zootot, Location != "Tule Red"))

PlotNMDS(zNMDS2, data = filter(zootot, Location != "Tule Red"), group = "SiteType2")
envfit(zNMDS2~SiteType2, data = filter(zootot, Location != "Tule Red"))

