#Quick look at the past years' nutrients

library(dplyr)
library(ggplot2)
library(reshape2)
library(lubridate)
library(RColorBrewer)
library(readxl)
library(visreg)
library(lme4)
library(lmerTest)
library(MuMIn)


#import our data
#First import all the anlyitical results
FRP2018 = read_xlsx("WQData2.xlsx", sheet = "WaterQuality")

#Now import some of the environmental information associated with th esamples
FRP2018samples = read_xlsx("WQData2.xlsx", sheet = "sampleinfo")

#Merge the datasets
FRP2018x = merge(FRP2018[,c(4:8)], FRP2018samples)


#put the YSI clorophyll and fDOM as their own constituants
YSI = group_by(FRP2018x, bryteID, site, sitetype, substa, Date) %>% 
  summarize(YSIchl = mean(YSIChl, na.rm = T),
            fDOM = mean(FDOM, na.rm = T))

YSIlong = melt(YSI, id.vars = c("bryteID", "site", "sitetype", "substa", "Date"), 
               variable.name = "analyte", value.name = "concentration")

#Get rid of the extraneous information
FRP = FRP2018x[,c(1,4,5,8,10:12)]
FRP = rbind(FRP, YSIlong)

#now import EMP's data and atach it
EMPWQ <- read_excel("EMPWQ.xlsx")
nuts2018 = rbind(FRP, EMPWQ)

#Import some more information about the sites
sites = read_xlsx("WQData2.xlsx", sheet = "sites")

#Do a little more data manipulation to clean it up
nuts2018 = merge(nuts2018, sites)
nuts2018$month = month(nuts2018$Date)
nuts2018$day = yday(nuts2018$Date)
nuts2018$substa = factor(nuts2018$substa, levels = c("inside", "breach", "outside", "EMP"))

###################################################################################################
#Now some exploritory plots
#I think it's going to be easiest to look at one analyte at a time

#set up a color palette
mypal = c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"))

#start with nitrate
nutsN = filter(nuts2018, analyte == "Dissolved Nitrate + Nitrite")

#Plot it over time, seperating by site
N = ggplot(nutsN, aes(x=Date, y = concentration))
N + geom_point(aes(col = substa)) +
  geom_line(aes(col = substa)) +
  scale_color_manual(values = mypal) +
  ylab("mg Nitrate + Nitrite") + facet_wrap(~site, scales = "free_y")

#graph just FRP samples
N = ggplot(filter(nutsN, substa != "EMP"), aes(x=Date, y = concentration))
N + geom_point(aes(col = substa)) +
  geom_line(aes(col = substa)) +
  scale_color_manual(values = mypal) +
  ylab("mg Nitrate + Nitrite") + facet_wrap(~site)

#zoom in
N + geom_point(aes(col = substa)) +
  geom_line(aes(col = substa)) +
  scale_color_manual(values = mypal) +
  ylab("mg Nitrate + Nitrite") + facet_wrap(~site)+
coord_cartesian(ylim = c(0, .7)) 

#facet by site type
N + geom_point(aes(col = site, pch = substa)) + 
  scale_color_manual(values = mypal) +
  ylab("mg Nitrate + Nitrite") + facet_grid(sitetype~., scales = "free_y")


#Let's look at a the summaries over time. First compute the averages and standard error
nuts2018$Date = date(nuts2018$Date)
nutsave = group_by(nuts2018, Date, site, sitetype, sitetype2, region, analyte) %>% 
  summarize(mcon = mean(concentration), seconcentration = sd(concentration)/length(concentration))
aveN = filter(nutsave, analyte == "Dissolved Nitrate + Nitrite")                     

#graph it
N2 = ggplot(aveN, aes(x=Date, y = mcon, col = site))
N2 + geom_point() + geom_line() + geom_errorbar(aes(ymin = mcon - seconcentration, ymax = mcon + seconcentration)) +
  scale_color_manual(values = mypal) + ylab("mg Nitrate + Nitrite")

#zoom in
N2 + geom_point() + geom_line() + coord_cartesian(ylim = c(0, 0.75))+ 
  geom_errorbar(aes(ymin = mcon - seconcentration, ymax = mcon + seconcentration)) +
  scale_color_manual(values = mypal)+ ylab("mg Nitrate + Nitrite")

#facet by site type
N2 + geom_point() + geom_line() +  
  geom_errorbar(aes(ymin = mcon - seconcentration, ymax = mcon + seconcentration)) +
  scale_color_manual(values = mypal)+ ylab("mg Nitrate + Nitrite") + 
  facet_grid(sitetype~., scales = "free_y") + 
  theme_bw()

#Try pairing restoration and reference sites
N3 = ggplot(aveN, aes(x=Date, y =mcon, color = site)) 
N3 + geom_line() + geom_point() + facet_grid(region~., scales = "free_y")

#Try pairing restoration and reference sites
nutsN$site = factor(nutsN$site, levels = c("Liberty", "Prospect", "Lindsey", "Winter", "Browns", "Stacys", 
                                           "Decker", "Blacklock", "Grizzly", "Wings"))
N3 = ggplot(filter(nutsN, substa != "EMP"), aes(x=Date, y =concentration, color = site)) 
N3  + geom_point(aes(pch = substa)) + facet_grid(region~., scales = "free_y")+
  scale_color_manual(values = mypal) +geom_smooth()

#just look at Decker versus Stacys
N4 = ggplot(filter(nutsN,  region =="SacSanJ"), aes(x=Date, y =concentration, color = site)) 
N4  + geom_point(aes(pch = substa)) + 
  scale_color_manual(values = mypal) + geom_smooth()

##############################################################################################
#NOw check out phosphorus

nutsP = filter(nuts2018, analyte == "Dissolved Ortho-phosphate")
aveP =  filter(nutsave, analyte == "Dissolved Ortho-phosphate")

P = ggplot(nutsP, aes(x=Date, y = concentration))
P + geom_point(aes(col = site)) +
  scale_color_manual(values = mypal) + ylab("dissolved ortho-phosphate")
P2 = ggplot(aveP, aes(x=Date, y = mcon, col = site))
P2 + geom_point() + geom_line() + geom_errorbar(aes(ymin = mcon - seconcentration, ymax = mcon + seconcentration))+
  scale_color_manual(values = mypal) + ylab("dissolved ortho-phosphate")

#just look at Decker versus Stacys
P4 = ggplot(filter(nutsP,  region =="SacSanJ"), aes(x=Date, y =concentration, color = site)) 
P4  + geom_point(aes(pch = substa)) + 
  scale_color_manual(values = mypal) + geom_smooth()

############################################################################################
#NOw check out organic nitrogen

nutsDON = filter(nuts2018, analyte == "Dissolved Organic Nitrogen")
aveDON =  filter(nutsave, analyte == "Dissolved Organic Nitrogen")

DON = ggplot(nutsDON, aes(x=Date, y = concentration))
DON + geom_point(aes(col = site)) +
  scale_color_manual(values = mypal) + ylab("dissolved organic N")
DON2 = ggplot(aveDON, aes(x=Date, y = mcon, col = site))
DON2 + geom_point() + geom_line() + geom_errorbar(aes(ymin = mcon - seconcentration, ymax = mcon + seconcentration))+
  scale_color_manual(values = mypal) + ylab("dissolved organic N")


#just look at Decker versus Stacys
NO4 = ggplot(filter(nutsDON,  region =="SacSanJ"), aes(x=Date, y =concentration, color = site)) 
NO4  + geom_point(aes(pch = substa)) + 
  scale_color_manual(values = mypal) + geom_smooth()

###########################################################################################

#Do a quick GLM to see if there are any statistically significant differences between site types or breach/outside/inside

g1 = glm(concentration~ site+substa + Date, data = nutsN)
summary(g1)
library(visreg)
visreg(g1)

#pairwise comparisons of restoration and reference sites
#Decker versus Stacys
g1ds = glm(concentration~ site+substa + Date, data = filter(nutsN, region == "SacSanJ"))
summary(g1ds)
visreg(g1ds)

#Prospect Versus Liberty
g1lp = glm(concentration~ site*substa + Date,
           data = filter(nutsN, site == "Liberty"| site == "Prospect"))
summary(g1lp)
visreg(g1lp)


g2 = glm(concentration~ sitetype + substa + Date, data = nutsN)
summary(g2)
visreg(g2)

g3 = glm(concentration~ sitetype*substa + Date, data = nutsN)
summary(g3)
visreg(g3)
visreg(g3, xvar = "substa", by = "sitetype")
#######################################################################################
#remove wings 'cause it's wierd

nutsN2 = filter(nutsN, site != "Wings")

#try it again
g3.1 = glm(concentration~ sitetype+substa + Date, data = nutsN2)
summary(g3.1)
visreg(g3.1, xvar = "substa", by = "sitetype")
visreg(g3.1)
visreg(g3.1, xvar = "Date")

g2.1 = glm(concentration~ site +substa + Date, data = nutsN2)
summary(g2.1)
visreg(g2.1)


#I think I need to have site or region in as an error term, or it won't work. 
#Prospect Versus Liberty

g2r = lmer(concentration~ sitetype + substa + day + (1|site), data = nutsN)
summary(g2r)

g3 = lmer(concentration~ sitetype*substa + (1|site) + day, data = nutsN)
summary(g3)
visreg(g3)
visreg(g3, xvar = "substa", by = "sitetype")

#remove wings 'cause it's wierd

nutsN2 = filter(nutsN, site != "Wings")

#try it again
g3.1r = lmer(concentration~ sitetype*substa + day +(1|site), data = nutsN2)
summary(g3.1r)
visreg(g3.1r)
visreg(g3.1r, xvar = "substa", by = "sitetype")


g2.1 = glm(concentration~ site +substa + Date, data = nutsN2)
summary(g2.1)
visreg(g2.1)

######################################################################################
#Compare just sites with EMP to without, remove months with no FRP

nuts2018$month = month(nuts2018$Date)
nutsPaired = filter(nuts2018, site != "Blacklock")
nutsPaired = filter(nutsPaired, month > 2 & month < 11 & month != 7 & month != 8)
nutsPairedN = filter(nutsPaired, analyte == "Dissolved Nitrate + Nitrite")
nutsPairedN$substa = factor(nutsPairedN$substa, levels = c("inside", "breach", "outside", "EMP"))
nutsPNave = group_by(nutsPairedN, month, site, sitetype, sitetype2, region, analyte) %>% 
  summarize(mcon = mean(concentration), seconcentration = sd(concentration)/length(concentration))
nutsPairedN$month = as.factor(nutsPairedN$month)

#try the GLMs again
gp1 = glm(concentration~ site*substa + Date, data = filter(nutsPairedN, site != "Wings"))
summary(gp1)
library(visreg)
visreg(gp1)
visreg(gp1, xvar = "substa", by = "site")

ga1 = aov(concentration~ substa*site + month, data = filter(nutsPairedN, site != "Wings"))
summary(ga1)
TukeyHSD(ga1)

N = ggplot(filter(nutsPairedN, site != "Wings"), aes(x=Date, y = concentration))
N + geom_point(aes(col = site, pch = substa)) +
  scale_color_manual(values = mypal) +
  ylab("mg Nitrate + Nitrite")

ggplot(filter(nutsPNave, site != "Wings"), aes(x=month, y = mcon)) +
geom_point(aes(col = site)) + geom_line(aes(col = site)) +
  scale_color_manual(values = mypal) +
  ylab("mg Nitrate + Nitrite") +
  geom_errorbar(aes(ymin = mcon - seconcentration, ymax = mcon + seconcentration), width = 0.1)+
  scale_color_manual(values = mypal) 

gp2 = glm(concentration~ sitetype + substa + Date, data = nutsPairedN)
summary(gp2)
visreg(gp2)

g3 = glm(concentration~ sitetype*substa + Date, data = nutsN)
summary(g3)
visreg(g3)
visreg(g3, xvar = "substa", by = "sitetype")

#remove wings 'cause it's wierd

nutsP2 = filter(nutsPairedN, site != "Wings")

#try it again
g3.1 = glm(concentration~ sitetype*substa + Date, data = nutsP2)
summary(g3.1)
visreg(g3.1, xvar = "substa", by = "sitetype")
visreg(g3.1)
visreg(g3.1, xvar = "Date")

############################################################################
#Try removing EMP, looking for differences between stations and site types.
NoEMP = filter(nuts2018, substa != "EMP")

NoEMPN = filter(NoEMP, analyte == "Dissolved Nitrate + Nitrite")



#try the GLMs again
gp1 = glm(concentration~ site+substa + Date, data = NoEMPN)
summary(gp1)
visreg(gp1)

gp2 = glm(concentration~ sitetype + substa + Date, data =NoEMPN)
summary(gp2)
visreg(gp2)

g3 = glm(concentration~ sitetype*substa + Date, data = nutsN)
summary(g3)
visreg(g3)
visreg(g3, xvar = "substa", by = "sitetype")

#remove wings 'cause it's wierd

NoEMP2 = filter(NoEMPN, site != "Wings")

#try it again
g3.1 = glm(concentration~ sitetype*substa + Date, data = NoEMP2)
summary(g3.1)
visreg(g3.1, xvar = "substa", by = "sitetype")
visreg(g3.1)
visreg(g3.1, xvar = "Date")

#try something silly
NoEMPN$month2 = as.factor(NoEMPN$month)

gN1 = glm(concentration~ site+substa + month2, data = NoEMPN)
summary(gN1)

gA1 = aov(concentration~ site+substa + month2, data = NoEMPN)
summary(gA1)
TukeyHSD(gA1)

gA2 = aov(concentration~ sitetype+substa + month2, data = NoEMPN)
summary(gA2)
TukeyHSD(gA2)
