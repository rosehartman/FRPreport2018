#Look at chlorophyll data. I think it will be a better indicator of
#productivity than the phytoplankton microscopy, since the 
#phytoplankton grabs were so variable.



library(ggplot2)
library(dplyr)
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


#set up a color pallete
mypal = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"), "black")

#set up a theme
mytheme = theme(strip.text.x = element_text(size=14),strip.text.y = element_text(size=18),
                axis.text.x  = element_text(angle=90, hjust=1, vjust = 0.5, size=12),
                axis.title.y  = element_text(size=16),
                axis.title.x  = element_text(size = 16), axis.text.y = element_text(size = 14),
                legend.text = element_text(size = 14),
                legend.title=element_blank(), legend.background =element_blank()) 



#import the data

chl <- read_excel("WaterQuality.xlsx", 
                             col_types = c("date", "text", "text", 
                                          "text", "text", "text", "numeric"))

#average the values and make a quick graph
chlave = group_by(chl, Location, Region, Sitetype) %>% 
  summarise(mchla = mean(concentration), sechla = sd(concentration)/length(concentration))

ch = ggplot(chlave, aes(x= Location, y = mchla))
ch + geom_bar(stat = "identity", aes(fill = Location)) +
  facet_grid(.~Region, scales = "free_x", space = "free") + 
  scale_fill_manual(values = mypal, guide = "none") +
  xlab("Site") + ylab("Chla concentration ug/L") +
  geom_errorbar(aes(ymin = mchla- sechla, ymax = mchla + sechla), group = "Region", width = 0.5) +
  mytheme

#GLMm of chlorohpyll concentration versus region, site type, and location
c1 = lmer(concentration~ Region + Sitetype + (1|Location), data = chl)
summary(c1)
#Well, tule red is higher than everything else. 

TRGZ = filter(chl, Location == "Grizzly"| Location == "Tule Red")

TRchl = ggplot(TRGZ, aes(x= Location, y = concentration))
TRchl + geom_boxplot( aes(fill = Location)) +
  scale_fill_manual(values = mypal, guide = "none") +
  xlab("Site") + ylab("Chla concentration ug/L") +
  mytheme

#################################################################################

#calculate coefficient of variation in CPUE
cCVs = summarize(group_by(chl, Location), mCPUE = mean(concentration), 
                 sdCPUE = sd(concentration), seCPUE = sd(concentration)/length(concentration), 
                 N = length(concentration), CV = sdCPUE/mCPUE)

#Calculate average withing-site CV and the CV of within-site means by sampling type.
cmCVs = summarize(cCVs, mCV = mean(CV, na.rm = T), mmCPUE = mean(mCPUE), 
                  sdmCPUE= sd(mCPUE), CV2 = sdmCPUE/mmCPUE)

