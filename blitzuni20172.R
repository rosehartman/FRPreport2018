#Linear MOdels and power analysis from teh 2017 and 2018 blitz
source("blitzclean.R")

################################################################################

#Can we do a total CPUE model of everythingat once?

#Because I have more than one site for some of the site types in some of the regions, I need "site" as an error 
#term in the model, or include it as a variable. I'm going to go for a mixed model
#since I'm more interetsed in the effect of Region2 and Site type than I am of site.

#if I try and do both years together, I probably want year in there too

#Linear mixed model of log total CPUE, as predicted by habitat type (targets2), Region2, and Site Type. 
#I'm going to take out the benthic samples since we changed the way we do those

blitz2 = lmer(logtot ~ targets2+ sitetype+ Region2  + Year2+ (1|site), 
              data = bugstotNoB, na.action = "na.fail")
summary(blitz2)
visreg(blitz2)

blitz2a = glm(logtot ~ targets2+ sitetype+ Region2+ Year  +site, 
             data = bugstotNoB, na.action = "na.fail")
summary(blitz2a)
visreg(blitz2a)


blitz2b = aov(log(tCPUE+1) ~ targets2 + sitetype+ Region2 +  Year2 +Error(site), 
              data = filter(bugstot, targets2 != "benthic"), na.action = "na.fail")
summary(blitz2b)
TukeyHSD(blitz2b)

#if I really want to pick out the differences between site types, maybe I should try getting
#rid of all the extra stuff
blitz2c = aov(log(tCPUE+1) ~ targets2 +sitetype + site, 
              data = filter(bugstot, targets2 != "benthic"))
summary(blitz2c)
TukeyHSD(blitz2c)

blitz2d = lmer(log(tCPUE +1) ~ sitetype+targets2+  (1|site), 
              data = filter(bugstot, targets2 != "benthic"))
summary(blitz2d)
visreg(blitz2d)

#That came out a lot better than I expected! 
#wetlands had higher abundance than channel, no other major differences

#but let's throw some AIC at it
dredge(blitz2b)
#the top four models all had identical AIC score
dredge(blitz2)


###########################################################################################
#Now I'll make a plot of total CPUE versus site and Region2

tot2 = ggplot(filter(bugstotave2, targets2 != "benthic"), aes(x = site, y = mCPUE))
tot2 + geom_bar(stat = "identity", aes(fill = site)) + 
  facet_grid(targets2~Region2, scales = "free", space = "free_x") +
  geom_errorbar(aes(ymin = mCPUE - seCPUE, ymax = mCPUE + seCPUE), width = 1) +
  geom_label(aes(label = paste("n = ", N), y = 0), label.padding = unit(0.1, "lines"), size = 3) +
  scale_fill_manual(values = mypal, guide = "none") + ylab("CPUE") 

#I should seperate by year...
bugstotave3 = summarize(group_by(bugstot, site, sitetype, targets2, 
                                 Region2, Year), mCPUE = mean(tCPUE, na.rm = T), 
                        sdCPUE = sd(tCPUE, na.rm = T), seCPUE = sd(tCPUE)/length(tCPUE), N = length(SampleID))


tot3 = ggplot(filter(bugstotave3, targets2 != "benthic"), aes(x = site, y = mCPUE))
tot3 + geom_bar(stat = "identity", aes(fill = site)) + 
  facet_grid(targets2~Year, scales = "free", space = "free_x") +
  geom_errorbar(aes(ymin = mCPUE - seCPUE, ymax = mCPUE + seCPUE), width = 0.7) +
  geom_label(aes(label = paste("n = ", N), y = mCPUE*1.1), label.padding = unit(0.1, "lines"), size = 3) +
  scale_fill_manual(values = mypal, guide = "none") + ylab("CPUE") 

#That's actually not super helpful, because we didn't see significant difference
#between years. Let's just do by site types

levels(bugstotave2$targets2) = c("benthic", "mysid", "neuston", "sweep net")
bugstotave2$Region2 = factor(bugstotave2$Region2, levels = c("Grizzly", "NurseDenverton", "Confluence", "SacSanJ", "Cache"))

tot4 = ggplot(filter(bugstotave2, targets2 != "benthic"), aes(x = site, y = mCPUE))
tot4 + geom_bar(stat = "identity", aes(fill = Region2)) + 
  facet_grid(targets2~sitetype, scales = "free", space = "free_x") +
  geom_errorbar(aes(ymin = mCPUE - seCPUE, ymax = mCPUE + seCPUE), width = 0.7) +
  geom_label(aes(label = paste("n = ", N), y = 0), label.padding = unit(0.1, "lines"), size = 3) +
  scale_fill_manual(values = mypal, name = "Region") + ylab("CPUE") +
  theme(strip.text = element_text(size = 16), axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16), legend.position = "bottom")



####################################################################################
#Look at the coefficient of variation within sites for each sample type

CVs = summarize(group_by(bugstot, site, SiteType, sitetype, targets2, 
                         Region2), mCPUE = mean(tCPUE), 
                sdCPUE = sd(tCPUE), seCPUE = sd(tCPUE)/length(tCPUE), N = length(SampleID), CV = sdCPUE/mCPUE)

#Calculate average withing-site CV and the CV of within-site means by sampling type.
mCVs = summarize(group_by(CVs, targets2), mCV = mean(CV, na.rm = T), mmCPUE = mean(mCPUE), 
                 sdmCPUE= sd(mCPUE), CV2 = sdmCPUE/mmCPUE)

#Repeat the analysis with the vegetation types sepearated
CVs2 = summarize(group_by(bugstot, site, SiteType, sitetype, targets, 
                         Region2), mCPUE = mean(tCPUE), 
                sdCPUE = sd(tCPUE), seCPUE = sd(tCPUE)/length(tCPUE), N = length(SampleID), CV = sdCPUE/mCPUE)

#Calculate average withing-site CV and the CV of within-site means by sampling type.
mCVs2 = summarize(group_by(CVs, targets), mCV = mean(CV, na.rm = T), mmCPUE = mean(mCPUE), 
                 sdmCPUE= sd(mCPUE), CV2 = sdmCPUE/mmCPUE)

