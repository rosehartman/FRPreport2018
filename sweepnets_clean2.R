#Look at just the sweep net samples
source("blitzclean.R")

#Note: for the sweep net analysis, i'm leaving mesozooplankton in, 
#because we have no other measure of mesozoops in vegetation
#subset teh samples from 2017 (Phase III)
sweep <- filter(inverts2, year(Date) == 2017)

sweepx = summarize(group_by(sweep, SampleID, Station, Location, Sampletype), total = sum(TotalCount))


#add some more info
sweep$Location = NULL
sweep$Station = NULL
sweep$Date = NULL
sweep = merge(sweep, targets, by = "SampleID")
sweep = filter(sweep, Geartype == "sweep net")


sweep.e = summarize(group_by(sweep, SampleID, Station, targets,VegWeight,
                                Region, Geartype, date, NetMeterSTart, NetMeterEnd, 
                                LatStart, LongStart, LatEnd, LongEnd, duration),
                       mean.en= mean(NetMeterSTart))


# Effort for sweepnets is vegwieght
sweep.e$effort = NA
sweep.e$effort[which(sweep.e$Geartype=="sweep net" & sweep.e$VegWeight>0)] =
  sweep.e$VegWeight[which(sweep.e$Geartype=="sweep net" & sweep.e$VegWeight>0)]


#everything else has an effort of "1"

sweep.e$effort[which(is.na(sweep.e$effort))] = 1

#attach the efforts to the origional dataset
sweep= merge(sweep, sweep.e)


# Adjust total count for subsampling
sweep$atotal = sweep$TotalCount*(1/(sweep$subsampled/100))

#Calculate CPUE
sweep$CPUE = sweep$atotal/sweep$effort

#add some analysis groups

sweep = merge(sweep, zoocodes)
sweep = sweep[order(sweep$SampleID),]

#attach info on the sites
#import some info on the sites
sitetypes = read_excel("sitetypes.xlsx")

sweep = merge(sweep, sitetypes)



#Summarize by sample and calculate the total CPUE of all the macroinvertebrates in the sample
sweep.1 = summarize(group_by(sweep, SampleID, Station, label, Location, targets, 
                                Region, Geartype, SiteType, SiteType2, date, effort),
                       tcount = sum(atotal), tCPUE = sum(CPUE), richness = length(unique(Analy)))

sweep.1$Region = as.factor(sweep.1$Region)
sweep.1$SiteType2 = as.factor(sweep.1$SiteType2)
sweep.1$targets = as.factor(sweep.1$targets)

#add the zeros back in so I can take an appropriate mean
ComMats = dcast(sweep, formula = SampleID~Analy, value.var="CPUE", 
               fun.aggregate = sum, fill = 0)

ComMat.1s = dcast(sweep, formula = SampleID~Analy2, value.var="CPUE", 
                 fun.aggregate = sum, fill = 0)


sweep.2x = melt(ComMat.1s, id.vars = "SampleID", variable.name = "Analy2", value.name = "CPUE")
sweep.2x = merge(sweep.1, sweep.2x)

#means by analysis group

sweepsum.1 = summarize(group_by(sweep.2x, Location, label, SiteType, SiteType2, targets, 
                               Region, Geartype, Analy2), mCPUE = mean(CPUE),
                      sdCPUE = sd(CPUE), seCPUE = sd(CPUE)/length(CPUE))


#total CPUE
sweeptot = summarize(group_by(sweep.2x, SampleID, Location, label, SiteType, SiteType2, targets, 
                             Region, Geartype), tCPUE = sum(CPUE))
sweeptotave = summarize(group_by(sweeptot, Location, SiteType, SiteType2, targets, label, 
                                Region), mCPUE = mean(tCPUE), 
                       sdCPUE = sd(tCPUE), seCPUE = sd(tCPUE)/length(tCPUE), N = length(SampleID))

#I exported "sweeptotave" and added the zeros back in excel, 'cause it was easier
#write.csv(sweeptotave, "sweeptotave.csv)
sweeptotave = read.csv("sweeptotave.csv")
######################################################################################################
#Now some plots...

#community composition
n2 = ggplot(sweep, aes(x= label, y = CPUE))
n2 + geom_bar(stat = "identity", aes(fill = Analy2), position = "fill") + 
  facet_grid(targets~Region, scales = "free", space = "free_x") + 
  scale_fill_manual(values = mypal) + mytheme + 
  ylab("Relative percent composition") + xlab("Site")

#CPUE
n3 = ggplot(sweeptotave, aes(x= LocLab, y = mCPUE))
n3 + geom_bar(stat = "identity") + 
  facet_grid(targets~SiteType2, scales = "free", space = "free_x") + 
  geom_errorbar(aes(ymin = mCPUE - seCPUE, ymax = mCPUE + seCPUE), width = 0.7) +
  geom_label(aes(label = paste("n = ", N), y = mCPUE*1.1), label.padding = unit(0.1, "lines")) +
  scale_fill_manual(values = mypal)



##################################################################################################

#community data matrix
ComMats2 = dcast(sweep.2x, formula = SampleID~Analy2, value.var="CPUE", 
                fun.aggregate = sum, fill = 0)
row.names(ComMats2) = ComMats2$SampleID
ComMat2s2 = ComMats2[,2:ncol(ComMats2)]
#make one based on proportion of total catch for each sample
ComMatps2 = ComMat2s2/rowSums(ComMat2s2)

#There was only one critter in one of the samples, and the outlier is majorly throwing off
#my analysis, so I'm going to remove it.
ComMatps2 = ComMatps2[which(sweep.1$SampleID != "EAV2-14FEB2017"),]
sweepdat = filter(sweep.1, SampleID != "EAV2-14FEB2017")


#Permanova

pb3x = adonis(ComMatps2~targets + SiteType2 + Region + Location, data = sweepdat)
pb3x


#NMDS
NMDS2b = metaMDS(ComMatps2, trymax = 2000)
#convered after a whole bunch of tries

source("plotNMDS.R")

#Plot it by site type
PlotNMDS(NMDS2b, data = sweepdat, group = "SiteType2")
#test for significance of the hulls
envfit(NMDS2b~SiteType2, data = sweepdat)

#plot it by region of hte estuary
PlotNMDS(NMDS2b, data = sweepdat, group = "Region")
#test for significance of the hulls
envfit(NMDS2b~Region, data = sweepdat)

#plot it by vegetation type
PlotNMDS(NMDS2b, data = sweepdat, group = "targets")
#test for significance of the hulls
envfit(NMDS2b~targets, data = sweepdat)


#######################################################################################
#some total CPUE models
sweeptots = group_by(sweep.1, targets, Region, SiteType, SiteType2) %>%
  summarize(mCPUE = mean(tCPUE), sdCPUE= sd(tCPUE), seCPUE = sd(tCPUE)/length(tCPUE))

s1 = lmer(log(tCPUE) ~ Region + SiteType2 + targets + (1|Location), data = sweep.1)
summary(s1)
#Wetalnds had higher catch (tule red), SAV had higher catch
visreg(s1)

s1.1 = glm(round(tCPUE) ~ Region + SiteType2 + targets, data = sweep.1, family = "poisson")
summary(s1.1)
#Meh.

sp = ggplot(sweeptots, aes(x=SiteType2, y = mCPUE))
sp + geom_bar(stat = "identity") + facet_grid(Region~targets, scale = "free") + 
  geom_errorbar(aes(ymin = mCPUE - seCPUE, ymax = mCPUE + seCPUE))

##########################################################################################
