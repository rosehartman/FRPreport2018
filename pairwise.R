#stacy wants some pairwise comparisons to include in each project report
source("blitzclean.R")

################################################################################

#Broad slough versus Browns Versus Winter
Conf = filter(bugsblitz, Region2 == "Confluence" & targets2 != "benthic")
Conf$Year = as.factor(year(Conf$Date))
Conf.1 = filter(bugsblitz.1, Region2 == "Confluence"& targets2 != "benthic")
Conf.1x = filter(bugsblitz.2x, Region2 == "Confluence"& targets2 != "benthic")
Confsum = filter(bugssum.1, Region2 == "Confluence"& targets2 != "benthic")
Conflogtot = filter(bugstotNoB, Region2 == "Confluence")

#First, are there differences in catch?
blitzConf = aov(logtot ~ site+ targets2 +Year2, 
              data = Conflogtot)
summary(blitzConf)
TukeyHSD(blitzConf)
visreg(blitzConf)
#Yes! Winter and browns have higher catch than broad slough!
#no differences 2017 v 2018

#now, are there differences in community composition

ComConf = dcast(Conf.1x, formula = SampleID~Analy2, value.var="CPUE", 
                 fun.aggregate = sum, fill = 0)
row.names(ComConf) = ComConf$SampleID
ComConf2 = ComConf[,2:ncol(ComConf)]
Conf2 = Conf.1[which(rowSums(ComConf2 )!=0),]
ComConf2 = ComConf2[which(rowSums(ComConf2)!=0),]
#make one based on proportion of total catch for each sample
ComConf2p = ComConf2/rowSums(ComConf2)


#PERMANOVA with relative abundance
Conf2$Year = as.factor(year(Conf2$Date))
pn1 = adonis(ComConf2p~site + Year + targets2, data = Conf2)
pn1
#NMDS

NMDSconf = metaMDS(ComConf2p, try = 50, trymax = 500)
NMDSconf

#plot it
PlotNMDS(NMDSconf, data = Conf2, group = "site")
#test significance
envfit(NMDSconf~site, data = Conf2)

#plot by year
PlotNMDS(NMDSconf, data = Conf2, group = "Year")
#test significance
envfit(NMDSconf~Year,data = Conf2)

#plot by year
PlotNMDS(NMDSconf, data = Conf2, group = "targets2")
#test significance
envfit(NMDSconf~targets2,data = Conf2)

#More differences by year and gear type than site

#stacked bar plot
ggplot(Conf, aes(x=site, y = CPUE, fill = Analy2)) + 
  geom_bar(stat = "identity", position = "fill")+
  facet_grid(targets2~Year, scales = "free", space = "free_x") +
  scale_fill_manual(values = mypal, name = NULL) + 
  xlab("Site")+ ylab("Relative percent composition") + mytheme

################################################################################################
#Prospect versus Miner Versus Liberty 

#I could add lindsey, but maybe just stick to the ones in the BA.

Cache = filter(bugsblitz, Region2 == "Cache" & site != "Lindsey" & targets2 != "benthic")
Cache.1 = filter(bugsblitz.1, Region2 == "Cache" &  site != "Lindsey"& targets2 != "benthic")
Cache.1x = filter(bugsblitz.2x, Region2 == "Cache"& targets2 != "benthic"&  site != "Lindsey")
Cachesum = filter(bugssum.1, Region2 == "Cache" & site != "Lindsey"& targets2 != "benthic")
Cachelogtot = filter(bugstotNoB, Region2 == "Cache" & site != "Lindsey"& targets2 != "benthic")


#First, are there differences in catch?
blitzCache = aov(logtot ~ site+ targets2 +Year2, 
                data = Cachelogtot)
summary(blitzCache)
TukeyHSD(blitzCache)
#Yes!

#now, are there differences in community composition

ComCache = dcast(Cache.1x, formula = SampleID~Analy2, value.var="CPUE", 
                fun.aggregate = sum, fill = 0)
row.names(ComCache) = ComCache$SampleID
ComCache2 = ComCache[,2:ncol(ComCache)]
Cache2 = Cache.1[which(rowSums(ComCache2 )!=0),]
ComCache2 = ComCache2[which(rowSums(ComCache2)!=0),]
#make one based on proportion of total catch for each sample
ComCache2p = ComCache2/rowSums(ComCache2)


#PERMANOVA with relative abundance
Cache2$Year = as.factor(year(Cache2$Date))
pn1 = adonis(ComCache2p~site + Year + targets2, data = Cache2)
pn1
#NMDS

NMDSCache = metaMDS(ComCache2p, try = 50, trymax = 1000)
NMDSCache
#it doesn't want to converge

#plot it
PlotNMDS(NMDSCache, data = Cache2, group = "site")
#test significance
envfit(NMDSCache~site, data = Cache2)

#plot by year
PlotNMDS(NMDSCache, data = Cache2, group = "Year")
#test significance
envfit(NMDSCache~Year,data = Cache2)

#plot by gear type
PlotNMDS(NMDSCache, data = Cache2, group = "targets2")
#test significance
envfit(NMDSCache~targets2,data = Cache2)


####################################################################################################
#Decker versus Horseshoe bend versus Stacy's

#Tule Red versus Ryer versus Grizzly

#Bradmoor versus Blacklock versus LHB

##################################################################################################
#can I write a function that makes all these graphs given the name of a region?

#start with macro inverts, just univariates
Regionalbugs = function(region, data) {
  require(tidyverse)
  require(reshape2)
  require(lubridate)
  require(visreg)
  #filter the whole blitz data set for the region of interest
  Reg = filter(data, Region2 == region & targets2 != "benthic")
  #turn tyear into a factor
  Reg$Year = as.factor(year(Reg$Date))
  
  #make a summary data set
  Reg.1 = summarize(group_by(Reg, SampleID, Station, Region2, Target, targets2, 
                                   Region, Sampletype, site, sitetype, Date),
                          tcount = sum(atotal, na.rm = T), tCPUE = sum(CPUE, na.rm = T), 
                    richness = length(unique(Analy)))
  
  #convert from strings to factors
  Reg.1$site = as.factor(Reg.1$site)
  Reg.1$sitetype = as.factor(Reg.1$sitetype)
  Reg.1$targets2 = as.factor(Reg.1$targets2)
  
  #I need to add the zeros in for taxa that were not present in a sample so I can take an appropriate mean
  #I usually do this by transfering it from "long" to "wide" and back again.
  
  #Make it wide
  RegMat.1 = dcast(Reg, formula = SampleID~Analy2, value.var="CPUE", 
                   fun.aggregate = sum, fill = 0)
  
  Reg.1x = melt(RegMat.1, id.vars = "SampleID", variable.name = "Analy2", value.name = "CPUE")
  Reg.1x = merge(Reg.1, Reg.1x)
  
  #Means by analysis group for each location
  Regsum = summarize(group_by(Reg.1x, site, Region2, sitetype, targets2, 
                                 Region, Sampletype, Analy2), mCPUE = mean(CPUE),
                        sdCPUE = sd(CPUE), seCPUE = sd(CPUE)/length(CPUE))
  
  #Now calculate total CPUE
  Regtot = summarize(group_by(Reg.1x, SampleID, Date, site, sitetype, Target, targets2, 
                               Region, Region2, Sampletype), tCPUE = sum(CPUE, na.rm = T))
  
  
  Regtot$Year = year(Regtot$Date)
  Regtot$Year2 = as.factor(Regtot$Year)
  Regtot$logtot = log(Regtot$tCPUE +1)
  RegtotNoB = filter(Regtot, targets2!= "benthic")
  
  #Mean CPUE, sample size, and standard error for each location and habitat so I can make a pretty graph.
  Regtotave = summarize(group_by(Regtot, site, sitetype, targets2, 
                                  Region2), mCPUE = mean(tCPUE, na.rm = T), 
                         sdCPUE = sd(tCPUE, na.rm = T), seCPUE = sd(tCPUE)/length(tCPUE), N = length(SampleID))
  
  #Some sites didn't have samples for a particular habitat type, so I added rows to make it easier to graph
  sitetarget = expand.grid(site=unique(Regtotave$site), targets2 = unique(Regtotave$targets2))
  foo = merge(sitetypes[,c(4,6,7)], sitetarget)
  foo = unique(foo)
  
  Regtotave2 = merge(foo, Regtotave, all.x = T)
  Regtotave2$mCPUE[which(is.na(Regtotave2$mCPUE))] = 0
  Regtotave2$N[which(is.na(Regtotave2$N))] = 0
  Regtotave2$site = factor(Regtotave2$site, levels = c("Ryer", "Grizzly", "Tule Red", "Blacklock", "Bradmoor", "LHB", "Browns",
                                                         "Winter", "Broad", "Horseshoe", "Decker", "Stacys", "Miner", "Prospect",
                                                         "Liberty", "Lindsey", "Flyway"))
  
  
  p1 = ggplot(Regtotave2, aes(x=site, y = mCPUE, fill = site))
  print(p1+geom_bar(stat = "identity") + 
          geom_errorbar(aes(ymin = mCPUE - seCPUE, ymax = mCPUE + seCPUE))+
      facet_grid(targets2~., scales = "free_y"))
  
  #First, are there differences in catch?
  blitzReg = aov(logtot ~ site+ targets2 +Year2, 
                  data = RegtotNoB)
  print(summary(blitzReg))
  print(TukeyHSD(blitzReg))
  visreg(blitzReg)
}


#tets
Regionalbugs("SacSanJ", bugsblitz)
Regionalbugs("Cache", bugsblitz)
Regionalbugs("Grizzly", bugsblitz)
Regionalbugs("Confluence", bugsblitz)

  #now, are there differences in community composition
  , scales = 
  ComConf = dcast(Conf.1x, formula = SampleID~Analy2, value.var="CPUE", 
                  fun.aggregate = sum, fill = 0)
  row.names(ComConf) = ComConf$SampleID
  ComConf2 = ComConf[,2:ncol(ComConf)]
  Conf2 = Conf.1[which(rowSums(ComConf2 )!=0),]
  ComConf2 = ComConf2[which(rowSums(ComConf2)!=0),]
  #make one based on proportion of total catch for each sample
  ComConf2p = ComConf2/rowSums(ComConf2)
  
  
  #PERMANOVA with relative abundance
  Conf2$Year = as.factor(year(Conf2$Date))
  pn1 = adonis(ComConf2p~site + Year + targets2, data = Conf2)
  pn1
  #NMDS
  
  NMDSconf = metaMDS(ComConf2p, try = 50, trymax = 500)
  NMDSconf
  
  #plot it
  PlotNMDS(NMDSconf, data = Conf2, group = "site")
  #test significance
  envfit(NMDSconf~site, data = Conf2)
  
  #plot by year
  PlotNMDS(NMDSconf, data = Conf2, group = "Year")
  #test significance
  envfit(NMDSconf~Year,data = Conf2)
  
  #plot by year
  PlotNMDS(NMDSconf, data = Conf2, group = "targets2")
  #test significance
  envfit(NMDSconf~targets2,data = Conf2)
  
  #More differences by year and gear type than site
  
  #stacked bar plot
  ggplot(Conf, aes(x=site, y = CPUE, fill = Analy2)) + 
    geom_bar(stat = "identity", position = "fill")+
    facet_grid(targets2~Year, scales = "free", space = "free_x") +
    scale_fill_manual(values = mypal, name = NULL) + 
    xlab("Site")+ ylab("Relative percent composition") + mytheme
  
  
  
  }