#stacy wants some pairwise comparisons to include in each project report
source("blitz2/blitzclean.R")

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
                                  Region2, Year2), mCPUE = mean(tCPUE, na.rm = T), 
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
      facet_grid(targets2~Year2, scales = "free_y"))
  
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


####################################################################################################
  #now, are there differences in community composition
RegMultibugs = function(region, data) {
  require(tidyverse)
  require(reshape2)
  require(lubridate)
  require(vegan)
  source("plotNMDS.R")
  Reg = filter(data, Region2 == region & targets2 != "benthic")
  #turn tyear into a factor
  Reg$Year = as.factor(year(Reg$Date))
  
  #make a summary data set for the environmental matrix
  Reg.1 = summarize(group_by(Reg, SampleID, Station, Region2, Target, targets2, 
                             Region, Sampletype, site, sitetype, Date),
                    tcount = sum(atotal, na.rm = T), tCPUE = sum(CPUE, na.rm = T), 
                    richness = length(unique(Analy)))
  
  #convert from strings to factors
  Reg.1$site = as.factor(Reg.1$site)
  Reg.1$sitetype = as.factor(Reg.1$sitetype)
  Reg.1$targets2 = as.factor(Reg.1$targets2)
  
  #set up the community matrix
  RegMat.1 = dcast(Reg, formula = SampleID~Analy2, value.var="CPUE", 
                   fun.aggregate = sum, fill = 0)
  

    row.names(RegMat.1) = RegMat.1$SampleID
  RegMat.12 = RegMat.1[,2:ncol(RegMat.1)]
  Reg2 = Reg.1[which(rowSums(RegMat.12 )!=0),]
  RegMat.12 = RegMat.12[which(rowSums(RegMat.12)!=0),]
  #make one based on proportion of total catch for each sample
  RegMat.12p = RegMat.12/rowSums(RegMat.12)
  
  
  #PERMANOVA with relative abundance
  Reg2$Year = as.factor(year(Reg2$Date))
  pn1 = adonis(RegMat.12p~site + Year + targets2, data = Reg2)
  print(pn1)
  #NMDS
  
  NMDSReg = metaMDS(RegMat.12p, try = 50, trymax = 500)
  print("NMDS of relative abundance"  )
 print( NMDSReg)
  
  #plot it
  PlotNMDS(NMDSReg, data = Reg2, group = "site")
  #test significance
print(  envfit(NMDSReg~site, data = Reg2))
  
  #plot by year
  PlotNMDS(NMDSReg, data = Reg2, group = "Year")
  #test significance
 print( envfit(NMDSReg~Year,data = Reg2))
  
  #plot by year
  PlotNMDS(NMDSReg, data = Reg2, group = "targets2")
  #test significance
  print(envfit(NMDSReg~targets2,data = Reg2))
  
  #More differences by year and gear type than site
  
  #stacked bar plot
 print( ggplot(Reg, aes(x=site, y = CPUE, fill = Analy2)) + 
    geom_bar(stat = "identity", position = "fill")+
    facet_grid(targets2~Year, scales = "free", space = "free_x") +
    scale_fill_manual(values = mypal, name = NULL) + 
    xlab("Site")+ ylab("Relative percent composition") + mytheme)
  
  
  
}

RegMultibugs("Confluence", bugsblitz)
