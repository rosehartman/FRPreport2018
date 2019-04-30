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
  Reg.1 = summarize(group_by(Reg, SampleID, Station, Region2, targets2, 
                                    Sampletype, site, sitetype, Date),
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
                                  Sampletype, Analy2), mCPUE = mean(CPUE),
                        sdCPUE = sd(CPUE), seCPUE = sd(CPUE)/length(CPUE))
  
  #Now calculate total CPUE
  Regtot = summarize(group_by(Reg.1x, SampleID, Date, site, sitetype, targets2, 
                               Region2, Sampletype), tCPUE = sum(CPUE, na.rm = T))
  
  
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
Regionalbugs("San-San Joaquin", bugsblitz)
Regionalbugs("Cache Slough", bugsblitz)
Regionalbugs("Suisun Bay", bugsblitz)
Regionalbugs("Confluence", bugsblitz)
Regionalbugs("Nurse-Denverton", bugsblitz)


####################################################################################################
  #now, are there differences in community composition
RegMultibugs = function(region, data) {
  require(tidyverse)
  require(reshape2)

  Reg = filter(data, Region2 == region & targets2 == "benthic")
  #turn tyear into a factor
  Reg$Year = as.factor(year(Reg$Date))
  
  #make a summary data set for the environmental matrix
  Reg.1 = summarize(group_by(Reg, SampleID, Station, Region2, Target, targets2, 
                             Region, Sampletype, site, sitetype, Date),
                    tcount = sum(atotal, na.rm = T), tCPUE = sum(CPUE, na.rm = T), 
                    richness = length(unique(Analy)), logtot = log(tCPUE + 1))
  
 
  #convert from strings to factors
  Reg.1$site = as.factor(Reg.1$site)
  Reg.1$sitetype = as.factor(Reg.1$sitetype)
  Reg.1$targets2 = as.factor(Reg.1$targets2)
  Reg.1$Year = as.factor(year(Reg.1$Date))
  
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
RegMultibugs("San-San Joaquin", bugsblitz)
RegMultibugs("Cache Slough", bugsblitz)
RegMultibugs("Suisun Bay", bugsblitz)
RegMultibugs("Nurse-Denverton", bugsblitz)

#########################################################################################################
#now let's do the clams
RegClam = function(region, data){
    require(tidyverse)
    require(reshape2)
    require(lubridate)
    require(visreg)
    #filter the whole blitz data set for the region of interest
    Reg = filter(data, Region2 == region & targets2 == "benthic")
    #turn tyear into a factor
    Reg$Year = as.factor(year(Reg$Date))
    
    #put in all the zeros for clams and remove the non-clams
   
    benMat.2 = dcast(Reg, formula = SampleID~CN, value.var="CPUE", 
                     fun.aggregate = sum, fill = 0)
    ben.2x = melt(benMat.2, id.vars = "SampleID", variable.name = "CN", value.name = "CPUE")
     ben.2x = filter(ben.2x, CN == "Corbicula" | CN == "Potamocorbula" | CN == "Clam Other" | CN == "No Catch")
    ben.2x = merge(Reg[,c(3,8,33,34,35,37)], ben.2x, by = "SampleID")
    ben.2x = ben.2x[(!duplicated(ben.2x)),]
    
    #Summary data set
    ben.1 = summarize(group_by(ben.2x, SampleID, 
                               Region2, site, sitetype, Date, Year),
                      tCPUE = sum(CPUE, na.rm = T))
    
    ben.1$site = as.factor(ben.1$site)
    ben.1$sitetype = as.factor(ben.1$sitetype)
    
    bensum = summarize(group_by(ben.2x, site, Region2,
                                sitetype, SampleID, Date, Year), 
                       tCPUE = sum(CPUE, na.rm = T), logCPUE = log(tCPUE + 1))
    
    
    #convert from strings to factors
    bensum$site = as.factor(bensum$site)
    bensum$sitetype = as.factor(bensum$sitetype)

    #Mean CPUE, sample size, and standard error for each location and habitat so I can make a pretty graph.
    Regtotave = summarize(group_by(bensum, site, sitetype,  
                                   Region2, Year), mCPUE = mean(tCPUE, na.rm = T), logmCPUE = mean(logCPUE, na.rm = T),
                          sdCPUE = sd(tCPUE, na.rm = T), seCPUE = sd(tCPUE)/length(tCPUE), N = length(SampleID))

    
    p1 = ggplot(Regtotave, aes(x=site, y = mCPUE, fill = site))
    print(p1+geom_bar(stat = "identity") + 
            geom_errorbar(aes(ymin = mCPUE - seCPUE, ymax = mCPUE + seCPUE))+
            facet_grid(.~Year, scales = "free_y"))
    
    #First, are there differences in catch?
    if (length(levels(bensum$Year)) >1 ) {
    
    blitzReg = aov(logCPUE ~ site+Year, 
                   data = bensum)
    } else  blitzReg = aov(logCPUE ~ site, 
                     data = bensum) 
    
    print(summary(blitzReg))
    print(TukeyHSD(blitzReg))
    visreg(blitzReg)
}

RegClam("Cache Slough", bugsblitz)
RegClam("Confluence", bugsblitz)
RegClam("San-San Joaquin", bugsblitz)
RegClam("Suisun Bay", bugsblitz)
RegClam("Nurse-Denverton", bugsblitz)

##########################################################################################################
#now let's do the zooplankton, first univariate total CPUE
RegZoop = function(region, data) {
  require(tidyverse)
  require(reshape2)
  require(lubridate)
  require(visreg)
  Reg = filter(data, Region2 == region)
  
  #total CPUE
  zootot = group_by(Reg, SampleID, Date, site, year, month, 
                    Region2, sitetype) %>% summarize(tCPUE = sum(CPUE, na.rm = T), logCPUE = log(tCPUE + 1))
  
  zoototave = group_by(zootot, site, Region2, sitetype, year) %>% 
    summarize(mCPUE = mean(tCPUE, na.rm = T), seCPUE = sd(tCPUE)/length(tCPUE), 
              N = length(tCPUE),  mlogCPUE = mean(logCPUE) )
  
  z3 = ggplot(zoototave, aes(x=site, y=mCPUE))
  z3 + geom_bar(stat = "identity", aes(fill = sitetype)) +
    facet_grid(year~Region2, scales = "free", space = "free_x") + 
    geom_errorbar(aes(ymin = mCPUE - seCPUE, ymax = mCPUE + seCPUE)) +
    scale_fill_manual(values = mypal)+
    geom_label(aes(label = paste("n = ", N), y = mCPUE+ 300), label.padding = unit(0.1, "lines"), size = 4) +
    ylab("CPUE") + xlab(label = NULL) 
  
  #log-transformed
  z3.1 = ggplot(zoototave, aes(x=site, y=mlogCPUE))
  z3.1 + geom_bar(stat = "identity", aes(fill = sitetype)) +
    facet_grid(year~Region2, scales = "free", space = "free_x") + 
    scale_fill_manual(values = mypal)+
    geom_label(aes(label = paste("n = ", N), y = 0), label.padding = unit(0.1, "lines"), size = 4) +
    ylab("mean log CPUE") + xlab(label = NULL)
  
  #GLMm of total CPUE ############################ This is probably the best one to use.
  zootot$sitetype = as.factor(zootot$sitetype)
  zootot$year = as.factor(zootot$year)
  zblitz = aov(logCPUE ~ site + year, data = zootot)
  print(summary(zblitz))
  visreg(zblitz)
  TukeyHSD(zblitz)
  
}

RegZoop("Cache Slough", zooblitz)
RegZoop("Confluence", zooblitz)
RegZoop("Suisun Bay", zooblitz)
RegZoop("San-San Joaquin", zooblitz)
RegZoop("Nurse-Denverton", zooblitz)

##############################################################################################################
#Multivariate zooplankton
RegMultiZoop = function(region, data){
  require(tidyverse)
  require(reshape2)
  require(vegan)
  source("plotNMDS.R")
  
  #subset the region we are interested in
  zooR = filter(data, Region2 == region)
  zoowide = dcast(zooR, SampleID+Region2+site+Date+month+year+sitetype~AnalyLS, 
                  fun.aggregate = sum, na.rm = T, value.var = "CPUE")
  
  #plot of community composition 
  z4 = ggplot(zooR, aes(x=site, y=CPUE))
  z4 = z4 + geom_bar(stat = "identity", aes(fill = AnalyLS), position = "fill") +
    facet_grid(year~., scales = "free", space = "free_x") + 
    scale_fill_manual(values = c(mypal, "grey"), name = NULL) +
    labs(y="Relative Abudance", x=NULL)
  print(z4)
  
 
  #Community matrix by CPUE
  zMat = zoowide[,8:ncol(zoowide)]
  #now by proprotions
  zMatp = zMat/rowSums(zMat)
  
  #PerMANOVA
  za1 = adonis(zMatp~ site + year, data = zoowide)
  print(za1)
  #region and site type are both significant, and year
  
  #NMDS
  zNMDS = metaMDS(zMatp, try = 200)
  print(zNMDS)
  
  #Make sure "site" is a factor
  zoowide$site = as.factor(zoowide$site)
  #plot the NMDS by site
  PlotNMDS(zNMDS, data = zoowide, group = "site")
  
  #make sure year is a factor
  zoowide$year = as.factor(zoowide$year)
  #Plot NMDS by year
  PlotNMDS(zNMDS, data = zoowide, group = "year")

}

RegMultiZoop("Cache Slough", zooblitz)
RegMultiZoop("San-San Joaquin", zooblitz)
