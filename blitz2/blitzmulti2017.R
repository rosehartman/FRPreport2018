#Multivariate analyses

source("blitz2/blitzclean.R")
source("plotNMDS.R")

#First we'll make a nice stacked bar graph of community composition by region and habitat type
#to see if anything jumps out at us.

ggplot(filter(bugsblitz, targets2!= "benthic"), aes(x=site, y = CPUE, fill = Analy2)) + geom_bar(stat = "identity", position = "fill")+
  facet_grid(targets2 + year~Region2, scales = "free", space = "free_x") +
  scale_fill_manual(values = mypal, name = NULL) + 
  xlab("Site")+ ylab("Relative percent composition") + mytheme

#clams = filter(bugs2017, Analy2 == "Mollusca" & targets2 == "benthic" & Region == "Suisun")

#do one by site type
ggplot(filter(bugsblitz, targets2!= "benthic"), aes(x=site, y = CPUE, fill = Analy2)) + 
  geom_bar(stat = "identity", position = "fill")+
  facet_grid(targets2~sitetype, scales = "free", space = "free_x") +
  scale_fill_manual(values = mypal, name = NULL) + 
  xlab("Site")+ ylab("Relative percent composition") + mytheme

#ow group everything within a site type together
ggplot(filter(bugsblitz, targets2!= "benthic"), aes(x=sitetype, y = CPUE, fill = Analy2)) + 
  geom_bar(stat = "identity", position = "fill")+
  facet_grid(targets2~., scales = "free", space = "free_x") +
  scale_fill_manual(values = mypal, name = NULL) + 
  xlab("Site Type")+ ylab("Relative percent composition") + mytheme

##############################################################################################
#community matricies for the different methods


#community data matrix

#Get rid of the first row (with sample IDs)
row.names(ComMat.1) = ComMat.1$SampleID
ComMat2.1 = ComMat.1[,2:ncol(ComMat.1)]
ComMat2.1 = ComMat2.1[which(bugstot$targets2 != "benthic"),]
#EAV2-14FEB2017 is giving us trouble. I'm going to take it out
ComMat2.1 = ComMat2.1[which(bugstot$SampleID!= "EAV2-14FEB2017"),]

#get rit of empty rows
ComMat2.2 = ComMat2.1[which(rowSums(ComMat2.1)!=0),]

#get rid of PVCs and ponars
bugstot2 = filter(bugstot, targets2 != "benthic" & SampleID != "EAV2-14FEB2017")[which(rowSums(ComMat2.1)!=0),]



#make a matrix based on proportion of total catch for each sample
ComMatp.1 = ComMat2.2/rowSums(ComMat2.2)

####################################################################################
#now do it for the smaller groups
#Get rid of the first row (with sample IDs)
row.names(ComMat.2) = ComMat.2$SampleID
ComMat.22 = ComMat.2[,2:ncol(ComMat.2)]
ComMat.22 = ComMat.22[which(bugstot$targets2 != "benthic"),]
#EAV2-14FEB2017 is giving us trouble. I'm going to take it out
ComMat.22 = ComMat.22[which(bugstot$SampleID!= "EAV2-14FEB2017"),]

#get rit of empty rows
ComMat.22 = ComMat.22[which(rowSums(ComMat.22)!=0),]


#make a matrix based on proportion of total catch for each sample
ComMat.22p = ComMat.22/rowSums(ComMat.22)

#now do it for the common names
#Get rid of the first row (with sample IDs)
row.names(ComMat.CN) = ComMat.CN$SampleID
Com.CN2 = ComMat.CN[,2:ncol(ComMat.CN)]
Com.CN2 = Com.CN2[which(bugstot$targets2 != "benthic"),]
#EAV2-14FEB2017 is giving us trouble. I'm going to take it out
Com.CN2 = Com.CN2[which(bugstot$SampleID!= "EAV2-14FEB2017"),]

#get rit of empty rows
Com.CN2 = Com.CN2[which(rowSums(Com.CN2)!=0),]


#make a matrix based on proportion of total catch for each sample
Com.CN2p = Com.CN2/rowSums(Com.CN2)


#now do it for the cleaned up common names
#Get rid of the first row (with sample IDs)
row.names(ComMat.CN2) = ComMat.CN2$SampleID
Com.CN3 = ComMat.CN2[,2:ncol(ComMat.CN2)]
Com.CN3 = Com.CN3[which(bugstot$targets2 != "benthic"),]
#EAV2-14FEB2017 is giving us trouble. I'm going to take it out
Com.CN3 = Com.CN3[which(bugstot$SampleID!= "EAV2-14FEB2017"),]

#get rit of empty rows
Com.CN3 = Com.CN3[which(rowSums(Com.CN3)!=0),]


#make a matrix based on proportion of total catch for each sample
Com.CN3p = Com.CN3/rowSums(Com.CN3)

##############################################################################
#look at multivariates
#First look at the "Analy2" version, with proportional catch
a1 = adonis(ComMatp.1 ~ Region2 + sitetype + site/targets2, strata = bugstot2$Region2, data = bugstot2)
a1
#The R2 is pretty small, but highly significant.

#now look at it with absolute catch
a2 = adonis(ComMat2.2 ~ sitetype +Region2 +  targets2 + site, data = bugstot2)
a2

#let's try an NMDS

mds3 = metaMDS(ComMat2.2, try = 50, trymax = 1000)
mds4 = metaMDS(ComMatp.1, try = 50, trymax = 1000)
#Nothing converged

#what if we used a different community matrix?
mds5 = metaMDS(ComMat.22, try = 50, trymax = 500)
mds6 = metaMDS(ComMat.22p, try = 50, trymax = 500)

mds7 = metaMDS(Com.CN2, try = 50, trymax = 500)
mds8 = metaMDS(Com.CN2p, try = 50, trymax = 500)

#Well, that didn't work. 

#Try summarizing by site and year


#Let's look at that "indicspecies" analysis I've done in the past, it would be cool if we can
#find something associated with managed wetlands

###########################################################################
#Package "indicspecies" can be used to find "indicator species" that are diagnostic
#of particular sites (or, in our case, habitat types). It doesn't show you the 
#statistical significance of the difference between every single speccies, but
#it pulls out the most important ones with the biggest differences, which I think
#is more useful.
library(indicspecies)
# This tests to see whether some groups are indicators of particular sites
indi1 = multipatt(ComMat2.2, bugstot2$sitetype, func = "r.g", control = how(nperm=999))
summary(indi1)

#Meh. Let's look at the finer-scale groupings
indi2 = multipatt(ComMat.22, bugstot2$sitetype, func = "r.g", control = how(nperm=999), max.order = 1)
summary(indi2)
#That's more interesting! 

#now try common names
indi3 = multipatt(Com.CN2, bugstot2$sitetype, func = "r.g", control = how(nperm=999), max.order = 1)
summary(indi3)

#now try cleaned up common names
indi4 = multipatt(Com.CN3, bugstot2$sitetype, func = "r.g", control = how(nperm=999), max.order = 1)
summary(indi4)
######################################################################3

#one habitat at a time

#######################################################################################

###############################################################################################
#mysids

#we only had one trawl in a diked wetland, so I am taking tha out
mys.2x = filter(mys.2x, site != "Bradmoor")

ComMatMx = dcast(mys.2x, formula = SampleID~Analy2, value.var="CPUE", 
                 fun.aggregate = sum, fill = 0)
row.names(ComMatMx) = ComMatMx$SampleID
ComMat2Mx = ComMatMx[,2:ncol(ComMatMx)]
mys.1x = mys.1[which(rowSums(ComMat2Mx )!=0),]
ComMat2Mx = ComMat2nx[which(rowSums(ComMat2Mx )!=0),]
#make one based on proportion of total catch for each sample
ComMatpMx = ComMat2nx/rowSums(ComMat2Mx)


#PERMANOVA with relative abundance
mys.1x$Year = year(mys.1x$Date)
pn1 = adonis(ComMatpMx~ sitetype+Year + Region2/site, strata = mys.1x$Region2, data = mys.1x)
pn1
#NMDS

NMDS2n = metaMDS(ComMatpnx, try = 50, trymax = 1000)
NMDS2n

#plot it
mys.1x$Region2 = as.factor(mys.1x$Region2)
PlotNMDS(NMDS2n, data = mys.1x, group = "Region2")
#test significance
envfit(NMDS2n~Region2,data = mys.1x)

#plot by site type
PlotNMDS(NMDS2n, data = mys.1x, group = "sitetype")
#test significance
envfit(NMDS2n~sitetype,data = mys.1x)

#Looks like wetland sites tend to be more highly variable than channel sites, but don't have
#all that different community compositions. Definitley some differences by region.


#########################################################################################
#neuston
#community data matrix
#I'm goig to try filtering out anything with less than three individuals
neu.2x = filter(neu.2x, tcount > 3)

ComMatn = dcast(neu.2x, formula = SampleID~Analy2, value.var="CPUE", 
                fun.aggregate = sum, fill = 0)
row.names(ComMatn) = ComMatn$SampleID
ComMat2n = ComMatn[,2:ncol(ComMatn)]
#make one based on proportion of total catch for each sample
ComMatpn = ComMat2n/rowSums(ComMat2n)

#PERMANOVA with relative abundance
neu.1$Year = year(neu.1$Date)
neu.1 = filter(neu.1, tcount >3)
pn1 = adonis(ComMatpn~sitetype +Year+ Region2/site, strata = neu.1$Region2, data = neu.1)
pn1
#NMDS

NMDS2nu = metaMDS(ComMatpn, trymax = 4999)
NMDS2nu

#Plot the NMDS by region
neu.1 = filter(neu.1, tcount >3)
neu.1$Region2 = as.factor(neu.1$Region2)
PlotNMDS(NMDS2nu, data = neu.1, group = "Region2")
#test the significance of the different groups
envfit(NMDS2nu~Region,data = neu.1)

#plot it based on site type
PlotNMDS(NMDS2nu, data = neu.1, group = "sitetype")
envfit(NMDS2nu~sitetype,data = neu.1)


#########################################################################################
#sweepnets
#community data matrix
sweeps.2xx = filter(sweeps.2x, SampleID != "EAV2-14FEB2017" & tcount >3)
ComMatsw = dcast(sweeps.2xx, formula = SampleID~Analy2, value.var="CPUE", 
                fun.aggregate = sum, fill = 0)
row.names(ComMatsw) = ComMatsw$SampleID
ComMat2sw = ComMatsw[,2:ncol(ComMatsw)]

#EAV2-14FEB2017 is giving us trouble. I'm going to take it ou

sw.1x = filter(sweeps.1[which(sweeps.1$tCPUE!=0),], SampleID != "EAV2-14FEB2017" & tcount>3)
ComMat2sw = ComMat2sw[which(rowSums(ComMat2sw, na.rm=T)!=0),]

#make one based on proportion of total catch for each sample
ComMatpsw = ComMat2sw/rowSums(ComMat2sw)

#PERMANOVA with relative abundance
sw.1x$Year = year(sw.1x$Date)
psw1 = adonis(ComMatpsw~sitetype+ Year+ Target+ Region2/site, data = sw.1x)
psw1
#NMDS

NMDS2sw = metaMDS(ComMatpsw, trymax = 2000)
NMDS2sw = metaMDS(ComMatpsw, trymax = 2000,  sratmax = 0.99999999, sfgrmin = 1e-8)
NMDS2sw = metaMDS(ComMatpsw, previous.best = NMDS2sw, trymax = 2000,  sratmax = 0.99999999, sfgrmin = 1e-8)

#Plot the NMDS by region
PlotNMDS(NMDS2sw, data = sw.1x, group = "Region2")
#test the significance of the different groups
envfit(NMDS2sw~Region2,data = sw.1x)



#Maybe if I just use 2018 I can get it to converge?
sw.1x$year = year(sw.1x$Date)
sw.1x$Region2 = as.factor(sw.1x$Region2)
SN2018 = filter(sw.1x, year == 2018)
Matsn2018 = ComMatpsw[which(sw.1x$year == 2018),]
sw2018 = metaMDS(Matsn2018, trymax = 1999)
sw2018


#Plot the NMDS by region
PlotNMDS(sw2018, data = SN2018, group = "Region2")

#plot it based on site type
PlotNMDS(sw2018, data = SN2018, group = "sitetype")


#plot it based on veg type
PlotNMDS(sw2018, data = SN2018, group = "Target")
