
#Look at the phytoplankton data for the blitz, comparing 2017 and 2018
library(readxl)
library(MuMIn)
library(lmerTest)
library(tidyverse)
library(lubridate)
library(reshape2)
library(vegan)
library(RColorBrewer)
library(indicspecies)

#load the function to get into the FRP database
#source("querydatabase.R")

#Now specify the path to the FRP database
#path = "U:/FRPA/MONITORING/Labs/Databases/FRPdata28DEC2018.accdb"

#Query the invertebrate data
#phytos <- GetFRPdata(path, type = "phytoplankton")
phytos <- read_excel("phytoplankton/phytoquery.xlsx")

stations = read_xlsx("blitz2/Stations2.xlsx")

phytos1 = merge(phytos, stations)
phytos1 = filter(phytos1, site != "Dow" & site != "Lindsey")

#see what taxa we got
FRPtax = group_by(phytos1, Taxon) %>% summarize(count = sum(CellspermL))

#lump the less common taxa together
FRPtax$taxon2 = as.character(FRPtax$Taxon)
FRPtax$prop = FRPtax$count/sum(FRPtax$count)
FRPtax$taxon2[which(FRPtax$prop < 0.0001)] = "other"
phytos1 = merge(phytos1, FRPtax[,c(1,3)])

############################################################

p1 = ggplot(phytos1, aes(x=SampleName, y = CellspermL, fill = taxon2))
p1+ geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))

#import tiffanys info
taxa2 <- read_excel("phytoplankton/FRPtaxa T Brown edits.xlsx", 
                    sheet = "Taxon Info")

#merge the new info onto the data set
phytos2 = merge(phytos1, taxa2, by = "Taxon")

phytos2$year = year(phytos2$Date)

#combine into one row per sample
physum = group_by(phytos2, SampleName, site, sitetype, Region2, Chlorophyll, year) %>% 
  summarize(CPUE = sum(CellspermL), BPUE = sum(BiovolumeperuL))



######################################################################################
#some exploritory plots
#set up a color pallete
mypal = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"), "black", "white", brewer.pal(8, "Set2"))

#set up a theme
mytheme = theme(strip.text.x = element_text(size=14),strip.text.y = element_text(size=18),
                axis.text.x  = element_text(angle=90, hjust=1, vjust = 0.5, size=12),
                axis.title.y  = element_text(size=16),
                axis.title.x  = element_text(size = 16), axis.text.y = element_text(size = 14),
                legend.text = element_text(size = 14),
                legend.title=element_blank(), legend.background =element_blank()) 

#graph by Region


#courser groups
p2.2 = ggplot(phytos2, aes(x=site, y = CellspermL, fill = Type2))
p2.2+ geom_bar(stat = "identity", position = "fill") + 
  facet_grid(year~Region2, scales = "free_x", space = "free_x")+
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = mypal, name = NULL) + ylab("Relative abundance")
p2.2 + geom_bar(stat = "identity", position = "fill") + facet_wrap(~sitetype, scales = "free")

####################################################################################################
#NMDS plot
Phy2 = group_by(phytos2, SampleName, `Habitat type`, Region2, Date, year, site, sitetype, taxon2) %>% summarize(cells = sum(CellspermL))
Mat = spread(Phy2, key = taxon2, value = cells, fill =0)
Mat$Region2 = as.factor(Mat$Region2)
Mat$sitetype = as.factor(Mat$sitetype)
Mat2 = Mat[, c(8:93)]
Mat2p = Mat2/rowSums(Mat2)

m1 = metaMDS(Mat2p, try = 99, trymax = 100)
m1
source("plotNMDS.R")
PlotNMDS2(m1, group = "Region2", data = Mat, xlimits = c(-1.5, 1.2), ylimits = c(-1.5, 1.5), textp = F)
PlotNMDS2(m1, group = "sitetype", data = Mat, xlimits = c(-1.5, 1.2), ylimits = c(-1.5, 1.5), textp = F)

a0 = adonis(Mat2p~ sitetype + year +Region2/site, strata = Mat$Region2, data = Mat)
a0


#Do it again with larger categories
Phy3 = group_by(phytos2, SampleName, `Habitat type`, Region2, Date, year, site, sitetype, Type2) %>% summarize(cells = sum(CellspermL))
Mat3 = spread(Phy3, key = Type2, value = cells, fill =0)
Mat3$Region2 = as.factor(Mat3$Region2)
Mat3$sitetype = as.factor(Mat3$sitetype)
Mat32 = Mat3[, c(8:21)]
Mat3p = Mat32/rowSums(Mat32)

m2 = metaMDS(Mat3p, trymax = 500)

PlotNMDS(m2, group = "Region2", data = Mat3,textp = T)
PlotNMDS(m2, group = "sitetype", data = Mat3,textp = T)

#Permanova
a1 = adonis(Mat3p~ Region2+ sitetype+ site, strata = Mat3$Region2, data = Mat3)
a1

###################################################################################################3
#look at the habitat-specific data from Liberty

PhyLib = filter(phytos2, site == "Liberty" & year == 2018)

pL = ggplot(PhyLib, aes(x=`Habitat type`, y = CellspermL, fill = Type2))
pL+ geom_bar(stat = "identity", position = "fill") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = mypal)

pL+ geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = mypal)

#Look at the individual species
pL2 = ggplot(PhyLib, aes(x=`Habitat type`, y = CellspermL, fill = taxon2))
pL2+ geom_bar(stat = "identity", position = "fill") + 
  theme(axis.text.x = element_text(angle = 90))

#We still have too many taxa. I'll lump some more into the "other" category

#lump the less common taxa together
Libtax = group_by(PhyLib, taxon2) %>% summarise(count = sum(CellspermL), prop = count/sum(PhyLib$CellspermL))
PhyLib = merge(PhyLib, Libtax)
PhyLib$taxon2[which(PhyLib$prop < 0.0001)] = "other"

pL3 = ggplot(PhyLib, aes(x=`Habitat type`, y = CellspermL, fill = taxon2))
pL3+ geom_bar(stat = "identity", position = "fill") + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_fill_manual(values = mypal)



####################################################################################################
#NMDS plot
PhyLib2 = group_by(PhyLib, SampleName, `Habitat type`, taxon2) %>% summarize(cells = sum(CellspermL))
LMat = spread(PhyLib2, key = taxon2, value = cells, fill =0)
LMat = mutate(LMat, Habitat = as.factor(`Habitat type`))
LMat$Habitat= as.factor(LMat$`Habitat type`)
LMat2 = LMat[, c(3:33)]
LMat2p = LMat2/rowSums(LMat2)

m1 = metaMDS(LMat2p)
source("plotNMDS.R")
PlotNMDS2(m1, group = "Habitat", data = LMat, xlimits = c(-1.5, 1.2), ylimits = c(-.7, 1.5))


#Do it again with larger categories
PhyLib3 = group_by(PhyLib, SampleName, `Habitat type`, Type2) %>% summarize(cells = sum(CellspermL))
LMat1 = spread(PhyLib3, key = Type2, value = cells, fill =0)
LMat1$Habitat= as.factor(LMat1$`Habitat type`)
LMat12 = LMat1[, c(3:11)]
LMat12p = LMat12/rowSums(LMat12)

m2 = metaMDS(LMat12p)

PlotNMDS2(m2, group = "Habitat", data = LMat, xlimits = c(-1, 1), ylimits = c(-1, 1))

#If I multiply the species scores from m2 by the speceis from a given sample, it might give me th
#coordinates for the NMDS plot

mscores = as.data.frame(m2$species)
mscores$Type2 = rownames(mscores)

#not all types were represented, but this is the best we got
Mat3x = cbind(Mat3[,1], Mat3p)
phytest = gather(Mat3x, key = "Type2", value = "prop", -SampleName)
phytest2 = merge(phytest, mscores)
phytest2$x = phytest2$MDS1*phytest2$prop
phytest2$y = phytest2$MDS2*phytest2$prop

#calculate coordinates
physcores = group_by(phytest2, SampleName) %>% summarize(xcoor = sum(x, na.rm=T), ycoor = sum(y, na.rm = T))

points(x=physcores$xcoor, y = physcores$ycoor)

#I'm not sure that was quite right, let's test it:
Mat3xl = cbind(LMat1[,1], LMat12p)
phytestl = gather(Mat3xl, key = "Type2", value = "prop", -SampleName)

phytest2l = merge(phytestl, mscores)
phytest2l$x = phytest2l$MDS1*phytest2l$prop
phytest2l$y = phytest2l$MDS2*phytest2l$prop

#calculate coordinates
physcoresl = group_by(phytest2l, SampleName) %>% summarize(xcoor = sum(x, na.rm=T), ycoor = sum(y, na.rm = T))

plot(m2, type="n", shrink = F, cex=1, xlim=c(-1, 1))
points(m2)
points(x=physcoresl$xcoor, y = physcoresl$ycoor, col = "red")

 #######################################################################################################
#Indicator species
LMat3 = droplevels(filter(LMat, Habitat != "6 = other"))
LMat3p = droplevels(LMat2p[which(LMat$Habitat != "6 = other"),])
m = multipatt(LMat3p, cluster = LMat3$Habitat, max.order = 1)
m
summary(m)
signassoc(LMat3p, cluster = LMat3$Habitat)
