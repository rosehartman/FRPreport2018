
#Look at the phytoplankton data EcoAnalysists sent us, just for the 2017 blitz

library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)
library(vegan)
library(RColorBrewer)

#read in the info EcoAnalysts send us
phytos <- read_excel("phytosR.xlsx", 
                     sheet = "phytosR", col_types = c("text", 
                                                      "date", "numeric", "numeric", "text", 
                                                      "text", "text", "text", "text", "text", 
                                                      "text", "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric", 
                                                      "numeric", "numeric", "numeric"))


phytob = filter(phytos, Date > "2017-3-1" & Date < "2017-5-1")

p1 = ggplot(phytob, aes(x=SiteID, y = CellspermL, fill = Division))
p1+ geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))

#Now lets get the site information:
phtyoquery <- read_excel("phtyoquery.xlsx", 
                         col_types = c("text", "text", 
                                       "date", "text", "text"))

phytos2 = merge(phytob, phtyoquery, all.x = T, by = "SiteID")

#import tiffanys info
taxa2 <- read_excel("FRPtaxa T Brown edits.xlsx", 
                    sheet = "Taxon Info")

#merge the new info onto the data set
phytos2 = merge(phytos2, taxa2, by = "Genus")



#combine into one row per sample
physum = group_by(phytos2, SampleName, Region2, site, sitetype, Date, year) %>% 
  summarize(CPUE = sum(CellspermL), BPUE = sum(BiovolumeperuL))

######################################################################################
#some exploritory plots
#set up a color pallete
mypal = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"), "black")

#set up a theme
mytheme = theme(strip.text.x = element_text(size=14),strip.text.y = element_text(size=18),
                axis.text.x  = element_text(angle=90, hjust=1, vjust = 0.5, size=12),
                axis.title.y  = element_text(size=16),
                axis.title.x  = element_text(size = 16), axis.text.y = element_text(size = 14),
                legend.text = element_text(size = 14),
                legend.title=element_blank(), legend.background =element_blank()) 

#graph by location

p2 = ggplot(phytos2, aes(x=SiteID, y = CellspermL, fill = Division))
p2+ geom_bar(stat = "identity") + facet_wrap(~Location, scales = "free")+
  theme(axis.text.x = element_text(angle = 90))

p21 = ggplot(phytos2, aes(x=Location, y = CellspermL, fill = Division))
p21+ geom_bar(stat = "identity", position = "fill") + facet_grid(.~Region, scales = "free_x", space = "free")
 


#finer taxonomic resolution:
p2.1 = ggplot(phytos2, aes(x=SiteID, y = CellspermL, fill = Order))
p2.1+ geom_bar(stat = "identity", position = "fill") + facet_wrap(~Location, scales = "free_x")+
  theme(axis.text.x = element_text(angle = 90))
p2.1 + geom_bar(stat = "identity") + facet_wrap(~Location, scales = "free")+
  theme(axis.text.x = element_text(angle = 90))


#finer taxonomic resolution:
p2.2 = ggplot(phytos2, aes(x=Location, y = CellspermL, fill = Order))
p2.2+ geom_bar(stat = "identity", position = "fill") + facet_grid(~Region, scales = "free_x", space = "free_x")+
  theme(axis.text.x = element_text(angle = 90))
p2.2 + geom_bar(stat = "identity", position = "fill") + facet_wrap(~SiteType, scales = "free")+
  theme(axis.text.x = element_text(angle = 90))


#Summarize by location and date
phytos3x = group_by(phytos2, Location, SiteID, Date.x, Division, SiteType) %>% summarise(Cells = sum(CellspermL))
phytos3 = group_by(phytos3x, Location, Division, SiteType) %>% summarize(Cells = mean(Cells), seCells = sd(Cells)/length(Cells))

p3 = ggplot(phytos3, aes(x=Division, y = Cells, fill = Division))
p3+ geom_bar(stat = "identity") + facet_wrap(~Location, scales = "free_x")

p3.1 = ggplot(phytos3, aes(x=Location, y = Cells, fill = Division))
p3.1+ geom_bar(stat = "identity") + facet_wrap(~SiteType, scales = "free_x")

p3.1+ geom_bar(stat = "identity", position = "fill") + facet_grid(.~SiteType, scales = "free_x")

p3.2 = ggplot(phytos2, aes(x=Location, y = CellspermL, fill = AlgalType))
p3.2+ geom_bar(stat = "identity", position = "fill") + 
  facet_grid(.~Region, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(mypal, "grey"), name = NULL) +
  mytheme + ylab("Relative Abudance") +xlab(label = NULL)

p3.2+ geom_bar(stat = "identity", position = "fill") + facet_grid(.~SiteType, scales = "free_x")


p3.3 = ggplot(phytos2, aes(x=Location, y = CellspermL, fill = Type2))
p3.3+ geom_bar(stat = "identity", position = "fill") + 
  facet_grid(.~Region, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(mypal, "grey"), name = NULL) +
  mytheme + ylab("Relative Abudance") +xlab(label = NULL)


p3.4 = ggplot(phytos2, aes(x=SiteType, y = CellspermL, fill = Type2))
p3.4+ geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = c(mypal, "grey"), name = NULL) +
  mytheme + ylab("Relative Abudance") +xlab(label = NULL)

p3.5 = ggplot(phytos2, aes(x=Region, y = CellspermL, fill = Type2))
p3.5+ geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = c(mypal, "grey"), name = NULL) +
  mytheme + ylab("Relative Abudance") +xlab(label = NULL)

########################################################################################################################
#I think I want to concentrate on community composition with this data, and get abundance with the chlorophyl a data



#Community matrix
PhyMat = dcast(phytos2, formula = SiteID~Division, value.var = "CellspermL", fun.aggregate = sum, fill = 0)
PhyMat = PhyMat[,2:8]

#PerMANOVA
a1 = adonis(PhyMat~Region + SiteType, data = physum)
a1

#different groupings

taxa2 <- read_excel("FRPtaxa T Brown edits.xlsx", 
                    sheet = "Taxon Info")

#merge the new info onto the data set
phytos2 = merge(phytos2, taxa2, by = "Genus")


#Community matrix
PhyMat2 = dcast(phytos2, formula = SiteID~AlgalType, value.var = "CellspermL", fun.aggregate = sum, fill = 0)
PhyMat2 = PhyMat2[,2:16]
PhyMat2p = PhyMat2/rowSums(PhyMat2)
#PerMANOVA
a2 = adonis(PhyMat2p~Region + SiteType + Location, data = physum)
a2

pNMDS = metaMDS(PhyMat2p)
pNMDS

source("plotNMDS.R")
physum$Region = as.factor(physum$Region)
physum$SiteType = as.factor(physum$SiteType)
PlotNMDS(pNMDS, data = physum, group = "Region")
envfit(pNMDS~Region, data = physum)
PlotNMDS(pNMDS, data = physum, group = "SiteType")
envfit(pNMDS~SiteType, data = physum)
