#Look at the phytoplankton data EcoAnalysists sent us, just for 2018
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


#see what taxa we got
FRPtax = group_by(phytos, Taxon, Division, Order, Family, Genus, Species) %>% summarize(count = sum(CellspermL))

############################################################
#Subset just 2018
phytob = filter(phytos, Date > "2018-1-1" )

p1 = ggplot(phytob, aes(x=SiteID, y = CellspermL, fill = Division))
p1+ geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))

#Now lets get the site information:
phtyoquery <- read_excel("phtyoquery.xlsx", 
                         col_types = c("text", "text", 
                                       "date", "text", "text"))

phytos2 = merge(phytob[,-2], phtyoquery, all.x = T, by = "SiteID")

#import tiffanys info
taxa2 <- read_excel("FRPtaxa T Brown edits.xlsx", 
                    sheet = "Taxon Info")

#merge the new info onto the data set
phytos2 = merge(phytos2, taxa2, by = "Genus")



#combine into one row per sample
physum = group_by(phytos2, SiteID, Region, Location, SiteType) %>% summarize(CPUE = sum(CellspermL))

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



#courser groups
p2.2 = ggplot(phytos2, aes(x=Location, y = CellspermL, fill = Type2))
p2.2+ geom_bar(stat = "identity", position = "fill") + 
  facet_grid(~Region, scales = "free_x", space = "free_x")+
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = mypal)
p2.2 + geom_bar(stat = "identity", position = "fill") + facet_wrap(~SiteType, scales = "free")+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = mypal)



#Summarize by location and date
phytos3x = group_by(phytos2, Location, SiteID, Date, Division, SiteType) %>% summarise(Cells = sum(CellspermL))
phytos3 = group_by(phytos3x, Location, Division, SiteType) %>% summarize(Cells = mean(Cells), seCells = sd(Cells)/length(Cells))

