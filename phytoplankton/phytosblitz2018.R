#Look at the phytoplankton data for the blitz, comparing 2017 and 2018
library(readxl)
library(MuMIn)
library(lmerTest)
library(tidyverse)
library(reshape2)
library(vegan)
library(RColorBrewer)

#load the function to get into the FRP database
source("querydatabase.R")

#Now specify the path to the FRP database
path = "U:/FRPA/MONITORING/Labs/Databases/FRPdata28DEC2018.accdb"

#Query the invertebrate data
phytos <- GetFRPdata(path, type = "phytoplankton")

stations = read_xlsx("blitz2/Stations2.xlsx")

phytos1 = merge(phytos, stations)

#see what taxa we got
FRPtax = group_by(phytos1, Taxon) %>% summarize(count = sum(CellspermL))

############################################################

p1 = ggplot(phytos1, aes(x=SampleName, y = CellspermL, fill = Taxon))
p1+ geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))

#import tiffanys info
taxa2 <- read_excel("phytoplankton/FRPtaxa T Brown edits.xlsx", 
                    sheet = "Taxon Info")

#merge the new info onto the data set
phytos2 = merge(phytos1, taxa2, by = "Taxon")

phytos2$year = year(phytos2$Date)

#combine into one row per sample
physum = group_by(phytos2, SampleName, site, sitetype, Region2) %>% summarize(CPUE = sum(CellspermL))

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

#graph by Region


#courser groups
p2.2 = ggplot(phytos2, aes(x=site, y = CellspermL, fill = Type2))
p2.2+ geom_bar(stat = "identity", position = "fill") + 
  facet_grid(year~Region2, scales = "free_x", space = "free_x")+
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = mypal)
p2.2 + geom_bar(stat = "identity", position = "fill") + facet_wrap(~sitetype, scales = "free")+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = mypal)



#Summarize by Region and date
phytos3x = group_by(phytos2, Region, SampleName, Date, Division, SiteType) %>% summarise(Cells = sum(CellspermL))
phytos3 = group_by(phytos3x, Region, Division, SiteType) %>% summarize(Cells = mean(Cells), seCells = sd(Cells)/length(Cells))

