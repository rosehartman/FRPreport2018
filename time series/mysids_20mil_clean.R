#Clean version of our mysid data that was taken concurrently with 20 mm in 2017 and 2018

#Load the required librarys


library(tidyverse)
library(lubridate)
library(reshape2)
library(vegan)
library(pwr)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(visreg)
library(readxl)
library(simr)
#set up a color pallete
mypal = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"), "black", "white", "grey", "pink")

#Upload the data
#First load the function I wrote to query the FRP database
source("querydatabase.R")

#Now specify the path to the FRP database
path = "U:/FRPA/MONITORING/Labs/Databases/FRPdata28DEC2018.accdb"

#Query the invertebrate data
inverts2 <- GetFRPdata(path, type = "inverts")

#Query the station information
stations = GetFRPdata(path, type = "stations")

#attatch station information
inverts2 = merge(inverts2, stations[,c(2,7)], all.x = T)

#assume that if "subsampled" is left blankn, it's 100%
inverts2$subsampled[which(is.na(inverts2$subsampled))]= 100

#convert the date from factor to date
inverts2$Date <- ymd(inverts2$Date)

#clean up the 'sample types' column
inverts2$Sampletype = as.character(inverts2$Sampletype)
inverts2$Sampletype[which(inverts2$Sampletype == "Mysid net (no sled)" | 
                            inverts2$Sampletype == "larval sled oblique"|
                            inverts2$Sampletype == "larval sled benthic")]=
  "Mysid net"

#

#subset teh mysid samples
mys <- filter(inverts2, substr(SampleID, 1, 3)=="MAC")

#import info on just the stations we are interested in
X20mil <- read_excel("time series/20mil.xlsx", sheet = "mysids", 
                     col_types = c( "text", "date", 
                                   "text",  "numeric"))



# Adjust total count for subsampling
mys$atotal = mys$TotalCount*(1/(mys$subsampled/100))

#Calculate CPUE
mys$CPUE = mys$atotal/mys$Volume

#Subset so it's just the time series data we took alongside 20mm, not the "blitz" samples

myss = merge(mys[,c(2, 4:26)], X20mil)


#add the zeros back in so we can correctly calculate averages
mys2c = dcast(myss, SampleID + Ggdist +  Station + Date ~CommonName, 
              value.var = "CPUE", fun.aggregate = sum)
mys2m = melt(mys2c, id.vars = c("SampleID",  "Station", "Date", "Ggdist"), 
             variable.name = "CommonName", value.name = "CPUE")


#import some other categories for the taxonomic groups
frptaxa <- read_excel("blitz2/zoocodes.xlsx")

mys2 = merge(mys2m, frptaxa[,c(3,10,13,14)])
mys2$month = month(mys2$Date)

#calculate averages
mysave = group_by(mys2, Station, month, CommonName, Analy2) %>% 
  summarize(CPUE = mean(CPUE))
mysave2 = group_by(mysave, month,  CommonName, Analy2) %>% 
  summarize(CPUE = mean(CPUE))

#################################################################################
#remove the mesozooplankton
mys3 = filter(mys2, Analy2 != "Copepoda" & Analy2 != "Cladocera" & Analy2 != "Ostracoda")
mys3$year = as.factor(year(mys3$Date))

#calculate averages
mysaveN = group_by(mys3, SampleID, Station, year, month, Analy2) %>% 
  summarize(CPUE = sum(CPUE))
mysave2N = group_by(mysaveN, year, month, Station, Analy2) %>% 
  summarize(CPUE = mean(CPUE))

#Fix the labels
myslab = c("Amphipoda", "Annelida", "Arachnida", "Cnidaria", "Collembola", 
           "Cumacea", "Decapoda", 
           "fish", "Insecta", "Isopoda", "Mollusca", "Mysidea", "Nematoda",
           "Ostracoda", "other", 
           "Platyhelminthes")

#Stacked bar plot of CPUE by location for each trawl.
N2 = ggplot(mysave2N, aes(x = month, y = CPUE))
N2 + geom_bar(stat = "identity", aes(fill = Analy2), position = "fill") + facet_grid(year~Station, scales = "free_y")+
  scale_fill_manual(values = mypal, labels = myslab, name = NULL) +
  xlab("Sampling Month") + ylab("CPUE (count per cubic meter)")+
  theme(legend.position = "bottom") 

##############################################################################
#For the June report, I will make some summeries of our data, but don't compare to any of IEP's data because
#we don't have the 20mm data yet.

#Calculate total CPUE for the entire sample (minus mesozooplankton)
myLitsum = group_by(mys3, SampleID, Station, year, month, Ggdist) %>% 
  summarize(tCPUE = sum(CPUE, na.rm = T), logCPUE =log(tCPUE+1) )

ggplot(myLitsum, aes(x=month, y = logCPUE, color = year)) + geom_point() +
  geom_line() +
  facet_wrap(.~Station)

#GLM of log CPUE as predicted by month of the year and distance from the Golden Gate
g1 = lm(log(CPUE+1)~month+Ggdist + year, data = myLitsum)
summary(g1)

#Check out the diagnostic plots
plot(g1)
#looks pretty good!

#check out the partial residual plots
visreg(g1)
#So there is a significant effect of distance from the Golden Gate, but not month of the year

##############################################################################
#multivariate stats

#set up a community data matrix
Commat = dcast(mys3, SampleID~Analy2, value.var = "CPUE", fun.aggregate = sum)
Commat = Commat[,2:16]

#get rid of empty rows
Commat2 = Commat[which(rowSums(Commat) !=0),]
myLitsum2 = myLitsum[which(rowSums(Commat) !=0),]

#Do a relative abundance matrix
Commatp = Commat2/rowSums(Commat2)

#PERMANOVA for community abundance:
adonis(Commat2~month + Ggdist + year, data = myLitsum2)

#Now for relative abundance:
am = adonis(Commatp~month + Ggdist + year, data = myLitsum2)
am
#Again, there is a significant effect of position on the ecocline (distance from Golden Gate)
#but not month of the year. R-squared is pretty low too. The relative abunance fits better than
#the absolute abundance.

#Calculate NMDS for both absoluate and proportinoal datasets. 
n1 = metaMDS(Commat2, trymax = 100)
n2 = metaMDS(Commatp, trymax = 100)
#Both NMDS work about the same, neither is great,
#I'm using the GLM of total catch anyway, so we will go with relative abundance from now on.

#plot the NMDS
source("plotNMDS.R")

PlotNMDS(n2, data = myLitsum2, group = "year", textp = "T") 
#Ignore the error messages

#plot GAMs of continuous variables on top of the plot
g = ordisurf(n2~Ggdist, data = myLitsum2, add = T, labcex = 1, col = "green")
g2 = ordisurf(n2~month, col = "yellow", data = myLitsum2, add = T, 
         levels = c(1,2,3,4,5, 6), labcex = 1)


legend(x="bottomright", 
       c("Distance from GG", "month"),
       col = c("green", "yellow"), lty = 1)

