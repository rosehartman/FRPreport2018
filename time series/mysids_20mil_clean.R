#Clean version of our mysid data that was taken concurrently with 20 mm in 2017

#Load the required librarys
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(lubridate)
library(reshape2)
library(vegan)
library(readxl)
library(visreg)

#set up a color pallete
mypal = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"), "black", "white", "grey", "pink")


#
#Upload the data
inverts2 <-  read_excel("inverts2.xlsx")

#convert the date from factor to date
inverts2$date <- ymd(inverts2$date)


#subset teh mysid samples from 2017 (Phase III)
mys2017 <- filter(inverts2, year(date) == 2017 & substr(SampleID, 1, 3)=="MAC")

#import info on just the stations we are interested in
X20mil <- read_excel("20mil.xlsx", sheet = "mysids", 
                     col_types = c("text", "text", "date", 
                                   "text", "text", "numeric"))

#get things organized to calculate effort
mys2017x = summarize(group_by(mys2017, SampleID, Station, Location, Sampletype), 
                     total = sum(TotalCount))

mys2017.e = summarize(group_by(mys2017, SampleID, Station, Location, 
                               date, NetMeterSTart, NetMeterEnd, LatStart, 
                               LongStart, LatEnd, LongEnd),
                       mean.en= mean(NetMeterSTart))

mys2017.e$volume = rep(NA, nrow(mys2017.e))


# Effort for benthic and oblique trawls is based on flowmeter readings and net mouth area

mys2017.e$volume= 
  (mys2017.e$NetMeterEnd- mys2017.e$NetMeterSTart)*0.026873027*0.16 #I've used the factory calibration here

# Fix any samples where the flowmeter ticked over to zero
mys2017.e$volume[which(mys2017.e$volume<0)] = ((1000000 - mys2017.e$NetMeterSTart[which(mys2017.e$volume<0)])+
                                                   mys2017.e$NetMeterEnd[which(mys2017.e$volume<0)])*0.026873027*0.16 #I've used the factory calibration here

#replace any "NA"s with the average volume for a 5-minute tow
mys2017.e$volume[which(is.na(mys2017.e$volume))] = mean(mys2017.e$volume, na.rm=T)

#something is wrong with a couple samples, so we will replace the origional
#flowmeter numbers with the average.
mys2017.e$volume[which(mys2017.e$SampleID=="MAC1-12APR2017" | mys2017.e$SampleID=="MAC1-8MAY2017" )] = 
  mean(mys2017.e$volume[which(mys2017.e$SampleID!="MAC1-12APR2017" & 
                                mys2017.e$SampleID!="MAC1-8MAY2017" )], na.rm=T)

#add the volumes to the origional data set
mys2017 = merge(mys2017, mys2017.e)

# Adjust total count for subsampling
mys2017$atotal = mys2017$TotalCount*(1/(mys2017$subsampled/100))

#Calculate CPUE
mys2017$CPUE = mys2017$atotal/mys2017$volume

#Subset so it's just the time series data we took alongside 20mm, not the "blitz" samples

mys2017s = merge(mys2017[,c(1, 5:22)], X20mil)


#add the zeros back in so we can correctly calculate averages
mys2c = dcast(mys2017s, SampleID + Ggdist + Location + Station + Date ~CommonName, 
              value.var = "CPUE", fun.aggregate = sum)
mys2m = melt(mys2c, id.vars = c("SampleID", "Location" , "Station", "Date", "Ggdist"), 
             variable.name = "CommonName", value.name = "CPUE")


#import some other categories for the taxonomic groups
frptaxa <- read_excel("FRP_EMPcrosswalk.xlsx", 
                      sheet = "allfrptaxa")

mys2 = merge(mys2m, frptaxa[,c(1,3,4,5,6)])
mys2$month = month(mys2$Date)

#calculate averages
mysave = group_by(mys2, Station, Location, month, CommonName, Analy) %>% 
  summarize(CPUE = mean(CPUE))
mysave2 = group_by(mysave, month, Location, CommonName, Analy) %>% 
  summarize(CPUE = mean(CPUE))

#################################################################################
#remove the mesozooplankton
mys3 = filter(mys2, Analy != "calanoid" & Analy != "cladocera" & Analy != "cyclopoid")

#calculate averages
mysaveN = group_by(mys3, SampleID, Station, Location, month, Analy) %>% 
  summarize(CPUE = sum(CPUE))
mysave2N = group_by(mysaveN, month, Station, Analy) %>% 
  summarize(CPUE = mean(CPUE))

#Fix the labels
myslab = c("Amphipoda", "Annelida", "Collembola", "Corophiidae", "Cumacea", "Decapoda", 
           "fish", "Gammaridea", "Insecta", "Isopoda", "Mollusca", "Mysidea", "Ostracoda", "other", 
           "Platyhelminthes", "Tanaidacea", "terrestrial")

#Stacked bar plot of CPUE by location for each trawl.
N2 = ggplot(mysave2N, aes(x = month, y = CPUE))
N2 + geom_bar(stat = "identity", aes(fill = Analy)) + facet_wrap(~Station, scales = "free_y")+
  scale_fill_manual(values = mypal, labels = myslab, name = NULL) +
  xlab("Sampling Month") + ylab("CPUE (count per cubic meter)")+
  theme(legend.position = "bottom")

##############################################################################
#For the June report, I will make some summeries of our data, but don't compare to any of IEP's data because
#we don't have the 20mm data yet.

#Calculate total CPUE for the entire sample (minus mesozooplankton)
myLitsum = group_by(mys3, SampleID, Station, Location, month, Ggdist) %>% 
  summarize(CPUE = sum(CPUE, na.rm = T))

#GLM of log CPUE as predicted by month of the year and distance from the Golden Gate
g1 = lm(log(CPUE)~month+Ggdist, data = myLitsum)
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
Commat = dcast(mys3, SampleID~Analy, value.var = "CPUE", fun.aggregate = sum)
Commat = Commat[,2:18]

#Do a relative abundance matrix
Commatp = Commat/rowSums(Commat)

#PERMANOVA for community abundance:
adonis(Commat~month + Ggdist, data = myLitsum)

#Now for relative abundance:
am = adonis(Commatp~month + Ggdist, data = myLitsum)
am
#Again, there is a significant effect of position on the ecocline (distance from Golden Gate)
#but not month of the year. R-squared is pretty low too. The relative abunance fits better than
#the absolute abundance.

#Calculate NMDS for both absoluate and proportinoal datasets. 
n1 = metaMDS(Commat, trymax = 100)
n2 = metaMDS(Commatp, trymax = 100)
#Both NMDS work about the same, neither is great,
#I'm using the GLM of total catch anyway, so we will go with relative abundance from now on.

#plot the NMDS
source("plotNMDS.R")

PlotNMDS(n2, data = myLitsum, group = "", textp = "T") 
#Ignore the error messages

#plot GAMs of continuous variables on top of the plot
g = ordisurf(n2~Ggdist, data = myLitsum, add = T, labcex = 1)
g2 = ordisurf(n2~month, col = "blue", data = myLitsum, add = T, 
         levels = c(1,2,3,4,5, 6), labcex = 1)
#The "month" GAM is so shitty it doesn't even show up. 

legend(par("usr")[1],par("usr")[3], 
       "Distance from GG",
       col = "red", lty = 1)

