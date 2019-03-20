# Eventually I should do something with all my lenghts

#Per Steve S:
#After generating wet mass per length, then convert to Carbon mass.   
#The most common conversion factors used are reported by Cushing (1958) with Wet mass to 
#Dry mass of * 0.20 then Dry Mass to Carbon mass conversion of * 0.60.   
#Or you can just do wet (mg) to carbon (mg C) conversion of * 0.12,
#as reporting in the attached NOAA plankton report Table 2 on page 4.  

# I don't have length-weight relationships for very many taxa. Mostly:
#gammarus (mabye)
#hyallela
#crangonyx
#Americorophium spinicorne
#coeagrionidae
#corixidae
#gnorimospereoma
#a few snails

#I have biomass estimates for the zooplankton too, average biomass per critter, not length-weight

#A few options:

#Analyze biomass differences with just amphipods (or somethign)
#look at lengns within a group, ignore biomass
#Comparisons between taxa only work if I have biomss...
#lump things into larger categories with really rough biomass estimates

#Can use equations from Benke et. al. 1999

#If I'm just going to use one taxon at a time, should I jsut use lengths?

#waht were the most common bugs I had to work with?
source("blitzclean.R")

bugranks = group_by(bugsblitz, CommonName) %>% summarize(tot = sum(CPUE))
bugranks = bugranks[order(bugranks$tot),]

#Upload the data with the lengths
bugslegnths <- read.csv("~/Documents/phaseIV/bugexportqrylegnths.txt")

#Extrapolate lengths to the plus-counts

buglengthsPC = group_by(bugslegnths, AdjCount, SampleID, CommonName) %>% 
  summarize(avelength = mean(length, na.rm = T), 
            medlength = median(length, na.rm=T),
            sdlength = sd(length, na.rm = T),
            nlength = length(length))

#remove the critters with no length meansurements
buglengthsPC = filter(buglengthsPC, is.nan(avelength)== F)
buglengthsPC$pluscount = ceiling(buglengthsPC$AdjCount) - buglengthsPC$nlength

#We need to make a new data frame with one row for each bug, plus-counted or not.
#let's start with the plus-counted
#get rid of the ones with no plus counts, and the ones with zeros for length
buglengthsPC2 = filter(buglengthsPC, pluscount >0)
buglengthsPC2 = filter(buglengthsPC2, avelength >0)

#now make a row for each bug
bugsPC = buglengthsPC2[rep(seq_len(nrow(buglengthsPC2)), buglengthsPC2$pluscount),]

#initialize a new "length" row
bugsPC$length = rep(NA, nrow(bugsPC))

#generate a random legnth for each bug
test3 = for (i in 1:nrow(bugsPC)) {
  av = unlist(bugsPC[i,4])
sdl = unlist(bugsPC[i,6])
if(is.na(sdl)) bugsPC$length[i] = av else
bugsPC$length[i] = unlist(abs(rnorm(1, av, sdl)))
i = i+1
}

#that took much too long. THere's got to be a better way, but I'm having trouble
#coding ti right now.

#Here's me trying to use the apply function
foo = function(m){
  avelength = unlist(m[4])
  sdlength = unlist(m[6])
  if(is.na(sdlength)) x = avelength else
  x = rnorm(1, avelength, sdlength)
  return(x)
}

test = sapply(bugsPC,  foo)
#that didn't work


#Well, now we can put this back with the bugs we did measure.
#########################################################################################################
#filter out all the measured bugs
bugleg2 = filter(bugslegnths, is.na(length)== F)[,c("SampleID", "CommonName", "length")]

#attatched the measured bugs to the extrapolated measurements
bugsall = rbind(bugleg2, bugsPC[,c("SampleID", "CommonName", "length")])

#Add the sample information back in
bugsamples = bugslegnths[!duplicated(bugslegnths$SampleID),1:19]
bugsall2 = merge(bugsamples, bugsall)

#add the analysis groups

#add some analysis groups
zoocodes <- read_excel("zoocodes.xlsx")
zoocodes <- zoocodes[,c(3,12,13, 14)]
bugsall2 = merge(bugsall2, zoocodes)
bugsall2 = bugsall2[order(bugsall2$SampleID),]

#we only want the blitz samples for nwo
blitzlength = inner_join(bugsall, bugsblitz)

############################################################################################################
#start with just the amphipods
amps = filter(blitzlength, Analy2 == "Amphipoda" & targets2 != "benthic")


ggplot(amps, aes(x=length)) + geom_density() + facet_grid(targets2~sitetype, scales = "free_y")
ggplot(amps, aes(x=length)) + geom_histogram() + facet_grid(targets2~sitetype, scales = "free_y")
ggplot(amps, aes(x=length)) + geom_histogram() + facet_grid(targets2~sitetype)

#I want to convert it to biomass. I don't have great equations for all the bugs, but I do have
#them for the major groups of amphipods

#Import the measurements to R to calculate regressions
LWreg <- read_excel("InvertLW_forR.xlsx", 
                    sheet = "LWreg")

#plot legth in mm versus weight in mg
ggplot(LWreg, aes(x= L.mm, y = DW.mg, color = CommonName)) + geom_point()  + 
  geom_smooth(method = "lm", formula = y~exp(x)) + xlab("length (mm)")+
  ylab("weight (mg)") 

#now log legth versus log weight
ggplot(LWreg, aes(x= Ln.L, y = logDWmg, color = CommonName)) + geom_point() + geom_smooth(method = "lm")+
  geom_abline(intercept = log(0.0058), slope = 3.015) + xlab("Log length (mm)") + ylab("log dry weight (mg)")

Hy = lm(logDWmg~Ln.L, data = filter(LWreg, CommonName == "Hyallela"))

LWregs = group_by(LWreg, CommonName) %>% summarize(a = lm(logDWmg~Ln.L)$coefficients[1], 
                                                   b = lm(logDWmg~Ln.L)$coefficients[2])
LWregs$exa = exp(LWregs$a)
LWregs$exb = exp(LWregs$b)
#Nicole's Gammarus and Corophium values look screwy, so I'll just use literature values for those

#I made a table with my Crangonyx and Hylella values, plus literature values (Benke et al. 1999)
#for the other amphipods
LWregs2 = read_excel("zoocodes.xlsx", sheet = "LWregs")

#first convert the length from microns to mm
amps$lengthmm = amps$length/10
amps = merge(amps, LWregs2)
amps$DW = amps$a*(amps$lengthmm)^(amps$b)

#Now sum by sample
amps$logDW = log(amps$DW)
ampssum = group_by(amps, SampleID, Region2, site, sitetype, targets2, Date, Sampletype) %>%
  summarize(ampW = sum(DW, na.rm = T), aveW = mean(logDW))

#average by site and sample type
ampsave = group_by(ampssum, site, sitetype, Region2, targets2) %>% summarize(ampave = mean(ampW, na.rm = T), 
                                                                    ampse = sd(ampW)/length(ampW))
ampsave$site = factor(ampsave$site, levels = c("Ryer", "Grizzly", "Tule Red", "Blacklock", "Bradmoor", "LHB", "Browns",
                                                 "Winter", "Broad", "Horseshoe", "Decker", "Stacys", "Miner", "Prospect",
                                                 "Liberty", "Lindsey", "Flyway"))

#plot average amphipod biomass by site
ggplot(ampsave, aes(x=site, y = ampave, fill = site)) + geom_bar(stat = "identity") + 
  facet_grid(targets2~sitetype, scales = "free", space = "free_x") + ylab("average mass of amphipods (mg)") +
  geom_errorbar(aes(ymin  = ampave + ampse, ymax = ampave- ampse), width = 0.25) +
  scale_fill_manual(values = mypal, guide = "none") + 
  theme(strip.text = element_text(size = 16), axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16)) 
 # coord_cartesian(ylim = c(0,0.25))

#break out the diked wetlands because Tule red was so big.
ggplot(filter(ampsave, sitetype == "diked"), aes(x=site, y = ampave, fill = site)) + geom_bar(stat = "identity") + 
  facet_grid(targets2~sitetype, scales = "free", space = "free_x") + ylab("average mass of amphipods (mg)") +
  geom_errorbar(aes(ymin  = ampave + ampse, ymax = ampave- ampse), width = 0.25) +
  scale_fill_manual(values = mypal, guide = "none") + 
  theme(strip.text = element_text(size = 16), axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16)) 

ggplot(filter(ampsave, site == "Tule Red"), aes(x=site, y = ampave, fill = site)) + geom_bar(stat = "identity") + 
  facet_grid(targets2~sitetype, scales = "free", space = "free_x") + ylab("average mass of amphipods (mg)") +
  geom_errorbar(aes(ymin  = ampave + ampse, ymax = ampave- ampse), width = 0.25) +
  scale_fill_manual(values = mypal, guide = "none") + 
  theme(strip.text = element_text(size = 16), axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16)) 


#average by site type
ampsave2 = group_by(ampsave, sitetype, Region2, targets2) %>% summarize(ampave2 = mean(ampave, na.rm = T), 
                                                               ampsd = sd(ampave)/length(ampave))

ggplot(ampsave2, aes(x=sitetype, y = ampave2)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = ampave2 + ampsd, ymax = ampave2 - ampsd), width = 0.25) +
  facet_grid(targets2~., scales = "free_x", space = "free_x") + ylab("average mass of amphipods (mg)")



#linear model of total amphipod mass
hist(ampssum$ampW)

LW1 = lmer(log(ampW)~sitetype +targets2+ (1|site), data = ampssum)
LW2 = aov(log(ampW)~sitetype +targets2+ site, data = ampssum)

summary(LW1)
visreg(LW1)
TukeyHSD(LW2)

#Diked wetlands had higher biomass than everything else, so did vegetation

#Maybe it would be better to look at size of ampihpod, rather than mean mass per sample

ggplot(amps, aes(x=log(DW))) + geom_density() + facet_grid(targets2~sitetype, scales = "free")
ggplot(amps, aes(x=DW)) + geom_histogram() + facet_grid(targets2~sitetype, scales = "free_y")
ggplot(amps, aes(x=DW)) + geom_histogram() + facet_grid(targets2~sitetype)

#I should calculate the average mass per amphipod, and model that

ampssum$Region2 = factor(ampssum$Region2, levels = c("Grizzly", "NurseDenverton", "Confluence", "SacSanJ", "Cache"))
LW3 = lmer(aveW~sitetype +targets2+ Region2 +(1|site), data = ampssum)
LW4 = aov(aveW~sitetype +targets2+ Region2 + site, data = ampssum)

summary(LW3)
visreg(LW3)
TukeyHSD(LW4)

#plot mean individual amphipod mass
ampsave2 = group_by(ampssum, site, sitetype, Region2, targets2) %>% summarize(ampave = mean(aveW, na.rm = T), 
                                                                    ampse = sd(aveW)/length(aveW))
ampsave2$ampDW = exp(ampsave2$ampave)

ggplot(ampsave2, aes(x=site, y = ampDW, fill = Region2)) + geom_bar(stat = "identity") + 
  facet_grid(targets2~sitetype, scales = "free_x", space = "free_x") + ylab("average amphipod mass (mg)") +
  geom_errorbar(aes(ymin  = exp(ampave + ampse), ymax = exp(ampave- ampse)), width = 0.25) +
  scale_fill_manual(values = mypal, name = "Region") + 
  theme(strip.text = element_text(size = 16), axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16), legend.position = "bottom") 
