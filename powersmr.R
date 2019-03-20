#Let's make a plot with power analyses for all the gear types on one graph
#Note: As written, this will take a very long time to run. To trouble shoot/test the code,
#replace "nsim = 999" with "nsim = 20"


#sweep nets
source("sweepnets_clean2.R")
fixef(s1)

#let's say we want to see a 0.5 or 1 log CPUE difference. I've 
#assumed things with higher differences in the origional model will probalby have more differences, so we
#can set them as larger in the power analysis

#first sweep nets
fixef(s1)

#what's the current power?
#spp2r = powerSim(s1, test = fixed("Region"), nsim = 999)
#spp2st = powerSim(s1, test = fixed("SiteType2"), nsim = 999)
#spp2t = powerSim(s1, test = fixed("targets"), nsim = 999)

#new model with set fixed effects
s2 = s1
fixef(s2)= c(0, -1, -1, -1, 1, 1, 1, 0.5, 1)

#let's look at up to twenty samples
s3 = extend(s2, within = "Region + SiteType2+ targets", n = 20)
spc2.1 = powerCurve(s3, nsim = 999, within = "Region + SiteType2 + targets", test = fixed("Region"), breaks = c(1,2,3,4,5,7,10, 15, 20))
spc2.1
plot(spc2.1)
#Hmm. It's not giving us more than 50% power no matter what.
snpwr = summary(spc2.1)
snpwr$predictor = rep("Region", 9)

#now by site type
spc2.2 = powerCurve(s3, nsim = 999, within = "Region + SiteType2 + targets", test = fixed("SiteType2"), breaks = c(1,2,3,4,5,7,10, 15, 20))
spc2.2
plot(spc2.2)

snpwr2 = summary(spc2.2)
snpwr2$predictor = rep("Site Type", 9)

#now by habitat type
spc2.3 = powerCurve(s3, nsim = 999, within = "Region + SiteType2 + targets", test = fixed("targets"), breaks = c(1,2,3,4,5,7,10, 15, 20))
spc2.3

snpwr3 = summary(spc2.3)
snpwr3$predictor = rep("Vegetation Type", 9)

#bind the results together
snpwr = rbind(snpwr, snpwr2, snpwr3)
snpwr$var = rep("sweep net", 27)
######################################################################################################################
#now for zoops

source("zoopsblitz2.R")
fixef(zblitz)

#what's the current power? of the whole model?
zpp2 = powerSim(zblitz,  nsim = 999)
print(zpp2)

#for differences between regions?
zpp2r = powerSim(zblitz, test = fixed("Region"), nsim = 299)
zpp2r

#for differences between site types?
zpp2st = powerSim(zblitz, test = fixed("SiteType2"), nsim = 999)
zpp2st

#Now a new model with fixed effects
zblitz2 = zblitz
fixef(zblitz2)= c(0, -1, -1, -1, 1, -1, -1)
zpc2 = powerSim(zblitz2,  test = fixed("Region"), nsim = 299)
zpc2

#So we are pretty good on power at that effect size.
zblitz2 = extend(zblitz, within = "Region + SiteType2", n = 20)
zpc2.1 = powerCurve(zblitz2, nsim = 999, within = "Region + SiteType2", test = fixed("Region"), breaks = c(1,2,3,4,5,7,10, 15, 20))
zpc2.1
#nice!
zpwr = summary(zpc2.1)
zpwr$predictor = rep("Region", 9)

zpc2.2 = powerCurve(zblitz2, nsim = 999, within = "Region + SiteType2", test = fixed("SiteType2"), breaks = c(1,2,3,4,5,7,10, 15, 20))
zpc2.2
#nice!
zpwr2 = summary(zpc2.2)
zpwr2$predictor = rep("Site Type", 9)

zpwr = rbind(zpwr, zpwr2)
zpwr$var = rep("zooplankton", 18)

###########################################################################################################
#chlorophyll
source("chloro_blitz2017.R")

#power with current sample size

#differences between regions
#pp1r = powerSim(c1, test = fixed("Region"), nsim = 999)
#pp1r

#differences between site types
#pp1st = powerSim(c1, test = fixed("SiteType"), nsim = 999)
#pp1st

#OK, it looks like we shouldn't use the observed effect size, because results can be misleading. Instead 
# we should decide what effect we are interested in. Let's say we want to be able to tell the
#difference between 3 ug/L chlorophyll
fixef(c1)

c2 = c1
fixef(c2) = c(0,-3,3,3,3,3,-3)

#So we are definitely low on power. But how many samples do we need?
c3 = extend(c2, within = "Region + Sitetype", n = 20)

#first look by site type
pc2.1 = powerCurve(c3, nsim = 999, within = "Region + Sitetype",test = fixed("Sitetype"), breaks = c(1,2,3,4,5,7,10,15,20))
summary(pc2.1)

cpwer = summary(pc2.1)
cpwer$predictor = rep("Site Type", 9)

#now by region
pc2.2 = powerCurve(c3, nsim = 999, within = "Region + Sitetype",test = fixed("Region"), breaks = c(1,2,3,4,5,7,10,15,20))
summary(pc2.2)

cpwer2 = summary(pc2.2)
cpwer2$predictor = rep("Region", 9)

cpwer = rbind(cpwer, cpwer2)
cpwer$var = rep("chlorophyll", 9)

#####################################################################################################
#other macroinvertebrate data
source("blitzuni20172.R")
#power with current sample size

#differences between regions
bpp1r = powerSim(blitz2, test = fixed("Region"), nsim = 99)
bpp1r

#differences between site types
#bpp1st = powerSim(blitz2, test = fixed("SiteType2"), nsim = 999)
#bpp1st

#differences between habitat types
#bpp1t = powerSim(blitz2, test = fixed("targets2"), nsim = 999)
#bpp1t



#OK, it looks like we shouldn't use the observed effect size, because results can be misleading. Instead 
# we should decide what effect we are interested in. Let's say we want to be able to tell the
#difference between 1log CPUE
fixef(blitz2)

blitz3 = blitz2
fixef(blitz3) = c(0,-1,-1,-1,-1,1,-1,1,1,1)

#So we are definitely low on power. But how many samples do we need?
blitz4 = extend(blitz3, within = "Region + SiteType2 + targets2", n = 20)

#look at region
bpc2.1 = powerCurve(blitz4, nsim = 999, within = "Region + SiteType2+targets2", test = fixed("Region"),breaks = c(1,2,3,4,5,7,10,15,20))
summary(bpc2.1)

bpwer = summary(bpc2.1)
bpwer$predictor = rep("Region", 9)

#look at site type
bpc2.2 = powerCurve(blitz4, nsim = 999, within = "Region + SiteType2+targets2", test = fixed("SiteType2"),breaks = c(1,2,3,4,5,7,10,15,20))
summary(bpc2.2)

bpwer2 = summary(bpc2.2)
bpwer2$predictor = rep("Site Type", 9)

#look at gear type
bpc2.3 = powerCurve(blitz4, nsim = 999, within = "Region + SiteType2+targets2", test = fixed("targets2"),breaks = c(1,2,3,4,5,7,10,15,20))
summary(bpc2.3)

bpwer3 = summary(bpc2.3)
bpwer3$predictor = rep("Gear Type", 9)

bpwer = rbind(bpwer, bpwer2, bpwer3)
bpwer$var = rep("macroinverts (all)", 27)

########################################################################################################
#benthic samples
benm = lmer(log(tCPUE)~ Region + SiteType2 + (1|Location), data = ben.1)
#power with current sample size and observed differences
#bp1r  = powerSim(benm, test = fixed("Region"), nsim = 999)
#bp1r
#bp1st  = powerSim(benm, test = fixed("SiteType2"), nsim = 999)
#bp1st

fixef(benm)
benm2 = benm
fixef(benm2) = c(0,-1,1,-1,1,1,1)

#expand to twenty samples per site, make power curves
benm3 = extend(benm2, within = "Region + SiteType2", n = 20)
bpc2 = powerCurve(benm3, nsim = 999, within = "Region + SiteType2",test = fixed("Region"), breaks = c(1,2,3,4,5,7,10,15,20))
summary(bpc2)

benpwer = summary(bpc2)
benpwer$predictor = rep("Region", 9)

#now look between site types
bpc3 = powerCurve(benm3, nsim = 999, within = "Region + SiteType2", test = fixed("SiteType2"), breaks = c(1,2,3,4,5,7,10,15,20))
summary(bpc3)

benpwer2 = summary(bpc3)
benpwer2$predictor = rep("Site Type", 9)

benpwer = rbind(benpwer, benpwer2)
benpwer$var = rep("benthic", 18)

###############################################################################################################

#mysids power analysis


my1 = lmer(log(tCPUE)~ Region + SiteType2 + (1|Location), data = mys.1)
#power with current sample size and observed differences
#myp1r  = powerSim(my1, test = fixed("Region"), nsim = 999)
#myp1r
#my1p1st  = powerSim(my1, test = fixed("SiteType2"), nsim = 999)
#my1p1st

fixef(my1)
my2 = my1
fixef(my2) = c(0,-1,-1,-1,1,1,1)

#expand to twenty samples per site, make power curves
my3 = extend(my2, within = "Region + SiteType2", n = 20)
mpc2 = powerCurve(my3, nsim = 999, within = "Region + SiteType2",test = fixed("Region"), breaks = c(1,2,3,4,5,7,10,15,20))
summary(mpc2)

mypwer = summary(mpc2)
mypwer$predictor = rep("Region", 9)

#now look between site types
mpc3 = powerCurve(my3, nsim = 999, within = "Region + SiteType2",test = fixed("SiteType2"), breaks = c(1,2,3,4,5,7,10,15,20))
summary(mpc3)

mypwer2 = summary(mpc3)
mypwer2$predictor = rep("Site Type", 9)

#combine results
mypwer = rbind(mypwer, mypwer2)
mypwer$var = rep("mysids", 18)

##########################################################################################################
#Neuston power analysis

#linear models
neu.1 = droplevels(ungroup(neu.1))
n1 = lmer(log(tCPUE)~ Region + SiteType2 + (1|Location), data = neu.1)
#power with current sample size and observed differences
#np1r  = powerSim(n1, test = fixed("Region"), nsim = 999)
#np1r
#np1st  = powerSim(n1, test = fixed("SiteType2"), nsim = 999)
#np1st

fixef(n1)
n2 = n1
fixef(n2) = c(0, 1,1,1,1,1)

#expand to twenty samples per site, make power curves
n3 = extend(n2, within = "Region + SiteType2", n = 20)

#first look at power to differentiate between regions
npc2 = powerCurve(n3, nsim = 100, within = "Region + SiteType2",test = fixed("Region"), breaks = c(1,2,3,4,5,7,10,15,20))
summary(npc2)
npwer = summary(npc2)
npwer$predictor = rep("Region", 9)

#now between site types
npc2 = powerCurve(n3, nsim = 100, within = "Region + SiteType2",test = fixed("SiteType2"), breaks = c(1,2,3,4,5,7,10,15,20))
summary(npc2)
npwer2 = summary(npc2)
npwer2$predictor = rep("Site Type", 9)

#combined data set of results
npwer = rbind(npwer, npwer2)
npwer$var = rep("neuston", 18)


###########################################################################################################
#bind the power analyses together
power.list = rbind(cpwer, zpwr, snpwr, bpwer, benpwer, mypwer, npwer)

pcp1 = ggplot(power.list, aes(x=nlevels, y = mean, color = var))
pcp1 + geom_point() + geom_line() + geom_errorbar(aes(ymin = lower, ymax = upper)) +
 geom_hline(aes(yintercept=0.8), lty = 2) +
  ylab("Predicted power") + xlab("Number of sampels per site") + 
  facet_wrap(~predictor, scales = "free_x") 

#Let's just plot region and site type
pcp2 = ggplot(filter(power.list, predictor != "Gear Type" & predictor != "Vegetation Type"), aes(x=nlevels, y = mean, color = var))
pcp2 + geom_point() + geom_line() + geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_hline(aes(yintercept=0.8), lty = 2) +
  ylab("Predicted power") + xlab("Number of samples per gear type per site") + 
  facet_wrap(~predictor, scales = "free_x") +
  scale_color_manual(values = mypal, name = NULL) +
  theme(legend.position = "bottom")

#Save the restults of the analysis so I don't have to re-run the whole thing to make the graph
write.csv(power.list, "powerlist.csv")
