#Look at just the neuston samples
source("2017samplesBlitz.R")

neu= filter(bugs2017, targets == "neuston")
neu.1 = filter(bugs2017.1, targets == "neuston")
neu.1$Region = as.factor(neu.1$Region)
neu.1$SiteType = as.factor(neu.1$SiteType)
neu.1$SiteType2 = as.factor(neu.1$SiteType2)
neu.2 = filter(bugs2017.2, targets == "neuston")
neusum = filter(bugssum, targets == "neuston")
neu.2x = filter(bugs2017.2x, targets == "neuston")
neusum.1 = filter(bugssum.1, targets == "neuston")


n1 = ggplot(neusum, aes(x= Location, y = mCPUE))
n1 + geom_bar(stat = "identity", aes(fill = Analy)) + facet_wrap(~Region, scales = "free_x")

n1.1 = ggplot(neusum.1, aes(x= Location, y = mCPUE))
n1.1 + geom_bar(stat = "identity", aes(fill = Analy2)) + 
  facet_wrap(~Region, scales = "free") + scale_fill_manual(values = mypal)

nGZ = ggplot(filter(neusum.1, Location == "Grizzly"), aes(x=Analy2, y = mCPUE, fill = Analy2))
nGZ + geom_bar(stat = "identity")+ scale_fill_manual(values = mypal) + 
  geom_errorbar(aes(ymin = mCPUE - seCPUE, ymax = mCPUE + seCPUE)) + 
  theme(axis.text.x = element_text(angle = 90))

#community data matrix
ComMatn = dcast(neu.2, formula = SampleID~Analy, value.var="CPUE", 
                fun.aggregate = sum, fill = 0)
row.names(ComMatn) = ComMatn$SampleID
ComMat2n = ComMatn[,2:ncol(ComMatn)]
#make one based on proportion of total catch for each sample
ComMatpn = ComMat2n/rowSums(ComMat2n)

#PERMANOVA with CPUE
ab1 = adonis(ComMat2n~Region+SiteType, data = neu.1)
ab1

#PERMANOVA with relative abundance
pn1 = adonis(ComMatpn~Region+SiteType2, data = neu.1)
pn1
#NMDS
NMDS1n = metaMDS(ComMat2n, trymax = 100)
NMDS2n = metaMDS(ComMatpn, trymax = 100)


source("plotNMDS.R")

PlotNMDS(NMDS1n, data = neu.1, group = "Region")
PlotNMDS(NMDS1n, data = neu.1, group = "SiteType2")

##########################################################################
#power analysis

#linear models
n1 = lm(log(tCPUE)~ Region + SiteType2, data = neu.1)
summary(n1)
n2 = lm(log(tCPUE) ~ Region , data = neu.1)
n4 = lm(log(tCPUE) ~ SiteType2 , data = neu.1)


#differences between regions
#power with given sample size
Rf = summary(n2)$r.squared[1]/(1-summary(n2)$r.squared[1])
pwr.f2.test(summary(n2)$df[1], v = n2$df.residual, f2 = Rf, sig.level = 0.05)
#sample size at 80% power
pwr.f2.test(summary(n2)$df[1],  f2 = Rf, sig.level = 0.05, power = 0.8)
#We need 29 samples for enough power to differentiate between regions


#differences between site types
#power with given sample size
stf = summary(n4)$r.squared[1]/(1-summary(n4)$r.squared[1])
pwr.f2.test(summary(n4)$df[1], v = n4$df.residual, f2 = stf, sig.level = 0.05)
#sample size at 80% power
pwr.f2.test(summary(n4)$df[1],  f2 = stf, sig.level = 0.05, power = 0.8)
#We need 70 samples for enough power to differentiate between site types (when there are 4 site types)

