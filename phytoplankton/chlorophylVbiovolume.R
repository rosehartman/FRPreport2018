#See whether the chlorophyll values match the phytoplankton counts

library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)
library(vegan)
library(RColorBrewer)
source("phytoplankton/phytosblitz2018.R")

#plot BPUE versus CPUE
p = ggplot(physum, aes(x=CPUE, y = BPUE))
p + geom_point()

#quick plot of chlorophyll versus biovolume from cell counts

p1 = ggplot(filter(physum, site != "Liberty"), aes(x=BPUE, y = Chlorophyll))
p1 + geom_point() + geom_smooth(method = "lm")+ xlab("biovolume in u3/mL") +
  ylab("Chlorophyll ug/L")

#I don't know why the model I fit with "lm" looks so much different than the one I fit with ggplot
l1 = glm(Chlorophyll~BPUE, data = filter(physum, !is.na(Chlorophyll)))
summary(l1)
plot(l1)

#now plot chlorophyll versus total cell count
p2 = ggplot(filter(physum, site != "Liberty"), aes(x=CPUE, y = Chlorophyll))
p2 + geom_point() + geom_smooth(method = "lm") + xlab("cells/mL") +
  ylab("Chlorophyll ug/L")


l2= lm(Chlorophyll~CPUE, data = physum)
summary(l2)

#try some log-transformations
p3 = ggplot(filter(physum, site != "Liberty"), aes(x=log(CPUE), y = Chlorophyll))
p3 + geom_point() + geom_smooth(method = "lm") + xlab("log(cells/mL)") +
  ylab("Chlorophyll ug/L")

#linear models
m1 = lm(Chlorophyll ~ log(CPUE), data = filter(physum, site != "Liberty"))
summary(m1)
#it is significant, but the R-squared is only 0.55

m2 = lm(Chlorophyll ~ vol, data = filter(physum, Location != "Liberty Is."))
summary(m2)

#let's try limiting it to chlorophyll values under 30 ug/L

#linear models
m3 = lm(Chlorophyll ~ log(CPUE), data = filter(physum, Chlorophyll < 30 &site != "Liberty"))
summary(m3)

m4 = lm(Chlorophyll ~vol, data = filter(physum, Chlorophyll < 30 &Location != "Liberty Is."))
summary(m4)
#so it's definitely breaking down at low chlorophyll values
#now plot chlorophyll versus total cell count


p1 + geom_point() + geom_smooth(method = "lm") + xlab("Biovolume u3/mL") +
  ylab("Chlorophyll ug/L") + scale_y_continuous(limits = c(0,30)) +
scale_x_continuous(limits = c(0, 300000000))

p2 + geom_point() + geom_smooth(method = "lm") + xlab("cells/mL") +
  ylab("Chlorophyll ug/L") + scale_y_continuous(limits = c(0,30))+
 scale_x_continuous(limits = c(0, 16000))

p3 + geom_point() + geom_smooth(method = "lm") + xlab("cells/mL") +
  ylab("Chlorophyll ug/L") + scale_y_continuous(limits = c(0,30))+
  scale_x_continuous(limits = c(4, 10))

