#GLM code for Dave
#By Rosemary Hartman, April 19th, 2019.

#load a few packages. 
library(tidyverse)
library(visreg)
library(lme4)
library(lmerTest)
library(MuMIn)

################################################################################
#GLMs are very easy. I have a data set of the CPUEs from all the blitz data (minus the bethic stuff)
BradSum = read.csv("Summer Bradmoor.csv")

#if I want to look at the structure of the dataset:
str(BradSum)

#to run a GLM on the data, all i need to do is

GLM1 = glm(CPUE~ Gear+Julian.Day, data = BradSum)
summary(GLM1)

#we can look at the effect of each term invididually using the "visreg" function
visreg(GLM1)

#the "default" for the "glm" function is to run a linear model assuming a normal distribution.
#however, if we look at the diagnostic plots for the model:
plot(GLM1)
#We definitely have a skewed distribution. Check out:
hist(BradSum$CPUE)

#Let's try a Poisson distribution, which is more common for count data.
#Note: Because we are using CPUEs insted of counts, we'll have to round to the nearest fish.

GLM2 = glm(round(CPUE)~Gear + Julian.Day, data = BradSum, family = "poisson")
summary(GLM2)
visreg(GLM2)
plot(GLM2)

#Alternatively, we can log-transform the CPUE and try the normal distribution
#Note: I added 1 to each catch because you can't take the log of 0
BradSum$logCPUE = log(BradSum$CPUE+1)
hist(BradSum$logCPUE)
#That looks more normal!
GLM3 = glm(logCPUE~Gear + Julian.Day, data = BradSum)
summary(GLM3)
visreg(GLM3)
plot(GLM3)

#Now we need to decide which predictor variables we want to include in the model. 
#We can use the "dredge" function to evaluate which predictor variables improve
#the model without making it too complicated.

#first I set up a model with everything
global = glm(logCPUE~ Gear + Temp + SpC + Turb + Tide + Julian.Day, data = BradSum, na.action = "na.fail")

#Now I see which of those terms is useful
dredge(global)

#Gear, tide, and turbidity came up as consistantly in all the top models. Julian day, temperature,
#and specific conductance, not so much. 
best.model = glm(logCPUE~Gear + Tide + Turb, data = BradSum)
summary(best.model)
visreg(best.model)

#That is how to do it for one site, but when we met before we discussed using all the sites together with 
#"site" as an error term. So let's load the other site and stick them together
Decker = read.csv("Summer Decker.csv")

#log-transform CPUE
Decker$logCPUE = log(Decker$CPUE + 1)

#add a "site" term to both data sets
Decker$site = "Decker"
BradSum$site = "Bradmoor"

#bind them together
twosites = rbind(Decker, BradSum)

#It might be interesting to compare "shallow" versus "deep" (combining seines and lamparas)
#Note: you could even add the FMWT data in to this data set and code it as "shallow" versus "deep" if you feel
#like living dangerously

#First I'll set up an empty column
twosites$geartype = rep(NA, nrow(twosites))

#now fill it in based on the "gear" collumn
twosites$geartype[which(twosites$Gear=="Townet")] = "deep"
twosites$geartype[which(is.na(twosites$geartype))]= "shallow"

#Now I"ll run a GLM of both sites together
GLM4 = glm(logCPUE~geartype + Turb + site, data = twosites)
summary(GLM4)
visreg(GLM4)

#In GLM4, I included "site" as a predictor variable. If I wan't all that interested in the effect
#of site, but I still wanted to account for differences in site, I could use a mixed model instead.
#This allows me to partition the variance. I don't know that it's needed for this case, but it's useful
#if you are comparing site types or regions and have multiple sites per region.

GLMM1 = lmer(logCPUE~geartype + Turb + (1|site), data = twosites)
summary(GLMM1)
