#look at whether nutrients are limiting by using the half-saturation constants for phytoplankton
#growth as cited by Cloern and Duffort (2005)
# K(N) = 0.7 micromolar
# K(P) = 0.1 micromolar
source("nutrient time series.R")

#The inorganic components are the ones that are most bioavailalbe.
totalN = filter(nuts2018, analyte == "Dissolved Nitrate + Nitrite" | analyte == "Dissolved Ammonia") %>%
  group_by(bryteID, site, sitetype, substa, region, Date, month) %>%
  summarize(DIN = sum(concentration))

#Convert mg/L to micromolar
totalN$micromolarN = totalN$DIN*1000/14.006

#calculate limitation index
totalN$limit = totalN$micromolarN/0.7


#Do a quick plot to see what the distribution of limitation indexes is.
ggplot(totalN, aes(x=limit, stat(density))) + geom_histogram() + 
  facet_grid(.~sitetype, scales = "free_x") + xlab("limitation index")

#How often is nitrogen limiting?
limiting = filter(totalN, limit <1)

#Let's look at phosphate
totalP = filter(nuts2018, analyte == "Dissolved Ortho-phosphate")

#Convert mg/L to micromolar
totalP$micromolarP = totalP$concentration*1000/30.973

#calculate limitation index
totalP$limit = totalP$micromolar/0.1
limiting = filter(totalP, limit <1)


#How often is phosphorus limiting?
ggplot(totalP, aes(x=limit, stat(density))) + geom_histogram() + facet_grid(.~sitetype)+
  xlab("limitation index")

#calculate N:P ratio
redfield = merge(totalN, totalP, by = "bryteID")
redfield$ratio = redfield$micromolarN/redfield$micromolarP
hist(redfield$ratio)

ggplot(redfield, aes(x=ratio, stat(density))) + geom_histogram() + 
  facet_grid(.~sitetype.x) + geom_vline(aes(xintercept = 16), color = "red")


####################################################################################################
#Let's see what percentage of nitrogen comes from different sources
allN = filter(nuts2018, analyte == "Total Kjeldahl Nitrogen" | analyte == "Dissolved Ammonia"|
            analyte ==   "Dissolved Nitrate + Nitrite"| analyte ==    "Dissolved Organic Nitrogen" )

#stacked bar graph of nitrogen types
ggplot(filter(allN, analyte != "Total Kjeldahl Nitrogen"), aes(x=site, y= concentration, fill = analyte)) +
  geom_bar(position = "fill", stat = "identity")


ggplot(filter(allN, anlyte != "Total Kjeldahl Nitrogen" & substa != "EMP"), aes(x=bryteID, y= concentration, fill = anlyte)) +
  geom_bar(position = "fill", stat = "identity") + facet_wrap(~site, scales = "free_x")


ggplot(filter(allN, analyte != "Total Kjeldahl Nitrogen" & substa != "EMP"), 
       aes(x=as.factor(month(Date)), y= concentration, fill = analyte)) +
  geom_bar(position = "fill", stat = "identity") + facet_wrap(~site, scales = "free_x") +
  xlab("Month")
