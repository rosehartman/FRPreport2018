#OK, let's work on the clams.
source("blitzclean.R")

#Linear mixed model of log total CPUE, as predicted by Region2, and Site Type and year. 
#Just the clams
ben.1$logtot = log(ben.1$tCPUE+1)
ben.1$Year2 = as.factor(year(ben.1$Date))

cl = lmer(logtot ~ sitetype+ Region2  + Year2+ (1|site), 
              data = ben.1)
summary(cl)
visreg(cl)


cl2 = aov(logtot ~ sitetype+ Region2 +  Year2 +Error(site), 
              data = ben.1)
summary(cl2)



###########################################################################################
#Now I'll make a plot of total CPUE versus site and Region2
clamsum = group_by(ben.1, sitetype, site, Region2) %>% 
  summarize(mCPUE = mean(tCPUE), sdCPUE = sd(tCPUE), 
            seCPUE = sd(tCPUE)/length(tCPUE), N = length(tCPUE))

clP = ggplot(clamsum, aes(x = site, y = mCPUE))
clP+ geom_bar(stat = "identity", aes(fill = site)) + 
  facet_grid(.~Region2, scales = "free", space = "free_x") +
  geom_errorbar(aes(ymin = mCPUE - seCPUE, ymax = mCPUE + seCPUE), width = 0.7) +
  geom_label(aes(label = paste("n = ", N), y = mCPUE*1.1), label.padding = unit(0.1, "lines"), size = 3) +
  scale_fill_manual(values = mypal, guide = "none") + ylab("CPUE") 

clP+ geom_bar(stat = "identity", aes(fill = site)) + 
  facet_grid(.~sitetype, scales = "free", space = "free_x") +
  geom_errorbar(aes(ymin = mCPUE - seCPUE, ymax = mCPUE + seCPUE), width = 0.7) +
  geom_label(aes(label = paste("n = ", N), y = mCPUE*1.1), label.padding = unit(0.1, "lines"), size = 3) +
  scale_fill_manual(values = mypal, guide = "none") + ylab("CPUE") 

