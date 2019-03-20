#Various ways to plot invertebrate data

# first load the ggplot2 and Rcolor Brewer libraries

library(ggplot2)
library(RColorBrewer)
library(dplyr)

#I made you some test data to play with
zoopstest = read.csv("zooptest.csv")
View(zoopstest)

######################################################################################################
#First, I'll show you how to make the stacked bar plot you wanted.

#First you tell ggplot what is related to what:
p1 = ggplot(zoopstest, aes(x=habtype, 
                               #Here I divide the data by habitat type, but you can put sampleID on the x axis,
                                   #or sampling site, or date, or whatever
                           y= count,  #We want count on the y axis. 
                                        #Right now it will give you the sum of the count in each sample, 
                                         #but you could calculate averages first if you need to
                           fill = taxa)) #Change the color of the bars based on the taxa

#Now we add the "geom_bar" to the plot
p1 + geom_bar(position = "stack", #this tells us to give you a stacked bar plot
              stat = "identity") + #we want the actual values, not the number of occurances
  ylab("total catch") + xlab("Habitat Type") # Some labels


######################################################################################################

#Anther thing I like to do:
p1 + geom_bar(position = "fill", #this tells us to give % composition instead of a regular stacked bar
              stat = "identity") +
  ylab("Percent Composition") + xlab("Habitat Type") # Some labels

#If you'd rather have the columns side-by-side
p1 + geom_bar(position = "dodge",
              stat = "identity") +
  ylab("total catch") + xlab("Habitat Type") 

######################################################################################################

#If you want to show differences between habitat types and sampling locations:
p1 + geom_bar(position = "stack",
              stat = "identity") +
  ylab("total catch") + xlab("Habitat Type") +
#facet_grid allows you to seperate by additional variables. write them in formula notation X~Y, with .~Y if you
  #only have one variable
  facet_grid(.~station) 



######################################################################################################

#It's not too bad if you only have 4 taxa, but if you have a lot you end up with a rainbow and it's
#hard to tell which is which the default color pallete. RColorBrewer gives you better contrasting colors


#set up a color pallete. I had 22 critters in my origional data set, so I combined a few of the brewer pallets
mypal = c(brewer.pal(9, "Set1"), brewer.pal(12, "Set3"), "black")

#If you want to see all the brewer pallettes that are available
display.brewer.all()
#if you are worried about color blind people
display.brewer.all(colorblindFriendly = T)

#Now add the pallette of your choice to your graph using the "scale_fill_manual" command
p1 + geom_bar(position = "stack",
              stat = "identity") +
  ylab("total catch") + xlab("Habitat Type") +
  facet_grid(.~station) +
  #add your pallet to the plot
  scale_fill_manual(values = mypal)

#or if you just want one of the default Brewer palletes
p1 + geom_bar(position = "stack",
              stat = "identity") +
  ylab("total catch") + xlab("Habitat Type") +
  facet_grid(.~station) +
  scale_fill_brewer(type = "qual", pallette = 1)

######################################################################################################

#if you really have a lot of taxa, you should probably lump some of the rare taxa to make
#the graph easier to read

#calculate % total abundance for all the taxa
zoopsTaxa = summarize(group_by(zoopstest, taxa),fract = sum(count)/sum(zoopstest$count))

#add the new data to the origional dataset
zoopstest2 = merge(zoopstest, zoopsTaxa)

#reassign the very small groups to "other"
zoopstest2$taxa[which(zoopstest2$fract<0.01)]="other"
zoopstest2 = droplevels(zoopstest2)

#Plot it again
p2 = ggplot(zoopstest2, aes(x=habtype, y = count, fill = taxa))
p2 + geom_bar(position = "stack", 
              stat = "identity") + 
  ylab("total catch") + xlab("Habitat Type") +
  scale_fill_manual(values = mypal)

######################################################################################################
#Plotting average catch with standard error bars is actually harder than stacked bar plots (in my humble opinion)

#First you have to summarize the data
zoopsSum = summarise(group_by(zoopstest, habtype, taxa), avgcount = mean(count), 
                     sdcount = sd(count), secount = sd(count)/length(count))

#now make the bar plot
p3 = ggplot(zoopsSum, aes(x=habtype, y=avgcount, fill = taxa))

p3 + geom_bar(position = "dodge", stat = "identity") +
  #putting the error bars in the right place gets complicated 
  geom_errorbar(aes(ymin=avgcount-secount, ymax=avgcount+secount, 
                    # it took me forever to figure out that I had to specify the "group" within this function
                    group = taxa),  position =position_dodge(width = 0.9), width = .5)

#if you want a stacked bar with averages instead of the totals:

p3 + geom_bar( stat = "identity")

#Change the labels to make them more informative
p3 + geom_bar( stat = "identity") + xlab("Habitat Type") + ylab("Average Catch")

               