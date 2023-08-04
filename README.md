# data
# Individual and composite life history trait data from water addition experiment
---

We conducted a water addition experiment in a common garden to determine if Brassica tournefortii exhibits phenotypic plasticity. 
Life history traits measurements that describe leaf architecture, plant architecture, plant size, and reproduction, are contained in this data set.  
We entered and cleaned the raw data collected from garden plants, and then performed principal component analysis to reduce trait variables 
to four composite trait variables. Thirteen (13) individual trait variables that we used for calculation of plasticity indices were  included 
in this data set. Relevant categorical variables that were used for our factorial experimental design are also included. The file is the final 
data set we used for all analyses in our study.


## Description of the data and file structure
This data set (.csv file) was collected from plant trait measurements from a garden experiment, conducted by Dr. Alfaro in 2017 
 the University of New Mexico. 
This data set can be opened in SAS and R, as well as other statistical programs that can analyze more than 10,000 data cells. The description of each field in this data set is decribed below:

population - source population and experimental population    
range    - source range of accession or population
treatment - experimental treatment, volume of water added per day
number    - database key for each data row
block    - experimental block
time    - number of days to planting relative to first day of experiment (last day of frost)
family    - seed family
leafnumber - number of leaves per plant    
rosdiam    - rosette diameter
flwtime    - number of days to flower
height - plant height
branchnumber - total number of branches 
branchlength - length of lateral branch
fruits - total number of fruit per plant    
lobes - number of lobes per leaf per plant
leaflength - leaf length 
fruitmass - mean individual fruit mass    
seednumber - mean number of seeds per fruit per plant    
vegbiomass - vegetative biomass    
repallocasin - reproductive allocation, arcsine squareroot transformed         
relfitasin - relative fitness, arcsine squareroot transformed    
vegpc1 - Vegetative Trait PC1    
vegpc2 - Vegetative Trait PC2
reppc1 - Reproductive Trait PC1
reppc2 - Reproductive Trait PC2

## Sharing/Access information

You can also find this data on: 

https://github.com/brianalf-hub/data/blob/main/a%26m2023ee

To retrieve the data on github, go to the link above and copy and paste the data set on to a text file or a spreadsheet. 

* 

Data was derived from the following sources:
Greenhouse experiment by Dr. Brian Alfaro


## Code/Software
The statistical tests and models are described in the manuscript, and can be run with statistical programs that can handle GLM-type ANOVA and ANCOVA models. 

Below are examples of scripts that were used to visualize results. 

#General code for boxcox-based transformations

box <- boxcox(leafnumber ~ 1,            
              lambda = seq(-6,6,0.1))      
cox <- data.frame(box$x, box$y)           
cox2 <- cox[with(cox, order(-cox$box.y)),] # Order the new data frame by decreasing y
cox2[1,]                        # Display the lambda with the greatest log likelihood
lambda <- cox2[1, "box.x"]                 # Extract that lambda
leafnumber.box <- (leafnumber ^ lambda - 1)/lambda   # Transform the original data
plotNormalHistogram(leafnumber.box)
shapiro.test(leafnumber.box)



#PCA procedure and combining files into one data set on a .csv

#prcomp generates an error with missing data, so I removed any rows of data with missing values

# Combining trait variables for Vegetative Trait PCs
veg.pca <- prcomp(~ leafnumber.box + rosdiam.box + height.box + lobes.box + leaflength.box 
                     + vegbiomass.box, center=TRUE, scale=TRUE)
summary(veg.pca)
loadings(veg.pca)
head(veg.pca$scores)
veg.pca

# prcomp codes; Use your pC axes as regression predictors
veg.axes <- predict(veg.pca, newdata = common)
veg.dat <- cbind(common, veg.axes[,1:3])
head(veg.axes, 4)
vegpc1 <- common
vegpc1[names(veg.dat)] <- veg.dat
vegpc1
write.table(vegpc1, "V:/CH1/a_reanalysis/vegpc1.csv", sep="\t") #file opened, formatted, renamed in Excel

# Combining trait variables for Reproductive Trait PCs
rep.pca <- prcomp(~  flwtime.box + branchlength.box
                     + branchnumber.box + fruits.box + fruitmass.box
                     + seednumber.box + repalloc
                   , center=TRUE, scale=TRUE)

summary(rep.pca)
loadings(rep.pca)
head(rep.pca$scores)
rep.pca

# prcomp codes; Use your pC axes as regression predictors
rep.axes <- predict(rep.pca, newdata = common)
rep.dat <- cbind(common, rep.axes[,1:3])
head(rep.axes, 4)
reppc1 <- common
reppc1[names(rep.dat)] <- rep.dat
reppc1
write.table(reppc1, "V:/CH1/a_reanalysis/reppc1.csv", sep="\t") #file opened,formatted, renamed in Excel

plasticity <- merge(
  lifepc, repropc,
  by="number")
write.table(plasticity, "K:/Phenotypic plasticy/reanalysis/gardendata.csv", sep="\t") 


################################################################################################################
#General prepwork for analysis and data visualization procedures

rm(list=ls(all=TRUE)) # Poof, all gone. Tabula rasa!
setwd("Z:/Ch 2 MS")

#Import data set 
library(readr)
rawdata <-read_csv("Z:/Ch 2 MS/gardendata.csv")

# call up the following installed packages-----------------
library(car)
library(Hmisc)
library(RColorBrewer)
library(wesanderson)
library(piecewiseSEM)
library(MASS)
library(splines)
library(extrafont)
library(ggplot2)
library(extrafont)
library(plyr)
loadfonts(device="win")
library(mi)
library(rcompanion)
#library(agricolae)
library(ggpubr)
theme_set(theme_pubr())
library(viridis)
library(mgcv)
library(ggborderline)

View(rawdata)

number <-rawdata$number
block <- rawdata$block
time <- rawdata$time
treatment <- rawdata$treatment
relfit <-rawdata$relfitasin
range <- rawdata$range
population <- rawdata$population
rawdata$block <- factor(rawdata$block)
rawdata$treatment <- factor(rawdata$treatment)
rawdata$family <- factor(rawdata$family)
rawdata$population <- factor(rawdata$population)
rawdata$treatment <- factor(rawdata$treatment)
rawdata$range <- factor(rawdata$range, levels=c("Native", "Invasive", "Landrace"))
rawdata$treatment <- factor(rawdata$treatment, levels=c("250", "450", "750"))

vegpc1 <-rawdata$vegpc1
vegpc2 <-rawdata$vegpc2
reppc1 <-rawdata$reppc1
reppc2 <-rawdata$reppc2


#### General codes for constructing box plots - Figure 3

p <- ggplot(rawdata, aes(x=range, y=reppc1, fill=range)) + 
      ylab("Reproductive Trait PC1\n(reproductive allocation,\n individual fruit mass,\n& number of fruits per plant)") +
     xlab("Range") 
     #guides(color=guide_legend(title="Range", alpha=0.5, title.position = "top", title.hjust = 0.5)) 
p <- p + geom_boxplot(notch="TRUE", size=0.25, width=0.5, outlier.size = 2, alpha=1, show.legend = FALSE) + 
      stat_summary(fun=mean, colour="grey30", geom="point", shape=18, size=5,show.legend = FALSE) 
p <- p + annotate("text", 
                  x = c(1,2,3),
                  y = c(-2.2, -3.6, -3.2),
                  label = c("a", "ab", "b"),
                  size=7, colour="grey30")
p <- p + annotate("text", 
         x = c(0.68, 0.68),
         y = c(-2, 2),
        label = c("High", "Low"),
         size=7, colour="darkgrey")
p <- p +  scale_fill_manual(values=c("darkorchid4", "forestgreen", "goldenrod3"))
p1 <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour="darkgrey", size=0.5))
p <- p + theme(axis.title.x = element_blank(),
               axis.text.x  = element_text(vjust=0.5, size=18, colour="darkgrey"))
p <- p+ theme(axis.title.y = element_text(size=18, colour="darkgrey"),
              axis.text.y  = element_text(size=16, colour="darkgrey"))
p <- p + scale_y_reverse()
p1 <- p + theme(legend.position="none")
print(p1) 

#### General codes for constructing reaction norms - Figure 4

trait <- ddply(rawdata,.(treatment,range),summarise, val = mean(reppc1))
p <- ggplot(rawdata, aes(x = treatment, y =reppc1, colour = range)) + 
  #geom_point(data = trait, aes(y = val)) +
  geom_line(data =trait, aes(y =val, group = range)) +
  theme_bw() +  ylab("Reproductive Trait PC1\n(reproductive allocation,\n individual fruit mass,\n& number of fruits per plant)") + 
  xlab("Water added to soil (ml/day)") 
p <- p + annotate("text", 
                  x = c(0.7, 0.7),
                  y = c(-.85, 1.15),
                  label = c("High", "Low"),
                  color="darkgrey", size=4)
p <- p + scale_y_reverse()               #We reverse the scale for RepPC1, RepPC2, and VegPC1 
p <- p +  scale_colour_manual(name  ="Range",values=c("darkorchid4", "forestgreen", "goldenrod3"))
#p <- p + stat_summary(fun = mean, geom = "point", size=2, aes(group=range))
p <- p + stat_summary(fun = mean, geom = "line", size=1.5, aes(group = range))
p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "darkgray", size=0.25))
p <- p + theme(axis.title.x = element_text(size=12, colour = "grey30"), axis.text.x  = element_text(vjust=0.5, size=12,  colour = "darkgray"))
p <- p + theme(axis.title.y = element_text(size=12, colour = "grey30"), axis.text.y  = element_text(vjust=0.5, size=12,  colour = "darkgray"))
p1<- p + theme(legend.position="none")
print(p1)

#### General codes for constructing gam-smoothed regression lines for planting date as a source of variation  - Figure 5

#Time x range
p <- ggplot(rawdata, aes(x = time, y = reppc1, colour = factor(range))) +
  ylab("Reproductive Trait PC1\n(reproductive allocation,\n individual fruit mass,\n& number of fruits per plant)") + xlab("Planting day\nfrom last frost") 
p <- p + stat_smooth(method = "gam", formula = y ~ s(x, k = 1),
                     se = F, size = 1) 
p <- p +  scale_colour_manual(name  ="Range",values=c("darkorchid4", "forestgreen", "goldenrod3"))
p <- p + annotate("text", 
                  x = c(40),
                  y = c(-1),
                  label = c("P = 0.021"),
                  size=3, colour = "grey30")
p <- p + annotate("text", 
                  x = c(2, 2),
                  y = c(-1.5, 1.15),
                  label = c("High", "Low"),
                  size=3, colour = "darkgray")
p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "darkgray", size=0.5))
p <- p + theme(axis.title.x = element_text(size=10, colour = "grey30"),
               axis.text.x  = element_text(vjust=0.5, size=10, colour = "darkgray"))
p <- p+ theme(axis.title.y = element_text(size=10, colour = "grey30"),
              axis.text.y  = element_text(size=10, colour = "darkgray"))
p <- p+ theme(legend.text = element_text(size = 7, colour = "grey30")) 
p <- p+ theme(legend.title = element_text(size = 7, colour = "grey30")) 
p <- p + scale_y_reverse()
p1 <- p + theme(legend.position = "top", legend.box = "horizontal")
print(p1) 

# General code for producing a panel of multiple graphs in a single figure
require(gridExtra)
tiff('samplefigure.tiff', units="in", width=12.5, height=6.5, res=1100)
samplefigure <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4, nrow = 2, labels = c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"))
samplefigure


# Population CV analyses
# Using the gardendata.csv data set, we used the PROC MEANS procedure in SAS 9.4 to obtain coefficient of variation values

# Individual trait Population CV box plot general code
p <- ggplot(fitness, aes(x=range, y=repalloc, fill=range)) + ylab("Reproductive allocation\n Population CV") + xlab("Range")
p <- p + geom_boxplot(notch=FALSE, size=1, width=0.5, outlier.size = 3,  alpha=0.75)
p1 <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray", size=0.5))
p <- p + theme(axis.title.x = element_blank(),
               axis.text.x  = element_text(vjust=0.5, size=16))
p <- p+ theme(axis.title.y = element_text(size=16),
              axis.text.y  = element_text(size=12))
p <- p +  scale_fill_manual(values=c("darkorchid4", "forestgreen", "goldenrod3"))
p1 <- p + theme(legend.position="none")
print(p1) 

# INdividual trait Population CV vs fitness gam plots general code
p <- ggplot(fitness, aes(x = repalloc, y = relfit, colour = range)) + xlab("Reproductive allocation\n Population CV") + ylab("Relative fitness\n Population mean") + geom_point(size=3)
p <- p +  scale_colour_manual(name  ="Range",values=c("darkorchid4", "forestgreen", "goldenrod3"))
p <- p + geom_smooth(method="lm", se = FALSE, size=1)
p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size=0.5))
p <- p + theme(axis.title.x = element_text(size=18),
               axis.text.x  = element_text(vjust=0.5, size=16))
p <- p+ theme(axis.title.y = element_text(size=18),
              axis.text.y  = element_text(size=16))
p1 <- p + theme(legend.position="none")
print(p1) 

