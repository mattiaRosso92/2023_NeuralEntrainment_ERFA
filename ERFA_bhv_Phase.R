# Import libraries
library(readxl)
library(ggplot2)
library(broom)
library(dplyr)
# Statistical package
library(lme4)
# Data arrangement
library(reshape2) # contains melt() function for long univariate format
library(ez)
library(psych)
### Normality tests ###
library(EnvStats)
### Colors in plot ###
library(RColorBrewer)
library(wesanderson)
### Tables for publication ###
library(sjPlot)
library(gridExtra)

# Import data (set user's path before file name)
dataPhase <- read_excel('/Users/mattiaipem/Documents/IPEM/Projects/ Lou _ Entrainment/ERFA_Metronomes/Data/Experiment/ERFA_bhv_Phase.xlsx')
# Convert variables to factors
dataPhase$Direction <- factor(dataPhase$Direction)
# Re-order factors, for interpretation
dataPhase$Direction <- relevel(dataPhase$Direction , 'Null')
contrasts(dataPhase$Direction) #set baseline for contrast across modalities

# Remove outlier: Subject 6 (missing data) and eventual outliers
dataPhase <- subset(dataPhase, Subject != 6)
# dataPhase <- subset(dataPhase, Subject != 1 & Subject != 6 & Subject != 17)


### Model fitting ###
model <- lmer(Integral ~ Direction + (1 | Subject) ,
              control = lmerControl(optimizer = "bobyqa") , 
              data = dataPhase , REML=FALSE)


#Get parameter estimates for all models
model.coefs <- data.frame(coef(summary(model)))

#...and estimate p-values
model.coefs$p <-  
  2*(1-pnorm(abs(model.coefs$t.value)))

# Compute residuals
model.res <- resid(model)

# Inspect residuals
hist(model.res)
qqnorm(model.res, pch = 1, frame = FALSE)
qqline(model.res, col = "steelblue", lwd = 2)



#Inspect results
model.coefs




