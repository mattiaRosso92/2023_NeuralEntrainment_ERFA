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
dataTempo <- read_excel('/Users/mattiaipem/Documents/IPEM/Projects/ Lou _ Entrainment/ERFA_Metronomes/Data/Experiment/ERFA_bhv_Tempo.xlsx')
# Convert variables to factors
dataTempo$Direction <- factor(dataTempo$Direction)
# Re-order factors, for interpretation
dataTempo$Direction <- relevel(dataTempo$Direction , 'Null')
contrasts(dataTempo$Direction) #set baseline for contrast across modalities
# Store number of steps
nsteps <- as.numeric(tail(names(dataTempo), n=1)); #convert last variable name to numeric, get last time value
stepn1 <- as.numeric(names(dataTempo))[3]; # get 1st time value (first index after factors)
# Compute downsampling-factor (if any)
dwn <- as.numeric(names(dataTempo)[4]) - as.numeric(names(dataTempo)[3])

# Remove outlier: Subject 6 (missing data) and eventual outliers
dataTempo <- subset(dataTempo, Subject != 6)
#dataTempo <- subset(dataTempo, Subject != 6 & Subject != 4 & Subject != 13 & Subject != 16 & Subject != 20)
 
# Convert to longitudinal data format: person, (other factors) , period
dataTempo = melt(dataTempo , id = c("Subject","Direction"))
# Rename serial position
dataTempo <- dataTempo %>% rename("Time" = "variable")
dataTempo <- dataTempo %>% rename("InstFrex" = "value")
# Check long format
summary(dataTempo)


#Create higher-order orthogonal polynomial
npols <- 2 #number of polynomials
t <- poly(seq(from = stepn1, to = nsteps, by = dwn),npols) #assign to 'time'; 'by' is downsampling factor, soft coded as interval between 2 datapoints
#Create time variable in data frame
dataTempo[,paste("ot",1:npols,sep="")] <-
  t[dataTempo$Time, 1:npols]


### Model fitting ###
#Fit quadratic model
model.quad <- lmer(InstFrex ~ (ot1+ot2)*Direction +
                   (ot1+ot2 | Subject) + #random effect of participant
                   (ot1+ot2 | Subject:Direction) , #...and interactions with factors
                 control = lmerControl(optimizer = "bobyqa") ,
                 data = dataTempo, REML=FALSE)
# final argument: if false, uses maximum likelihood estimation to fit the model 
# (vs restricted maximum likelihood estimation)



#Get parameter estimates for all models
model.coefs <- data.frame(coef(summary(model.quad)))

#...and estimate p-values
model.coefs$p <-  
  2*(1-pnorm(abs(model.coefs$t.value)))

# Compute residuals
model.res <- resid(model.quad)

# Inspect residuals
hist(model.res)
qqnorm(model.res, pch = 1, frame = FALSE)
qqline(model.res, col = "steelblue", lwd = 2)


# Visualize grand-average per conditions
ggplot(data = dataTempo , aes(x = Time, y = InstFrex, color = Direction)) +
  facet_grid(cols = vars(Direction)) +
  geom_point() +
  scale_color_brewer(palette="Set1") + #manually set color mapping
  scale_shape_manual(values=c(16, 1))+  #manually set shape mapping
  stat_summary(fun = mean , geom = "line", size = 2, color = 'black') + #mean timecourse
  stat_summary(fun.data = mean_se, geom = "pointrange", color = 'black', alpha = 1) + #st error time course
  xlab("Time (ms)") +
  ylab("Frequency (Hz)") # +
  #scale_x_discrete(breaks=c("1","32","64"),
       #            labels=c("0", expression(pi), expression("2"~pi~"  ")))

# Visualize fit per conditions
ggplot(data = dataTempo , aes(x = Time, y = fitted(model.quad), color = Direction)) +
  facet_grid(cols = vars(Direction)) +
  geom_point() +
  scale_color_brewer(palette="Set1") + #manually set color mapping
  scale_shape_manual(values=c(16, 1))+  #manually set shape mapping
  stat_summary(fun = mean , geom = "point", size = 2, color = 'black') + #mean timecourse
  stat_summary(fun.data = mean_se, geom = "pointrange", color = 'black', alpha = 1) + #st error time course
  xlab("Time (ms)") +
  ylab("Frequency (Hz)") # +
  #scale_x_discrete(breaks=c("1","32","64"),
  #            labels=c("0", expression(pi), expression("2"~pi~"  ")))




#Inspect results
model.coefs




# # Generate table (see 'ERPA_eeg_Tempo.R')


