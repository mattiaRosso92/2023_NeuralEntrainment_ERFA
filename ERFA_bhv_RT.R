library(readxl)

data <- read_excel('/Users/mattiaipem/Documents/IPEM/Projects/ Lou _ Entrainment/ERFA_Metronomes/Data/Experiment/ERFA_rt_Phase.xlsx')

# Remove outliers
data <- subset(data, Subject != 14)

model <- lm(data = data , RT ~ Direction)
summary(model)

hist(residuals(model))
