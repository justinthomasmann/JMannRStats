setwd("C:/Users/Karin/Desktop/BIOSTATS")
library(nlme)
cod <- read.table("ParasiteCod.txt",header = T)
cod
cod <- na.omit(cod)
cod$fArea <- factor(cod$Area)
cod$fYear <- factor(cod$Year)
FM <- glm(Prevalence~fArea*fYear+Length, family = binomial, data = cod)
FM
summary(FM)
