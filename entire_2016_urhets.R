setwd("C:/Users/Karin/Desktop/Data")
library(vegan)
library(ggplot2)

#Full season het urs for 2016 without focal species changes at e sites
ent.urhets.16 <- read.csv("urs_diversity_2016.csv", h = T)
ent.urhets.16[is.na(ent.urhets.16)] <- 0
urhets.df <- ent.urhets.16[,4:15] 
urhets.df[is.na(urhets.df)] <- 0

div.hets.e <- diversity(urhets.df[ent.urhets.16$treat=="e"], index = "shannon")
