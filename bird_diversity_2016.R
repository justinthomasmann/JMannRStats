setwd("C:/Users/Karin/Desktop/Data")
library(vegan)
library(ggplot2)

bdiv <- read.csv("bird_diversity_2016.csv", h=T)
bdiv

colnames(bdiv)
bdiv[is.na(bdiv)] <- 0

tot.div <- diversity(bdiv[1:48,4:71], index = "shannon")
tot.div

barplot(tot.div)

plot(bdiv$site, tot.div,
     ylim = c(1.5,3.5),
     col = c("orange","green", "blue", "yellow"),
     main="Shannon Diversity by Site", 
     xlab = "Site", 
     ylab = "H'")

bdiv.no.sglb <- data.frame(subset(bdiv, bdiv$site != "sglb"))
bdiv.no.sglb$site <- factor(bdiv.no.sglb$site)
bdiv.no.sglb
bdiv.no.sglb$site


plot(tot.div[bdiv$site != "sglb"],bdiv.no.sglb$uv15,
     pch = 19,
     main = "Shannon Diversity ~ Understory Vegetation Density",
     xlab = "H'",
     ylab = "Understory vegetation density (%)",
     col = ifelse(bdiv.no.sglb$site== "sgl", "blue", 
                  ifelse(bdiv.no.sglb$site=="fss", "green", "orange")))
     abline(lm(bdiv.no.sglb$uv15~tot.div[bdiv$site != "sglb"]))
        legend("topleft", c("sgl", "fss", "bunp"), 
        col = c("blue", "green", "orange"), pch = 19)
        text(2,75, labels = "adjusted r^2 = 0.3")

     
summary(lm(bdiv.no.sglb$uv15~tot.div[bdiv$site != "sglb"]))

#diversity ~ treatment
e.div <- tot.div[bdiv.no.sglb$treat=="e"]
e.div
mean(e.div)

c.div <- tot.div[bdiv.no.sglb$treat=="c"]
c.div
mean(c.div)

boxplot(e.div,c.div,
        names = c("Experimental", "Control"),
        main = "Comparison of Shannon diversity by site treatment type",
        xlab = "Treatment",
        ylab = "H'",
        col = c("magenta", "grey"))

var.test(e.div,c.div)
t.test(e.div,c.div)

summary(lm(tot.div[bdiv$site != "sglb"]~bdiv.no.sglb$treat== "e"))








sgl.a.div <- diversity(bdiv[1:12,4:71], index = "shannon")
sgl.a.div

sgl.a.factor <- rep("sgl.a", length(sgl.a.div))
sgl.a.df <- data.frame(sgl.a.factor, sgl.a.div)
## sgla_df$sgla_factor <- factor(sgla_df$sgla_factor)
plot(sgl.a.df)

sgl.b.div <- diversity(bdiv[13:22,4:71], index = "shannon")
sgl.b.div

sgl.b.factor <- rep("sgl.b", length(sgl.b.div))
sgl.b.df <- data.frame(sgl.b.factor, sgl.b.div)
plot(sgl.b.df)


fss.div <- diversity(bdiv[23:34,4:71], index = "shannon")
fss.div

fss.factor <- rep("fss", length(fss.div))
fss.df <- data.frame(fss.factor, fss.div)
plot(fss.df)

bunp.div <- diversity(bdiv[35:48,4:71], index = "shannon")
bunp.div

bunp.factor <- rep("bunp", length(bunp.div))
bunp.df <- data.frame(bunp.factor, bunp.div)
plot(bunp.df)



mean(sgla_div)
mean(sglb_div)
mean(fss_div)
mean(bunp_div)

var(sgla_div)
var(sglb_div)
var(fss_div)
var(bunp_div)

div_nosglb <- data.frame(subset(bdiv, bdiv$site != "sglb"))
uv_nosglb <- data.frame(subset(bdiv, bdiv))
total_div <- c(tot_div)
total_div

lm(bdiv$uv15~total_div)
uv_div <- lm(bdiv$uv15~total_div)

x1 <- factor(bdiv$site)
x1

bdiv$site
anova()

t.test(fss_div,bunp_div)
t.test(fss_div, sgla_div)

bdiv$uv15
tot_div
