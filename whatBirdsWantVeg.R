setwd("D:/data")

library(plyr) #summary statistics
library(ggplot2) ; theme_set(theme_classic())
library(lme4)  
library(sjPlot)
library(ggpubr)
library(glmmTMB)
library(MCMCglmm)
library(multcomp)
library(ggbiplot)
library(lmerTest)
library(ggeffects)
library(bayesplot)
library(rstanarm)
library(shinystan)
launch_shinystan_demo()
citation()

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt)
}
## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

wbw.df <- data.frame(read.csv("birdsEyeViewOnVeg.csv",h=T))

wbw.c.df <- data.frame(subset(wbw.df, wbw.df$treat != "e"))

colnames(wbw.c.df)

c.with1.df <- data.frame(subset(wbw.c.df,wbw.c.df$baww.cumm > 0))

cSite <- wbw.c.df$site
cStop <- wbw.c.df$stop
cYear <- wbw.c.df$year
cBawwCumm <- wbw.c.df$baww.cumm
cBawwB <- wbw.c.df$baww.b
cTotUrs <- wbw.c.df$tot.urs.c
cUvP <- wbw.c.df$uv*100
cGvP <- wbw.c.df$gv*100
cCvP <- wbw.c.df$cv*100
cCvbP <- wbw.c.df$cvb*100
cUvVar <- wbw.c.df$uvVar
cUvRange <- wbw.c.df$uvRange
cUvIqr75 <- wbw.c.df$iqr75
cUvIqr50 <- wbw.c.df$iqr50
cUvMax <- wbw.c.df$max
cUvMin <- wbw.c.df$min
cUvMed <- wbw.c.df$median

pairs(~cBawwCumm + cUvIqr75,
      diag.panel = panel.hist,
      lower.panel = panel.smooth,
      upper.panel = panel.cor)

pairs(~ c.with1.df$baww.cumm + c.with1.df$median,
      diag.panel = panel.hist,
      lower.panel = panel.smooth,
      upper.panel = panel.cor)



with0.cumm.mdl <- glmmTMB(wbw.c.df$baww.cumm~wbw.c.df$median, family = nbinom1)
summary(with0.cumm.mdl)



with0.cumm.75.mdl <- glmmTMB(wbw.c.df$baww.cumm~wbw.c.df$iqr75, family = nbinom1)
summary(with0.cumm.75.mdl)

anova(with0.cumm.75.mdl,with0.cumm.mdl)#median model fits significantly better


##########no experimental stop data: plot relationship between median uv and cummulative baww obs#########
ggplot(data = wbw.c.df, aes(x= wbw.c.df$median, y= wbw.c.df$baww.cumm, color = wbw.c.df$site))+
  geom_point(shape = 1, size = 2)+
  stat_smooth(method="glm", method.args=list(family="poisson"))+
  xlab("Median understory vegetation cover (%)")+
  ylab("Cummulative BAWW observations")+
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18, face = "bold"))+
  annotate("text", x = 60, y = 5, label = "Bayes R^2 = 0.34", size = 5)

stan.c.median.mdl <- stan_glm(cBawwCumm ~ cUvMed,
         algorithm = "sampling",
         data = wbw.c.df,
         family = poisson,
         prior = student_t(df = 7), 
         prior_intercept = student_t(df = 7))

tab_model(stan.c.median.mdl)

#######relationship between cummulative obs and median UV########
cumm.75.mdl <- glm(c.with1.df$baww.cumm~c.with1.df$iqr75, family = poisson)
summary(cumm.75.mdl)
rsq(cumm.75.mdl, adj = T)#adj R^2 = 0.18

cumm.mdl <- glm(c.with1.df$baww.cumm~c.with1.df$median, family = poisson)
summary(cumm.mdl)
rsq(cumm.mdl, adj = T)#adj R^2 = 0.26

anova(cumm.75.mdl,cumm.mdl)#median model fits data slightly better


ggplot(data = c.with1.df, aes(x= c.with1.df$median, y= c.with1.df$baww.cumm))+
  geom_point()+
  stat_smooth(method="glm", method.args=list(family="poisson"))+
  xlab("Median understory vegetation cover (%)")+
  ylab("Breeding-season pressence/absense")+
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18, face = "bold"))





#######relationship between breeding season obs and median UV########
bre.mdl <- glm(wbw.c.df$baww.b~wbw.c.df$median, family = binomial)
summary(bre.mdl)


#goodness of fit chi squared test indicates model fits the data well because the residual deviance is not significantly 
#different from the ideal "perfect fit" model
with(bre.mdl, cbind(res.deviance = deviance, df = df.residual,
               p = pchisq(deviance, df.residual, lower.tail=FALSE)))#p = 0.2

plot_model(bre.mdl,
           show.values = T)

ggplot(data = wbw.c.df, aes(x= wbw.c.df$median, y= wbw.c.df$baww.b))+
  geom_point(shape = 1, size = 2)+
  stat_smooth(method="glm", method.args=list(family="binomial"))+
  xlab("Median understory vegetation cover (%)")+
  ylab("Breeding-season BAWW pressence/absense")+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1.00), labels = c("0", "","","", "1"))+
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18, face = "bold"))+
  annotate("text", x = 15, y = .85, label = "Bayes R^2 = 0.22", size = 5)

stan.c.median.bre.mdl <- stan_glm(wbw.c.df$baww.b ~ cUvMed,
                              algorithm = "sampling",
                              data = wbw.c.df,
                              family = binomial,
                              prior = student_t(df = 7), 
                              prior_intercept = student_t(df = 7))

tab_model(stan.c.median.bre.mdl)


#####Only SGL non-experimental data####

sglVeg.df <- data.frame(read.csv("sglBawwVeg.csv", h=T))
sglVegYear.df <- subset(sglVeg.df, sglVeg.df$year != "2018")


#####BAWW cumm vs uv Median PP: One year of data from sgla and b with no experimental treatment#####
ggplot(data = sglVegYear.df, aes(x= sglVegYear.df$median, y= sglVegYear.df$baww.cumm, shape = sglVegYear.df$site))+
  geom_point(size = 3)+
  stat_smooth(method="glm", method.args=list(family="poisson"))+
  xlab("Median understory vegetation cover (%)")+
  ylab("Cummulative BAWW observations")+
  scale_shape_manual(values = c(8,6),
                     name = "SGL Survey Routes",
                     breaks = c("sgla","sglb"),
                     labels = c("A", "B"))+
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = c(0.15,0.9),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black"))

sglVeg.mdl <- stan_glm(sglVeg.df$baww.cumm~log(sglVeg.df$median), 
         data = sglVeg.df,
         algorithm = "sampling",
         family = poisson,
         prior = student_t(df = 7), 
         prior_intercept = student_t(df = 7),
         chains = 10,
         seed = 12345)

tab_model(sglVeg.mdl,transform = NULL,
          show.intercept = FALSE)


hist(log(sglVeg.df$median))

ggplot(data = sglVegYear.df, aes(x= sglVegYear.df$median, y= sglVegYear.df$baww.b))+
  geom_point(shape = 1, size = 2)+
  stat_smooth(method="glm", method.args=list(family="binomial"))+
  xlab("Median understory vegetation cover (%)")+
  ylab("Breeding season BAWW presence/absence")+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1.00), labels = c("0", "","","", "1"))+
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18, face = "bold"))

range(log(sglbVeg.df$median))
log(sglbVeg.df$median)
#####Just sglb plots#####
sglbVeg.df <- subset(sglVegYear.df, sglVegYear.df$site == "sglb")


ggplot(data = sglbVeg.df, aes(x= log(sglbVeg.df$median), y= sglbVeg.df$baww.cumm))+
  geom_point(shape = 1, size = 3)+
  stat_smooth(method="glm", method.args=list(family="poisson"))+
  xlab("Log-transformed median understory vegetation cover")+
  ylab("Cummulative BAWW observations")+
  scale_x_continuous(breaks = seq(from = 2.4, to = 4.6, by = 0.2))+
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16, face = "bold"))+
  annotate("text", x = 3.5, y = 6, label = "Bayes R^2 = 0.62", size = 6)

sglbVeg.mdl <- stan_glm(sglbVeg.df$baww.cumm~log(sglbVeg.df$median), 
                       data = sglbVeg.df,
                       algorithm = "sampling",
                       family = poisson,
                       prior = student_t(df = 7), 
                       prior_intercept = student_t(df = 7),
                       chains = 10,
                       seed = 12345)

tab_model(sglbVeg.mdl,transform = NULL,
          show.intercept = FALSE)
posterior_interval(sglbVeg.mdl)
summary(sglbVeg.mdl)
plot(sglbVeg.mdl)

ggplot(data = sglbVeg.df, aes(x= sglbVeg.df$median, y= sglbVeg.df$baww.b))+
  geom_point(shape = 1, size = 2)+
  stat_smooth(method="glm", method.args=list(family="binomial"))+
  xlab("Median understory vegetation cover (%)")+
  ylab("Breeding season BAWW presence/absence")+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1.00), labels = c("0", "","","", "1"))+
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18, face = "bold"))

sglb.bre.Veg.mdl <- stan_glm(sglbVeg.df$baww.b~log(sglbVeg.df$median), 
                        data = sglbVeg.df,
                        algorithm = "sampling",
                        family = binomial(link = "logit"),
                        prior = student_t(df = 7), 
                        prior_intercept = student_t(df = 7),
                        chains = 10,
                        seed = 12345)

tab_model(sglb.bre.Veg.mdl,transform = NULL,
          show.intercept = FALSE)


##not a strong relationship between uv median and baww.cumm @ sgla
sglaVeg.df <- subset(sglVegYear.df, sglVegYear.df$site == "sgla")
sglaVeg.mdl <- stan_glm(sglaVeg.df$baww.cumm~log(sglaVeg.df$median), 
                        data = sglaVeg.df,
                        algorithm = "sampling",
                        family = poisson,
                        prior = student_t(df = 7), 
                        prior_intercept = student_t(df = 7),
                        chains = 10,
                        seed = 12345)

tab_model(sglaVeg.mdl,transform = NULL,
          show.intercept = FALSE)

#####More stuff in the works#####
pairs(~ cBawwCumm + cBawwB + cTotUrs + cGvP + cCvbP + cUvRange + cUvIqr75 + cUvMax + cUvMin + cUvMed + pc1, 
            diag.panel = panel.hist,
            lower.panel = panel.smooth,
            upper.panel = panel.cor)

veg.mdl1 <- stan_glm(cBawwCumm ~ cUvP + cGvP + cCvbP + cUvRange + cUvIqr75,
         data = wbw.c.df,
         family = poisson)
summary(veg.mdl1)

veg.mdl2 <- stan_glm(cBawwCumm ~ cGvP + cUvIqr75,
                     data = wbw.c.df)
summary(veg.mdl2)

plot_model(veg.mdl2)

launch_shinystan(veg.mdl2)

veg.mdl.b <- stan_glmer(cBawwB ~ cUvMed + cGvP + (1|cSite),
                      data = wbw.c.df,
                      family = binomial(link = "logit"))
summary(veg.mdl.b)
plot_model(veg.mdl.b)

freq.veg.mdl.b <- glm(cBawwB ~ cUvP + cGvP + cCvbP + cUvRange + cUvIqr75,
                      data = wbw.c.df,
                      family = binomial(link = "logit"))
summary(freq.veg.mdl.b)

lmerTest::step(freq.veg.mdl.b)

freq.veg.mdl.BF.b <- glm(cBawwB ~ cGvP + cUvIqr75,
                         data = wbw.c.df,
                         family = binomial(link = "logit"))
summary(freq.veg.mdl.BF.b)


freq.veg.mdl1 <- glm(cBawwCumm ~ cUvP + cGvP + cCvbP + cUvRange + cUvIqr75,
    data = wbw.c.df,
    family = poisson)

lmerTest::step(freq.veg.mdl1)


##########PCA for veg###########
wbw.c.df$gvT <- wbw.c.df$gv*100
wbw.c.df$bawwCummT <- wbw.c.df$baww.cumm*10
wbw.c.df$cvbT <- wbw.c.df$cvb*100

vegPca <- prcomp(wbw.c.df[,c(14,18,19)])
summary(vegPca)
vegPca
pc1 <- vegPca$x[,1]
#wbw.c.df$pc1T <- pc1*(-1)#Reverse signs so that values increase with increasing vegetation metrics
#pc1T <- wbw.c.df$pc1T


ggbiplot(vegPca, ellipse = T, groups = wbw.c.df$site, labels = wbw.c.df$stop)

totUrs.freq.mdl <- glm(cTotUrs ~ pc1-1,
                         family = poisson)
summary(totUrs.freq.mdl)
plot_model(totUrs.freq.mdl)

totUrs.freq.mdl1 <- glm(cTotUrs ~ cUvMed + cGvP -1,
                        family = poisson)
summary(totUrs.freq.mdl1)
plot_model(totUrs.freq.mdl1,
           show.values = T)



 

freq.veg.BF.mdl1 <- glm(cBawwCumm ~ cGvP + cUvIqr75,
                        data = wbw.c.df,
                        family = poisson)
summary(freq.veg.BF.mdl1)
plot_model(freq.veg.BF.mdl1)

freq.veg.pca.mdl <- glm(cBawwCumm ~ pc1T,
                        data = wbw.c.df,
                        family = poisson)
summary(freq.veg.pca.mdl)

anova(freq.veg.pca.mdl,freq.veg.BF.mdl1)

bay.pca.mdl <- stan_glmer(cBawwCumm ~ pc1T + (1|cSite),
                          data = wbw.c.df,
                          family = poisson)
summary(bay.pca.mdl)
plot_model(bay.pca.mdl)
plot_model(bay.pca.mdl,
           type = "re")

wbw.noSglB.df <- data.frame(subset(wbw.c.df, wbw.c.df$site!="sglb"))

noSglB.stan.mdl <- stan_glm(wbw.noSglB.df$baww.cumm ~ wbw.noSglB.df$gvT + wbw.noSglB.df$uvRange + wbw.noSglB.df$iqr75,
                              data = wbw.noSglB.df,
                              family = poisson)
plot_model(noSglB.stan.mdl,
           transform = NULL)

noSglB.freq.mdl <- glm(wbw.noSglB.df$baww.cumm ~ wbw.noSglB.df$gvT + wbw.noSglB.df$iqr75,
                            data = wbw.noSglB.df,
                            family = poisson)
summary(noSglB.freq.mdl)

pairs(~cBawwCumm + pc1T,
      diag.panel = panel.hist,
      lower.panel = panel.smooth,
      upper.panel = panel.cor)


glm(cBawwCumm ~ cGvP, family = poisson)
ggplot(data = wbw.c.df, aes(x=cGvP, y = cBawwCumm))+
  geom_point()+
  geom_abline(intercept = -0.27,slope = 0.027)


glm(cBawwCumm ~ cUvIqr75, family = poisson)
ggplot(data = wbw.c.df, aes(x=cUvIqr75, y = cBawwCumm))+
  geom_point()+
  geom_abline(intercept = -0.72,slope = 0.022)

glm(cBawwCumm ~ cUvMed, family = poisson)
