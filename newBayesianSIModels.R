setwd("D:/data")

library(plyr) #summary statistics
library(dplyr)
library(ggplot2) ; theme_set(theme_classic())
library(lme4) #frequentist models  
library(sjPlot) #plot_model & tab_model (blue = #377EB8)
library(ggpubr) #ggarrange
library(glmmTMB) #zi & overdispersed models
library(ggbiplot) #pca plots
library(lmerTest) #lmerTest::step
library(ggeffects) #ggpredict
library(rstanarm) #stan models
library(shinystan) #stan model evaluation
library(DescTools) #Dunn's rank sum test
library(loo) #use loo() to compare fits between bayesian models
library(BayesianFirstAid) #bayes.prop.test for comparison of response to treatment between sites
launch_shinystan_demo()


####FULL DATAFRAME####

#The following csv has a dummy level for site (adummy) and treatment (a). This avoids R using 
#BUNP and control stops as reference levels, resulting in BUNP and control stops estimates being 
#included in all model outputs. For all vegetation metrics, the dummy level is set to 0. 

# year
# x=dummy
# y=2016
# z=2017

df <- data.frame(read.csv("soc_info_expt_complete_NoNAs.csv", h=T))
colnames(df)

levels(df$site)

df$uvRangeT <- df$uvRange*0.1
df$uvVarT <- df$uvVar*0.01
df$uvT <- df$uv*10
df$gvT <- df$gv*10
df$cvT <- df$cv*10
df$uvM <- df$uvMetric*0.1
df$iqr75T <- df$iqr75*0.1
df$medianT <- df$median*0.1

#####FUll DATAFRAME VARIABLES#####

site <- as.factor(df$site)
stop <- as.factor(df$stop)
year <- as.factor(df$year)
treat <- as.factor(df$treat)
priorFspT <- as.factor(df$prior.fsp.t)
priorFspB <- as.factor(df$prior.fsp.b)
bawwCumm <- df$baww.cumm
fsiT <- as.factor(df$fsi.t)
addit <- as.factor(df$add)
fsiP <- df$fsi.p
fsiB <- as.factor(df$fsi.b)
uvR <- df$uvRange
uvV <- df$uvVar
gv <- df$gv
cvb <- df$cvb
uvT <- df$uvT
gvT <- df$gvT
cvT <- df$cvT
uvMed <- df$median
uvVarT <- df$uvVarT
uvRangeT <- df$uvRangeT
uvMedT <- df$medianT
uv75T <- df$iqr75T

summary(df)

#response at FSS
length(df$fsi.t[df$site=="fss"])
length(df$fsi.t[df$fsi.t == 1 & df$site == "fss" & df$treat=="e"])#10
length(df$fsi.t[df$fsi.t == 0 & df$site == "fss" & df$treat=="e"])#3
10/13 #77%
#increase at control
length(df$fsi.t[df$fsi.t == 1 & df$site == "fss" & df$treat=="c"])#3 
length(df$fsi.t[df$site == "fss" & df$treat=="c"])#out of 11
3/11 #27%

#response at SGL
length(df$fsi.t[df$site=="sgla"])
length(df$fsi.t[df$fsi.t == 1 & df$site == "sgla" & df$treat=="e"])#3
length(df$fsi.t[df$fsi.t == 0 & df$site == "sgla" & df$treat=="e"])#3
#50%
#increase at control
length(df$fsi.t[df$fsi.t == 1 & df$site == "sgla" & df$treat=="c"])#0
length(df$fsi.t[df$site == "sgla" & df$treat=="c"])#out of 6 c
#0%

#response at BUNP
length(df$fsi.t[df$site=="bunp"])
length(df$fsi.t[df$fsi.t == 1 & df$site == "bunp" & df$treat=="e"])#4
length(df$fsi.t[df$fsi.t == 0 & df$site == "bunp" & df$treat=="e"])#9
4/13 #31%
#increase at control
length(df$fsi.t[df$fsi.t == 1 & df$site == "bunp" & df$treat=="c"])#0
length(df$fsi.t[df$site == "bunp" & df$treat=="c"])#out of 11 c 


#about column "add"
#if add = 3, prior baww = 1 & fsi = 0
#if add = 2, prior baww = 1 & fsi = 1
#if add = 1, prior baww = 0 & fsi = 1
#if add = 0, prior baww = 0 & fsi = 0 
df$add[df$prior.fsp.t == "2" & df$fsi.t == "0"] <- 0
levels(addit)

length(which(df$add == "2"))#only 4 stops had prior & fsi

length(which(df$add == "2" & df$treat == "e"))#two were experimental @ sgl 
add2 <- rbind(df[which(df$add =="2"),])#sgl 4e,12e & fss 3c & 6c
length(which(df$add == "1"))
length(which(df$add == "3"))

#Baseline BAWW FSS | SGL
fs.base <- c(2,11)
fs.bn <- c(12,12)
bayes.prop.test(fs.base,fs.bn)
plot(bayes.prop.test(fs.base,fs.bn))

#FSS | SGL Bayesian proportion test
fs.fsi <- c(10,3)
fs.n <- c(13, 6)
bayes.prop.test(fs.fsi,fs.n)
plot(bayes.prop.test(fs.fsi,fs.n))


#FSS | BUNP Bayesian proportion test
length(fsiT[site=="bunp" & treat=="e"])
fb.fsi <- c(10,4)
fb.n <- c(13,13)
bayes.prop.test(fb.fsi,fb.n)
plot(bayes.prop.test(fb.fsi,fb.n))


#SGL | BUNP Bayestian proportion test 
sb.fsi <- c(3,4)
sb.n <- c(6,13)
bayes.prop.test(sb.fsi,sb.n)
plot(bayes.prop.test(sb.fsi,sb.n))
citation("BayesianFirstAid")

?prop.test()




#####NO DUMMY LEVEL DATAFRAME#####
nodum.df <- df[1:60,]
length(nodum.df$treat[nodum.df$treat == "e"])#32 e across both years
length(nodum.df$treat[nodum.df$treat == "c"])#28 c "    "    "

levels(nodum.df$treat)

#####no dummy variables#####

ndSite <- as.factor(nodum.df$site)
ndStop <- as.factor(nodum.df$stop)
ndYear <- as.factor(nodum.df$year)
ndTreat <- as.factor(nodum.df$treat)
ndBawwCumm <- nodum.df$baww.cumm
ndFsiT <- as.factor(nodum.df$fsi.t)
ndFsiP <- as.factor(nodum.df$fsi.p)
ndFsiB <- as.factor(nodum.df$fsi.b)
ndUvR <- nodum.df$uvRange
ndUvV <- nodum.df$uvVar
ndGv <- nodum.df$gv
ndCvb <- nodum.df$cvb
ndUvT <- nodum.df$uvT
ndGvT <- nodum.df$gvT
ndCvT <- nodum.df$cvT
ndUvVarT <- nodum.df$uvVarT
ndUvRangeT <- nodum.df$uvRangeT
ndUvMedT <- nodum.df$medianT
ndUv75T <- nodum.df$iqr75T



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#####models with dummy level#####

#full-season focal species increase: 
#year is n.s. 
#interactions between treat, uvMedT, and site are all n.s.
mdl.fsiT <- stan_glm(fsiT  ~ site + treat + uvMedT, 
                          data = df,
                          algorithm = "sampling",
                          family = binomial(link = "logit"),
                          prior = student_t(df = 7), 
                          prior_intercept = student_t(df = 7),
                     chains = 10,
                     seed = 12345)
summary(mdl.fsiT)

#the mean of fsiT matches the mean_PPD
length(fsiT[df$fsi.t=="1"])/length(fsiT)#rough fit diagnostic is good
mdl.fsiT.rsq <- bayes_R2(mdl.fsiT)
print(median(mdl.fsiT.rsq))

plot_model(mdl.fsiT,
           transform = NULL,
           colors = "#377EB8",
           bpe = "median",
           vline.color = "black",
           prob.inner = .50,
           prob.outer = .95,
           show.values = TRUE, 
           value.offset = .3,
           line.size = 1.5)

tab_model(mdl.fsiT,
          transform = NULL,
          show.intercept = FALSE)
posterior_interval(mdl.fsiT, prob = .95)

#Get 500 draws from the posterior predictive distribution
?posterior_predict()
pred.full <- posterior_predict(mdl.fsiT, draws = 500)
pred.full

dim(pred.full)

#Calculate the proportion of those 500 draws in which the model predicts 
#a focal species increase
df$fullPostPred <- colMeans(pred.full)

#Calculate the mean and 95% CI for predicted focal species increase at c & e stops
fullMeanC <- MeanCI(df$fullPostPred[df$treat=="c"])
fullMeanE <- MeanCI(df$fullPostPred[df$treat=="e"])
fullMeanC#10% chance of increase at c
fullMeanE#51% chance of increase at e 
#Dataframe to plot the treatment means with 95% CI
fullPlot.df <- data.frame(rbind(fullMeanC,fullMeanE))
fullPlot.df
fullPlot.df$treat <- c("c","e")

propFullFsiPEB <- ggplot(fullPlot.df, aes(x=treat, y = mean))+
  geom_point(size = 2)+
  geom_errorbar(fullPlot.df, mapping=aes(x=treat, ymin = lwr.ci, ymax=upr.ci), width=0.05, size=1)+
  labs(title = "B",
       x = "",
       y = "Probability of BAWW increase")+
  scale_x_discrete(breaks = c("c","e"),
                   labels = c("Control", "Playback"))+
  scale_y_continuous(limits = c(0.0,1.0))+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.margin = unit(c(3,3,3,3), "lines"))
propFullFsiPEB  

propFullFsiPEB + geom_point(mapping = aes(x=newdat$treat, y=newdat$fsi.t), 
                            shape = 1, size = 2, 
                            position = position_jitter(width = 0.05, 
                                                       height = 0.05))


####fsiTByTreatmentPP with PPD probabilities####
newdat <- data.frame(read.csv("soc_info_expt_complete_NoDummy.csv", h=T))

ggplot()+
  geom_point(data = newdat, mapping = aes(x=treat, y=fsi.t),
             shape = 1, size = 3, 
             position = position_jitter(width = 0.05, height = 0.05))+
  geom_pointrange(data = fullPlot.df, 
                  mapping = aes(x = treat, y=mean, ymin = lwr.ci, ymax = upr.ci),
                  size = 1.3,
                  shape = 20)+
  scale_x_discrete(breaks = c("c", "e"), labels = c("Control","Playback"))+
  xlab("")+
  scale_y_continuous(name = expression("Probability of BAWW increase"), 
                     sec.axis = sec_axis(~ . * 1, name = "Observed increase",
                                         breaks = c(0.00,0.25,0.50,0.75,1.00),
                                         labels = c("No", "","","","Yes")))+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        text = element_text(family = "Arial"))



#####cummulative baww counts model#####
#year is n.s.

mdl.cumm <- stan_glm(bawwCumm ~ site + treat + uvMedT, 
                     data = df,
                     algorithm = "sampling",
                     family = poisson,
                     prior = student_t(df = 7), 
                     prior_intercept = student_t(df = 7))

plot_model(mdl.cumm,
           transform = NULL,
          # colors = "#377EB8",
          bpe = "median",
          vline.color = "black",
          prob.inner = .50,
          prob.outer = .95,
          show.values = TRUE, 
          value.offset = .3,
          line.size = 1.5)

tab_model(mdl.cumm,
          transform = NULL,
          show.intercept = FALSE)


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#####year has no predictive value in the model. I also tested interaction between year*treat
mdl.pre.year <- stan_glm(fsiP  ~ site + treat + uvMedT + year, 
                    data = df,
                    algorithm = "sampling",
                    family = binomial(link = "logit"),
                    prior = student_t(df = 7), 
                    prior_intercept = student_t(df = 7),
                    chains = 10,
                    seed = 12345)

summary(mdl.pre.year)
posterior_interval(mdl.pre.year)
plot(mdl.pre.year, prob = 0.8)

####prebreeding season model####
#Student t distribution with 7 df 
#This model includes dummy levels for site and treatment.
#At the dummy levels, fsi and vegetation metrics are all zero
#$...i.e. nothing happened at the reference site and treatment
mdl.pre <- stan_glm(fsiP  ~ site + treat + uvMedT, 
                    data = df,
                    algorithm = "sampling",
                    family = binomial(link = "logit"),
                    prior = student_t(df = 7), 
                    prior_intercept = student_t(df = 7),
                    chains = 10,
                    seed = 12345)

summary(mdl.pre)
posterior_interval(mdl.pre)
plot(mdl.pre, prob = 0.8)
posterior_vs_prior(mdl.pre, prob = 0.5)

launch_shinystan(mdl.pre)
loo.pre <- loo(mdl.pre)

loo_compare(mdl.pre, mdl.pre.norm)

#Using default normal distribution argument...no appreciable difference...see loo_compare
mdl.pre.norm <- stan_glm(fsiP  ~ site + treat + uvMedT, 
                         data = df,
                         algorithm = "sampling",
                         family = binomial(link = "logit"),
                         prior = normal(), 
                         prior_intercept = normal(),
                         chains = 10,
                         seed = 12345)

plot(mdl.pre.norm, prob = 0.8)
posterior_vs_prior(mdl.pre.norm, prob = 0.5)

loo.pre.norm <- loo(mdl.pre.norm)
loo_compare(loo.pre, loo.pre.norm)

#the no-dummy model...no directional differences...fss becomes stronger, but I'm more comfortable 
#interpreting the dummy model...without the dummy level, R uses control stops at BUNP
#as the reference level, which allow interpretation of effects at BUNP
mdl.pre.nd <- stan_glm(ndFsiP  ~ ndSite + ndTreat + ndUvMedT, 
                    data = df,
                    algorithm = "sampling",
                    family = binomial(link = "logit"),
                    prior = student_t(df = 7), 
                    prior_intercept = student_t(df = 7),
                    chains = 10,
                    seed = 12345)

plot(mdl.pre.nd, prob = 0.8)


####prePM####
prePM <-plot_model(mdl.pre,
           transform = NULL,
           #colors = "#377EB8",
           bpe = "median",
           vline.color = "black",
           prob.inner = .50,
           prob.outer = .95,
           show.values = TRUE, 
           value.offset = .3,
           line.size = 1,
           title = "A", 
           #axis.title = "",
           axis.labels = c("Median understory vegetation cover", "Playback", "Control", "SGL","FSS","BUNP"))+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.margin = unit(c(3,3,3,3), "lines"))

prePM

#Get 500 draws from the posterior predictive distribution
?posterior_predict()
pred.pre <- posterior_predict(mdl.pre, draws = 500)
pred.pre

dim(pred.pre)

#Calculate the proportion of those 500 draws in which the model predicts 
#a focal species increase
df$PrePostPred <- colMeans(pred.pre)

#Calculate the mean and 95% CI for predicted focal species increase at c & e stops
preMeanC <- MeanCI(df$PrePostPred[df$treat=="c"])
preMeanE <- MeanCI(df$PrePostPred[df$treat=="e"])
preMeanC
preMeanE
#Dataframe to plot the treatment means with 95% CI
prePlot.df <- data.frame(rbind(preMeanC,preMeanE))
prePlot.df$treat <- c("c","e")

propPreFsiPEB <- ggplot(prePlot.df, aes(x=treat, y = mean))+
  geom_point(size = 2)+
  geom_errorbar(prePlot.df, mapping=aes(x=treat, ymin = lwr.ci, ymax=upr.ci), width=0.05, size=1)+
  labs(title = "B",
       x = "",
       y = "Probability of BAWW increase")+
  scale_x_discrete(breaks = c("c","e"),
                   labels = c("Control", "Playback"))+
  scale_y_continuous(limits = c(0.0,0.6))+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.margin = unit(c(3,3,3,3), "lines"))
propPreFsiPEB  

preFsiPlotWrap <- ggarrange(prePM,propPreFsiPEB, widths = c(1.5,1))
preFsiPlotWrap

preAnnotated <- annotate_figure(preFsiPlotWrap, top = text_grob("Pre-breeding season", size = 20, face = "bold"))


#In this wrap... 
#Left-panel inner bars are 50% Uncertainty Intervals, i.e. High Density Intervals (HDI) 
#and outer bars are 95% HDI.
#Right-panel represents the predicted probability of BAWW increase at c and e stops
#based on 500 draws from the posterior predictive distribution. 
#Thus, simulating the outcome at each stop 500 times, I calculated the proportion of simulations 
#in which the model predicted BAWW increase and averaged those proportions across treatment types.
#Error bars indicate the 95% CI around that mean.
?plot_model() #see arguments for prob.inner & prob.outer
?rstanarm::posterior_predict()



#breeding season model
mdl.bre <- stan_glm(fsiB  ~ site + treat + uvMedT, 
                    data = df,
                    algorithm = "sampling",
                    family = binomial(link = "logit"),
                    prior = student_t(df = 7), 
                    prior_intercept = student_t(df = 7),
                    chains = 10,
                    seed = 12345)
posterior_interval(mdl.bre)
summary(mdl.bre)
launch_shinystan(mdl.bre)

brePM <- plot_model(mdl.bre,
           transform = NULL,
           #colors = "#377EB8",
           bpe = "median",
           vline.color = "black",
           prob.inner = .50,
           prob.outer = .95,
           show.values = TRUE, 
           value.offset = .3,
           line.size = 1,
           title = "C", 
           axis.labels = c("Median understory vegetation cover", "Playback", "Control", "SGL","FSS","BUNP"))+  
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.margin = unit(c(3,3,3,3), "lines"))

brePM

tab_model(mdl.bre,
          transform = NULL,
          show.intercept = FALSE)

pred.bre <- posterior_predict(mdl.bre, draws = 500)
pred.bre

dim(pred.bre)

#Calculate the proportion of those 500 draws in which the model predicts 
#a focal species increase
df$BrePostPred <- colMeans(pred.bre)

#Calculate the mean and 95% CI for predicted focal species increase at c & e stops
breMeanC <- MeanCI(df$BrePostPred[df$treat=="c"])
breMeanE <- MeanCI(df$BrePostPred[df$treat=="e"])
breMeanC
breMeanE
#Dataframe to plot the treatment means with 95% CI
brePlot.df <- data.frame(rbind(breMeanC,breMeanE))
brePlot.df$treat <- c("c","e")

propBreFsiPEB <- ggplot(brePlot.df, aes(x=treat, y = mean))+
  geom_point(size = 2)+
  geom_errorbar(brePlot.df, mapping=aes(x=treat, ymin = lwr.ci, ymax=upr.ci), width=0.05, size=1)+
  labs(title = "D",
       x = "",
       y = "Probability of BAWW increase")+
  scale_x_discrete(breaks = c("c","e"),
                   labels = c("Control", "Playback"))+
  scale_y_continuous(limits = c(0.0,0.6))+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.margin = unit(c(3,3,3,3), "lines"))


breFsiPlotWrap <- ggarrange(brePM,propBreFsiPEB, widths = c(1.5,1))
breFsiPlotWrap 

breAnnotated <- annotate_figure(breFsiPlotWrap, top = text_grob("Breeding season", size = 20, face = "bold"))



####PreBrePlots####
ggarrange(preAnnotated,breAnnotated, ncol = 1)








#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#Prior BAWW models
#If there was already a BAWW at a stop at during the previous year, 
#how does that affect probability of FSI?

#full season
mdl.priorT <- stan_glm(fsiT  ~ site + treat + priorFspT + uvMedT, 
         data = df,
         algorithm = "sampling",
         family = binomial(link = "logit"),
         prior = student_t(df = 7), 
         prior_intercept = student_t(df = 7))

plot_model(mdl.priorT,
           transform = NULL,
           #colors = "#377EB8",
           bpe = "median",
           vline.color = "black",
           prob.inner = .75,
           prob.outer = .95,
           show.values = TRUE, 
           value.offset = .3,
           line.size = 1,
           title = "Pre-breeding Season", 
           axis.title = ""#,
           #axis.labels = c("Median understory vegetation cover", "Experimental", "Control", "SGL","FSS","BUNP"))
           )+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16))


#breeding season prior vs. fsiT

mdl.priorB <- stan_glm(fsiT  ~ site + treat + priorFspB + uvMedT, 
                       data = df,
                       algorithm = "sampling",
                       family = binomial(link = "logit"),
                       prior = student_t(df = 7), 
                       prior_intercept = student_t(df = 7))

plot_model(mdl.priorB,
           transform = NULL,
           #colors = "#377EB8",
           bpe = "median",
           vline.color = "black",
           prob.inner = .75,
           prob.outer = .95,
           show.values = TRUE, 
           value.offset = .3,
           line.size = 1,
           title = "Pre-breeding Season", 
           axis.title = ""#,
           #axis.labels = c("Median understory vegetation cover", "Experimental", "Control", "SGL","FSS","BUNP"))
)+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16))

df$ifprior <- df$prior.fsp.t - df$fsi.t
ifPrior <- df$ifprior

add.df <- data.frame(cbind(df$site,df$stop,df$treat,df$prior.fsp.t,df$fsi.t,df$ifprior))




#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





#####INDIVIDUAL SITE MODELS#####


#####fss model#####
fss.df <- data.frame(subset(df, df$site == "fss"))

mdl.fss.fsiT <- stan_glm(fss.df$fsi.t  ~ fss.df$treat + fss.df$medianT, 
                         data = fss.df,
                         algorithm = "sampling",
                         family = binomial(link = "logit"),
                         prior = student_t(df = 7), 
                         prior_intercept = student_t(df = 7),
                         chains = 10,
                         seed = 12345)
summary(mdl.fss.fsiT)
launch_shinystan(mdl.fss.fsiT)

fssPM <- plot_model(mdl.fss.fsiT,
           transform = NULL,
           colors = "#377EB8",
           bpe = "median",
           vline.color = "black",
           prob.inner = .50,
           prob.outer = .95,
           show.values = TRUE, 
           value.offset = .3,
           line.size = 1,
           title = "C",
           axis.labels = c("Median understory vegetation cover","Playback" ))+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.margin = unit(c(3,3,3,3), "lines"))   
fssPM

pred.fss <- posterior_predict(mdl.fss.fsiT, draws = 500)
pred.fss

mean(pred.fss[,13])
mean(pred.fss[,1])

fss.flip <- t(pred.fss)
fss.flip

fss.pred.df <- as.data.frame(fss.flip)
fss.pred.df$treat <- fss.df$treat
fss.pred.df$treat
fss.pred.df$fsi.t <- fss.df$fsi.t

length(fss.pred.df[1,])

dim(pred.fss)

#Calculate the proportion of those 500 draws in which the model predicts 
#a focal species increase
fss.df$PostPred <- colMeans(pred.fss)

#Calculate the mean and 95% CI for predicted focal species increase at c & e stops
fssMeanC <- MeanCI(fss.df$PostPred[fss.df$treat=="c"])
fssMeanE <- MeanCI(fss.df$PostPred[fss.df$treat=="e"])
fssMeanC
fssMeanE

length(fss.df$PostPred[fss.df$treat=="e"])

#Dataframe to plot the treatment means with 95% CI
fssPlot.df <- data.frame(rbind(fssMeanC,fssMeanE))
fssPlot.df$treat <- c("c","e")

propFssFsiPEB <- ggplot(fssPlot.df, aes(x=treat, y = mean))+
  geom_point(size = 2)+
  geom_errorbar(fssPlot.df, mapping=aes(x=treat, ymin = lwr.ci, ymax=upr.ci), width=0.05, size=1)+
  labs(title = "D",
       x = "",
       y = "Probability of BAWW increase")+
  scale_x_discrete(breaks = c("c","e"),
                   labels = c("Control", "Playback"))+
  scale_y_continuous(limits = c(0.0,0.8))+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.margin = unit(c(3,3,3,3), "lines"))
propFssFsiPEB

fssFsiPlotWrap <- ggarrange(fssPM,propFssFsiPEB, widths = c(1.5,1))
fssFsiPlotWrap

fssAnnotated <- annotate_figure(fssFsiPlotWrap, top = text_grob("FSS", size = 20, face = "bold"))
fssAnnotated


#####sgl model#####
sgl.df <- data.frame(subset(df, df$site == "sgla"))
mdl.sgl.fsiT <- stan_glm(sgl.df$fsi.t ~ sgl.df$treat + sgl.df$medianT,
                         data = sgl.df,
                         algorithm = "sampling",
                         family = binomial(link = "logit"),
                         prior = student_t(df = 7), 
                         prior_intercept = student_t(df = 7),
                         chains = 10,
                         seed = 12345)
summary(mdl.sgl.fsiT)
launch_shinystan(mdl.sgl.fsiT)

sglPM <- plot_model(mdl.sgl.fsiT,
                    transform = NULL,
                    #colors = "#377EB8",
                    bpe = "median",
                    vline.color = "black",
                    prob.inner = .50,
                    prob.outer = .95,
                    show.values = TRUE, 
                    value.offset = .3,
                    line.size = 1,
                    title = "E",
                    axis.labels = c("Playback", "Median understory vegetation cover"))+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.margin = unit(c(3,3,3,3), "lines"))
sglPM

pred.sgl <- posterior_predict(mdl.sgl.fsiT, draws = 500)
pred.sgl

dim(pred.sgl)

#Calculate the proportion of those 500 draws in which the model predicts 
#a focal species increase
sgl.df$PostPred <- colMeans(pred.sgl)

#Calculate the mean and 95% CI for predicted focal species increase at c & e stops
sglMeanC <- MeanCI(sgl.df$PostPred[sgl.df$treat=="c"])
sglMeanE <- MeanCI(sgl.df$PostPred[sgl.df$treat=="e"])
sglMeanE
#Dataframe to plot the treatment means with 95% CI
sglPlot.df <- data.frame(rbind(sglMeanC,sglMeanE))
sglPlot.df$treat <- c("c","e")

propSglFsiPEB <- ggplot(sglPlot.df, aes(x=treat, y = mean))+
  geom_point(size = 2)+
  geom_errorbar(sglPlot.df, mapping=aes(x=treat, ymin = lwr.ci, ymax=upr.ci), width=0.05, size=1)+
  labs(title = "F",
       x = "",
       y = "Probability of BAWW increase")+
  scale_x_discrete(breaks = c("c","e"),
                   labels = c("Control", "Playback"))+
  scale_y_continuous(limits = c(0.0,0.8))+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.margin = unit(c(3,3,3,3), "lines"))
propSglFsiPEB

sglFsiPlotWrap <- ggarrange(sglPM,propSglFsiPEB, widths = c(1.5,1))
sglFsiPlotWrap
sglAnnotated <- annotate_figure(sglFsiPlotWrap, top = text_grob("SGL", size = 20, face = "bold"))
sglAnnotated



#####bunp model#####
bunp.df <- data.frame(subset(df, df$site == "bunp"))
mdl.bunp.fsiT <- stan_glm(bunp.df$fsi.t ~ bunp.df$treat + bunp.df$medianT,
                         data = bunp.df,
                         algorithm = "sampling",
                         family = binomial(link = "logit"),
                         prior = student_t(df = 7), 
                         prior_intercept = student_t(df = 7),
                         chains = 10,
                         seed = 12345)
summary(mdl.bunp.fsiT)
launch_shinystan(mdl.bunp.fsiT)
posterior_interval(mdl.bunp.fsiT)

bunpPM <- plot_model(mdl.bunp.fsiT,
                     transform = NULL,
                     #colors = "#377EB8",
                     bpe = "median",
                     vline.color = "black",
                     prob.inner = .50,
                     prob.outer = .95,
                     show.values = TRUE, 
                     value.offset = .3,
                     line.size = 1,
                     title = "A",
                     axis.labels = c("Playback", "Median understory vegetation cover"))+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.margin = unit(c(3,3,3,3), "lines"))   
bunpPM

pred.bunp <- posterior_predict(mdl.bunp.fsiT, draws = 500)
pred.bunp

dim(pred.bunp)

#Calculate the proportion of those 500 draws in which the model predicts 
#a focal species increase
bunp.df$PostPred <- colMeans(pred.bunp)

#Calculate the mean and 95% CI for predicted focal species increase at c & e stops
bunpMeanC <- MeanCI(bunp.df$PostPred[bunp.df$treat=="c"])
bunpMeanE <- MeanCI(bunp.df$PostPred[bunp.df$treat=="e"])
bunpMeanE

#Dataframe to plot the treatment means with 95% CI
bunpPlot.df <- data.frame(rbind(bunpMeanC,bunpMeanE))
bunpPlot.df$treat <- c("c","e")

propBunpFsiPEB <- ggplot(bunpPlot.df, aes(x=treat, y = mean))+
  geom_point(size = 2)+
  geom_errorbar(bunpPlot.df, mapping=aes(x=treat, ymin = lwr.ci, ymax=upr.ci), width=0.05, size=1)+
  labs(title = "B",
       x = "",
       y = "Probability of BAWW increase")+
  scale_x_discrete(breaks = c("c","e"),
                   labels = c("Control", "Playback"))+
  scale_y_continuous(limits = c(0.0,0.8))+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.margin = unit(c(3,3,3,3), "lines"))
propBunpFsiPEB

bunpFsiPlotWrap <- ggarrange(bunpPM,propBunpFsiPEB, widths = c(1.5,1))
bunpFsiPlotWrap
bunpAnnotated <- annotate_figure(bunpFsiPlotWrap, top = text_grob("BUNP", size = 20, face = "bold"))
bunpAnnotated



#####SitesWrap#####
ggarrange(bunpAnnotated,fssAnnotated,sglAnnotated, ncol = 1)


#####odds ratios for site models##### 
#BUNP
exp(2.96) #19 

#SGL
exp(3.36) #29

#FSS
exp(2.11) #8
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



#####VEGETATION COMPARISONS BY SITE#####
just2017.df <- data.frame(subset(nodum.df, df$year=="z"))
#uv percent
just2017.df$uvPercent <- just2017.df$uv*100

#arcsine transform percentages 
just2017.df$uvArc <- asin(sqrt(just2017.df$uvPercent/100))
just2017.df$medianArcT <- asin(sqrt(just2017.df$median/100))
just2017.df$iqr75ArcT <- asin(sqrt(just2017.df$iqr75/100))
just2017.df$maxArcT <- asin(sqrt(just2017.df$max/100))

medianArcT <- just2017.df$iqr75ArcT

TukeyHSD(aov(just2017.df$medianArcT~just2017.df$site))
kruskal.test(medianArcT~just2017.df$site, data = just2017.df)
pairwise.wilcox.test(medianArcT, just2017.df$site, p.adjust.method = "none")
medianArcT

#####TREE, SAPLING, SEEDLING DATAFRAME#####
tree.df <- data.frame(subset(just2017.df, just2017.df$treeN > 0))

vegSite <- factor(tree.df$site)
seedN <- tree.df$seedN

seedH <- tree.df$seedH
seedL50 <- tree.df$seedL50
seedG50 <- tree.df$seedG50

sapL1 <- tree.df$sapL1
sap1.3 <- tree.df$sap1.3
sap3.6 <- tree.df$sap3.6
sap6.10 <- tree.df$sap6.10
sapN <- tree.df$sapN
sapH <- tree.df$sapH

treeN <- tree.df$treeN
treeH <- tree.df$treeH
treeEq <- tree.df$treeEq
treeBA <- tree.df$treeTotBA

summary(tree.df)


#####uv median#####
#for colors add scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
medianUvBP <- ggplot(data=just2017.df, aes(x=just2017.df$site, y=just2017.df$median, fill = just2017.df$site))+
  geom_violin()+
  geom_boxplot(fill = "white", width = 0.1)+
  scale_fill_grey(name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  labs(title = "Median")+
  scale_y_continuous(name = "",
                     limits=c(0,80))+
  scale_x_discrete(name = "",
                   breaks = c("bunp","fss","sgla"),
                   labels = c("BUNP", "FSS","SGL"))+
  theme(axis.text.x = element_text(size = 16, color = "black"))+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(axis.title.y = element_text(size = 18,face = "bold"))+
  theme(axis.text.y = element_text(size = 16))+
  theme(legend.position = "none")+
  theme(text = element_text(family = "Arial"))+
  annotate("text", x = as.factor(unique(just2017.df$site)), y = c(52, 48, 19),
           label = c("B", "B", "A"), size = 6)
medianUvBP

#Pairwise comparisons
kruskal.test(just2017.df$medianArcT~just2017.df$site)
DunnTest(just2017.df$medianArcT~just2017.df$site, method = "none")

#####uv maximum#####
# for color add scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
maxUvBP <- ggplot(data=just2017.df, aes(x=just2017.df$site, y=just2017.df$max, fill = just2017.df$site))+
  geom_violin()+
  geom_boxplot(fill = "white", width = 0.1)+
  scale_fill_grey(name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  labs(title = "Maximum")+
  scale_y_continuous(name = "",
                     limits=c(0,100))+
  scale_x_discrete(name = "",
                   breaks = c("bunp","fss","sgla"),
                   labels = c("BUNP", "FSS","SGL"))+
  theme(axis.text.x = element_text(size = 16, color = "black"))+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(axis.title.y = element_text(size = 18,face = "bold"))+
  theme(axis.text.y = element_text(size = 16))+
  theme(legend.position = "none")+
  theme(text = element_text(family = "Arial"))+
  annotate("text", x = as.factor(unique(just2017.df$site)), y = c(100, 90, 66),
           label = c("B", "B", "A"), size = 6)
maxUvBP

#Pairwise comparisons
kruskal.test(just2017.df$maxArcT~just2017.df$site)
DunnTest(just2017.df$maxArcT~just2017.df$site, method = "none")

#####uv minimum#####
minUvBP <- ggplot(data=just2017.df, aes(x=just2017.df$site, y=just2017.df$min, fill = just2017.df$site))+
  geom_violin()+
  geom_boxplot(fill = "white", width = 0.1)+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  labs(title = "")+
  scale_y_continuous(name = "Minimum understory vegetation cover (%)",
                     limits=c(0,30))+
  scale_x_discrete(name = "",
                   breaks = c("bunp","fss","sgla"),
                   labels = c("BUNP", "FSS","SGL"))+
  theme(axis.text.x = element_text(size = 16, color = "black", face = "bold"))+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(axis.title.y = element_text(size = 18,face = "bold"))+
  theme(axis.text.y = element_text(size = 16))+
  theme(legend.position = "none")+
  annotate("text", x = as.factor(unique(just2017.df$site)), y = c(52, 48, 20),
           label = c("B", "B", "A"), size = 6)
minUvBP



#####uv mean#####
meanUvBP <- ggplot(data=just2017.df, aes(x=just2017.df$site, y=just2017.df$uvPercent, fill = just2017.df$site))+
  geom_boxplot()+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  labs(title = "Mean")+
  scale_y_continuous(name = "",
                     limits=c(0,80))+
  scale_x_discrete(name = "",
                   breaks = c("bunp","fss","sgla"),
                   labels = c("BUNP", "FSS","SGL"))+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(axis.text.x = element_text(size = 16, color = "black"))+
  theme(axis.title.y = element_text(size = 18,face = "bold"))+
  theme(axis.text.y = element_text(size = 16))+
  annotate("text", x = as.factor(unique(just2017.df$site)), y = c(52, 48, 22),
           label = c("B", "B", "A"), size = 6)
meanUvBP

#use tukey test for multiple comparisons 
TukeyHSD(aov(just2017.df$uvArc~just2017.df$site))


#####uv iqr75#####
#for colors add scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
iqr75UvBP <- ggplot(data=just2017.df, aes(x=just2017.df$site, y=just2017.df$iqr75, fill = just2017.df$site))+
  geom_violin()+
  geom_boxplot(fill = "white", width = .1)+
  scale_fill_grey(name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  labs(title = "3rd IQR")+
  scale_y_continuous(name = "",
                     limits=c(0,80))+
  scale_x_discrete(name = "",
                   breaks = c("bunp","fss","sgla"),
                   labels = c("BUNP", "FSS","SGL"))+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(axis.text.x = element_text(size = 16, color = "black"))+
  theme(axis.title.y = element_text(size = 18,face = "bold"))+
  theme(axis.text.y = element_text(size = 16))+
  theme(text = element_text(family = "Arial"))+
  annotate("text", x = as.factor(unique(just2017.df$site)), y = c(68, 60, 27),
           label = c("B", "B", "A"), size = 6)
iqr75UvBP

#Pairwise comparisons
kruskal.test(just2017.df$iqr75ArcT~just2017.df$site)
DunnTest(just2017.df$iqr75ArcT~just2017.df$site, method = "none")


#####UV METRICS WRAP####

uvMetricsWrap <- ggarrange(medianUvBP,iqr75UvBP,maxUvBP,
                           ncol = 3, legend = "none")

uvMetricsWrap

annotate_figure(uvMetricsWrap,
                top = text_grob("",
                                face = "bold",
                                size = 18),
                left = text_grob("% Understory vegetation cover",
                                 size = 18,
                                 face = "bold",
                                 rot = 90))



ggplot(data = tree.df, aes(x = vegSite, y = sapN, fill = vegSite))+
  geom_violin()+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  geom_boxplot(fill = "white", width = 0.2)+
  ylab("Number of saplings")+
  scale_x_discrete(breaks = c("bunp", "fss", "sgla"),
                   labels = c("BUNP", "FSS", "SGL"))+
  xlab("")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))

#Pairwise comparisons
kruskal.test(sapN ~ vegSite, data = tree.df)
pairwise.wilcox.test(tree.df$sapN, tree.df$site, p.adjust.method = "none")
#No significant differences



#####SapL1 BVP####
ggplot(data = tree.df, aes(x = vegSite, y = sapL1, fill = vegSite))+
  geom_violin()+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  geom_boxplot(fill = "white", width = 0.2)+
  ylab("Number of saplings > 1 cm")+
  scale_x_discrete(breaks = c("bunp", "fss", "sgla"),
                   labels = c("BUNP", "FSS", "SGL"))+
  xlab("")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))

#Pairwise comparisons
wilcox_test(tree.df$sapL1~tree.df$site)
?pairwise.wilcox.test()
kruskal.test(sapL1 ~ vegSite, data = tree.df)
pairwise.wilcox.test(tree.df$sapL1, tree.df$site, p.adjust.method = "none", exact = FALSE, paired = FALSE)
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  tree.df$sapL1 and tree.df$site 
# 
#         bunp   fss   
#   fss  0.0030 -     
#   sgla 0.0047 0.1735
# 
# P value adjustment method: none 
DunnTest(tree.df$sapL1, tree.df$site, method = "none")

####SapH BVP####
#for colors add scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),

ggplot(data = tree.df, aes(x = vegSite, y = sapH, fill = vegSite))+
  geom_violin()+
  scale_fill_grey(name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  geom_boxplot(fill = "white", width = 0.2)+
  ylab("Sapling species diversity (H')")+
  scale_x_discrete(breaks = c("bunp", "fss", "sgla"),
                   labels = c("BUNP", "FSS", "SGL"))+
  xlab("")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))

#Pairwise comparisons
kruskal.test(sapH ~ vegSite, data = tree.df)
DunnTest(tree.df$sapH, tree.df$site, method = "none")
#No significant differences



####SeedN BVP####
tree.df$seedN
#for colors add scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),

seed.n.bp <- ggplot(data = tree.df, aes(x = vegSite, y = seedN, fill = vegSite))+
  geom_violin()+
  scale_fill_grey(name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  geom_boxplot(fill = "white", width = 0.2)+
  scale_x_discrete(breaks = c("bunp", "fss", "sgla"),
                   labels = c("BUNP", "FSS", "SGL"))+
  scale_y_continuous(name = "Number of seedlings",
                     limits=c(0,1800))+
  xlab("")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        text = element_text(family = "Arial"))+
  annotate("text", x = as.factor(unique(tree.df$site)), y = c(1795, 1200,1075), 
           label = c("B","B","A"), size = 6)
seed.n.bp

#Pairwise comparisons
kruskal.test(seedN ~ vegSite, data = tree.df)
DunnTest(tree.df$seedN, tree.df$site, method = "none")


####SeedL50 BVP####
# compute lower and upper whiskers for each group


seed.l50.bp <- ggplot(data = tree.df, aes(x = vegSite, y = seedL50, fill = vegSite))+
  geom_violin()+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  geom_boxplot(fill = "white", width = 0.2)+
  scale_y_continuous(name = "Number of seedlings < 50cm",
                     limits=c(0,1600))+
  scale_x_discrete(breaks = c("bunp", "fss", "sgla"),
                   labels = c("BUNP", "FSS", "SGL"))+
  xlab("")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))+
  annotate("text", x = as.factor(unique(tree.df$site)), y = c(1550, 1000, 1050),
           label = c("B", "AB", "A"), size = 6)
seed.l50.bp

#Pairwise comparisons
kruskal.test(seedL50 ~ vegSite, data = tree.df)
pairwise.wilcox.test(tree.df$seedL50, tree.df$site, p.adjust.method = "none")
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  tree.df$seedL50 and tree.df$site 
# 
#        bunp  fss  
# fss  0.073 -    
# sgla 0.022 0.485
# 
# P value adjustment method: none 



####SeedG50 BVP####
length(seedG50)
seed.g50.bp <- ggplot(data = tree.df, aes(x = vegSite, y = seedG50, fill = vegSite))+
  geom_violin()+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  geom_boxplot(fill = "white", width = 0.2)+
  scale_y_continuous(name = "Number of seedlings > 50cm",
                     limits=c(0,350))+
  scale_x_discrete(breaks = c("bunp", "fss", "sgla"),
                   labels = c("BUNP", "FSS", "SGL"))+
  xlab("")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))+
  annotate("text", x = as.factor(unique(tree.df$site)), y = c(315, 300, 65),
           label = c("B", "B", "A"), size = 6)
seed.g50.bp

#Pairwise comparisons
kruskal.test(seedG50 ~ vegSite, data = tree.df)
pairwise.wilcox.test(tree.df$seedG50, tree.df$site, p.adjust.method = "none")
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  tree.df$seedG50 and tree.df$site 
# 
#       bunp   fss   
# fss  0.0030 -     
# sgla 0.0047 0.3939
# 
# P value adjustment method: none 




####SeedH BVP####
#for colors add scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),

seed.h.bp <- ggplot(data = tree.df, aes(x = vegSite, y = seedH, fill = vegSite))+
  geom_violin()+
  scale_fill_grey(name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  geom_boxplot(fill = "white", width = 0.2)+
  scale_y_continuous(name = "Seedlings diversity (H')",
                     limits=c(0,2.3))+
  scale_x_discrete(breaks = c("bunp", "fss", "sgla"),
                   labels = c("BUNP", "FSS", "SGL"))+
  xlab("")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        text = element_text(family = "Arial"))+
  annotate("text", x = as.factor(unique(tree.df$site)), y = c(2.1, 1.95, 1.75),
           label = c("B", "A", "A"), size = 6)
seed.h.bp

#Pairwise comparisons
kruskal.test(seedH ~ vegSite, data = tree.df)
DunnTest(tree.df$seedH, tree.df$site, method = "none")




####TreeN BVP####
tree.n.bp <- ggplot(data = tree.df, aes(x = vegSite, y = treeN, fill = vegSite))+
  geom_violin()+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  geom_boxplot(fill = "white", width = 0.2)+
  scale_y_continuous(name = "Number of trees",
                     limits=c(10,55))+
  scale_x_discrete(breaks = c("bunp", "fss", "sgla"),
                   labels = c("BUNP", "FSS", "SGL"))+
  xlab("")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))+
  annotate("text", x = as.factor(unique(tree.df$site)), y = c(41, 44, 51),
           label = c("A", "B", "A"), size = 6)
tree.n.bp

#Pairwise comparisons
kruskal.test(treeN ~ vegSite, data = tree.df)
DunnTest(tree.df$treeN, tree.df$site, method = "none")




####TreeH BVP####
tree.h.bp <- ggplot(data = tree.df, aes(x = vegSite, y = treeH, fill = vegSite))+
  geom_violin()+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  geom_boxplot(fill = "white", width = 0.2)+
  ylab("Diversity (H')")+
  scale_x_discrete(breaks = c("bunp", "fss", "sgla"),
                   labels = c("BUNP", "FSS", "SGL"))+
  xlab("")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
tree.h.bp

#Pairwise comparisons
kruskal.test(treeH ~ vegSite, data = tree.df)
DunnTest(tree.df$treeH, tree.df$site, method = "none")
#No significant differences




####TreeBA BVP####
tree.ba.bp <- ggplot(data = tree.df, aes(x = vegSite, y = treeBA, fill = vegSite))+
  geom_violin()+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  geom_boxplot(fill = "white", width = 0.2)+
  ylab("Total basal area (dm^2)")+
  scale_x_discrete(breaks = c("bunp", "fss", "sgla"),
                   labels = c("BUNP", "FSS", "SGL"))+
  xlab("")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
tree.ba.bp

#Pairwise comparisons
kruskal.test(treeBA ~ vegSite, data = tree.df)

#No significant differences


tree.df$MeanBA <- tree.df$treeTotBA/tree.df$treeN

tree.meanBa.bp <- ggplot(data = tree.df, aes(x = vegSite, y = tree.df$MeanBA, fill = vegSite))+
  geom_violin()+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
                    breaks = c("bunp","fss", "sgla"),
                    labels = c("BUNP","FSS", "SGL"))+
  geom_boxplot(fill = "white", width = 0.2)+
  ylab("Mean basal area (dm^2)")+
  scale_x_discrete(breaks = c("bunp", "fss", "sgla"),
                   labels = c("BUNP", "FSS", "SGL"))+
  xlab("")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))
tree.meanBa.bp
pairwise.wilcox.test(tree.df$MeanBA, tree.df$site, p.adjust.method = "none")








siteDiffWrap <- ggarrange(seed.l50.bp,seed.g50.bp,seed.h.bp,tree.n.bp, legend = "none")
siteDiffWrap

annotate_figure(treeWrap,
                top = text_grob("Tree metrics", size = 20, face = "bold"))

#####Seedling Wrap BVP#####

seedWrap <- ggarrange(seed.n.bp,seed.h.bp, legend = "none")
seedWrap



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#####SAPLING SIZE CLASS DATAFRAME#####

sap.df <- data.frame(site = c("BUNP", "FSS", "SGL"),
                     sapL1Mean = c(mean(sapL1[tree.df$site=="bunp"]),
                                   mean(sapL1[tree.df$site=="fss"]),
                                   mean(sapL1[tree.df$site=="sgla"])),
                     sapL1SD = c(sd(sapL1[tree.df$site=="bunp"]),
                                 sd(sapL1[tree.df$site=="fss"]),
                                 sd(sapL1[tree.df$site=="sgla"])),
                     sap1.3Mean = c(mean(sap1.3[tree.df$site=="bunp"]),
                                    mean(sap1.3[tree.df$site=="fss"]),
                                    mean(sap1.3[tree.df$site=="sgla"])),
                     sap1.3SD = c(sd(sap1.3[tree.df$site=="bunp"]),
                                  sd(sap1.3[tree.df$site=="fss"]),
                                  sd(sap1.3[tree.df$site=="sgla"])),
                     sap3.6Mean = c(mean(sap3.6[tree.df$site=="bunp"]),
                                    mean(sap3.6[tree.df$site=="fss"]),
                                    mean(sap3.6[tree.df$site=="sgla"])),
                     sap3.6SD = c(sd(sap3.6[tree.df$site=="bunp"]),
                                  sd(sap3.6[tree.df$site=="fss"]),
                                  sd(sap3.6[tree.df$site=="sgla"])),
                     sap6.10Mean = c(mean(sap6.10[tree.df$site=="bunp"]),
                                     mean(sap6.10[tree.df$site=="fss"]),
                                     mean(sap6.10[tree.df$site=="sgla"])),
                     sap6.10SD = c(sd(sap6.10[tree.df$site=="bunp"]),
                                   sd(sap6.10[tree.df$site=="fss"]),
                                   sd(sap6.10[tree.df$site=="sgla"]))
                     
)

sap2.df <- data.frame(read.csv("sapSiteCompare.csv", h=T))
sap2.df$size <- factor(sap2.df$size, levels = c("sapL1", "sap1.3", "sap3.6", "sap6.10"))
max(sap1.3)


#####Sapling size class comparison LP#####
# scale_color_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
#                    name = "Site")+
#   scale_x_discrete(name = "Sapling diameter classes (cm)",
#                    breaks = c("sapL1", "sap1.3", "sap3.6", "sap6.10"),
#                    labels = c("<1", "1-3", "3-6", "6-10"))+
#   scale_y_continuous(name = "Number of saplings",
#                      limits = c(-5,150),
#                      breaks = seq(0,150,25))

sapLP <- ggplot(sap2.df, aes(x = sap2.df$size, 
                             y = sap2.df$mean,
                             group = sap2.df$site,
                             color = sap2.df$site))+
  
  geom_line(position = position_dodge(0.3), size = 1.3)+
  geom_point(position = position_dodge(0.3), size = 3)+
  geom_errorbar(aes(ymin = sap2.df$mean-sap2.df$sd, ymax = sap2.df$mean+sap2.df$sd), width = .2,
                position = position_dodge(0.3))+
  scale_color_grey(name = "Site")+
  scale_x_discrete(name = "Sapling diameter classes (cm)",
                   breaks = c("sapL1", "sap1.3", "sap3.6", "sap6.10"),
                   labels = c("<1", "1-3", "3-6", "6-10"))+
  scale_y_continuous(name = "Number of saplings",
                     limits = c(-5,150),
                     breaks = seq(0,150,25))

sapLP + theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
              axis.text.x = element_text(size = 14, face = "bold"),
              axis.title.x = element_text(size = 18, face = "bold"),
              axis.text.y = element_text(size = 14, face = "bold"),
              axis.title.y = element_text(size = 18, face = "bold"),
              legend.title = element_text(size = 14, face = "bold"),
              legend.text = element_text(size = 12),
              legend.background = element_blank(),
              legend.box.background = element_rect(color = "black"),
              legend.position = c(0.8,0.8),
              text = element_text(family = "Arial"))+
  annotate("text", x = c(0.9,1,1.1,1.9,2,2.1,2.9,3,3.1,3.9,4,4.1), y = c(12,127,97,68,145,95,80,63,41,30,23,21),
           label = c("A","B","B","A","A","A","A","A","A","A","B","B"), size = 4)+
  annotate("rect", xmin = 0.8, xmax = 1.2, ymin=-5, ymax = 150, alpha = 0.05)+
  annotate("rect", xmin = 1.8, xmax = 2.2, ymin=-5, ymax = 150, alpha = 0.05)+
  annotate("rect", xmin = 2.8, xmax = 3.2, ymin=-5, ymax = 150, alpha = 0.05)+
  annotate("rect", xmin = 3.8, xmax = 4.2, ymin=-5, ymax = 150, alpha = 0.05)

mean(sapL1[tree.df$site=="fss"])
mean(sapL1[tree.df$site=="sgla"])


#Pairwise comparisons
kruskal.test(sapL1 ~ vegSite, data = tree.df)
DunnTest(tree.df$sapL1, tree.df$site, method = "none")

#Pairwise comparisons
kruskal.test(sap1.3 ~ vegSite, data = tree.df)
DunnTest(tree.df$sap1.3, tree.df$site, method = "none")

#Pairwise comparisons
kruskal.test(sap3.6 ~ vegSite, data = tree.df)
DunnTest(tree.df$sap3.6, tree.df$site, method = "none")

#Pairwise comparisons
kruskal.test(sap6.10 ~ vegSite, data = tree.df)
DunnTest(tree.df$sap6.10, tree.df$site, method = "none")


#####Only SGL non-experimental data####
sglVeg.df <- data.frame(read.csv("sglBawwVeg.csv", h=T))
sglVegYear.df <- subset(sglVeg.df, sglVeg.df$year != "2018")

#####BAWW cumm vs uv Median PP: One year of data from sgla and b with no experimental treatment#####
ggplot(data = sglVegYear.df, aes(x= sglVegYear.df$median, y= sglVegYear.df$baww.cumm, shape = sglVegYear.df$site))+
  geom_point(size = 3)+
  stat_smooth(method="glm", method.args=list(family="poisson"), color = c("black"), linetype = c("solid"))+
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

summary(sglVeg.mdl)
posterior_interval(sglVeg.mdl)
sglVeg.mdl.rsq <- bayes_R2(sglVeg.mdl)

#####Just sglb plots#####
sglbVeg.df <- subset(sglVegYear.df, sglVegYear.df$site == "sglb")
sglbVeg.df$medianArcT <- asin(sqrt(sglbVeg.df$median/100))

sglaVeg.df <- subset(sglVegYear.df, sglVegYear.df$site == "sgla")
sglaVeg.df$medianArcT <- asin(sqrt(sglaVeg.df$median/100))



ggplot(data = sglbVeg.df, aes(x= sglbVeg.df$medianArcT, y= sglbVeg.df$baww.cumm))+
  geom_point(shape = 1, size = 3)+
  stat_smooth(method="glm", method.args=list(family="poisson"), color = "black")+
  xlab("Arcsine-transformed median understory vegetation cover")+
  ylab("Cummulative BAWW observations")+
  scale_x_continuous(breaks = seq(from = 0.2, to = 1.3, by = 0.2))+
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18, face = "bold"),
        text = element_text(family = "Arial"))

sglbVeg.mdl <- stan_glm(sglbVeg.df$baww.cumm~sglbVeg.df$medianArcT, 
                        data = sglbVeg.df,
                        algorithm = "sampling",
                        family = poisson,
                        prior = student_t(df = 7), 
                        prior_intercept = student_t(df = 7),
                        chains = 10,
                        seed = 12345)
launch_shinystan(sglbVeg.mdl)
posterior_interval(sglbVeg.mdl, prob = 0.95)
summary(sglbVeg.mdl)
plot(sglbVeg.mdl)


