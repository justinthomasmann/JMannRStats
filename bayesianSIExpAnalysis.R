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



####FULL DATASET####

data.comp <- data.frame(read.csv("soc_info_expt_complete_NoNAs.csv", h=T))
colnames(data.comp)

levels(data.comp$site)

dummystan.mdl <- stan_glm(data.comp$fsi.t  ~ factor(data.comp$site) + data.comp$treat + data.comp$median, 
                          data = data.comp,
                          algorithm = "sampling",
                          family = binomial(link = "logit"),
                          prior = student_t(df = 7), 
                          prior_intercept = student_t(df = 7))

plot_model(dummystan.mdl,
           colors = "#377EB8",
           show.values = T,
           sort.est = T)
tab_model(dummystan.mdl)

SI.df <- data.frame(subset(data.comp, data.comp$site!= "adummy"))
colnames(SI.df)

SI.df$uvRangeT <- SI.df$uvRange*0.1
SI.df$uvVarT <- SI.df$uvVar*0.01
SI.df$uvT <- SI.df$uv*10
SI.df$gvT <- SI.df$gv*10
SI.df$cvT <- SI.df$cv*10
SI.df$uvM <- SI.df$uvMetric*0.1
SI.df$iqr75T <- SI.df$iqr75*0.1
SI.df$iqr50T <- SI.df$iqr50*0.1
SI.df$medianT <- SI.df$median*0.1


######no significant differences between median uv at control and experimental sites#####
#obviously for fss and bunp since c & e were switched year to year, but importantly also not at SGL, 
#which did not switch

t.test(SI.df$medianT[SI.df$site=="fss" & SI.df$treat=="c"], SI.df$medianT[SI.df$site=="fss" & SI.df$treat=="e"])
#t = -0.14705, df = 20.674, p-value = 0.8845

t.test(SI.df$medianT[SI.df$site=="sgla" & SI.df$treat=="c"], SI.df$medianT[SI.df$site=="sgla" & SI.df$treat=="e"])
#t = -0.63692, df = 8.9928, p-value = 0.54

t.test(SI.df$medianT[SI.df$site=="bunp" & SI.df$treat=="c"], SI.df$medianT[SI.df$site=="bunp" & SI.df$treat=="e"])
#t = 0.25551, df = 20.953, p-value = 0.8008



###########full dataset variables##############

site <- as.factor(SI.df$site)
stop <- as.factor(SI.df$stop)
year <- as.factor(SI.df$year)
treat <- as.factor(SI.df$treat)
bawwCumm <- SI.df$baww.cumm
fsiT <- as.factor(SI.df$fsi.t)
fsiP <- as.factor(SI.df$fsi.p)
fsiB <- as.factor(SI.df$fsi.b)
uvR <- SI.df$uvRange
uvV <- SI.df$uvVar
gv <- SI.df$gv
cvb <- SI.df$cvb
uvT <- SI.df$uvT
gvT <- SI.df$gvT
cvT <- SI.df$cvT
uvVarT <- SI.df$uvVarT
uvRangeT <- SI.df$uvRangeT
uvMedT <- SI.df$medianT
uv75T <- SI.df$iqr75T
uv50T <- SI.df$iqr50T
summary(SI.df)


pairs(~bawwCumm + fsiT + fsiP + fsiB + uvRangeT + uvMedT + uv75T + gvT, 
      diag.panel = panel.hist,
      lower.panel = panel.smooth,
      upper.panel = panel.cor)

vegPca <- prcomp(SI.df[,c(24,26,28)])
summary(vegPca)
vegPca
ggbiplot(vegPca)

pc1 <- vegPca$x[,1]
SI.df$pc1T <- pc1*(-1)#Reverse signs so that values increase with increasing vegetation metrics
pc1T <- SI.df$pc1T

######bayesian models with rstanarm#########
freq.test.mdl <- glm(fsiT  ~ factor(site) + year + treat + uv75T, 
                     data = SI.df,
                     family = binomial(link = "logit"))
summary(freq.test.mdl)
lmerTest::step(freq.test.mdl)#drops year first, then uv75

freq.test.mdl2 <- glm(fsiT  ~ factor(site) + treat + pc1T, 
                      data = SI.df,
                      family = binomial(link = "logit"))
summary(freq.test.mdl2)


######Uv75 model fits the data slightly better (non-significantly)########
#However when modeling the data in "whatBirdsWantVeg.R" there is a stronger better fit between median
#and cummulative BAWW obs. 
med.mdl <- glm(fsiT  ~ factor(site) + treat + uvMedT, 
               data = SI.df,
               family = binomial(link = "logit"))

uv75.mdl <- glm(fsiT  ~ factor(site) + treat + uv75T, 
                data = SI.df,
                family = binomial(link = "logit"))
summary(med.mdl)
summary(uv75.mdl)
anova(med.mdl,uv75.mdl)


stan.test.mdl <- stan_glm(fsiT  ~ factor(site) + treat + uvMedT, 
                          data = SI.df,
                          family = binomial(link = "logit"),
                          prior = student_t(df = 7), 
                          prior_intercept = student_t(df = 7))
####Use the PPcheck tab on this webpage to validate models!!!####
launch_shinystan(stan.test.mdl)
stan.test.mdl
summary(stan.test.mdl)


########Bayesian FsiT plot_model#########
plot_model(stan.test.mdl,
           show.intercept = T,
           transform = NULL,
           bpe = "median",
           bpe.style = "dot",
           vline.color = "black",
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Full-season focal species increase")+
  theme_classic2()+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

####
##########plot the predicted curve of treatment effects across all three sites
par(mar = c(4, 4, 1, 1)) # Reduce some of the margins so that the plot fits better
dat <- as.data.frame(cbind(site,treat,fsiT))
dat
newdat <- dat-1
newdat


ggplot(data = newdat, aes(x=treat, y=fsiT)) +
  geom_point(shape = 1, size = 2, position = position_jitter(width = 0.05, height = 0.05)) +
  stat_smooth(method="glm", method.args=list(family="binomial"))+
  scale_x_continuous(breaks = c(1.00,1.2,1.50,1.8,2.00), labels = c("Control", "", "","", "Experimental"))+
  scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00), labels = c("0", "", "","", "1"))+
  ylab("Full-season focal species increase")+
  xlab("Treatment")+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))+
  annotate("text", x = 1.5, y = .85, label = "Bayes R^2 = 0.34", size = 5)

stan.fsiT.treat.mdl <- stan_glm(fsiT  ~ factor(site) + treat + uvMedT, 
                                data = SI.df,
                                family = binomial(link = "logit"),
                                prior = student_t(df = 7), 
                                prior_intercept = student_t(df = 7))
tab_model(stan.fsiT.treat.mdl)


stan.pre.mdl <- stan_glm(fsiP  ~ factor(site) + treat + uvMedT, 
                          data = SI.df,
                          family = binomial(link = "logit"),
                          prior = student_t(df = 7), 
                          prior_intercept = student_t(df = 7))
get_model_data(stan.pre.mdl,
               transform = NULL,
               show.intercept = TRUE,
               ci.lvl = 0.95)

tab_model(stan.pre.mdl, stan.breed0.mdl,
          transform = NULL, 
          pred.labels = c("BUNP", "FSS", "SGL", "Experimental Treatment", "Median understory vegetation cover"),
          dv.labels = c("Prebreeding Season", "Breeding Season"))

ci95 <- posterior_interval(stan.pre.mdl, prob = 0.95, pars = "treate")
ci95

stanPrePM <- plot_model(stan.pre.mdl,
           show.intercept = T,
           transform = NULL,
           order.terms =  c(1,2,3,4,5),
           bpe = "median",
           bpe.style = "dot",
           vline.color = "black",
           prob.inner = .5,
           prob.outer = .95,
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Prebreeding-season", 
           axis.title = "",
           axis.labels = c("Median understory vegetation cover", "Experimental treatment", "SGL","FSS","BUNP"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16))
stanPrePM

########breeding season data frame (if baww.cumm > 0)#########
breed.df <- data.frame(subset(SI.df, SI.df$baww.cumm > 0))

bSite <- as.factor(breed.df$site)
bStop <- as.factor(breed.df$stop)
bYear <- as.factor(breed.df$year)
bTreat <- as.factor(breed.df$treat)
bBawwCumm <- breed.df$baww.cumm
bFsiT <- as.factor(breed.df$fsi.t)
bFsiP <- as.factor(breed.df$fsi.p)
bFsiB <- as.factor(breed.df$fsi.b)
bUvR <- breed.df$uvRange
bUvV <- breed.df$uvVar
bGv <- breed.df$gv
bCvb <- breed.df$cvb
uvT <- breed.df$uvT
gvT <- breed.df$gvT
cvT <- breed.df$cvT
bUvVarT <- breed.df$uvVarT
bUvRangeT <- breed.df$uvRangeT
bUvMedT <- breed.df$medianT
bUv75T <- breed.df$iqr75T
bUv50T <- breed.df$iqr50T


#########stan breeding model without zeros, so only data from sites with at least one BAWW#########
stan.breed.mdl <- stan_glm(bFsiB ~  factor(bSite) + bTreat + bUvMedT,
                      algorithm = "sampling",
                      data = breed.df,
                      family = binomial(link = "logit"),
                      prior = student_t(df = 7), 
                      prior_intercept = student_t(df = 7))
summary(stan.breed.mdl)

get_model_data(stan.breed.mdl,
               transform = NULL,
               show.intercept = TRUE)

stanBreedPM <- plot_model(stan.breed.mdl,
           show.intercept = T,
           transform = NULL,
           order.terms = c(1,2,3,4,5),
           bpe = "median",
           bpe.style = "dot",
           prob.inner = .5,
           prob.outer = .95,
           vline.color = "black",
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Breeding-season", 
           axis.title = "",
           axis.labels = c("", "", "","",""))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_blank())
stanBreedPM



########stan breeding model with zeros included############
stan.breed0.mdl <- stan_glm(fsiB ~  factor(site) + treat + uvMedT,
                            algorithm = "sampling",
                            data = SI.df,
                            family = binomial(link = "logit"),
                            prior = student_t(df = 7), 
                            prior_intercept = student_t(df = 7))

stanBreed0PM <- plot_model(stan.breed0.mdl,
                          show.intercept = T,
                          transform = NULL,
                          order.terms = c(1,2,3,4,5),
                          bpe = "median",
                          bpe.style = "dot",
                          prob.inner = .5,
                          prob.outer = .95,
                          vline.color = "black",
                          show.values = TRUE, 
                          value.offset = .3,
                          dot.size = 3,
                          line.size = 1,
                          title = "Breeding-season", 
                          axis.title = "",
                          axis.labels = c("", "", "","",""))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_blank())
stanBreed0PM


#########preVsBreFsiPM wrap up##########
stanFsiWrap <- ggarrange(stanPrePM,stanBreed0PM, widths = c(.6,.5),
          ncol = 2)
stanFsiWrap

annotate_figure(stanFsiWrap,
                top = text_grob("",
                                face = "bold",
                                size = 20),
                bottom = text_grob("Log-Odds",
                                   size = 18,
                                   face = "bold")
)



#########cummulative observations model##########

stan.cumm.mdl <- stan_glm(bawwCumm ~  factor(site) + treat + uvMedT,
                          algorithm = "sampling",
                          data = SI.df,
                          family = poisson(),
                          prior = student_t(df = 7), 
                          prior_intercept = student_t(df = 7))
summary(stan.cumm.mdl)

get_model_data(stan.cumm.mdl,
               transform = NULL,
               show.intercept = TRUE)

stanCummPM <- plot_model(stan.cumm.mdl,
                        show.intercept = T,
                        transform = NULL,
                        order.terms =  c(1,2,3,4,5),
                        bpe = "median",
                        bpe.style = "dot",
                        vline.color = "black",
                        prob.inner = .5,
                        prob.outer = .9,
                        show.values = TRUE, 
                        value.offset = .3,
                        dot.size = 3,
                        line.size = 1,
                        title = "Full-season cummulative observations", 
                        axis.title = "",
                        axis.labels = c("Median understory vegetation cover", "Experimental treatment", "SGL","FSS","BUNP"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16))
stanCummPM

noBu.df <- data.frame(subset(SI.df, SI.df$site!="bunp"))

noBu.stan.cumm.mdl <- stan_glm(noBu.df$baww.cumm ~ noBu.df$treat + noBu.df$iqr75T,
                               algorithm = "sampling",
                               data = noBu.df,
                               family = poisson(),
                               prior = student_t(df = 7), 
                               prior_intercept = student_t(df = 7))
plot_model(noBu.stan.cumm.mdl,
           transform = NULL,
           show.intercept = T)

noBu.c.df <- data.frame(subset(noBu.df, noBu.df$treat == "c"))

noBu.c.stan.cumm.mdl <- stan_glm(noBu.c.df$baww.cumm ~ noBu.c.df$medianT,
                               algorithm = "sampling",
                               data = noBu.c.df,
                               family = poisson(),
                               prior = student_t(df = 7), 
                               prior_intercept = student_t(df = 7))
plot_model(noBu.c.stan.cumm.mdl,
           transform = NULL)


#######FSS data frame#######

fss.df <- data.frame(subset(SI.df, SI.df$site == "fss"))
fssTreat <- fss.df$treat
fssFsiT <- fss.df$fsi.t
fssFsiP <- fss.df$fsi.p
fssFsiB <- fss.df$fsi.b
fssCumm <- fss.df$baww.cumm
fssMedT <- fss.df$medianT
fssYear <- fss.df$year

stan.fss.mdl <- stan_glm(fssFsiT ~ fssTreat + fssMedT,
                         algorithm = "sampling",
                         data = fss.df,
                         family = binomial,
                         prior = student_t(df = 7), 
                         prior_intercept = student_t(df = 7))
stanFssPM <- plot_model(stan.fss.mdl,
           transform = NULL,
           order.terms =  c(1,2),
           colors = "#377EB8",
           bpe = "median",
           bpe.style = "dot",
           vline.color = "black",
           prob.inner = .5,
           prob.outer = .95,
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           axis.labels = c("Median understory vegetation cover", "Experimental treatment"),
           axis.title = "",
           title = "FSS")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_blank())
stanFssPM

library(RColorBrewer)
display.brewer.pal(n=8, name = "Set1")
brewer.pal(n = 8, name = "Set1")

######SGL data frame#######

sgl.df <- data.frame(subset(SI.df, SI.df$site == "sgla"))
sglTreat <- sgl.df$treat
sglFsiT <- sgl.df$fsi.t
sglFsiP <- sgl.df$fsi.p
sglFsiB <- sgl.df$fsi.b
sglCumm <- sgl.df$baww.cumm
sglMedT <- sgl.df$medianT
sglYear <- sgl.df$year

stan.sgl.mdl <- stan_glm(sglFsiT ~ sglTreat + sglMedT,
                         algorithm = "sampling",
                         data = sgl.df,
                         family = binomial,
                         prior = student_t(df = 7), 
                         prior_intercept = student_t(df = 7))
stanSglPM <- plot_model(stan.sgl.mdl,
           transform = NULL,
           order.terms =  c(1,2),
           #colors = "#4182dd",
           bpe = "median",
           bpe.style = "dot",
           vline.color = "black",
           prob.inner = .5,
           prob.outer = .95,
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           axis.labels = c("Median understory vegetation cover", "Experimental treatment"),
           axis.title = "",
           title = "SGL")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_blank())
stanSglPM





stan.sglcumm.mdl <- stan_glm(sglCumm ~ sglTreat + sglMedT,
                         algorithm = "sampling",
                         data = sgl.df,
                         family = poisson,
                         prior = student_t(df = 7), 
                         prior_intercept = student_t(df = 7))
plot_model(stan.sglcumm.mdl,
           transform = NULL,
           order.terms =  c(1,2),
           colors = "#4182dd",
           bpe = "median",
           bpe.style = "dot",
           vline.color = "black",
           prob.inner = .5,
           prob.outer = .95,
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1)


######BUNP data frame#######

bunp.df <- data.frame(subset(SI.df, SI.df$site == "bunp"))
bunpTreat <- bunp.df$treat
bunpFsiT <- bunp.df$fsi.t
bunpFsiP <- bunp.df$fsi.p
bunpFsiB <- bunp.df$fsi.b
bunpCumm <- bunp.df$baww.cumm
bunpMedT <- bunp.df$medianT
bunpYear <- bunp.df$year

stan.bunp.mdl <- stan_glm(bunpFsiT ~ bunpTreat + bunpMedT,
                         algorithm = "sampling",
                         data = bunp.df,
                         family = binomial,
                         prior = student_t(df = 7), 
                         prior_intercept = student_t(df = 7))
stanBunpPM <- plot_model(stan.bunp.mdl,
           transform = NULL,
           order.terms =  c(1,2),
           #colors = "#4182dd",
           bpe = "median",
           bpe.style = "dot",
           vline.color = "black",
           prob.inner = .5,
           prob.outer = .95,
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           axis.labels = c("Median understory vegetation cover", "Experimental treatment"),
           axis.title = "",
           title = "BUNP")+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"))
stanBunpPM



stanSitesWrap <- ggarrange(stanBunpPM,stanFssPM,stanSglPM, widths = c(.38,.18,.18),
          ncol = 3)

annotate_figure(stanSitesWrap,
                top = text_grob("",
                                face = "bold",
                                size = 20),
                bottom = text_grob("Log-Odds",
                                   size = 18,
                                   face = "bold")
)
