setwd("D:/data")


library(nlme) #contains logistic regression 
library(AICcmodavg) #contains AICc
library(plyr) #summary statistics
library(ggplot2) ; theme_set(theme_classic())
library(ROCR) #ROC curve
library(lme4)  
library(sjPlot)
library(ggpubr)
library(glmmTMB)
library(MCMCglmm)
library(multcomp)
library(ggbiplot)
library(lmerTest)
library(ggeffects)
library(effsize)
library(rstanarm)
library(shinystan)
library(effects)
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

data.comp <- read.csv("soc_info_expt_complete_NoNAs.csv", h=T)
colnames(data.comp)

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


#########only sgl dataframe and models#################

sgl.df <- data.frame(subset(SI.df, SI.df$site=="sgla"))
sglTreat <- sgl.df$treat
sglBawwCumm <- sgl.df$baww.cumm
sglUvRangeT <- sgl.df$uvRangeT
sglFsiT <- sgl.df$fsi.t
sglFsiP <- sgl.df$fsi.p
sglFsiB <- sgl.df$fsi.b
sglIqr75 <- sgl.df$iqr75T

length(sgl.df$fsi.t[sgl.df$fsi.t==1 & sgl.df$treat=="e"])#3
length(sgl.df$fsi.t[sgl.df$fsi.t==1 & sgl.df$treat=="c"])#0
#In 2017 at SGL, fsi at 3 e and 0 c (50% vs 0%)

sgl.mdl <- glm(sglFsiT ~ sglTreat + sglUvRangeT, family = binomial(link = "logit"))
summary(sgl.mdl)#No effect of treatment or uv Range. 

sgl.bawwCumm.mdl <- glm(sglBawwCumm ~ sglTreat + sglUvRangeT, family = poisson)
summary(sgl.bawwCumm.mdl)#At sgl, only uvRangeT is predictive of increased cummulative counts

sgl.Cumm.iqr.mdl <- glm(sglBawwCumm ~ sglTreat + sglIqr75, family = poisson)
summary(sgl.Cumm.iqr.mdl)

sgl.bawwCumm.PM <- plot_model(sgl.bawwCumm.mdl,
           type = "est",
           transform = NULL,
           vline.color = "black", 
           colors = "#4182dd",
           sort.est = T, 
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Cummulative observations of black-and-white warblers at SGL", 
           axis.title = "",
           axis.labels = c("Range of understory vegetation cover","Experimental treatment"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

sgl.bawwCumm.PM


sgl.bawwCumm.UvR <- plot_model(sgl.bawwCumm.mdl,
           type = "eff",
           terms = "sglUvRangeT",
           title = "SGL")+
  labs(x = "Range of understory vegetation cover",
       y = "Cummulative observations of black-and-white warblers")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = .5),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"))

sgl.bawwCumm.UvR


# sgl.bawwCumm.wrap <- ggarrange(sgl.bawwCumm.PM,sgl.bawwCumm.UvR,
#                           ncol = 1,
#                           nrow = 2,
#                           common.legend = TRUE)
# sgl.bawwCumm.wrap
# 
# 
# annotate_figure(sgl.bawwCumm.wrap,
#                 top = text_grob("",
#                                 face = "bold",
#                                 size = 18),
#                 bottom = text_grob("Effect Size",
#                                    size = 18,
#                                    face = "bold")
# )





sgl.fsiB.mdl <- glm(sglFsiP ~ sglTreat + sglUvRangeT, family = binomial)
summary(sgl.fsiB.mdl)#No effects





#############only fss dataframe and models################

fss.df <- data.frame(subset(SI.df, SI.df$site=="fss"))

fssTreat <- fss.df$treat
fssFsiT <- fss.df$fsi.t
fssUvR <- fss.df$uvRangeT
fssMedT <- fss.df$medianT

t.test(fssMedT[fssTreat=="e"], fssMedT[fssTreat=="c"])

fss.mdl <- glm(fssFsiT ~ fssTreat + fssUvR, family = binomial())
summary(fss.mdl)

#if I adjust the p value for multiple comparisons, 
p.adjust(0.211, method = "holm", n = 2)

fss.logi.PM <- plot_model(fss.mdl,
           type = "est",
           transform = NULL,
           vline.color = "black", 
           colors = "#4182dd",
           sort.est = T, 
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "", 
           axis.title = "Effect size",
           axis.labels = c("Range of understory vegetation cover","Experimental treatment"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

fss.logi.treatEFF <- plot_model(fss.mdl,
           type = "eff",
           terms = "fssTreat")+
  labs(title = "",
       x = "Treatment",
       y = "Full-season focal species increase")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

fss.logi.wrap <- ggarrange(fss.logi.PM,fss.logi.treatEFF,
          ncol = 2)

annotate_figure(fss.logi.wrap,
                top = text_grob("Full-season focal species increase at FSS",
                                face = "bold",
                                size = 18),
                bottom = text_grob("",
                                   size = 18,
                                   face = "bold")
)


fss.fsiB.mdl <- glm(fss.df$fsi.b ~ fss.df$treat + fss.df$uvRangeT, family = binomial)
summary(fss.fsiB.mdl)# no effect of treatment or uv range

fss.bawwCumm.mdl <- glm(fss.df$baww.cumm ~ fss.df$treat + fss.df$uvRangeT, family = poisson)
summary(fss.bawwCumm.mdl) #significant positive effect of treatment but not uv range


length(fss.df$fsi.t[fss.df$fsi.t==1 & fss.df$treat=="e"])#10/13
length(fss.df$fsi.t[fss.df$fsi.t==1 & fss.df$treat=="c"])#3/11
#Across 2016 & 2017 at FSS, fsi at 10 e and 3 c (77% vs 27%)





###########only bunp dataframe and models#############################

bunp.df <- data.frame(subset(SI.df, SI.df$site=="bunp"))
bunp.mdl <- glm(bunp.df$fsi.t ~ bunp.df$treat + bunp.df$uvRangeT, family = binomial())
summary(bunp.mdl)

length(bunp.df$treat[bunp.df$treat=="e"])#13 e and 11 c total
length(bunp.df$fsi.t[bunp.df$fsi.t==1 & bunp.df$treat=="e"])#4
length(bunp.df$fsi.t[bunp.df$fsi.t==1 & bunp.df$treat=="c"])#0
#Across 2016 & 2017 at FSS, fsi at 11 e and 3 c (31% vs 0%)

noBunp.df <- data.frame(subset(SI.df, SI.df$site!="bunp"))




###########just 2017 dataframe and models##############
just2017dum.df <- data.frame(subset(SI.df, SI.df$year=="2017"))

just2017.df <- data.frame(subset(SI.df, SI.df$year=="2017" & SI.df$site != "adummy"))

just2017.mdl <- glmer(just2017.df$fsi.t ~ just2017.df$treat + just2017.df$uvRangeT + (1|just2017.df$site), family = binomial(link = "logit"))
summary(just2017.mdl)

#########plot_model for 2017 data with site as random effect########### 
just2017.mdl.PM <- plot_model(just2017.mdl,
                            type = "est",
                            transform = NULL,
                            colors = "#4182dd",
                            vline.color = "black", 
                            sort.est = T, 
                            show.values = TRUE, 
                            value.offset = .3,
                            dot.size = 3,
                            line.size = 1,
                            title = "Fixed Effects", 
                            axis.title = "",
                            axis.labels = c("Range of understory vegetation cover","Experimental treatment"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

just2017.mdl.PM

just2017.REmdl.PM <- plot_model(just2017.mdl,
                              type = "re",
                              transform = NULL,
                              vline.color = "black",
                              show.values = TRUE,
                              value.offset = .3,
                              dot.size = 3,
                              line.size = 1,
                              title = "Random Effects", 
                              axis.title = "",
                              axis.labels = c( "BUNP","FSS","SGL"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

just2017.REmdl.PM

just2017Wrap <- ggarrange(just2017.mdl.PM,just2017.REmdl.PM,
                        ncol = 2,
                        nrow = 1,
                        common.legend = TRUE)
just2017Wrap


annotate_figure(just2017Wrap,
                top = text_grob("",
                                face = "bold",
                                size = 18),
                bottom = text_grob("Effect Size",
                                   size = 18,
                                   face = "bold")
)






###########full dataset variables#################

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
uvM <- SI.df$uvM
uv75T <- SI.df$iqr75T
uv50T <- SI.df$iqr50T
summary(SI.df)


pairs(~bawwCumm + fsiT + fsiP + fsiB + uvRangeT + uvT + uv75T + uv50T, 
      diag.panel = panel.hist,
      lower.panel = panel.smooth,
      upper.panel = panel.cor)

uvR.mdl <- glm(fsiT ~ factor(site) + treat + uvRangeT + year, 
               data = SI.df, 
               family = binomial(link = "logit"))
summary(uvR.mdl)

lmerTest::step(uvR.mdl)
uvR.BF.mdl <- glm(formula = fsiT ~ factor(site) + treat, family = binomial(link = "logit"))

stan.test.mdl <- stan_glm(fsiT  ~ factor(site) + treat, 
         data = SI.df,
         family = binomial(link = "logit"),
         prior = student_t(df = 7), 
         prior_intercept = student_t(df = 7))

stan.test.mdl
summary(stan.test.mdl)
plot_model(stan.test.mdl,
           show.intercept = T,
           transform = NULL,
           bpe = "mean",
           bpe.style = "dot",
           vline.color = "black",
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Full-season focal species increase", 
           axis.title = "",
           axis.labels = c("Experimental treatment", "SGL","FSS","BUNP"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

dat <- ggpredict(stan.test.mdl, terms = "treat")
dat
dimnames(post)
plot(dat,rawdata = T)

launch_shinystan(stan.test.mdl)

library(bayesplot)
post <- as.array(stan.test.mdl)
post
mcmc_areas(post,
           pars = c("(Intercept)", "factor(site)fss", "factor(site)sgla", "treate"),
           point_est = "mean")+
  vline_0(size = 1, color = "black")+
  xaxis_text()

?xaxis_text()

plot_model(uvR.BF.mdl,
           type = "est",
           show.intercept = T,
           transform = NULL,
           vline.color = "black",
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Full-season focal species increase", 
           axis.title = "",
           axis.labels = c("Experimental treatment", "SGL","FSS","BUNP"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

plot_model(dum.uvR.mdl,
           type = "eff",
           terms = "uvRangeT")


dum.uv.mdl <- glm(fsiT ~ factor(site)~1 + factor(treat) + uvT, family = binomial(link = "logit"))
summary(dum.uv.mdl)

plot_model(dum.uv.mdl,
           type = "est",
           rm.terms = "adummy")

dum.uvM.mdl <- glm(fsiT ~ factor(site) -1 + factor(treat) + uvM, family = binomial(link = "logit"))
summary(dum.uvM.mdl)

dum.uv50.mdl <- glm(fsiT ~ factor(site) -1 + treat + uv50T, family = binomial(link = "logit"))
summary(dum.uv50.mdl)

dum.uv75.mdl <- glm(fsiT ~ factor(site) -1 + treat + uv75T, family = binomial(link = "logit"))
summary(dum.uv75.mdl)


anova(dum.uvR.mdl,dum.uvM.mdl,dum.uv.mdl,dum.uv75.mdl,dum.uv50.mdl)

plot(dum.uvR.mdl$residuals)


######What's going on with uvRange!!!!!!!!!!###############
####test if reducing the demensions of my vegetation variables improves model fit####

SI.noDum.df <- data.frame(subset(SI.df, SI.df$site != "adummy"))

test.mdl <- glm(SI.noDum.df$fsi.t ~ SI.noDum.df$site +  SI.noDum.df$treat, family = binomial)
summary(test.mdl)

plot_model(test.mdl, 
           show.intercept = T,
           axis.labels = c("Experimental treatment", "SGL", "FSS", "BUNP"),
           show.values = T)

test <- glht(test.mdl,linfct=mcp(SI.noDum.df$site ="Tukey"))

test2.mdl <- glm(SI.noDum.df$fsi.t ~ SI.noDum.df$treat, family = binomial)
summary(test2.mdl)

test1.mdl <- glmer(SI.noDum.df$fsi.t ~ SI.noDum.df$treat + SI.noDum.df$site + (1|SI.noDum.df$stop), family = binomial)
summary(test1.mdl)

vegPcaFull <- prcomp(SI.noDum.df[,c(18:25)])
ggbiplot(vegPcaFull, ellipse = T, groups = SI.noDum.df$site, labels = SI.noDum.df$stop)

vegPca1 <- prcomp(SI.noDum.df[,c(18:20,24,25)])
ggbiplot(vegPca1, ellipse = T, groups = SI.noDum.df$site, labels = SI.noDum.df$stop)
summary(vegPca1)
vegPca1
pc1 <- vegPca1$x[,1]
SI.noDum.df$pc1T <- pc1*(-1)#Reverse signs so that values increase with increasing vegetation metrics
pc1T <- SI.noDum.df$pc1T

pc.mdl <- glmer(SI.noDum.df$fsi.t ~ factor(SI.noDum.df$treat) + factor(SI.noDum.df$year) + pc1T + (1|SI.noDum.df$site), family = binomial(link = "logit"))
summary(pc.mdl)

range.mdl <- glmer(SI.noDum.df$fsi.t ~ factor(SI.noDum.df$treat) + factor(SI.noDum.df$year) + SI.noDum.df$uvRangeT + (1|SI.noDum.df$site), family = binomial(link = "logit"))
summary(range.mdl)

mean.mdl <- glmer(SI.noDum.df$fsi.t ~ factor(SI.noDum.df$treat) + factor(SI.noDum.df$year) + SI.noDum.df$uvT + (1|SI.noDum.df$site), family = binomial(link = "logit"))
summary(mean.mdl)

iqr75.mdl <- glmer(SI.noDum.df$fsi.t ~ factor(SI.noDum.df$treat) + factor(SI.noDum.df$year) + SI.noDum.df$iqr75T + (1|SI.noDum.df$site), family = binomial(link = "logit"))
summary(iqr75.mdl)

anova(pc.mdl,range.mdl, mean.mdl, iqr75.mdl)

anova(mean.mdl, iqr75.mdl)

anova(range.mdl, iqr75.mdl)






###########uv range provides the best fit of all vegetation metrics##########
totFsi.uvR.mdl <- glmer(fsiT~treat + uvRangeT + year + (1|site), family = binomial)
summary(totFsi.uvR.mdl)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
par(mar = c(4, 4, 1, 1)) # Reduce some of the margins so that the plot fits better


dat <- as.data.frame(cbind(treat,fsiT))
na.omit(dat)
newdat <- dat-1
newdat
plot(jitter(newdat$treat, .2),jitter(newdat$fsiT, .2) )
plot(treat,fsiT,xlab="treat",ylab="fsiT") 
g=glm(fsiT~treat,family=binomial,data = newdat) 

curve(predict(g,data.frame(treat=x),type="resp"),add=TRUE) 

ggplot(data = newdat, aes(x=treat, y=fsiT)) +
  geom_point(shape = 1, size = 2, position = position_jitter(width = 0.05, height = 0.05)) +
  stat_smooth(method="glm", method.args=list(family="binomial"))+
  scale_x_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00), labels = c("Control", "", "","", "Experimental"))+
  scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00), labels = c("0", "", "","", "1"))+
  ylab("Full-season focal species increase")+
  xlab("Treatment")+
  theme(plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))


dat.uv <- as.data.frame(cbind(uvRangeT,fsiT,treat))
dat.uv$fsiT-1
ggplot(data = dat.uv, aes(x=uvRangeT,y=(fsiT-1)))+
  geom_point(aes(color = factor(treat)), shape = 1, size = 2, position = position_jitter(width = 0.05, height = 0.05)) +
  stat_smooth(method="glm", method.args=list(family="binomial"))
  # scale_x_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00), labels = c("Control", "", "","", "Experimental"))+
  # scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00), labels = c("0", "", "","", "1"))+
  # ylab("Full-season focal species increase")+
  # xlab("Treatment")+
  # theme(plot.title = element_text(size = 18, face = "bold"),
  #       axis.text.x = element_text(size = 12, face = "bold"),
  #       axis.title.x = element_text(size = 14, face = "bold"),
  #       axis.text.y = element_text(size = 12, face = "bold"),
  #       axis.title.y = element_text(size = 14, face = "bold"))







totFsi.varT.mdl <- glmer(fsiT~treat + uvVarT + year + (1|site), family = binomial)
summary(totFsi.varT.mdl)

anova(totFsi.mdl, totFsi.uvR.mdl)#pc1 model fits better than uvT model
anova(totFsi.pc1.mdl, totFsi.uvR.mdl)#but...uvRangeT model fits better that pc1 and uvT
anova(totFsi.uvR.mdl, totFsi.varT.mdl)#uvVarT is not a better fit than uvRangeT


#########plot_model for best fit total fsi model########### 
totFsi.mdl.PM <- plot_model(totFsi.uvR.mdl,
           type = "est",
           transform = NULL,
           vline.color = "black", 
           sort.est = T, 
           show.values = TRUE, 
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Fixed Effects", 
           axis.title = "",
           axis.labels = c( "Year","Range of understory vegetation cover","Experimental treatment"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

totFsi.mdl.PM

totFsi.REmdl.PM <- plot_model(totFsi.uvR.mdl,
           type = "re",
           transform = NULL,
           vline.color = "black",
           show.values = TRUE,
           value.offset = .3,
           dot.size = 3,
           line.size = 1,
           title = "Random Effects", 
           axis.title = "",
           axis.labels = c( "BUNP","FSS","SGL"))+
  theme_classic2()+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

totFsi.REmdl.PM

totFsiWrap <- ggarrange(totFsi.mdl.PM,totFsi.REmdl.PM,
                                 ncol = 2,
                        nrow = 1,
                                 common.legend = TRUE)
totFsiWrap


annotate_figure(totFsiWrap,
                top = text_grob("",
                                face = "bold",
                                size = 18),
                bottom = text_grob("Effect Size",
                                   size = 18,
                                   face = "bold")
)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

plot_model(totFsi.uvR.mdl,
           type = "eff",
           terms = "treat")

plot_model(totFsi.uvR.mdl,
           type = "eff",
           terms = "uvRangeT")



confint(totFsi.mdl)

contrasts(site)

lmerTest::step(totFsi.mdl)
#Best fit model 
totFsi.BF.mdl <- glm(fsiT~ site + treat, family = binomial)
summary(totFsi.BF.mdl)

model.matrix(totFsi.BF.mdl)

plot_model(totFsi.BF.mdl)
tab_model(totFsi.BF.mdl)
print(summary(glht(totFsi.BF.mdl, linfct=mcp(site="Tukey", treat="Tukey"))))

confint(glht(totFsi.BF.mdl))

contrasts(SI.df$site)

totFsi2017.mdl <- glm(just2017.df$fsi.t~just2017.df$treat + just2017.df$site + just2017.df$uv, family = binomial)
summary(totFsi2017.mdl)
tab_model(totFsi2017.mdl)

preFsi2017.mdl <- glm(just2017.df$fsi.p~just2017.df$treat + just2017.df$site + just2017.df$uv, family = binomial)
summary(preFsi2017.mdl)
tab_model(preFsi2017.mdl)

breedFsi2017.mdl <- glm(just2017.df$fsi.b~just2017.df$treat + just2017.df$site + just2017.df$uv, family = binomial)
summary(breedFsi2017.mdl)
tab_model(breedFsi2017.mdl)
plot(breedFsi2017.mdl)

pre.mdl <- glm(SI.df$fsi.p ~ SI.df$treat + SI.df$uv, data = SI.df, family = binomial)
summary(pre.mdl)

breed.mdl <- glm(fsiB ~ site + treat + year, family = binomial)
summary(breed.mdl)

nobu.breed.mdl <- glm(noBunp.df$fsi.b ~ noBunp.df$treat,
                        family = binomial)
summary(nobu.breed.mdl)

fss.breed.mdl <- glmer(fss.df$fsi.b ~ fss.df$treat + fss.df$uvRangeT + (1|fss.df$year), family = binomial)
summary(fss.breed.mdl)

lmerTest::step(breed.mdl)


FM <- glmer(fsiB~treat + uvT + year + (1|site), family = binomial)
summary(FM)

multi.mdl <- glm(cbind(c(fsiT,fsiB,fsiP)) ~ SI.df$treat + site + uvT + year, family = binomial)

uv.mdl <- glmer(fsiT~uv + (1|site), family = binomial)
summary(uv.mdl)

summary(glm(fsiT~gv, family = binomial))


#bawwCumm isn't significantly zero-inflated
summary(glmmTMB(bawwCumm~treat + uvT + (1|site), family = poisson,
                ziformula = ~1))




ggplot(data = data.comp, aes(x=treat, y=bawwCumm))+
  geom_point()+
  geom_abline()

plot_model(FM)

ggplot(data = data.comp, aes(x = uv, y = bawwCumm, color = site))+
  geom_point()+
  geom_abline()




#fsi pre

pre.mdl <- glmmTMB(fsiP ~ treat + uv + site + year, family = gaussian)
summary(pre.mdl)

length(data.comp$fsi.p[data.comp$treat == "e"])
length(data.comp$fsi.p[data.comp$treat=="c"])



length(fsiP)#78
length(fsiP[fsiP==1 & treat== "e"])#30
length(fsiP[fsiP==0 & treat== "e"])#28
length(fsiP[fsiP==1 & treat== "c"])#7
length(fsiP[fsiP==0 & treat== "c"])

length(fsiP[treat =="c"])
length(fsiP[treat =="e"])

ggplot(data = data.comp, aes(x = treat, y = fsiP, color = site))+
  geom_point()


breeding.mdl <- glm(fsiB ~ treat + uv + year + site, family = binomial)
summary(breeding.mdl)

ggplot(data = data.comp, aes(x = treat, y = fsiB, color = site))+
  geom_point()



#Likelihood Ratio Test:
anova(FM, test="Chisq")


####just 2017 dataframe####

just2017.df <- data.frame(subset(SI.df, SI.df$year=="2017"))
na.omit(just2017.df)

####Multiple pairwise comparisons of mean understory vegetation cover####

#uv percent
just2017.df$uvPercent <- just2017.df$uv*100

#arcsine transform percentages 
just2017.df$uvArc <- asin(sqrt(just2017.df$uvPercent/100))

#use tukey test for multiple comparisons 

TukeyHSD(aov(just2017.df$uvArc~just2017.df$site))

#####plot uv per site####


ggplot(data=just2017.df, aes(x=just2017.df$site, y=just2017.df$uvPercent))+
  geom_boxplot()+
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5)+
  scale_y_continuous(name = "Understory vegetation cover (%)",
                     limits=c(0,60))+
  scale_x_discrete(name = "",
                   breaks = c("bunp","fss","sgla"),
                   labels = c("BUNP", "FSS","SGL"))+
  theme(axis.text.x = element_text(size = 18, face = "bold", color = "black"))+
  theme(axis.title.y = element_text(size = 18,face = "bold"))+
  theme(axis.text.y = element_text(size = 16))+
  annotate("text", x = as.factor(unique(just2017.df$site)), y = c(52, 48, 22),
           label = c("B", "B", "A"), size = 6)




#Bayesian stuff####

library(MCMCpack)
library(mcmcplots)
library(coda)

help(package= "MCMCpack")

MCMC.test <- MCMClogit(fsiT~treat + uv + site, data = SI.df)
summary(MCMC.test)

plot(MCMC.test)

codamenu()

MCMC.pre <- MCMClogit(fsiP~treat + uv + site, data = SI.df)
summary(MCMC.pre)
plot(MCMC.pre)

MCMC.breed <- MCMClogit(fsiB~treat + uv + site, data = SI.df)
summary(MCMC.breed)
plot(MCMC.breed)

#Plots

dat <- as.data.frame(cbind(treat,fsiT))
na.omit(dat)

plot(treat,fsiT,xlab="treat",ylab="fsiT") 
g=glm(fsiT~treat,family=binomial,dat) 

curve(predict(g,data.frame(treat=x),type="resp"),add=TRUE) 
fsiT

dat.uv <- as.data.frame(cbind(uv,fsi))
na.omit(dat.uv)

plot(uv,fsi,xlab="uv",ylab="fsi") 
g2=glm(fsi~uv,family=binomial,dat) 
curve(predict(g2,data.frame(uv=x),type="resp"),add=TRUE) 



#ROC curve (receiver operating characteristic)
#Measures the predictive performance of the model
#The closer the auc is to one, 
#the better the predictive ability of the model, "auc = 0.90460453"

p <- predict(FM, newdata=subset(data.comp,select=c(2,3,4,5,6,7,8)), type="response")
pr <- prediction(p, y)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

SI.NoDum.df <- data.frame(subset(SI.df, SI.df$site != "adummy"))
set.seed(364)
sample <- sample(nrow(SI.df),floor(nrow(SI.df)*0.8))
train <- SI.df[sample,]
test <- SI.df[-sample,]

test$pred <- predict(dum.uvR.mdl, test, type="response")

test$good_pred <- ifelse(test$pred > 0.80, "good", "bad")
confusionMatrix(test$good_pred, test$good)

