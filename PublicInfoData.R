setwd("D:/data")
library(unmarked)
library(AICcmodavg) #aictab() for AICc model comparisons
library(ggplot2) ; theme_set(theme_classic())
library(MuMIn) #dredge()
library(DescTools) #MeanCI()
library(sjPlot)

((15^2)*3.14)*2

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

veg.df <- read.csv("PublicInfoUnmarkedVeg.csv", h=T) #2017 data, first four surveys before pb treatment
                                                     #use for baseline habitat relations prior to expt

#transformations rescale metrics to counts
veg.count <- as.matrix(veg.df[,4:7])
veg.treat <- veg.df$treat
veg.site <- veg.df$site
veg.uv <- veg.df$uv*10  
veg.medUv <- veg.df$median*0.1 
veg.cvb <- veg.df$cvb*10
veg.gv <- veg.df$gv*10
veg.varUv <- veg.df$uvVar*0.01
veg.ranUv <- veg.df$uvRange*0.1
veg.uv75 <- veg.df$iqr75*0.1
veg.uv50 <- veg.df$iqr50*0.1
veg.maxUv <- veg.df$max*0.1
veg.minUv <- veg.df$min*0.1

ufPC.veg <- unmarkedFramePCount(y = veg.count, siteCovs = data.frame(veg.treat = veg.treat, 
                                                                 veg.site = veg.site, 
                                                                 veg.uv = veg.uv,
                                                                 veg.medUv = veg.medUv,
                                                                 veg.cvb = veg.cvb,  
                                                                 veg.gv = veg.gv,
                                                                 veg.varUv = veg.varUv,
                                                                 veg.ranUv = veg.ranUv,
                                                                 veg.uv75 = veg.uv75,
                                                                 veg.uv50 = veg.uv50,
                                                                 veg.maxUv = veg.maxUv,
                                                                 veg.minUv = veg.minUv
))
head(ufPC.veg)

pairs(~veg.df$totDet + veg.cvb + veg.gv + veg.uv + veg.medUv + veg.maxUv + veg.uv50 + veg.uv75,
      diag.panel = panel.hist,
      lower.panel = panel.smooth,
      upper.panel = panel.cor)


veg1 <- pcount(~1 ~veg.uv,
               data = ufPC.veg, mixture = "P", K = 102, se = TRUE)
veg2 <- pcount(~1 ~veg.medUv,
              data = ufPC.veg, mixture = "P", K = 102, se = TRUE)
veg3 <- pcount(~1 ~veg.uv75,
               data = ufPC.veg, mixture = "P", K = 102, se = TRUE)
veg4 <- pcount(~1 ~veg.ranUv,
               data = ufPC.veg, mixture = "P", K = 102, se = TRUE)
veg5 <- pcount(~1 ~veg.maxUv,
               data = ufPC.veg, mixture = "P", K = 102, se = TRUE)
veg6 <- pcount(~1 ~veg.varUv,
               data = ufPC.veg, mixture = "P", K = 102, se = TRUE)

vegModels <- list(veg1, veg2, veg3, veg4, veg5, veg6)
vegModnames <- c("veg1", "veg2", "veg3", "veg4", "veg5", "veg6")
aictab(cand.set = vegModels,modnames = vegModnames, sort = TRUE)

veg3 #75IQR provides best fit 

veg3quad <-pcount(~1 ~I(veg.uv75^2),
                          data = ufPC.veg, mixture = "P", K = 102, se = TRUE)
veg3quad #quadratic term produces higher AICc

veg3det <- pcount(~veg.uv75 ~1,
               data = ufPC.veg, mixture = "P", K = 102, se = TRUE)
veg3det #75IQR in detection submodel also produces higher AICc



?Nmix.gof.test()
gof.veg3 <- Nmix.gof.test(veg3, nsim = 1000)
gof.veg3 #c-hat is 1.35, but NB model indicates no significant overdispersion.
#gof test indicates no significant difference between observed and simulated statistic (P = 0.07)

veg3NB <- pcount(~1 ~veg.uv75,
               data = ufPC.veg, mixture = "NB", K = 102, se = TRUE)
veg3NB #not significantly overdispersed

veg3ZIP <- pcount(~1 ~veg.uv75,
               data = ufPC.veg, mixture = "ZIP", K = 102, se = TRUE)
veg3ZIP #not significantly zero-inflated


predVeg <- predict(veg3, type="state")
predVeg

veg.df$pred <- predVeg$Predicted
veg.df$up <- predVeg$upper
veg.df$lo <- predVeg$lower

abd_uv <- ggplot(data = veg.df)+
  geom_line(aes(x=iqr75 ,y=pred), size = 1)+
  geom_ribbon(aes(x=iqr75,y=pred,ymin=lo, ymax=up), alpha=0.1)+
  xlab("Understory vegetation cover (75% IQR)")+
  ylab("Predicted abundance")+
  theme(plot.title = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 22, face = "bold"),
        text = element_text(family = "Arial"),
        plot.margin = margin(15,15,15,15))
suppressWarnings(print(abd_uv))



df <- read.csv("PublicInfoUnmarked.csv", h=T) #2018 data following pb treatment in 2017
head(df)

#transformations rescale metrics to counts
count <- as.matrix(df[,4:8])
treat <- df$treat
site <- df$site
uv <- df$uv*10  
medUv <- df$median*0.1 
cvb <- df$cvb*10
gv <- df$gv*10
varUv <- df$uvVar*0.01
ranUv <- df$uvRange*0.1
uv75 <- df$iqr75*0.1
uv50 <- df$iqr50*0.1
maxUv <- df$max*0.1
minUv <- df$min*0.1
#?unmarkedFramePCount() 
ufPC <- unmarkedFramePCount(y = count, siteCovs = data.frame(treat = treat, 
                                                             site = site, 
                                                             uv = uv,
                                                             medUv = medUv,
                                                             varUv = varUv,
                                                             ranUv = ranUv,
                                                             uv75 = uv75,
                                                             uv50 = uv50,
                                                             maxUv = maxUv,
                                                             minUv = minUv
                                                             ))
head(ufPC)

#?pcount(): pcount(~detection ~abundance, data = df, mixture = P, NB, or ZIP)
m0 <- pcount(~1 ~1, data = ufPC, mixture = "P", K = 102, se = TRUE) #Null model
m0.5 <- pcount(~1 ~site, data = ufPC, mixture = "P", K = 102, se = TRUE)
m0.7 <- pcount(~site ~1, data = ufPC, mixture = "P", K = 102, se = TRUE)

#Model selection using aictab() for AICc (corrected for small sample sizes) 
#Null model (intercept only) vs. Study area detection model 
nullModels <- list(m0, m0.5, m0.7)
nullModnames <- c("m0","m0.5","m0.7")
aictab(cand.set = nullModels,modnames = nullModnames, sort = TRUE) 
#site in abundance submodel is only nominally better fit than null model,
#therefore differences between sgla and sglb are not explaining variance in data.
#sglb is nonsignficantly negatively associated with abundance (Est = -1.26; P = 0.11)

m1 <- pcount(~1 ~treat * uv75, data = ufPC, mixture = "P", K = 102, se = TRUE)
m1

m2 <- pcount(~1 ~treat + uv75, data = ufPC, mixture = "P", K = 102, se = TRUE)
m2

m3 <- pcount(~1 ~uv75, data = ufPC, mixture = "P", K = 102, se = TRUE)
m3


treatModels <- list(m1,m2,m3)
treatModNames <- c("m1","m2","m3")
aictab(cand.set = treatModels, modnames = treatModNames, sort = TRUE)

#interaction model is the best fit
m1

gof.m1 <- Nmix.gof.test(m1, nsim = 1000)
gof.m1 #m1 fits the data well 
#c-hat = 1.04
#P = 0.41

m1NB <- pcount(~1 ~treat * uv75, data = ufPC, mixture = "NB", K = 102, se = TRUE)
m1ZIP <- pcount(~1 ~treat * uv75, data = ufPC, mixture = "ZIP", K = 102, se = TRUE)

distModels <- list(m1,m1NB,m1ZIP)
distModNames <- c("m1","m1NB","m1ZIP")
aictab(cand.set = distModels, modnames = distModNames, sort = TRUE)
#Poisson distribution fits best

confint(m1, type = "state" )

#m1 predictions####
predTreat <- predict(m1, type = "state")
predTreat

df$pred <- predTreat$Predicted
df$Lo <- predTreat$lower
df$Up <- predTreat$upper

#PLOT--interaction between UV and treat####
abd_uv <- ggplot(data = df)+
  geom_line(aes(x=uv75*10,y=pred,color = treat), size = 1)+
  geom_ribbon(aes(x=uv75*10,y=pred,ymin=Lo,ymax=Up,group = treat), alpha=0.1)+
  xlab("Understory vegetation cover (75% IQR)")+
  ylab("Predicted abundance")+
  scale_color_manual(breaks = c("c","e"),
                       labels = c("Control", "Playback"),
                       values = c("#F0E442","#0072B2"))+
  theme(plot.title = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 22, face = "bold"),
        text = element_text(family = "Arial"),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = c(0.8,0.9),
        plot.margin = margin(15,15,15,15))
suppressWarnings(print(abd_uv))




#PLOT--Predicted abundance by treatment ####
PbMean <- MeanCI(df$pred[df$treat=="e"])
PbMean
ContMean <- MeanCI(df$pred[df$treat=="c"])
ContMean
treatPlot <- data.frame(rbind(PbMean,ContMean))
treatPlot$treat <- c("Playback", "Control")

abd_treatment <- ggplot(data = treatPlot)+
  geom_pointrange(aes(x = treat, y=mean, ymin = lwr.ci, ymax = upr.ci),
                  size = 1.2,
                  shape = 20)+
  xlab("")+
  scale_x_discrete(labels = c("c"= "Control", "e" = "Playback"))+
  ylab("Predicted abundance")+
  ggtitle("")+
  theme(plot.title = element_text(size = 24),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 22, face = "bold"),
        text = element_text(family = "Arial"),
        plot.margin = margin(15,15,15,15))
suppressWarnings(print(abd_treatment))

plotModel.df <- data.frame(params = c("Playback","UVC75","Playback*UVC75"),
                           esti = c(-6.634,-0.597,1.175),
                           se = c(3.195,0.345,0.510),
                           z = c(-2.08,-1.73,2.30),
                           p = c(0.038,0.084,0.021),
                           ci.lo = c(-12.89,-1.27,0.18),
                           ci.up = c(-0.37,0.079,2.17))


interModelPlot <- ggplot(data=plotModel.df)+
  geom_pointrange(aes(x=params, y=esti,
                      ymin=ci.lo, ymax=ci.up), size=1.3, shape=20)+
  coord_flip()+
  xlab("")+
  ylab("Model Estimates")+
  ggtitle("")+
  geom_hline(yintercept = 0, color = "red", size = 1)+
  theme(plot.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 22, face = "bold"),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        text = element_text(family = "Arial"),
        plot.margin = margin(15,15,15,15))+
  annotate("text", x=1.2, y=-6.63, label = "-2.08 (0.04)", size = 5)+
  annotate("text", x=2.2, y=1.175, label = "2.30 (0.02)", size = 5)+
  annotate("text", x=3.2, y=-0.597, label = "-1.73 (0.0081)", size = 5)

interModelPlot
