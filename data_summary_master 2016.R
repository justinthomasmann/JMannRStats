if(!require(installr)) {
  install.packages("installr"); require(installr)}
updateR()


setwd("C:/Users/Karin/Desktop/Data")
datamaster <- read.csv("data_summary_master.csv", h=T)

colnames(datamaster)

#contains logistic regression 
library(nlme)
#contains AICc 
library(AICcmodavg)
#Summary statistics
library(plyr)
library(ggplot2)

barplot(datamaster$uv15,
        ylim = c(0,80),
        main = "Understory Vegetation Density",
        ylab = "%",
        col = ifelse(datamaster$site== "sgl", "blue", 
                     ifelse(datamaster$site=="fss", "green", "orange")))
legend("topright", c("sgl", "fss", "bunp"), col = c("blue", "green", "orange"), pch = 15)


plot(datamaster$site,datamaster$uv15, 
     ylim = c(0,80),
     col = c("orange","green", "blue"),
     main="Understory Vegetation Density", 
     xlab = "Site", 
     ylab = "%")
     
plot(datamaster$site,datamaster$turs15, 
     col = c("orange","green", "blue"), 
     main = "Understory-reliant Birds", 
     xlab = "Site", 
     ylab = "Number of understory-reliant individuals")


plot(datamaster$uv15,datamaster$turs15,
     pch = 19,
     main = "Understory-reliant Birds ~ Understory Vegetation Density",
     xlab = "Understory vegetation density (%)",
     ylab = "Total number of understory-reliant individuals",
     col = ifelse(datamaster$site== "sgl", "blue", 
                  ifelse(datamaster$site=="fss", "green", "orange")))
abline(lm(datamaster$turs15~datamaster$uv15))
legend(60, 5, c("bunp", "fss", "sgl"), col = c("orange", "green", "blue", "black"), pch = 19)
text(60,10, labels = "r^2 = 0.5255") 

uvturs <- lm(datamaster$turs15~datamaster$uv15)
summary(uvturs)

#model variables: y=binary response variable (focal species increase)
# xi=independent variables (1,5 = factors)
y <- datamaster$fsi
y
length(y)

x1 <- factor(datamaster$treat)
x1

x2 <- datamaster$uv15
x2

x3 <- datamaster$gv15
x3

x4 <- datamaster$tot_invert16
x4

x5 <- factor(datamaster$site)
x5

#logistic regression uses Bernoulli distribution
#full model - focal species increase across all sites

#Hypotheses: 
#Null: Beta = 0
#H1: Beta does not = 0


FM <- glm(y~x1+x2+x3+x4+x5,data=datamaster,family=quasipoisson)
summary(FM)

#Test for dispersion of data

FM.disp.test <- glm(y~x1+x2+x3+x4+x5,data=datamaster,family=quasibinomial)
summary(FM) 

pchisq(summary(FM.disp.test)$dispersion * FM$df.residual, 
       FM$df.residual, lower = F)

#Dispersion ratio is < 1, indicating underdispesion.

#Likelihood Ratio Test:
anova(FM, test="Chisq")

step(FM,direction = "backward")

drop1(FM,test="Chisq")


###Predictive ability of the model


#ROC curve (receiver operating characteristic)
#Measures the predictive performance of the model
#The closer the auc is to one, the better the predictive ability of the model.
library(ROCR)
p <- predict(FM, newdata=subset(datamaster,select=c(2,3,4,5,6,7,8)), type="response")
pr <- prediction(p, y)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc




###Model validations!!!
library(caret)

#Likelihood ratio test compares the performance of models,
#one with the variable of interest and one without it.
FM <- glm(y~x1+x2+x3+x4+x5,data=datamaster,family=binomial)
summary(FM)
RM <- glm(y~x2+x3+x4+x5,data=datamaster,family=binomial)
summary(RM)

anova(FM,RM, test = "Chisq")
# The addition of x1 to the model reduces the residual deviance from 
#39.6 to 27.8; therefore x1 explains a significant proportion of the variation

#Alone, x1 explains 30% of the variation in the data.
(1-27.8/39.6)*100

#variable importance 

varImp(FM)




#k-fold cross validation
Train <- createDataPartition(datamaster$fsi, p=0.6, list=FALSE)
training <- datamaster[ Train, ]
testing <- datamaster[ -Train, ]

ctrl <- trainControl(method = "repeatedcv", 
                     number = 10, 
                     savePredictions = TRUE)

mod_fit <- train(y~x1+x2+x3+x4+x5,data=datamaster, 
                 method="glm", 
                 family="binomial",
                 trControl = ctrl, 
                 tuneLength = 5)

pred = predict(mod_fit, newdata=testing)
confusionMatrix(data=pred, testing$Class)




#AICc candidate models
cand.set<-list()
cand.set[[1]]<- glm(y~x1,data = datamaster,family = binomial)
cand.set[[2]]<- glm(y~x1+x2+x3+x4+x5,data=datamaster,family=binomial)
cand.set[[3]]<- glm(y~x1+x5,data=datamaster,family=binomial)
cand.set[[4]]<- glm(y~x2+x3,data=datamaster,family=binomial)
cand.set[[5]]<- glm(y~x2+x3+x4, data = datamaster,family = binomial)
cand.set[[6]]<- glm(y~x1+x2+x5, data = datamaster,family = binomial)


#assign model names for table
Modnames.b <- c("treat", "full", "treat+site", "uv+gv", "uv+gv+totinvert","treat+UV+site")

#create AICc table
aicc_table <- aictab(cand.set, modnames = Modnames.b, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1)
aicc_table

evidence(aic.table = aicc_table)

#Evidence ratio between models 'treat+UV+site' and 'treat+site':
#  1.28 which means the top model is 1.28 times more parsimonious than the second ranked model

#candidate model #2: includes interaction between site and treat
cand.set2<-list()
cand.set2[[1]]<- lm(y~x1,data = datamaster)
cand.set2[[2]]<- lm(y~x1+x2+x3+x4+x5,data=datamaster)
cand.set2[[3]]<- lm(y~x1*x5,data=datamaster)
cand.set2[[4]]<- lm(y~x2+x3,data=datamaster)
cand.set2[[5]]<- lm(y~x2+x3+x4, data = datamaster)
cand.set2[[6]]<- lm(y~x1+x2, data = datamaster)

Modnames.2 <- c("treat", "full", "treat*site", "uv+gv", "uv+gv+totinvert","treat+uv")


aicc_table2 <- aictab(cand.set2, modnames = Modnames.2, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1)
aicc_table2

evidence(aic.table = aicc_table2)







#######summary statistics 
#assign variables 

fsi_bunp<-sum(datamaster$fsi[datamaster$site=="bunp"])
fsi_fss<-sum(datamaster$fsi[datamaster$site=="fss"])
fsi_sgl<-sum(datamaster$fsi[datamaster$site=="sgl"])-1

fsi_bunp
fsi_fss
fsi_sgl

fsi_bunp/7
fsi_fss/7
fsi_sgl/9

mean(datamaster$uv15[datamaster$site=="bunp"])
sqrt(var(datamaster$uv15[datamaster$site=="bunp"]))
mean(datamaster$uv15[datamaster$site=="fss"])
sqrt(var(datamaster$uv15[datamaster$site=="fss"]))
mean(datamaster$uv15[datamaster$site=="sgl"])
sqrt(var(datamaster$uv15[datamaster$site=="sgl"]))









########### Generalized Linear Mixed Effects Models

#experiment here with model changes!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#PREEBREEDING

#new y, includes continuous data on any (+,-) PREBREEDING change 
#in focal species presence
y.fscp <- datamaster$fscpt
y.fscp

x1 <- factor(datamaster$treat)
x1

x2 <- log(datamaster$uv15+1)
x2

x3 <- log(datamaster$gv15+1)
x3

x4 <- log(datamaster$tot_invert16+1)
x4

x5 <- factor(datamaster$site)
x5

require(scatterplot3d)

threeD <- scatterplot3d(y.fscp, x2, x1,
                        pch = 16,
                        type = "h",
                        highlight.3d = TRUE)
lm3d <- lm(y.fscp~ x1 + x2)
threeD$plane3d(lm3d, lty.box = "solid")

plot(y.fscp, x2)

#no 1*2 interaction

lm.noint <- lm(y.fscp~x1+x2+x3+x4+x5, data = datamaster)
summary(lm.noint)
anova(lm.noint)

resid_lm.noint <- resid(lm.noint)
hist(resid_lm.noint, xlab = "Residuals")

# 1*2 interaction

lm.int <- lm(y.fscp~x1*x2+x3+x4+x5, data = datamaster)
summary(lm.int)
anova(lm.int)

resid_lm.int <- resid(lm.int)
hist(resid_lm.int, xlab = "Residuals")

drop1(lm.int, test = "Chisq")
step(lm.int, direction = "backward")

cand.set4<-list()
cand.set4[[1]]<- lm(y.fscp~x1,data = datamaster)
cand.set4[[2]]<- lm(y.fscp~x1+x2+x3+x4+x5,data=datamaster)
cand.set4[[3]]<- lm(y.fscp~x1*x2,data=datamaster)
cand.set4[[4]]<- lm(y.fscp~x2+x3,data=datamaster)
cand.set4[[5]]<- lm(y.fscp~x2+x3+x4, data = datamaster)
Modnames.4 <- c("treat", "full", "treat*uv", "uv+gv", "uv+gv+totinvert")
aicc_table4 <- aictab(cand.set4, modnames = Modnames.4, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1)
aicc_table4

#interaction between x1 and x5
lm.int.mod <- lm(y.fscp~x1*x5+x2+x3+x4, data = datamaster)
summary(lm.int.mod)
anova(lm.int.mod, test = "chisq")

drop1(lm.int.mod, test = "Chisq")
step(lm.int.mod, direction = "backward")

#only treatment and treatment by site interaction prebreeding change
treatmod <- lm(y.fscp~x1*x5)
summary(treatmod)
anova(treatmod)

#new y, includes continuous data on any (+,-) BREEDING change 
#in focal species presence
y.fscb <- datamaster$fscb
y.fscb

lm.breeding.mod <- lm(y.fscb~x1+x2+x3+x4+x5, data = datamaster)
summary(lm.breeding.mod)
anova(lm.breeding.mod)

#breeding model, just treatment by site interaction
lm.breeding.int.mod <- lm(y.fscb~x1*x5)
summary(lm.breeding.int.mod)
anova(lm.breeding.int.mod)

#No BU
datamaster_nobu <- data.frame(subset(datamaster, datamaster$site != "bunp"))
datamaster_nobu$site <- factor(datamaster_nobu$site)
y.fscp_nobu <- datamaster_nobu$fscpt
x1_nobu <- factor(datamaster_nobu$treat)
x2_nobu <- datamaster_nobu$uv15
x3_nobu <- datamaster_nobu$gv15
x4_nobu <- datamaster_nobu$tot_invert16
x5_nobu <- factor(datamaster_nobu$site)

lm.nobu <- lm(y.fscp_nobu~x1_nobu*x2_nobu+x3_nobu+x4_nobu+x5_nobu, data = datamaster_nobu)
summary(lm.nobu)
anova(lm.nobu)

resid(lm.nobu)
hist(resid(lm.nobu))

#################################################


#Just BAWW
datamaster_baww <- data.frame(subset(datamaster, datamaster$site != "sgl"))
datamaster_baww$site <- factor(datamaster_baww$site)
datamaster_baww$site
y.fscp_baww <- (datamaster_baww$fscpt)
x1_baww <- factor(datamaster_baww$treat)
x2_baww <- log(datamaster_baww$uv15+1)
x3_baww <- log(datamaster_baww$gv15+1)
x4_baww <- log(datamaster_baww$tot_invert16+1)
x5_baww <- factor(datamaster_baww$site)

#just baww no interaction
lm.baww <- lm(y.fscp_baww~x1_baww+x2_baww+x3_baww+x4_baww+x5_baww, data = datamaster_baww)
summary(lm.baww)
anova(lm.baww)
hist(resid(lm.baww))
model.matrix(lm.baww)

#Likelihood Ratio Test:
anova(lm.baww)
step(lm.baww, direction = "backward")


#just baww interaction between site and treat
lm.baww.int <- lm(y.fscp_baww~x1_baww*x2_baww+x3_baww+x4_baww+x5_baww, data = datamaster_baww)
summary(lm.baww.int)
#Likelihood Ratio Test:
anova(lm.baww.int, test="Chisq")
step(lm.baww.int, direction = "backward")

hist(resid(lm.baww.int))
#baww aic with interaction

cand.set3<-list()
cand.set3[[1]]<- lm(y.fscp_baww~x1_baww,data = datamaster_baww)
cand.set3[[2]]<- lm(y.fscp_baww~x1_baww+x2_baww+x3_baww+x4_baww+x5_baww,data=datamaster_baww)
cand.set3[[3]]<- lm(y.fscp_baww~x1_baww*x5_baww,data=datamaster_baww)
cand.set3[[4]]<- lm(y.fscp_baww~x2_baww+x3_baww,data=datamaster_baww)
cand.set3[[5]]<- lm(y.fscp_baww~x2_baww+x3_baww+x4_baww, data = datamaster_baww)
cand.set3[[6]]<- lm(y.fscp_baww~x1_baww+x2_baww, data = datamaster_baww)

Modnames.3 <- c("treat", "full", "treat*site", "uv+gv", "uv+gv+totinvert","treat+uv")


aicc_table3 <- aictab(cand.set3, modnames = Modnames.3, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1)
aicc_table3

evidence(aic.table = aicc_table3)

#baww aic without interaction

cand.set4<-list()
cand.set4[[1]]<- lm(y.fscp_baww~x1_baww,data = datamaster_baww)
cand.set4[[2]]<- lm(y.fscp_baww~x1_baww+x2_baww+x3_baww+x4_baww+x5_baww,data=datamaster_baww)
cand.set4[[3]]<- lm(y.fscp_baww~x2_baww+x3_baww,data=datamaster_baww)
cand.set4[[4]]<- lm(y.fscp_baww~x2_baww+x3_baww+x4_baww, data = datamaster_baww)
cand.set4[[5]]<- lm(y.fscp_baww~x1_baww+x2_baww, data = datamaster_baww)

Modnames.4 <- c("treat", "full", "uv+gv", "uv+gv+totinvert","treat+uv")


aicc_table4 <- aictab(cand.set4, modnames = Modnames.4, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1)
aicc_table4

evidence(aic.table = aicc_table4)

#just BAWW logistic
y_fsi <- datamaster_baww$fsi
logistic_baww<- glm(y_fsi~x1_baww+x2_baww+x3_baww+x4_baww+x5_baww,data=datamaster_baww,family=binomial)
summary(logistic_baww)

anova(logistic_baww, test="Chisq")
step(logistic_baww, direction = "backward")

cand.set5<-list()
cand.set5[[1]]<- lm(y_fsi~x1_baww,data = datamaster_baww)
cand.set5[[2]]<- lm(y_fsi~x1_baww+x2_baww+x3_baww+x4_baww+x5_baww,data=datamaster_baww)
cand.set5[[3]]<- lm(y_fsi~x1_baww*x5_baww,data=datamaster_baww)
cand.set5[[4]]<- lm(y_fsi~x2_baww+x3_baww,data=datamaster_baww)
cand.set5[[5]]<- lm(y_fsi~x2_baww+x3_baww+x4_baww, data = datamaster_baww)
cand.set5[[6]]<- lm(y_fsi~x1_baww+x2_baww, data = datamaster_baww)

Modnames.5 <- c("treat", "full", "treat*site", "uv+gv", "uv+gv+totinvert","treat+uv")


aicc_table5 <- aictab(cand.set5, modnames = Modnames.5, second.ord = TRUE, nobs = NULL, sort = TRUE, c.hat = 1)
aicc_table5

evidence(aic.table = aicc_table5)


plot(x1_baww,y.fscp_baww)

#box plot of UV at BUNP and FSS
plot(datamaster_baww$site, datamaster_baww$uv15, ylab = "Understory Vegetation Density (%)", xlab = "Site")
#box plot of fscp comparison by site (BUNP and FSS)
plot(datamaster_baww$site, datamaster_baww$fscp)

fscp.fss <- datamaster$fscp[datamaster$site=="fss"]
fscp.fss
fss.e.change <- fscp.fss[datamaster$treat=="e"]
fss.e.change
mean(fss.e.change)

#Non-parametric pairwise comparisons

pairwise.wilcox.test(datamaster$uv15,datamaster$site,alternative="two.sided")

pairwise.wilcox.test(datamaster$turs15,datamaster$site,alternative="two.sided")

pairwise.wilcox.test(datamaster$cv15,datamaster$site,alternative="two.sided")


uv.turs.reg <- lm(datamaster$uv15~datamaster$turs15)
summary(uv.turs.reg)

lm1 <- lm(datamaster$uv15~datamaster$site)
summary(lm1)
lm2<- lm(datamaster$fscp~datamaster$treat)
summary(lm2)
lm3<-lm(datamaster$fscp~datamaster$treat*datamaster$uv15)
summary(lm3)

plot(datamaster$treat, datamaster$fscp, xlab = "Treatment", ylab = "Focal Species Change (prebreeding)")

plot(x5,y.fscp)
#Ideas and questions:

anova(lm(y.fscp~x1*x5))
anova(lm(y.fscp~x1%in%x5))
TukeyHSD(aov(y.fscp~x5*x1))
pairwise.t.test(y.fscp, x5)

t.test(datamaster$fsi[datamaster$site=="fss"], datamaster$fsi[datamaster$site=="bunp"])

###### Test Models ######

lm.baww.int1 <- lm(y.fscp_baww~x1_baww+x1_baww*x5_baww, data = datamaster_baww)
summary(lm.baww.int1)
#Likelihood Ratio Test:
anova(lm.baww.int1, test="Chisq")
step(lm.baww.int1, direction = "backward")


#### Model validation 
resid_baww <- resid(lm.baww)
resid_baww
hist(resid_baww, xlab = "Residuals",
     col = ifelse(datamaster$site=="fss", "green", "orange"))
plot(resid_baww,
     col = ifelse(datamaster$site=="fss", "green", "orange"))

resid_FM <- resid(FM)
hist(resid_FM, xlab = "Residuals")

resid_full.lm <- resid(lm.mod)
hist(resid_full.lm, xlab = "Residuals")

#Total binary response
#e
17/23
#c
2/15
#FSS.p.response
6/7
#BUNP.p.response
2/7
#HOWA.p.response
4/6
#cawa.p.response
3/3
