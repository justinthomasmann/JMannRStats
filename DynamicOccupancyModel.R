setwd("D:/data")
library(unmarked)

SI.df <- read.csv("bawwPCountOpen.csv", h=T)

treat <- as.matrix(SI.df[,17:19])
season <- as.matrix(SI.df[,23:37])

?unmarkedFramePCO()

y.baww <- SI.df[2:16] #count of male black-and-white warblers within 100m
                      #of survey point [0,1,2] (zero inflated?)
TreatByYear <- as.matrix(SI.df[17:19]) #Site and Year specific treatment (baseline=initial occupancy, control, playback)
CovBySite <- SI.df[38:39] #Site ID and Understory Vegetation Cover 
SeasonByObs <- as.matrix(SI.df[23:37]) #season (prospect or breeding)
prim<- as.matrix(SI.df[20:22]) #primary period, year: 1,2,3 M*T

umfPCO <- unmarkedFramePCO(y = y.baww, siteCovs = CovBySite, 
                           yearlySiteCovs = list(Treat=TreatByYear),
                           obsCovs = list(Season=SeasonByObs),
                           numPrimary = 3, primaryPeriod = prim)
summary(umfPCO)
?pcountOpen

#####Null model#####
m0 <- pcountOpen(lambdaformula = ~1, #initial abundance
                 gammaformula = ~1, #recruitment rate
                 omegaformula = ~1, #survival prob
                 pformula = ~1, #detection prob
                 data = umfPCO,
                 mixture = "ZIP", #zeroinflation?
                 K = 25,
                 dynamics = "constant", #constant (no expected relationship between omega & gamma)
                 se = TRUE)
summary(m0)

m1 <- pcountOpen(lambdaformula = ~1, 
                 gammaformula = ~Treat, 
                 omegaformula = ~1, 
                 pformula = ~1, 
                 data = umfPCO,
                 mixture = "ZIP", 
                 K = 25,
                 dynamics = "trend", 
                 se = TRUE)

summary(m1)
coef(m0)

m2 <- colext(psiformula= ~Treat*site, 
             gammaformula = ~1, 
             epsilonformula = ~ 1, 
             pformula = ~site, 
             data = umf)
summary(m1)

?colext()
?unmarkedMultFrame

crossbill




data("crossbill")
colnames(crossbill)
dates <- select(crossbill, starts_with("da"))
detects <- select(crossbill, starts_with("det"))
DATE <- as.matrix(crossbill[,32:58]) 
y.cross <- as.matrix(crossbill[,5:31]) 
y.cross[is.na(DATE) != is.na(y.cross)] <- NA
sd.DATE <- sd(c(DATE), na.rm=TRUE) 
mean.DATE <- mean(DATE, na.rm=TRUE) 
DATE <- (DATE - mean.DATE) / sd.DATE
years <- as.character(1999:2007) 
years <- matrix(years, nrow(crossbill), 9, byrow=TRUE)
umf <- unmarkedMultFrame(y=y.cross, siteCovs=crossbill[,2:3], 
                         yearlySiteCovs=list(year=years), 
                         obsCovs=list(date=DATE), 
                         numPrimary=9)
summary(umf)


#####pcountOpen documentation#####
# Repeated count data with 4 primary periods and
# no 2 secondary sampling periods (ie J=2)
y2 <- matrix(c(
  0,0,  2,2,  3,2,  2,2,
  2,2,  2,1,  3,2,  1,1,
  1,0,  1,1,  0,0,  0,0,
  0,0,  0,0,  0,0,  0,0), nrow=4, ncol=8, byrow=TRUE)


# Site-specific covariates
sc2 <- data.frame(x1 = 1:4, x2 = c('A','A','B','B'))

# Observation-specific covariates
#oc2 <- list(
  x3 = matrix(1:8, nrow=4, ncol=8, byrow=TRUE)#,
  x4 = matrix(letters[1:8], nrow=4, ncol=8, byrow=TRUE)#)

# Yearly-site covariates
#ysc <- list(
  x5 = matrix(c(
    1,2,3,4,
    1,2,3,4,
    1,2,3,4,
    1,2,3,4), nrow=4, ncol=4, byrow=TRUE)#)

# Primary periods of surveys
primaryPeriod2 <- matrix(as.integer(c(
  1,2,5,7,
  1,2,3,4,
  1,2,4,5,
  1,3,5,6)), nrow=4, ncol=4, byrow=TRUE)

# Create the unmarkedFrame
umf2 <- unmarkedFramePCO(y=y2, siteCovs=sc2, obsCovs=oc2,
                         yearlySiteCovs=ysc,
                         numPrimary=4, primaryPeriod=primaryPeriod2)

# Take a look
umf2
summary(umf2)

