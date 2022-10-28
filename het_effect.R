setwd("C:/Users/Karin/Desktop/Data")
library(vegan)
library(ggplot2)

?boxplot()

#het effects without focal species
hets.15 <- read.csv("bird_hets_2015.csv", h=T)
hets.16 <- read.csv("bird_hets_2016.csv", h=T)
colnames(hets.16)
colnames(hets.15)

#het effects with focal species
hets.16.wfs <- read.csv("bird_hets_2016_wfs.csv", h=T)
colnames(hets.16.wfs)


#data frames to analyze het effects in fss and sgl only without focal species
hets.16.fss.sgl <- data.frame(subset(hets.16, hets.16$site!= "bunp"))
hets.15.fss.sgl <- data.frame(subset(hets.15, hets.15$site!= "bunp"))
#with focal species
hets.16.fss.sgl.wfs <- data.frame(subset(hets.16.wfs, hets.16.wfs$site!= "bunp"))

#het effects: parulidae only @ sgl and fss
hets.15.par <- hets.15.fss.sgl[,9:15]
hets.16.par <- hets.16.fss.sgl[,9:15]

hets.15.par[is.na(hets.15.par)] <- 0
hets.16.par[is.na(hets.16.par)] <- 0
diversity(hets.15.par, index = "shannon")
diversity(hets.16.par, index = "shannon")

t.test(diversity(hets.15.par, index = "shannon"), diversity(hets.16.par, index = "shannon"))

#without focal species change, experimental site comparison by year of total understory reliant individuals 
hets.16.fss.sgl$total[hets.16.fss.sgl$treat=="e"]
hets.15.fss.sgl$total[hets.15.fss.sgl$treat=="e"]
#with focal species change
hets.16.fss.sgl.wfs$total[hets.16.fss.sgl.wfs$treat=="e"]



#t test between yearly turs WITHOUT focal species change at experimental sites
var.test(hets.15.fss.sgl$total[hets.15.fss.sgl$treat=="e"], 
         hets.16.fss.sgl$total[hets.16.fss.sgl$treat=="e"], 
         alternative = "two.sided", conf.level = 0.95)

t.test(hets.15.fss.sgl$total[hets.15.fss.sgl$treat=="e"], 
       hets.16.fss.sgl$total[hets.16.fss.sgl$treat=="e"],
       alternatve="two.sided",var.equal=FALSE, conf.level=0.95)

boxplot(hets.15.fss.sgl$total[hets.15.fss.sgl$treat=="e"], 
        hets.16.fss.sgl$total[hets.16.fss.sgl$treat=="e"],
        names = c("Pretreatment", "Treatment"),
        col = c("grey", "purple"),
        main = "Experimental sites: total understory-reliant individuals",
        ylab = "Number of indivduals",
        xlab = "Year")

#t test between yearly turs WITH focal species change at experimental sites
var.test(hets.15.fss.sgl$total[hets.15.fss.sgl$treat=="e"], 
         hets.16.fss.sgl.wfs$total[hets.16.fss.sgl.wfs$treat=="e"], 
         alternative = "two.sided", conf.level = 0.95)

t.test(hets.15.fss.sgl$total[hets.15.fss.sgl$treat=="e"], 
       hets.16.fss.sgl.wfs$total[hets.16.fss.sgl.wfs$treat=="e"], 
       alternative = "two.sided", conf.level = 0.95)

boxplot(hets.15.fss.sgl$total[hets.15.fss.sgl$treat=="e"], 
        hets.16.fss.sgl.wfs$total[hets.16.fss.sgl.wfs$treat=="e"],
        names = c("Pretreatment", "Treatment"),
        col = c("grey", "purple"),
        main = "Experimental sites: total understory-reliant individuals (including fsc) ",
        ylab = "Number of indivduals",
        xlab = "Year")


#control site comparison by year of total understory reliant individuals
hets.16.fss.sgl$total[hets.16.fss.sgl$treat=="c"]
hets.15.fss.sgl$total[hets.15.fss.sgl$treat=="c"]

var.test(hets.15.fss.sgl$total[hets.15.fss.sgl$treat=="c"], 
         hets.16.fss.sgl$total[hets.16.fss.sgl$treat=="c"],
         alternative = "two.sided", conf.level = 0.95)
t.test(hets.15.fss.sgl$total[hets.15.fss.sgl$treat=="c"], 
       hets.16.fss.sgl$total[hets.16.fss.sgl$treat=="c"])

boxplot(hets.15.fss.sgl$total[hets.15.fss.sgl$treat=="c"], hets.16.fss.sgl$total[hets.16.fss.sgl$treat=="c"],
        names = c("Pretreatment", "Treatment"),
        col = c("grey", "purple"),
        main = "Control sites: total understory-reliant individuals",
        ylab = "Number of indivduals",
        xlab = "Year")



#comparison by year of urs richness at experimental sites
hets.16.fss.sgl$richness[hets.16.fss.sgl$treat=="e"]
hets.15.fss.sgl$richness[hets.15.fss.sgl$treat=="e"]

#t test between yearly urs richness
var.test(hets.15.fss.sgl$richness[hets.15.fss.sgl$treat=="e"], 
         hets.16.fss.sgl$richness[hets.16.fss.sgl$treat=="e"],
         alternative = "two.sided", conf.level = 0.95)
t.test(hets.15.fss.sgl$richness[hets.15.fss.sgl$treat=="e"], 
       hets.16.fss.sgl$richness[hets.16.fss.sgl$treat=="e"])

boxplot(hets.15.fss.sgl$richness[hets.15.fss.sgl$treat=="e"], 
        hets.16.fss.sgl$richness[hets.16.fss.sgl$treat=="e"],
        names = c("2015", "2016"),
        col = c("grey", "purple"),
        main = "Understory-reliant species richness at experimental sites",
        ylab = "Number of understory-relaint species",
        xlab = "Year")

#t test between yearly urs richness
var.test(hets.15.fss.sgl$richness[hets.15.fss.sgl$treat=="c"], 
         hets.16.fss.sgl$richness[hets.16.fss.sgl$treat=="c"],
         alternative = "two.sided", conf.level = 0.95)
t.test(hets.15.fss.sgl$richness[hets.15.fss.sgl$treat=="c"], 
       hets.16.fss.sgl$richness[hets.16.fss.sgl$treat=="c"])

boxplot(hets.15.fss.sgl$richness[hets.15.fss.sgl$treat=="c"], 
        hets.16.fss.sgl$richness[hets.16.fss.sgl$treat=="c"],
        names = c("2015", "2016"),
        col = c("grey", "purple"),
        main = "Understory-reliant species richness at control sites",
        ylab = "Number of snderstory-reliant species",
        xlab = "Year")


#Year comparison of Shannon diversity in URS 
hets.15[is.na(hets.15)] <- 0
div.15 <- diversity(hets.15[1:38,4:18], index = "shannon")
div.15

hets.16[is.na(hets.16)] <- 0 
div.16 <- diversity(hets.16[1:38,4:18], index = "shannon")
div.16


var.test(div.15,div.16, alternative = "two.sided", conf.level = 0.95)
t.test(div.15,div.16)

boxplot(div.15, div.16,
        names = c("Pretreatment", "Treatment"),
        col = c("grey", "magenta"),
        main = "Comparison by year of Shannon diversity",
        ylab = "H'",
        xlab = "Year")

div.15.fss <- div.15[hets.15$site == "fss"]  
div.15.fss
div.16.fss <- div.16[hets.16$site=="fss"]
div.16.fss


boxplot(div.15.fss, div.16.fss)


