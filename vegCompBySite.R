setwd("D:/data")

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
library(rstanarm)
library(shinystan)

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



####FULL SI EXPT DATASET####

data.comp <- read.csv("soc_info_expt_complete_NoNAs.csv", h=T)
colnames(data.comp)

SI.df <- data.frame(subset(data.comp, data.comp$site!= "adummy"))

SI.df$uvRangeT <- SI.df$uvRange*0.1
SI.df$uvVarT <- SI.df$uvVar*0.01
SI.df$uvT <- SI.df$uv*10
SI.df$gvT <- SI.df$gv*10
SI.df$cvT <- SI.df$cv*10
SI.df$uvM <- SI.df$uvMetric*0.1
SI.df$iqr75T <- SI.df$iqr75*0.1
SI.df$iqr50T <- SI.df$iqr50*0.1
SI.df$medianT <- SI.df$median*0.1
colnames(SI.df)



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



####Multiple pairwise comparisons of mean understory vegetation cover####
just2017.df <- data.frame(subset(SI.df, SI.df$year=="2017"))
#uv percent
just2017.df$uvPercent <- just2017.df$uv*100

#arcsine transform percentages 
just2017.df$uvArc <- asin(sqrt(just2017.df$uvPercent/100))
just2017.df$medianArcT <- asin(sqrt(just2017.df$median/100))
just2017.df$iqr75ArcT <- asin(sqrt(just2017.df$iqr75/100))
just2017.df$maxArcT <- asin(sqrt(just2017.df$max/100))

medianArcT <- just2017.df$iqr75ArcT

TukeyHSD(aov(just2017.df$medianArcT~just2017.df$site))
kruskal.test(medianArcT~vegSite, data = just2017.df)
pairwise.wilcox.test(medianArcT, just2017.df$site, p.adjust.method = "none")
medianArcT

pairwise.wilcox.test(just2017.df$cvb, just2017.df$site, p.adjust.method = "none")
mean(just2017.df$cvb[just2017.df$site=="fss"])
mean(just2017.df$cvb[just2017.df$site=="bunp"])
mean(just2017.df$cvb[just2017.df$site=="sgla"])

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  medianArcT and just2017.df$site 
# 
#       bunp    fss 
# fss  3.7e-05 -   
# sgla 2.7e-06 0.38

#####plot uv per site####

#median
medianUvBP <- ggplot(data=just2017.df, aes(x=just2017.df$site, y=just2017.df$median, fill = just2017.df$site))+
  geom_violin()+
  geom_boxplot(fill = "white", width = 0.1)+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
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
  annotate("text", x = as.factor(unique(just2017.df$site)), y = c(52, 48, 19),
           label = c("B", "B", "A"), size = 6)
medianUvBP
pairwise.wilcox.test(just2017.df$medianArcT, just2017.df$site, p.adjust.method = "none")

#Maximum UV
maxUvBP <- ggplot(data=just2017.df, aes(x=just2017.df$site, y=just2017.df$max, fill = just2017.df$site))+
  geom_violin()+
  geom_boxplot(fill = "white", width = 0.1)+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
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
  annotate("text", x = as.factor(unique(just2017.df$site)), y = c(100, 90, 66),
           label = c("B", "B", "A"), size = 6)
maxUvBP
pairwise.wilcox.test(just2017.df$maxArcT, just2017.df$site, p.adjust.method = "none")

#Minimum UV
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



#mean
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


#iqr75
iqr75UvBP <- ggplot(data=just2017.df, aes(x=just2017.df$site, y=just2017.df$iqr75, fill = just2017.df$site))+
  geom_violin()+
  geom_boxplot(fill = "white", width = .1)+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
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
  annotate("text", x = as.factor(unique(just2017.df$site)), y = c(68, 60, 27),
           label = c("B", "B", "A"), size = 6)
iqr75UvBP
pairwise.wilcox.test(just2017.df$iqr75ArcT, just2017.df$site, p.adjust.method = "none")

#use tukey test for multiple comparisons 
TukeyHSD(aov(just2017.df$iqr75ArcT~just2017.df$site))

#####Uv Metrics Wrap####

uvMetricsWrap <- ggarrange(medianUvBP,iqr75UvBP,maxUvBP,
          ncol = 3, legend = "none")

uvMetricsWrap

annotate_figure(uvMetricsWrap,
                top = text_grob("",
                                face = "bold",
                                size = 18),
                left = text_grob("%",
                                   size = 18,
                                   face = "bold",
                                 rot = 90))



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



#####Veg var data####
veg.df <- data.frame(read.csv("vegVarComplete.csv", h=T))
tree.sglb.df <- data.frame(subset(veg.df, veg.df$treeN > 0 & veg.df$year=="2017"))
tree.df <- data.frame(subset(tree.sglb.df, tree.sglb.df$site!="sglb"))

nobuVeg.df <- data.frame(subset(tree.sglb.df, tree.sglb.df$site != "bunp"))

#pairs sgl 
pairs(~nobuVeg.df$bawwC + nobuVeg.df$bawwB + nobuVeg.df$totUrsSpp + nobuVeg.df$totUrsB + nobuVeg.df$uv + nobuVeg.df$gv + 
        nobuVeg.df$seedN + nobuVeg.df$seedH + nobuVeg.df$treeTotBA,
      diag.panel = panel.hist,
      lower.panel = panel.smooth,
      upper.panel = panel.cor)

#pairs with sglb
pairs(~tree.sglb.df$totUrsB + tree.sglb.df$totUrsSpp + tree.sglb.df$totUrsC + tree.sglb.df$bawwC + tree.sglb.df$seedL50 + 
        tree.sglb.df$seedH + tree.sglb.df$sapN + tree.sglb.df$treeTotBA,
      diag.panel = panel.hist,
      lower.panel = panel.smooth,
      upper.panel = panel.cor)

#pairs no sglb
pairs(~tree.df$totUrsC + seedN + seedH + sapN + treeN + treeH + treeBA,
      diag.panel = panel.hist,
      lower.panel = panel.smooth,
      upper.panel = panel.cor)

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

tree.df$site <- ordered(tree.df$site, levels = c("bunp", "fss", "sgla"))

library(dplyr)
library(RColorBrewer)
group_by(tree.df, tree.df$site) %>%
  summarise(
    sapCount = n(),
    sapMean = mean(sapN),
    sapSd = sd(sapN),
    sapMedian = median(sapN),
    sapIqr = IQR(sapN)
  )

######tree/seedling/sapling comparisons by site########

#####SapN BVP####
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
kruskal.test(sapL1 ~ vegSite, data = tree.df)
pairwise.wilcox.test(tree.df$sapL1, tree.df$site, p.adjust.method = "none")

####SapH BVP####
ggplot(data = tree.df, aes(x = vegSite, y = sapH, fill = vegSite))+
  geom_violin()+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
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
pairwise.wilcox.test(tree.df$sapH, tree.df$site, p.adjust.method = "none")
#No significant differences



####SeedN BVP####
tree.df$seedN

seed.n.bp <- ggplot(data = tree.df, aes(x = vegSite, y = seedN, fill = vegSite))+
  geom_violin()+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
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
        legend.text = element_text(size = 12))+
  annotate("text", x = as.factor(unique(tree.df$site)), y = c(1795, 1200,1075), 
           label = c("B","B","A"), size = 6)
seed.n.bp

#Pairwise comparisons
kruskal.test(seedN ~ vegSite, data = tree.df)
pairwise.wilcox.test(tree.df$seedN, tree.df$site, p.adjust.method = "none")
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  tree.df$seedN and tree.df$site 
# 
#      bunp  fss  
# fss  0.022 -    
# sgla 0.022 0.937
# 
# P value adjustment method: none 

is.numeric(tree.df$site)
is.integer(tree.df$site)
####SeedL50 BVP####
# compute lower and upper whiskers for each group
ylims <- seedL50 %>%
  group_by(tree.df$site) %>%
  summarise(Q1 = quantile(Result, 1/4), Q3 = quantile(Result, 3/4)) %>%
  ungroup() %>%
  #get lowest Q1 and highest Q3
  summarise(lowQ1 = min(Q1), highQ3 = max(Q3))

plot + coord_cartesian(ylim = as.numeric(ylims)*1.05)

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
seed.h.bp <- ggplot(data = tree.df, aes(x = vegSite, y = seedH, fill = vegSite))+
  geom_violin()+
  scale_fill_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                    name = "Site",
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
        legend.text = element_text(size = 12))+
  annotate("text", x = as.factor(unique(tree.df$site)), y = c(2.1, 1.95, 1.75),
           label = c("B", "A", "A"), size = 6)
seed.h.bp

#Pairwise comparisons
kruskal.test(seedH ~ vegSite, data = tree.df)
pairwise.wilcox.test(tree.df$seedH, tree.df$site, p.adjust.method = "none")
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  tree.df$seedH and tree.df$site 
# 
#        bunp   fss   
# fss  0.4452 -     
# sgla 0.0047 0.0260
# 
# P value adjustment method: none 




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
pairwise.wilcox.test(tree.df$treeN, tree.df$site, p.adjust.method = "none")
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  tree.df$treeN and tree.df$site 
# 
#       bunp  fss  
# fss  0.045 -    
# sgla 0.617 0.024
# 
# P value adjustment method: none 




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
pairwise.wilcox.test(tree.df$treeH, tree.df$site, p.adjust.method = "none")
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
pairwise.wilcox.test(tree.df$treeBA, tree.df$site, p.adjust.method = "none")
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

sapLP <- ggplot(sap2.df, aes(x = sap2.df$size, 
                y = sap2.df$mean,
                group = sap2.df$site,
                color = sap2.df$site))+
  
  geom_line(position = position_dodge(0.3), size = 1.3)+
  geom_point(position = position_dodge(0.3), size = 3)+
  geom_errorbar(aes(ymin = sap2.df$mean-sap2.df$sd, ymax = sap2.df$mean+sap2.df$sd), width = .2,
                position = position_dodge(0.3))+
  scale_color_manual(values = c("#FF7F00","#377EB8","#4DAF4A"),
                     name = "Site")+
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
              legend.position = c(0.8,0.8))+
  annotate("text", x = c(0.9,1,1.1,1.9,2,2.1,2.9,3,3.1,3.9,4,4.1), y = c(12,127,97,68,145,95,80,63,41,30,23,21),
           label = c("A","B","B","A","A","A","A","AB","B","A","AB","B"), size = 4)+
  annotate("rect", xmin = 0.8, xmax = 1.2, ymin=-5, ymax = 150, alpha = 0.15)+
  annotate("rect", xmin = 1.8, xmax = 2.2, ymin=-5, ymax = 150, alpha = 0.15)+
  annotate("rect", xmin = 2.8, xmax = 3.2, ymin=-5, ymax = 150, alpha = 0.15)+
  annotate("rect", xmin = 3.8, xmax = 4.2, ymin=-5, ymax = 150, alpha = 0.15)


#Pairwise comparisons
kruskal.test(sapL1 ~ vegSite, data = tree.df)
pairwise.wilcox.test(tree.df$sapL1, tree.df$site, p.adjust.method = "none")

#Pairwise comparisons
kruskal.test(sap1.3 ~ vegSite, data = tree.df)
pairwise.wilcox.test(tree.df$sap1.3, tree.df$site, p.adjust.method = "none")

#Pairwise comparisons
kruskal.test(sap3.6 ~ vegSite, data = tree.df)
pairwise.wilcox.test(tree.df$sap3.6, tree.df$site, p.adjust.method = "none")

#Pairwise comparisons
kruskal.test(sap6.10 ~ vegSite, data = tree.df)
pairwise.wilcox.test(tree.df$sap6.10, tree.df$site, p.adjust.method = "none")
