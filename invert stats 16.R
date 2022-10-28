setwd("C:/Users/Karin/Desktop/Research/Mann Data")
inverts16 <- read.csv("inverts 2016.csv", h=T)
inverts16

Mbunp <- mean(inverts16$totinvert[inverts16$site=="bunp"])
Mbunp

Mfss <- mean(inverts16$totinvert[inverts16$site=="fss"])
Mfss

Msgl <- mean(inverts16$totinvert[inverts16$site=="sgl"])
Msgl

Sbunp <- var(inverts16$totinvert[inverts16$site=="bunp"])
Sbunp

sqrt(Sbunp)

Sbunp <- var(inverts16$totinvert[inverts16$site=="bunp"])
Sbunp

sqrt(Sbunp)

Sbunp <- var(inverts16$totinvert[inverts16$site=="bunp"])
Ssgl

sqrt(Ssgl)

