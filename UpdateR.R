#Run in R gui 
library(installr)
updateR()

#transfer libraries from old versions to new version
#change the end of the next line to the version of R 
#from which you want to move packages
lib_loc <- "C:/Users/Justin/Documents/R/win-library/4.0.2"
to_install <- unname(installed.packages(lib.loc = lib_loc)[, "Package"])
to_install
install.packages(pkgs = to_install)
