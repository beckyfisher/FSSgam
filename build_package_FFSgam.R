

# This script includes the basic steps required to create and
#    build the FFSgam package

library(devtools)
library(roxygen2)
library(knitr)
library(R.rsp)

setwd("C:/Users/rfisher/OneDrive - Australian Institute of Marine Science/Documents/AIMS/EcologicalRiskModelling/Full_subsets_methods/FSSgam_package")
# copy function files over to FFSgam/R

devtools::document()

devtools::use_package("doSNOW")
devtools::use_package("MuMIn")
devtools::use_package("gamm4")
devtools::use_package("mgcv")
devtools::use_package("nnet")
devtools::build()

#devtools::install_github("beckyfisher/FSSgam_package")
