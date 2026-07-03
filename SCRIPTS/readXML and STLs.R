
# initialise --------------------------------------------------------------


rm(list = ls())


library(morphomap)
require(rgl)
require(Rvcg)
require(Morpho)


# readfiles ---------------------------------------------------------------

source("./SCRIPTS/analysis functions.R")


extfiles <- list.files("./STLs", pattern = "ext")
intfiles <- list.files("./STLs", pattern = "int")
