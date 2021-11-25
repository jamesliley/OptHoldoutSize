################################################################################
## R script for functions in OptHoldoutSize package                           ##
################################################################################
##
## Sami Haidar-Wehbe, Sam Emerson, James Liley
## October 2021
##

################################################################################
## Packages and scripts                                                       ##
################################################################################

#' @import matrixStats
#' @import mnormt
#' @import mvtnorm
#' @import mle.tools

require("matrixStats")  # Matrix row- and column-wise operations
require("mnormt")       # Multivariate normal distributions
require("mvtnorm")      # Multivariate normal and t distributions
require("ranger")       # Random forests
require("mle.tools")    # Gaussian processes

################################################################################
## Small auxiliary functions                                                  ##
################################################################################

## These functions are not annotated.

logit=function(x) 1/(1+exp(-x))
logistic=function(x) -log((1/x)-1)


#### Run if manually using package
if (FALSE) {
  lx=list.files("R",full.names=T)
  for (i in 1:length(lx)) source(lx[i])
}

