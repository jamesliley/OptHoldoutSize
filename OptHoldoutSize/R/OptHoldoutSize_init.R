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

##' Logit
##'
##' @export
##' @name logit
##' @description Logit function: 1/(1+exp(-x))
##' @param x argument
##' @return value of logit(x)
##' @examples
##'
##' # Plot
##' x=seq(-5,5,length=1000)
##' plot(x,logit(x),type="l")
logit=function(x) 1/(1+exp(-x))

##' Logistic
##'
##' @export
##' @name logistic
##' @description Logistic function: -log((1/x)-1)
##' @param x argument
##' @return value of logit(x); na if x is outside (0,1)
##' @examples
##'
##' # Plot
##' x=seq(0,1,length=100)
##' plot(x,logistic(x),type="l)
##'
##' # Logit and logistic are inverses
##' x=seq(-5,5,length=1000)
##' plot(x,logistic(logit(x)),type="l)
logistic=function(x) -log((1/x)-1)


