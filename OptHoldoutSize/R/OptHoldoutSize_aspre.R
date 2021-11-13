################################################################################
## R script for functions in OptHoldoutSize package                           ##
################################################################################
##
## Sami Haidar-Wehbe, Sam Emerson, James Liley
## October 2021
##

## TODO: update with Bayesian emulation for curve estimation
## TODO: tidy up confidence interval stuff
## TODO: complete documentation of return etc, write examples for some of these


################################################################################
## Functions for ASPRE-related analyses                                       ##
################################################################################

##' Computes sensitivity of a risk score at a threshold at which
##'  10% of samples (or some proportion pi_int) are above it.
##'
##' @name sens10
##' @description Sensitivity at theshold quantile 10%
##' @keywords aspre
##' @param Y True labels
##' @param Ypred Predictions
##' @param pi Compute sensitivity when a proportion pi_int of samples exceed threshold, default 0.1
##' @return Sensitivity at this threshold
##' @examples
##'
##' # Lorem Ipsum
sens10=function(Y,Ypred,pi_int=pi_intervention) {
  qa=quantile(Ypred,1-pi_int)
  return(mean(Y*(Ypred >= qa))/pi_int)
}



##' Fit power law curve
##' Find least-squares solution: MLE of (a,b,c) under model
##'  y= ax^-b + c + e, e~N(0,s^2)
##'
##' @name powersolve
##' @description Fit power law curve
##' @keywords aspre
##' @param x x values
##' @param y y values
##' @param init initial values of (a,b,c) to start
##' @param ... further parameters passed to optim
##' @return MLE values of (a,b,c)
##' @examples
##'
##' # Lorem Ipsum
powersolve=function(x,y,init=c(1,-0.8,0.0165),...) {
  fabc=function(abc) sum( (y- (abc[1]*(x^(abc[2])) + abc[3]))^2 )
  out=optim(par=init,fn=fabc,...)
  return(out)
}


##' Generate random population of parameters where properties match that of params
##' Assume independence of parameter variation for simplicity.
##'
##' @name sim_random_aspre
##' @description Simulate random dataset similar to ASPRE training data
##' @keywords aspre
##' @param n size of population
##' @param pars list of parameters
##' @return
##' @examples
##'
##' # Load ASPRE related data
##' data(params_aspre)
##'
##' # etc
sim_random_aspre=function(n,params=params_aspre) {
  out=list()
  ind=1
  for (i in 1:length(params)) {
    if (all(names(params[[i]])==c("median","IQR"))) {
      vals=rnorm(n,mean=params[[i]]$median,sd=(params[[i]]$IQR[2]-params[[i]]$median)/qnorm(0.75))
      out[[i]]=vals
    }
    if (all(names(params[[i]])==c("name","freq"))) {
      vals=factor(sample(params[[i]]$name,n,replace=TRUE,prob=params[[i]]$freq)) # note - strings not factors
      out[[i]]=vals
    }
    if (all(names(params[[i]])=="freq")) {
      vals=sample(0:1,n,replace=TRUE,prob=c(1-params[[i]]$freq,params[[i]]$freq))
      out[[i]]=vals
    }
  }
  out=as.data.frame(out)
  colnames(out)=names(params)
  return(out)
}


##' Add various interaction terms to X. Interaction terms correspond to those in ASPRE,
##'
##' @name add_aspre_interactions
##' @keywords aspre
##' @description Adds interaction terms corresponding to ASPRE model
##' @param X data frame
##' @return New data frame containing interaction terms.
##' @examples
##'
##' # Lorem Ipsum
add_aspre_interactions=function(X) {

  # Pre-processing
  X$age_x=pmax(X$age-35,0) # age-35 if age>35, else 0
  X$height_x=X$height-164
  X$weight_ch=(X$weight-69)*(X$htn_hist==0) # If history of hypertension, do not count weight
  X$fam_hx_ch=(X$fam_hx)*(X$htn_hist==0) # If history of hypertension, do not count family history
  X$dm_hist_ch=(X$dm_hist)*(X$htn_hist==0) # If history of hypertension, do not count history of diabetes

  X$ga_prev_PE2=(X$last_delivery_gestational_age-24)^2 *(X$parity=="Par_prev_PE") # 0 if no previous PE
  X$years_no_prev_PE_1=((X$years_since_last_pregnancy)^(-1))*(X$parity=="Par_no_prev_PE")# etc
  X$years_no_prev_PE_sqrt=(pmax(X$years_since_last_pregnancy,1)^(-0.5))*(X$parity=="Par_no_prev_PE")
  X$ga_no_prev_PE2=(X$last_delivery_gestational_age-24)^2 *(X$parity=="Par_no_prev_PE") # 0 if no previous PE

  return(X)
}


##' Computes ASPRE model prediction on a matrix X of covariates
##'
##' Full ASPRE model, from https://www.nejm.org/doi/suppl/10.1056/NEJMoa1704559/suppl_file/nejmoa1704559_appendix.pdf
##' Model is to predict gestational age at PE, eg higher=lower PE risk, so coefficients are negated for model to predict PE risk
##'
##' @name aspre
##' @description Computes ASPRE score give matrix of covariates
##' @keywords aspre
##' @param X matrix, assumed to be output of sim_random_aspre with parameter params=params_aspre and transformed using add_aspre_interactions
##' @return vector of scores.
##' @examples
##'
##' # Lorem Ipsum
aspre=function(X) {

  X1pred=data.frame(
    age=X$age_x,
    height=X$height_x,
    afc=(X$race=="Afro_Caribbean"),
    sa=(X$race=="S_Asian"),
    ch=X$htn_hist,
    sle_apl=X$sle_apl_hist,
    invitro=(X$conception=="In_vitro"),
    parous_prev_PE=(X$parity=="Par_prev_PE"),
    ga_prev_PE2=X$ga_prev_PE2,
    parous_no_prev_PE=(X$parity=="Par_no_prev_PE"),
    years_no_prev_PE_1=X$years_no_prev_PE_1,
    years_no_prev_PE_sqrt=X$years_no_prev_PE_sqrt,
    ga_no_prev_PE2=X$ga_no_prev_PE2,
    weight_ch=X$weight_ch,
    fam_hx_ch=X$fam_hx_ch,
    dm_hist_ch=X$dm_hist_ch
  )

  Xxpred=data.frame(
    map=pmax(X$map,0.1),
    uterine_artery_PI=pmax(X$uterine_artery_PI,0.1),
    serum_PPA=pmax(X$serum_PPA,0.1),
    serum_PPG=pmax(X$serum_PPG,0.1)
  )

  coefs=c(-0.206886, 0.11711, -2.6786, -1.129, -7.2897,
    -3.0519, -1.6327, -8.1667, 0.0271988, -4.335, -4.15137651,
    9.21473572, 0.01549673, -0.0694096, -1.7154, -3.3899
  )

  intercept_g=54.3637

  mu_g=(as.matrix(X1pred) %*% coefs) + intercept_g
  sd_g=6.8833

  beta0_x=c(0.09564,0.54453,-0.62165,-0.93687)
  beta1_x=c(-0.001824,-0.013143,0.014692,0.021930)


  cor_x=rbind(
    c(1,        -0.05133, -0.00497, -0.02791),
    c(-0.05133,        1, -0.15992, -0.15084),
    c(-0.00497, -0.15992,        1,  0.32085),
    c(-0.02791, -0.15084,  0.32085,        1))

  sd_x=c(0.03724,  0.12894,  0.23539,  0.17723)

  var_x=cor_x*outer(sd_x,sd_x)


  # Compute density at g
  dg=function(g) dmvnorm(log10(as.matrix(Xxpred)),mean=pmax(beta0_x+ beta1_x*g,0),sigma=var_x)*
    dnorm(g,mu_g,sd_g)

  # Integrate to find result
  numerator=rep(0,dim(X)[1]); denominator=rep(0,dim(X)[1])
  for (i in 24:37) numerator=numerator + dg(i)
  for (i in 38:100) denominator=denominator + dg(i)

  score=pmin(pmax(numerator/denominator,1e-3),0.95)

  # quick imputation for missing values, unimportant
  mscore=mean(score[which(!is.na(score))])
  score[which(is.na(score))]=mscore

  return(score)

}
















# Possibly include confidence interval for OHS and cost given point estimates of k1, N, theta and an estimate
#  of sigma, assuming normality. Option of empirical or asymptotic.

##' To include
##'   Functions to find asymptotic confidence interval (given sigma estimate)
##'   Function to find bootstrap confidence interval; maybe given sampler
##'   Function to generate simulated data
##'   Functions to call simulated ASPRE data and model

