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

##' Sensitivity at theshold quantile 10%
##'
##'
##' @name sens10
##' @description Computes sensitivity of a risk score at a threshold at which 10% of samples (or some proportion pi_int) are above the threshold.
##' @keywords aspre
##' @param Y True labels (1 or 0)
##' @param Ypred Predictions (univariate; real numbers)
##' @param pi_int Compute sensitivity when a proportion pi_int of samples exceed threshold, default 0.1
##' @return Sensitivity at this threshold
##' @examples
##'
##' # Simulate
##' set.seed(32142)
##'
##' N=1000
##' X=rnorm(N); Y=rbinom(N,1,prob=logit(X/2))
##'
##' pi_int=0.1
##' q10=quantile(X,1-pi_int) # 10% of X values are above this threshold
##'
##' print(length(which(Y==1 & X>q10))/length(which(X>q10)))
##' print(sens10(Y,X,pi_int))
##'
sens10=function(Y,Ypred,pi_int=0.1) {
  qa=quantile(Ypred,1-pi_int)
  return(mean(Y*(Ypred >= qa))/pi_int)
}



##' Fit power law curve
##'
##' @name powersolve
##' @description  Find least-squares solution: MLE of (a,b,c) under model
##'  ``y_i= a x_i^-b + c + e_i; e_i~N(0,s^2 y_var_i^2)``
##' @keywords aspre
##' @param x X values
##' @param y Y values
##' @param init Initial values of (a,b,c) to start. Default c(20000,2,0.1)
##' @param y_var Optional parameter giving sampling variance of each y value. Defaults to 1.
##' @param estimate_s Parameter specifying whether to also estimate s (as above). Defaults to FALSE (no).
##' @param ... further parameters passed to optim. We suggest specifying lower and upper bounds for (a,b,c); e.g. lower=c(1,0,0),upper=c(10000,3,1)
##' @return MLE values of (a,b,c)
##' @examples
##'
##' # Retrieval of original values
##' A_true=2000; B_true=1.5; C_true=0.3; sigma=0.002
##'
##' X=1000*abs(rnorm(10000,mean=4))
##' Y=A_true*(X^(-B_true)) + C_true + rnorm(length(X),sd=sigma)
##'
##' c(A_true,B_true,C_true)
##' powersolve(X[1:10],Y[1:10])$par
##' powersolve(X[1:100],Y[1:100])$par
##' powersolve(X[1:1000],Y[1:1000])$par
##' powersolve(X[1:10000],Y[1:10000])$par
powersolve=function(x,y,init=c(20000,2,0.1),y_var=rep(1,length(y)),estimate_s=FALSE,...) {
  sc=mean(x)^2 # scale by this to avoid overflow errors
  if (!estimate_s) {
    fabc=function(abc) {
      out=sum( ((y- (abc[1]*(x^(-abc[2])) + abc[3]))^2)/(y_var))
      if (is.finite(out)) return(out) else return(1e10)
    }
    out=suppressWarnings(optim(par=init,fn=fabc,...))
  } else {
    fabcs=function(abcs) {
      out=-(sum( -((y- (abcs[1]*(x^(-abcs[2])) + abcs[3]))^2 / (2*y_var*(abcs[4]^2))) - log(sqrt(2*3.1415*y_var)*abcs[4])))
      if (is.finite(out)) return(out) else return(1e10)
    }
    out=suppressWarnings(optim(par=c(init,0.05),fn=fabcs,...))
  }
  return(out)
}



##' Standard error matrix for learning curve parameters (power law)
##'
##'
##' @name powersolve_se
##' @description Find approximate standard error matrix for ``(a,b,c)`` under power law model for learning curve.
##'
##' Assumes that
##'
##'   ``y_i= a x_i^-b + c + e, e~N(0,s^2 y_var_i^2)``
##'
##' Standard error can be computed either asymptotically using Fisher information (`method='fisher'`) or boostrapped (`method='bootstrap'`)
##'
##' These estimate different quantities: the asymptotic method estimates
##'
##' ``Var[MLE(a,b,c)|X,y_var]``
##'
##' and the boostrap method estimates
##'
##' ``Var[MLE(a,b,c)]``.
##'
##' @keywords aspre
##' @param x X values (typically training set sizes)
##' @param y Y values (typically observed cost per individual/sample)
##' @param method One of 'fisher' (for asymptotic variance via Fisher Information) or 'bootstrap' (for Bootstrap)
##' @param init Initial values of (a,b,c) to start when computing MLE. Default c(20000,2,0.1)
##' @param y_var Optional parameter giving sampling variance of each y value. Defaults to 1.
##' @param n_boot Number of bootstrap resamples. Only used if method='bootstrap'. Defaults to 1000
##' @param seed Random seed for bootstrap resamples. Defaults to NULL.
##' @param ... further parameters passed to optim. We suggest specifying lower and upper bounds; since optim is called on (a*1000^-b,b,c), bounds should be relative to this; for instance, lower=c(0,0,0),upper=c(100,3,1)
##' @return Standard error matrix; approximate covariance matrix of MLE(a,b,c)
##' @examples
##'
##' A_true=10; B_true=1.5; C_true=0.3; sigma=0.1
##'
##' set.seed(31525)
##'
##' X=1+3*rchisq(10000,df=5)
##' Y=A_true*(X^(-B_true)) + C_true + rnorm(length(X),sd=sigma)
##'
##' # 'Observations' - 100 samples
##' obs=sample(length(X),100,rep=F)
##' Xobs=X[obs]; Yobs=Y[obs]
##'
##' # True covariance matrix of MLE of a,b,c on these x values
##' ntest=1000
##' abc_mat_xfix=matrix(0,ntest,3)
##' abc_mat_xvar=matrix(0,ntest,3)
##' E1=A_true*(Xobs^(-B_true)) + C_true
##' for (i in 1:ntest) {
##'   Y1=E1 + rnorm(length(Xobs),sd=sigma)
##'   abc_mat_xfix[i,]=powersolve(Xobs,Y1)$par # Estimate (a,b,c) with same X
##'
##'   X2=1+3*rchisq(length(Xobs),df=5)
##'   Y2=A_true*(X2^(-B_true)) + C_true + rnorm(length(Xobs),sd=sigma)
##'   abc_mat_xvar[i,]=powersolve(X2,Y2)$par # Estimate (a,b,c) with variable X
##' }
##'
##' Ve1=var(abc_mat_xfix) # empirical variance of MLE(a,b,c)|X
##' Vf=powersolve_se(Xobs,Yobs,method='fisher') # estimated SE matrix, asymptotic
##'
##' Ve2=var(abc_mat_xvar) # empirical variance of MLE(a,b,c)
##' Vb=powersolve_se(Xobs,Yobs,method='bootstrap') # estimated SE matrix, bootstrap
##'
##' cat("Empirical variance of MLE(a,b,c)|X\n")
##' print(Ve1)
##' cat("\n")
##' cat("Asymptotic variance of MLE(a,b,c)|X\n")
##' print(Vf)
##' cat("\n\n")
##' cat("Empirical variance of MLE(a,b,c)\n")
##' print(Ve2)
##' cat("\n")
##' cat("Bootstrap-estimated variance of MLE(a,b,c)\n")
##' print(Vb)
##' cat("\n\n")
##'
powersolve_se=function(x,y,method='fisher',init=c(20000,2,0.1),y_var=rep(1,length(y)),n_boot=1000,seed=NULL,...) {
  if (method=="fisher") {
    ## Fisher Information Matrix - straightforward to compute analytically
    FI_mat=function(x,v,a,b,c,s)
      -(1/(s^2 * v))*cbind(
        c(-x^(-2*b),a*(x^(-2*b))*log(x),-x^(-b),0),
        c(a*(x^(-2*b))*log(x),-(a^2)*(x^(-2*b))*(log(x)^2),a*(x^(-b))*log(x),0),
        c(-x^(-b),a*(x^(-b))*log(x),-1,0),
        c(0,0,0,-2*v))

    ## Need to estimate s as well; call powersolve
    abcs=powersolve(x,y,y_var=y_var,init=init,estimate_s=TRUE,...)$par

    fmat=matrix(0,4,4); for (i in 1:length(x)) fmat=fmat + FI_mat(x[i],y_var[i],abcs[1],abcs[2],abcs[3],abcs[4])

    # Check if fmat is invertible
    if ("matrix" %in% class(try(solve(fmat),silent=T))) {
      vmat=solve(fmat)[1:3,1:3]
      return(vmat)
    } else {
      return(matrix(NA,3,3))
    }
  }
  if (method=="bootstrap") {

    # set seed
    if (!is.null(seed)) set.seed(seed)

    # Compute bootstrap resamples
    xboot=matrix(0,n_boot,3)
    for (i in 1:n_boot) {
      sub=sample(length(x),length(x),replace=TRUE)
      xboot[i,]=powersolve(x[sub],y[sub],y_var=y_var[sub],init=init,...)$par
    }

    vmat=var(xboot)
    return(vmat)
  }
}


##' Simulate random dataset similar to ASPRE training data
##'
##'
##' @name sim_random_aspre
##' @description Generate random population of individuals (e.g., newly pregnant women) with given population parameters
##'
##' Assumes independence of parameter variation. This is not a realistic assumption, but is satisfactory for our purposes.
##' @keywords aspre
##' @param n size of population
##' @param pars list of parameters
##' @return
##' @examples
##'
##' # Load ASPRE related data
##' data(params_aspre)
##'
##' X=sim_random_aspre(1000,params_aspre)
##'
##' print(c(median(X$age),params_aspre$age$median))
##'
##' print(rbind(table(X$parity)/1000,params_aspre$parity$freq))
##'
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


##' Add interaction terms corresponding to ASPRE model
##'
##'
##' @name add_aspre_interactions
##' @keywords aspre
##' @description Add various interaction terms to X. Interaction terms correspond to those in ASPRE.
##' @param X data frame
##' @return New data frame containing interaction terms.
##' @examples
##'
##' # Load ASPRE related data
##' data(params_aspre)
##'
##' X=sim_random_aspre(1000,params_aspre)
##' Xnew=add_aspre_interactions(X)
##'
##' print(colnames(X))
##' print(colnames(Xnew))
##'
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


##' Computes ASPRE score
##'
##'
##' @name aspre
##' @description Computes ASPRE model prediction on a matrix `X` of covariates
##'
##' Full ASPRE model from https://www.nejm.org/doi/suppl/10.1056/NEJMoa1704559/suppl_file/nejmoa1704559_appendix.pdf
##'
##' Model is to predict gestational age at PE; that is, a higher score indicates a lower PE risk, so coefficients are negated for model to predict PE risk.
##' @keywords aspre
##' @param X matrix, assumed to be output of sim_random_aspre with parameter params=params_aspre and transformed using add_aspre_interactions
##' @return vector of scores.
##' @examples
##'
##' # Load ASPRE related data
##' data(params_aspre)
##'
##' X=sim_random_aspre(1000,params_aspre)
##' Xnew=add_aspre_interactions(X)
##'
##' aspre_score=aspre(Xnew)
##'
##' plot(density(aspre_score))
##'
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


