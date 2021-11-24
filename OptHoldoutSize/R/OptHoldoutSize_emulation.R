################################################################################
## R script for functions in OptHoldoutSize package                           ##
################################################################################
##
## Sami Haidar-Wehbe, Sam Emerson, James Liley
## October 2021
##


######################################################################
## Functions for Bayesian emulation of the learning curve           ##
######################################################################


##' Covariance function for Gaussian process
##'
##'
##' @export
##' @name cov_fn
##' @description Radial kernel covariance function for Gaussian process.
##'
##' Used for a Gaussian process `GP(m,k(.,.))`, an instance X of which has covariance k(n,n') between X(n) and X(n').
##'
##' Covariance function is parametrised by `var_u` (general variance) and `k_width` (kernel width)
##'
##' @keywords emulation
##' @param n Argument 1: kernel is a function of `ndash-n`
##' @param ndash Argument 2: kernel is a function of `ndash-n`
##' @param var_u Global variance
##' @param k_width Kernel width
##' @return Covariance value
##' @examples
##'
##' ##' # We will sample from Gaussian processes GP(0,k(.,.)=cov_fn(.,.;var_u,theta)) at these values of n
##' nvals=seq(1,300,length=100)
##'
##' # We will consider two theta values
##' kw1=10; kw2=30
##'
##' # We will consider two var_u values
##' var1=1; var2=10
##'
##' # Covariance matrices
##' cov11=outer(nvals,nvals,function(n,ndash) cov_fn(n,ndash,var_u=var1,k_width=kw1))
##' cov12=outer(nvals,nvals,function(n,ndash) cov_fn(n,ndash,var_u=var1,k_width=kw2))
##' cov21=outer(nvals,nvals,function(n,ndash) cov_fn(n,ndash,var_u=var2,k_width=kw1))
##' cov22=outer(nvals,nvals,function(n,ndash) cov_fn(n,ndash,var_u=var2,k_width=kw2))
##'
##' # Dampen slightly to ensure positive definiteness
##' damp=1e-5
##' cov11=(1-damp)*(1-diag(length(nvals)))*cov11 + diag(length(nvals))*cov11
##' cov12=(1-damp)*(1-diag(length(nvals)))*cov12 + diag(length(nvals))*cov12
##' cov21=(1-damp)*(1-diag(length(nvals)))*cov21 + diag(length(nvals))*cov21
##' cov22=(1-damp)*(1-diag(length(nvals)))*cov22 + diag(length(nvals))*cov22
##'
##' # Sample
##' set.seed(35243)
##' y11=rmnorm(1,mean=0,varcov=cov11)
##' y12=rmnorm(1,mean=0,varcov=cov12)
##' y21=rmnorm(1,mean=0,varcov=cov21)
##' y22=rmnorm(1,mean=0,varcov=cov22)
##'
##' # Plot
##' rr=max(abs(c(y11,y12,y21,y22)))
##' plot(0,xlim=range(nvals),ylim=c(-rr,rr+10),xlab="n",ylab=expression("GP(0,cov_fn(.,.;var_u,theta))"))
##' lines(nvals,y11,lty=1,col="black")
##' lines(nvals,y12,lty=2,col="black")
##' lines(nvals,y21,lty=1,col="red")
##' lines(nvals,y22,lty=2,col="red")
##' legend("topright",c("k_width=10, var_u=1", "k_width=30, var_u=1", "k_width=10, var_u=10","k_width=30, var_u=10"),
##'   lty=c(1,2,1,2),col=c("black","black","red","red"))
##'
cov_fn = function(n,ndash,var_u,k_width){
  return(var_u*exp(-(((n-ndash)/k_width)^2)))
}







##' Power law function
##'
##' @export
##' @name mu_fn
##' @description Power law function for modelling learning curve (taken to mean change in expected loss per sample with training set size)
##'
##' Recommended in [review of learning curve forms](https://arxiv.org/abs/2103.10948)
##'
##' If `theta=c(a,b,c)` then models as `a n^(-b) + c`. Note `b` is negated.
##'
##' Note that `powerlaw(n,c(a,b,c))` has limit `c` as `n` tends to infinity, if `a,b > 0`
##'
##' @param n Set of training set sizes to evaluate
##' @param theta Parameter of values
##' @return Vector of values of same length as `n`
##' @examples
##'
##' ncheck=seq(1000,10000)
##' plot(ncheck, powerlaw(ncheck, c(5e3,1.2,0.3)),type="l",xlab="n",ylab="powerlaw(n)")
##'
powerlaw=function(n,theta)  (theta[1] *n^(-theta[2]) + theta[3])


##' Updating function for mean.
##'
##' @export
##' @name mu_fn
##' @description Posterior mean for emulator given points `n`.
##' @param n Set of training set sizes to evaluate
##' @param nset Training set sizes for which a loss has been evaluated
##' @param d Loss at training set sizes `nset`
##' @param var_w Variance of error in loss estimate at each training set size.
##' @param N Total number of samples on which the model will be fitted/used
##' @param k1 Mean loss per sample with no predictive score in place
##' @param var_u Marginal variance for Gaussian process kernel. Defaults to 1e7
##' @param k_width Kernel width for Gaussian process kernel. Defaults to 5000
##' @param mean_fn Functional form governing expected loss per sample given sample size. Should take two parameters: n (sample size) and theta (parameters). Defaults to function `powerlaw`.
##' @param theta Current estimates of parameter values for mean_fn. Defaults to the MLE power-law solution corresponding to n,d, and var_w.
##' @return
##' @examples
##'
##' # Suppose we have population size and cost-per-sample without a risk score as follows
##' N=100000
##' k1=0.4
##'
##' # Suppose we begin with loss estimates at n-values
##' nset0=c(10000,20000,30000)
##'
##' # with cost-per-individual estimates
##' d0=c(0.35,0.26,0.31)
##'
##' # and associated error on those estimates
##' var_w0=c(0.1^2,0.08^2,0.12^2)
##'
##' # We estimate theta from these three points
##' theta0=powersolve(nset0,d0,y_var=var_w0)$par
##'
##' # We will estimate the posterior at these values of n
##' n=seq(1000,50000,length=1000)
##'
##' # Mean and variance
##' p_mu=mu_fn(n,nset=nset0,d=d0,var_w = var_w0, N=N,k1=k1,theta=theta0,k_width=5000,var_u=3000000)
##' p_var=psi_fn(n,nset=nset0,N=N,var_w = var_w0,k_width=5000,var_u=3000000)
##' p_var=p_var[cbind(1:length(n),1:length(n))] # just the variances
##'
##' # Plot
##' plot(0,xlim=range(n),ylim=range(c(p_mu - 3*sqrt(p_var),p_mu + 3*sqrt(p_var))),type="n",
##'   xlab="Training/holdout set size",
##'   ylab="Total cost (= num. cases)")
##' lines(n,p_mu,col="blue")
##' lines(n,p_mu - 3*sqrt(p_var),col="red")
##' lines(n,p_mu + 3*sqrt(p_var),col="red")
##' points(nset0,k1*nset0 + d0*(N-nset0),pch=16,col="purple")
##' lines(n,k1*n + powerlaw(n,theta0)*(N-n),lty=2)
##' legend("topright",
##'   c(expression(mu(n)),
##'     expression(mu(n) %+-% 3*sqrt(psi(n))),
##'     "prior(n)",
##'     "d"),
##'   lty=c(1,1,2,NA),lwd=c(1,1,1,NA),pch=c(NA,NA,NA,16),pt.cex=c(NA,NA,NA,1),
##'   col=c("blue","red","black","purple"),bg="white")
##'
mu_fn= function(n,nset,d,var_w,N,k1,
  var_u=1e7,
  k_width=5000,
  mean_fn=powerlaw,
  theta=powersolve(nset,d,y_var=var_w)$par){

  loss_fn=function(n,theta) k1*n + mean_fn(n,theta)*(N-n)

  mean_n <- loss_fn(n,theta) # Expected loss at values n according to current parameters theta
  cov_n_nset <- outer(n,nset,cov_fn,var_u=var_u,k_width=k_width) # GP covariance between nstar values and current values n
  cov_nset_nset <- outer(nset,nset,cov_fn,var_u=var_u,k_width=k_width) # Internal covariance matrix for existing values n
  var_w_mat <- diag(var_w*(N-nset)^2) # variance on the scale of loss for whole population

  xd=k1*nset + d*(N-nset) # outputs on scale of whole-dataset loss

  mean_nset <- loss_fn(nset,theta)
  #browser()
  return(as.numeric(mean_n + cov_n_nset%*%solve(cov_nset_nset + var_w_mat)%*%(xd - mean_nset)))
}


##' Updating function for variance.
##'
##' @export
##' @name psi_fn
##' @description Posterior variance for emulator given points `n`.
##' @param n Set of training set sizes to evaluate at
##' @param nset Training set sizes for which a loss has been evaluated
##' @param var_w Variance of error in loss estimate at each training set size.
##' @param N Total number of samples on which the model will be fitted/used. Only used to rescale var_w
##' @param var_u Marginal variance for Gaussian process kernel. Defaults to 1e7
##' @param k_width Kernel width for Gaussian process kernel. Defaults to 5000
##' @return Vector Psi of same length of n where Psi[i]=var(posterior(cost(n[i])))
##' @examples
##'
##' # See examples for `mu_fn`
##'
psi_fn =function(n,nset,var_w,N,
  var_u=1e7,
  k_width=5000){

  cov_n_n <- cov_fn(n,n,var_u=var_u,k_width=k_width)

  cov_n_nset <- outer(n,nset,cov_fn,var_u=var_u,k_width=k_width)
  cov_nset_nset <- outer(nset,nset,cov_fn,var_u=var_u,k_width=k_width)
  var_w_mat <- diag(var_w*(N-nset)^2)
  xsol=solve(cov_nset_nset + var_w_mat)
  return(cov_n_n - apply(cov_n_nset,1,function(x) t(x)%*%xsol%*%x))
}


##' Expected improvement
##'
##' @export
##' @name exp_imp_fn
##' @description Expected improvement
##'
##' Essentially chooses the `next point` to add to `n`, called `n*`, in order to minimise the expectation of `loss(n*)`.
##'
##' @param n Set of training set sizes to evaluate
##' @param nset Training set sizes for which a loss has been evaluated
##' @param d Loss at training set sizes `nset`
##' @param var_w Variance of error in loss estimate at each training set size.
##' @param N Total number of samples on which the model will be fitted/used
##' @param k1 Mean loss per sample with no predictive score in place
##' @param var_u Marginal variance for Gaussian process kernel. Defaults to 1e7
##' @param k_width Kernel width for Gaussian process kernel. Defaults to 5000
##' @param mean_fn Functional form governing expected loss per sample given sample size. Should take two parameters: n (sample size) and theta (parameters). Defaults to function `powerlaw`.
##' @param theta Current estimates of parameter values for mean_fn. Defaults to the MLE power-law solution corresponding to n,d, and var_w.
##' @return
##' @examples
##'
##' # Set seed.
##' set.seed(24015)
##'
##' # Kernel width and Gaussian process variance
##' kw0=5000
##' vu0=1e7
##'
##' # Include legend on plots or not; inclusion can obscure plot elements on small figures
##' inc_legend=FALSE
##'
##' # Suppose we have population size and cost-per-sample without a risk score as follows
##' N=100000
##' k1=0.4
##'
##' # Suppose that true values of a,b,c are given by
##' theta_true=c(10000,1.2,0.2)
##' theta_lower=c(1,0.5,0.1) # lower bounds for estimating theta
##' theta_upper=c(20000,2,0.5) # upper bounds for estimating theta
##'
##'
##'
##' # We start with five random holdout set sizes (nset0),
##' #  with corresponding cost-per-individual estimates d0 derived
##' #  with various errors var_w0
##' nstart=4
##' vwmin=0.001; vwmax=0.005
##' nset0=round(runif(nstart,1000,N/2))
##' var_w0=runif(nstart,vwmin,vwmax)
##' d0=rnorm(nstart,mean=fn(nset0,theta_true),sd=sqrt(var_w0))
##'
##' # We estimate theta from these three points
##' theta0=powersolve(nset0,d0,y_var=var_w0,lower=theta_lower,upper=theta_upper,init=theta_true,control=list(parscale=theta_true))$par
##'
##' # We will estimate the posterior at these values of n
##' n=seq(1000,N,length=1000)
##'
##' # Mean and variance
##' p_mu=mu_fn(n,nset=nset0,d=d0,var_w = var_w0, N=N,k1=k1,theta=theta0,k_width=kw0,var_u=vu0)
##' p_var=psi_fn(n,nset=nset0,N=N,var_w = var_w0,k_width=kw0,var_u=vu0)
##'
##' # Plot
##' yrange=c(-30000,100000)
##' plot(0,xlim=range(n),ylim=yrange,type="n",
##'   xlab="Training/holdout set size",
##'   ylab="Total cost (= num. cases)")
##' lines(n,p_mu,col="blue")
##' lines(n,p_mu - 3*sqrt(p_var),col="red")
##' lines(n,p_mu + 3*sqrt(p_var),col="red")
##' points(nset0,k1*nset0 + d0*(N-nset0),pch=16,col="purple")
##' lines(n,k1*n + powerlaw(n,theta0)*(N-n),lty=2)
##' lines(n,k1*n + powerlaw(n,theta_true)*(N-n),lty=3,lwd=3)
##' if (inc_legend) {
##'   legend("topright",
##'     c(expression(mu(n)),
##'       expression(mu(n) %+-% 3*sqrt(psi(n))),
##'       "prior(n)",
##'       "True",
##'       "d"),
##'     lty=c(1,1,2,3,NA),lwd=c(1,1,1,3,NA),pch=c(NA,NA,NA,NA,16),pt.cex=c(NA,NA,NA,NA,1),
##'     col=c("blue","red","black","purple"),bg="white")
##' }
##'
##' ## Add line corresponding to recommended new point
##' exp_imp_em <- exp_imp_fn(n,nset=nset0,d=d0,var_w = var_w0, N=N,k1=k1,theta=theta0,k_width=kw0,var_u=vu0)
##' abline(v=n[which.max(exp_imp_em)])
##'
exp_imp_fn = function(n,nset,d,var_w,N,k1,
  var_u=1e7,
  k_width=5000,
  mean_fn=powerlaw,
  theta=powersolve(nset,d,y_var=var_w)$par){
  mu_val=mu_fn(n,nset,d,var_w,N,k1,var_u,k_width,mean_fn,theta)
  d_min <- min(d)
  psi_val <- psi_fn(n,nset,var_w,N,var_u,k_width)
  Z <- (-mu_val+d_min)/sqrt(psi_val)
  EI <- ((-mu_val+d_min)*pnorm(Z)) + (sqrt(psi_val)*dnorm(Z))
  return(EI)
}



