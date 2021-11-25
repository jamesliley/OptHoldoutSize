################################################################################
## R script for functions in OptHoldoutSize package                           ##
################################################################################
##
## Sami Haidar-Wehbe, Sam Emerson, James Liley
## October 2021
##

## TODO: Possibly include function for confidence interval for OHS and cost given
#  point estimates of k1, N, theta and an estimate of sigma, assuming normality.




################################################################################
## Optimal holdout size, confidence interval, and gradient                    ##
################################################################################

##' Estimate optimal holdout size under parametric assumptions
##'
##'
##' @export
##' @name optimal_holdout_size
##' @description Compute optimal holdout size for updating a predictive score given appropriate parameters of cost function
##'
##' Evaluates empirical minimisation of cost function ``l(n;k1,N,theta) = k1 n + k2(n;theta) (N-n)``.
##'
##' The function will return `Inf` if no minimum exists. It does not check if the minimum is unique, but this can be guaranteed using the assumptions for theorem 1 in the manuscript.
##'
##' This calls the function `optimize` from package `stats`.
##' @keywords estimation
##' @param N Total number of samples on which the predictive score will be used/fitted. Can be a vector.
##' @param k1 Cost value in the absence of a predictive score. Can be a vector.
##' @param theta Parameters for function k2(n) governing expected cost to an individual sample given a predictive score fitted to n samples. Can be a matrix of dimension n x n_par, where n_par is the number of parameters of k2.
##' @param k2 Function governing expected cost to an individual sample given a predictive score fitted to n samples. Must take two arguments: n (number of samples) and theta (parameters). Defaults to a power-law form ``powerlaw(n,c(a,b,c))=a n^(-b) + c``.
##' @param round_result Set to TRUE to solve over integral sizes
##' @param ... Passed to function `optimize`
##' @return List/data frame of dimension (number of evaluations) x (4 + n_par) containing input data and results. Columns `size` and `cost` are optimal holdout size and cost at this size respectively. Parameters N, k1, theta.1, theta.2,...,theta.{n_par} are input data.
##' @examples
##'
##' # Evaluate optimal holdout set size for a range of values of k1 and two values of N, some of which lead to infinite values
##' N1=10000; N2=12000
##' k1=seq(0.1,0.5,length=20)
##' A=3; B=1.5; C=0.15; theta=c(A,B,C)
##'
##' res1=optimal_holdout_size(N1,k1,theta)
##' res2=optimal_holdout_size(N2,k1,theta)
##'
##' par(mfrow=c(1,2))
##' plot(0,type="n",ylim=c(0,500),xlim=range(res1$k1),xlab=expression("k"[1]),ylab="Optimal holdout set size")
##'   lines(res1$k1,res1$size,col="black")
##'   lines(res2$k1,res2$size,col="red")
##'   legend("topright",as.character(c(N1,N2)),title="N:",col=c("black","red"),lty=1)
##' plot(0,type="n",ylim=c(1500,1600),xlim=range(res1$k1),xlab=expression("k"[1]),ylab="Minimum cost")
##'   lines(res1$k1,res1$cost,col="black")
##'   lines(res2$k1,res2$cost,col="red")
##'   legend("topleft",as.character(c(N1,N2)),title="N:",col=c("black","red"),lty=1)
optimal_holdout_size=function(
  N,
  k1,
  theta,
  k2 = powerlaw,
  round_result=FALSE,
  ...
) {

  ## Error handlers
  if (!is.numeric(c(N,k1))) stop("Parameters N and k1 must be numeric")
  if (!(is.numeric(theta)|is.matrix(theta))) stop("Parameter theta must be a vector or matrix")
  if (!is.function(k2)) stop("Parameter k2 must be a function taking two arguments: n and theta")
  if (length(as.list(args(k2)))!=3) stop("Parameter k2 must be a function taking two arguments: n and theta")
  if ((length(N)>1 | length(k1)>1)|!is.null(dim(theta))) {
    n=max(length(N),length(k1),dim(theta)[1])
    if (is.null(dim(theta))) theta=t(matrix(theta,length(theta),n))
    if (!all(c(length(k1),length(N),dim(theta)[1]) %in% c(1,n))) stop("If vectors, N and k1 must have the same length, or have length one. Parameter theta must either be a vector or a matrix of dimension length(n) x n_par, where n_par is the number of parameters of k2.")
    if (length(N)==1) N=rep(N,n)
    if (length(k1)==1) k1=rep(k1,n)
  } else theta=matrix(theta,1,length(theta))

  par_mat=cbind(N,k1,theta)
  np=dim(par_mat)[2]

    ## Solve
  out=apply(par_mat,1,function(x) {
    unlist(optimize(function(n) x[2]*n + k2(n,x[3:np])*(x[1]-n),c(1,N),...))
  })

  ## Post-processing
  out=cbind(t(out),par_mat)
  out=as.data.frame(out)
  colnames(out)=c("size","cost","N","k1",paste0("theta.",1:dim(theta)[2]))

  ## Replace errors with infinity
  no_min=which((out$size> (N-1))|(out$size < 2))
  out$size[no_min]=Inf
  out$cost[no_min]=Inf

  # Round result if required. Minimum will be at either floor(size) or ceiling(size)
  if (round_result) {
    size1=floor(out$size)
    size2=ceiling(out$size)
    k2_1=apply(cbind(size1,theta),1,function(x) k2(x[1],x[2:(1+dim(theta)[2])]))
    k2_2=apply(cbind(size2,theta),1,function(x) k2(x[1],x[2:(1+dim(theta)[2])]))
    cost1=out$k1*size1 + k2_1*(out$N-size1)
    cost2=out$k1*size2 + k2_1*(out$N-size2)

    w=which(cost2<cost1)
    out_size=size1; out_cost=cost1
    out_size[w]=size2[w]; out_cost[w]=cost2[w]

    out$size=out_size
    out$cost=out_cost
  }

  class(out)=c("optholdoutsize",class(out))
  return(out)
}

##' Plot estimated cost function
##'
##' @export
##' @name plot.optholdoutsize
##' @keywords estimation
##' @description Plot estimated cost function, when parametric method is used for estimation.
##'
##' Draws cost function as a line and indicates minimum. Assumes a power-law form of k2 unless parameter k2 is set otherwise.
##'
##' @param X Object of type `optholdoutsize`
##' @param k2 Function governing expected cost to an individual sample given a predictive score fitted to n samples. Must take two arguments: n (number of samples) and theta (parameters). Defaults to a power-law form ``powerlaw(n,c(a,b,c))=a n^(-b) + c``.
##' @param ... Other arguments passed to `plot()` and `lines()`
##' @examples
##'
##' # Simple example
##'
##' N=100000;
##' k1=0.3
##' A=8000; B=1.5; C=0.15; theta=c(A,B,C)
##'
##' res1=optimal_holdout_size(N,k1,theta)
##'
##' plot(res1)
##'
plot.optholdoutsize=function(X,k2=powerlaw,...) {

  N=X$N; k1=X$k1; size=X$size; cost=X$cost;
  if (length(N)==1) {
    theta=unlist(X[5:length(X)])
    xx=seq(1,N,length.out=1000)
    yy= k1*xx + k2(xx,theta)*(N-xx)
    costN=yy[length(xx)]

    plot(xx,yy,type="l",ylim=c(2*cost-costN,costN + cost),xlab="Holdout size",ylab="Total cost",...)
    abline(v=size,lty=2); abline(h=cost,lty=2)
    points(size,cost,pch=16,col="black")
    legend("bottomright",c("Est. cost","OHS"),lty=c(1,NA),pch=c(NA,16),col=c("black","black"),bty="n")
  } else {
    Nmax=max(N)
    theta=as.data.frame(X[5:length(X)])

    xx=seq(1,Nmax,length.out=1000)
    yy=matrix(0,length(N),length(xx))
    for (i in 1:length(N)) yy[i,]=k1[i]*xx + k2(xx,as.numeric(theta[i,]))*(N[i]-xx)

    costN=yy[,length(xx)]

    plot(0,type="n",xlim=c(1,Nmax),ylim=c(2*min(cost)-mean(costN),max(costN) + min(cost)),
      xlab="Holdout size",ylab="Total cost",...)
    for (i in 1:length(N)) lines(xx,yy[i,],...)
    points(size,cost,pch=16)
    legend("bottomright",c("Est. cost","OHS"),lty=c(1,NA),pch=c(NA,16),col=c("black","black"),bty="n")
  }
}





##' Estimate optimal holdout size under semi-parametric assumptions
##'
##'
##' @export
##' @name optimal_holdout_size_emulation
##' @keywords estimation,emulation
##' @description Compute optimal holdout size for updating a predictive score given a set of training set sizes and estimates of mean cost per sample at those training set sizes.
##'
##' This is essentially a wrapper for function `mu_fn()`.
##'
##'
##' @param nset Training set sizes for which a loss has been evaluated
##' @param d Loss at training set sizes `nset`
##' @param var_w Variance of error in loss estimate at each training set size.
##' @param N Total number of samples on which the model will be fitted/used
##' @param k1 Mean loss per sample with no predictive score in place
##' @param var_u Marginal variance for Gaussian process kernel. Defaults to 1e7
##' @param k_width Kernel width for Gaussian process kernel. Defaults to 5000
##' @param mean_fn Functional form governing expected loss per sample given sample size. Should take two parameters: n (sample size) and theta (parameters). Defaults to function `powerlaw`.
##' @param theta Current estimates of parameter values for mean_fn. Defaults to the MLE power-law solution corresponding to n,d, and var_w.
##' @param npoll Check npoll equally spaced values between 1 and N for minimum. If NULL, check all values (this can be slow). Defaults to 1000
##' @param ... Passed to function `optimise()`
##' @return
##' @examples
##'
##'
optimal_holdout_size_emulation= function(nset,d,var_w,N,k1,
  var_u=1e7,
  k_width=5000,
  mean_fn=powerlaw,
  theta=powersolve(nset,d,y_var=var_w)$par,
  npoll=1000,
  ...){

  if (!is.null(npoll)) n=seq(1,N,length=npoll) else n=1:N
  xv=mu_fn(n,nset,d,var_w,N,k1,var_u,k_width,mean_fn,theta)
  w=which.min(xv)

  out=list(xv[w],n[w],nset,d,var_w,N,k1,var_u,k_width,theta)
  ## The 'optimise' below doesn't tend to work bimodally
  #minf=function(n) mu_fn(n,nset,d,var_w,N,k1,var_u,k_width,mean_fn,theta)
  #out=optimise(minf,c(1,N))

  names(out)=c("cost","size","nset","d","var_w","N","k1","var_u","k_width","theta")
  class(out)=c("optholdoutsize_emul",class(out))
  return(out)
}

##' Plot estimated cost function using emulation (semiparametric)
##'
##' @export
##' @name plot.optholdoutsize_emul
##' @keywords estimation
##' @description Plot estimated cost function, when semiparametric (emulation) method is used for estimation.
##'
##' Draws posterior mean of cost function as a line and indicates minimum. Also draws mean +/- 3 SE.
##'
##' Assumes a power-law form of k2 unless parameter k2 is set otherwise.
##'
##' @param X Object of type `optholdoutsize_emul`
##' @param k2 Function governing expected cost to an individual sample given a predictive score fitted to n samples. Must take two arguments: n (number of samples) and theta (parameters). Defaults to a power-law form ``powerlaw(n,c(a,b,c))=a n^(-b) + c``.
##' @param ... Other arguments passed to `plot()`
##' @examples
##'
##' # Simple example
##'
##' # Parameters
##' N=100000;
##' k1=0.3
##' A=8000; B=1.5; C=0.15; theta=c(A,B,C)
##'
##' # True mean function
##' k2_true=function(n) powerlaw(n,theta)
##'
##' # Values of n for which cost has been estimated
##' np=50 # this many points
##' nset=round(runif(np,1,N))
##' var_w=runif(np,0.001,0.002)
##' d=rnorm(np,mean=k2_true(nset),sd=sqrt(var_w))
##'
##' # Compute OHS
##' res1=optimal_holdout_size_emulation(nset,d,var_w,N,k1)
##'
##' # Plot
##' plot(res1)
plot.optholdoutsize_emul=function(X,k2=powerlaw,...) {

  # Plot at these values of n
  nn=seq(1,N,length.out=1000)

  # Mean function
  xv=mu_fn(nn,X$nset,X$d,X$var_w,X$N,X$k1,X$var_u,X$k_width,mean_fn=k2,X$theta)

  # Variance function
  psiv=pmax(0,psi_fn(nn, X$nset, X$var_w, X$N, X$var_u, X$k_width))

  # 3 SD bounds
  xlo=xv - 3*sqrt(psiv)
  xhi=xv + 3*sqrt(psiv)

  # m(n)
  xm=X$k1*nn + k2(nn,X$theta)*(X$N-nn)

  xq=quantile(c(xlo,xhi),c(0.1,0.9)); xmid=median(xv)
  yr=c(xmid - 2*(xmid-xq[1]),xmid + 2*(xq[2]-xmid))

  plot(0,type="n",xlab="Holdout size",ylab="Total cost",xlim=c(1,N),ylim=yr,...)
  lines(nn,xm,lty=2,col="black")
  lines(nn,xv,col="blue")
  lines(nn,xlo,col="red")
  lines(nn,xhi,col="red")
  points(X$nset,X$k1*X$nset + X$d*(X$N - X$nset),col="purple",pch=16)
  abline(h=X$cost,lty=2); abline(v=X$size,lty=2)
  legend("bottomright",
    c(expression(mu(n)),
      expression(mu(n) %+-% 3*sqrt(psi(n))),
      "m(n)",
      "d"),
    lty=c(1,1,2,NA),lwd=c(1,1,1,NA),pch=c(NA,NA,NA,16),pt.cex=c(NA,NA,NA,1),
    col=c("blue","red","black","purple"),bg="white",bty="n")
}




##' Confidence interval for optimal holdout size, when estimated using parametric method
##'
##' @export
##' @name ci_ohs
##' @description Compute confidence interval for optimal holdout size given either a standard error covariance matrix or a set of n_e estimates of parameters.
##'
##' This can be done either asymptotically, using a method analogous to the Fisher information matrix, or empirically (using bootstrap resampling)
##'
##' If sigma (covariance matrix) is specified and method='bootstrap', a confidence interval is generated assuming a Gaussian distribution of (N,k1,theta). To estimate a confidence interval assuming a non-Gaussian distribution, simulate values under the requisite distribution and use then as parameters N,k1, theta, with sigma set to NULL.
##'
##' @keywords estimation
##' @param N Vector of estimates of total number of samples on which the predictive score will be used/fitted. Can be a vector.
##' @param k1 Vector of estimates of cost value in the absence of a predictive score. Can be a vector.
##' @param theta Matrix of estimates of parameters for function k2(n) governing expected cost to an individual sample given a predictive score fitted to n samples. Can be a matrix of dimension n x n_par, where n_par is the number of parameters of k2.
##' @param alpha Construct 1-alpha confidence interval. Defaults to 0.05
##' @param k2 Function governing expected cost to an individual sample given a predictive score fitted to n samples. Must take two arguments: n (number of samples) and theta (parameters). Defaults to a power-law form ``k2(n,c(a,b,c))=a n^(-b) + c``.
##' @param grad_nstar Function giving partial derivatives of optimal holdout set, taking three arguments: N, k1, and theta. Only used for asymptotic confidence intervals. F NULL, estimated empirically
##' @param sigma Standard error covariance matrix for (N,k1,theta), in that order. If NULL, will derive as sample covariance matrix of parameters. Must be of the correct size and positive definite.
##' @param n_boot Number of bootstrap resamples for empirical estimate.
##' @param seed Random seed for bootstrap resamples. Defaults to NULL.
##' @param mode One of 'asymptotic' or 'empirical'. Defaults to 'empirical'
##' @param ... Passed to function `optimize`
##' @return A vector of length two containing lower and upper limits of confidence interval.
##' @examples
##' ## We will assume that our observations of N, k1, and theta=(a,b,c) are distributed with mean mu_par and variance sigma_par
##' mu_par=c(N=10000,k1=0.35,A=3,B=1.5,C=0.1)
##' sigma_par=cbind(
##'   c(100^2,       1,      0,       0,       0),
##'   c(    1,  0.07^2,      0,       0,       0),
##'   c(    0,       0,  0.5^2,    0.05,  -0.001),
##'   c(    0,       0,   0.05,   0.4^2,  -0.002),
##'   c(    0,       0, -0.001,  -0.002,  0.02^2)
##' )
##'
##' # Firstly, we make 500 observations
##' par_obs=rmnorm(500,mean=mu_par,varcov=sigma_par)
##'
##' # Optimal holdout size and asymptotic and empirical confidence intervals
##' ohs=optimal_holdout_size(N=mean(par_obs[,1]),k1=mean(par_obs[,2]),theta=colMeans(par_obs[,3:5]))$size
##' ci_a=ci_ohs(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],alpha=0.05,seed=12345,mode="asymptotic")
##' ci_e=ci_ohs(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],alpha=0.05,seed=12345,mode="empirical")
##'
##'
##' # Assess cover at various n_e
##' n_e_values=c(20,30,50,100,150,200,300,500,750,1000,1500)
##' ntrial=5000
##' alpha_trial=0.1 # use 90% confidence intervals
##' nstar_true=optimal_holdout_size(N=mu_par[1],k1=mu_par[2],theta=mu_par[3:5])$size
##'
##' ## The matrices indicating cover take are included in this package but take around 30 minutes to generate. They are generated using the code below (including random seeds).
##' data(ci_cover_a_yn)
##' data(ci_cover_e_yn)
##'
##' if (!exists("ci_cover_a_yn")) {
##'   ci_cover_a_yn=matrix(NA,length(n_e_values),ntrial) # Entry [i,j] is 1 if ith asymptotic CI for jth value of n_e covers true nstar
##'   ci_cover_e_yn=matrix(NA,length(n_e_values),ntrial) # Entry [i,j] is 1 if ith empirical CI for jth value of n_e covers true nstar
##'
##'   for (i in 1:length(n_e_values)) {
##'     n_e=n_e_values[i]
##'     for (j in 1:ntrial) {
##'       # Set seed
##'       set.seed(j*ntrial + i + 12345)
##'
##'       # Make n_e observations
##'       par_obs=rmnorm(n_e,mean=mu_par,varcov=sigma_par)
##'       ci_a=ci_ohs(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],alpha=alpha_trial,mode="asymptotic")
##'       ci_e=ci_ohs(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],alpha=alpha_trial,mode="empirical",n_boot=500)
##'
##'       if (nstar_true>ci_a[1] & nstar_true<ci_a[2]) ci_cover_a_yn[i,j]=1 else ci_cover_a_yn[i,j]=0
##'       if (nstar_true>ci_e[1] & nstar_true<ci_e[2]) ci_cover_e_yn[i,j]=1 else ci_cover_e_yn[i,j]=0
##'     }
##'     print(paste0("Completed for n_e = ",n_e))
##'   }
##'
##' }
##'
##' # Cover at each n_e value and standard error
##' cover_a=rowMeans(ci_cover_a_yn)
##' cover_e=rowMeans(ci_cover_e_yn)
##' zse_a=2*sqrt(cover_a*(1-cover_a)/ntrial)
##' zse_e=2*sqrt(cover_e*(1-cover_e)/ntrial)
##'
##'
##' # Draw plot. Convergence to 1-alpha cover is evident. Cover is not far from alpha even at small n_e.
##'
##' plot(0,type="n",xlim=range(n_e_values),ylim=c(0.7,1),xlab=expression("n"[e]),ylab="Cover")
##'
##' # Asymptotic cover and 2*SE pointwise envelope
##' polygon(c(n_e_values,rev(n_e_values)),c(cover_a+zse_a,rev(cover_a-zse_a)),
##'   col=rgb(1,1,1,alpha=0.3),border=NA)
##' lines(n_e_values,cover_a,col="black")
##'
##' # Empirical cover and 2*SE pointwiseenvelope
##' polygon(c(n_e_values,rev(n_e_values)),c(cover_e+zse_e,rev(cover_e-zse_e)),
##'   col=rgb(0,0,1,alpha=0.3),border=NA)
##' lines(n_e_values,cover_e,col="blue")
##'
##' abline(h=1-alpha_trial,col="red")
##' legend("bottomright",c("Asym.","Emp.",expression(paste("1-",alpha))),lty=1,col=c("black","blue","red"))
##'
ci_ohs=function(
  N,
  k1,
  theta,
  alpha=0.05,
  k2 = powerlaw,
  grad_nstar=NULL,
  sigma=NULL,
  n_boot=10000,
  seed=NULL,
  mode="empirical",
  ...
) {

  ## Error handlers
  if (length(N)<5 & is.null(sigma)) stop("At least five samples necessary for estimation if parameter sigma not specified")
  if (!(is.numeric(N) & is.numeric(k1) & is.numeric(theta))) stop("Parameters N, k1 and theta must be numeric")
  if (!is.function(k2)) stop("Parameter k2 must be a function taking two arguments: n and theta")
  if (length(as.list(args(k2)))!=3) stop("Parameter k2 must be a function taking two arguments: n and theta")
  if (!is.null(grad_nstar)) {
    if (!is.function(grad_nstar)) stop("Parameter grad_nstar must be a function taking three arguments: N, k1, and theta")
    if (length(as.list(args(grad_nstar)))!=4) stop("Parameter grad_nstar must be a function taking three arguments: N, k1, and theta")
  }
  if (!(mode %in% c("asymptotic","empirical"))) stop("Parameter mode must be either 'asymptotic' or 'empirical'")
  if (mode=="empirical") if (n_boot*alpha < 10) stop("Parameter nboot too low for this level of alpha")
  if (!is.numeric(alpha)) stop("Parameter alpha must be numeric")
  if (alpha <= 0 | alpha>=1) stop("Parameter alpha must be >0 and <1")

  if (mode=="asymptotic") {
    if (is.null(sigma)) {
      par_mat=cbind(N,k1,theta)
      mu=colMeans(par_mat)
      sigma_hat=var(par_mat)
    } else {
      mu=c(N,k1,theta)
      sigma_hat=sigma
    }

    ## Determine function for grad_nstar if not specified
    if (is.null(grad_nstar)) {
      grad_nstar=function(N,k1,theta) {
        dx=1e-5 # estimate gradient by using this difference
        par2=outer(rep(1,2+length(theta)),c(N,k1,theta))
        par2=par2 + dx*diag(dim(par2)[1]) # parameters, shifted by dx one-at-a-time
        ohs2=optimal_holdout_size(N=par2[,1],k1=par2[,2],theta=par2[,3:dim(par2)[2]],k2=k2)$size
        ohs1=optimal_holdout_size(N,k1,theta,k2=k2)$size
        return(t((ohs2-ohs1)/dx))
      }
    }

    # Parameters for asymptotic confidence interval
    nstar=optimal_holdout_size(mu[1],mu[2],mu[3:length(mu)],k2=k2)$size
    beta_est=grad_nstar(mu[1],mu[2],mu[3:length(mu)])
    z_a=-qnorm(alpha/2)
    n_e=length(N)
    # Variation from estimated value in asymptotic confidence interval
    di=z_a*sqrt(((beta_est) %*% sigma_hat %*% t(beta_est))/n_e)

    cx=nstar+c(-di,di)
    names(cx)=c("lower","upper")
    return(cx)
  }

  if (mode=="empirical") {

    if (is.null(sigma)) {
      par_mat=cbind(N,k1,theta)
      mu=colMeans(par_mat)
      sigma_hat=var(par_mat)
      n_e=dim(par_mat)[1]
    } else {
      mu=c(N,k1,theta)
      sigma_hat=sigma
      n_e=n_boot

      # Sample, allowing that some parameters may have zero variance
      w=which(sigma_hat[cbind(1:length(mu),1:length(mu))]<1e-20)
      wc=setdiff(1:length(mu),w)
      if (length(w)>0) {
        psub=outer(rep(1,n_boot),mu[w])
        psubc=rmnorm(n_boot,mean=mu[wc],varcov=sigma_hat[wc,wc])
        par_mat=outer(rep(1,n_boot),rep(0,length(mu)))
        par_mat[,w]=psub; par_mat[,wc]=psubc
      } else par_mat=rmnorm(n_boot,mean=mu,varcov=sigma_hat)
    }

    # Random seed
    if (!is.null(seed)) set.seed(seed)

    # Populate with mean parameters from bootstrap resamples
    ci_mat=matrix(0,n_boot,dim(par_mat)[2])
    for (i in 1:n_boot) {
      sboot=sample(n_e,n_e,replace=TRUE)
      ci_mat[i,]=colMeans(par_mat[sboot,])
    }

    # Compute OHSs for bootstrap resamples
    ohs_boot=optimal_holdout_size(ci_mat[,1],ci_mat[,2],ci_mat[,3:dim(ci_mat)[2]],k2=k2)

    # Estimate confidence interval as quantile
    cx=quantile(ohs_boot$size,c(alpha/2,1-(alpha/2)))
    names(cx)=c("lower","upper")
    return(cx)
  }
}



##' Gradient of optimal holdout size (power law)
##'
##'
##' @export
##' @name grad_nstar_powerlaw
##' @description Compute gradient of optimal holdout size assuming a power-law form of k2
##'
##' Assumes cost function is ``l(n;k1,N,theta) = k1 n + k2(n;theta) (N-n)`` with ``k2(n;theta)=k2(n;a,b,c)= a n^(-b) + c``
##' @keywords estimation
##' @param N Total number of samples on which the predictive score will be used/fitted. Can be a vector.
##' @param k1 Cost value in the absence of a predictive score. Can be a vector.
##' @param theta Parameters for function k2(n) governing expected cost to an individual sample given a predictive score fitted to n samples. Can be a matrix of dimension n x n_par, where n_par is the number of parameters of k2.
##' @return List/data frame of dimension (number of evaluations) x 5 containing partial derivatives of nstar (optimal holdout size) with respect to N, k1, a, b, c respectively.
##' @examples
##'
##' # Evaluate optimal holdout set size for a range of values of k1, and compute derivative
##' N=10000;
##' k1=seq(0.1,0.5,length=20)
##' A=3; B=1.5; C=0.15; theta=c(A,B,C)
##'
##' nstar=optimal_holdout_size(N,k1,theta)
##' grad_nstar=grad_nstar_powerlaw(N,k1,theta)
##'
##' plot(0,type="n",ylim=c(-2000,500),xlim=range(res1$k1),xlab=expression("k"[1]),ylab="Optimal holdout set size")
##' lines(res$k1,res$size,col="black")
##' lines(res$k1,grad_nstar[,2],col="red")
##' legend("bottomright",c(expression("n"["*"]),expression(paste(partialdiff[k1],"n"["*"]))),
##'     col=c("black","red"),lty=1)
##'
grad_nstar_powerlaw=function(
  N,
  k1,
  theta
) {

  ## Error handlers
  if (!is.numeric(c(N,k1))) stop("Parameters N and k1 must be numeric")
  if (!(is.numeric(theta)|is.matrix(theta))) stop("Parameter theta must be a vector or matrix")
  if ((length(N)>1 | length(k1)>1)|!is.null(dim(theta))) {
    n=max(length(N),length(k1),dim(theta)[1])
    if (is.null(dim(theta))) theta=t(matrix(theta,length(theta),n))
    if (!all(c(length(k1),length(N),dim(theta)[1]) %in% c(1,n))) stop("If vectors, N and k1 must have the same length, or have length one. Parameter theta must either be a vector or a matrix of dimension length(n) x n_par, where n_par is the number of parameters of k2.")
    if (length(N)==1) N=rep(N,n)
    if (length(k1)==1) k1=rep(k1,n)
  } else theta=matrix(theta,1,length(theta))
  if (dim(theta)[2]!=3) stop("Parameter theta must have length 3 or dimension n x 3")


  ns=optimal_holdout_size(N,k1,theta)$size

  A=theta[,1]; B=theta[,2]; C=theta[,3]


  # partial derivatives of ns*
  dnstarda = (B*N*ns - (B-1)*(ns^2))/(B*(B+1)*A*N - B*(B-1)*A*ns)
  dnstardb = -(ns*N*( -1 + B*log(ns))-(ns^2)*(- 1 + (B-1)*log(ns)))/(B*(B+1)*N - B*(B-1)*ns)
  dnstardc = 1/(B*(B+1)*A*N*ns^(-B-2) - B*(B-1)*A*ns^(-B-1))
  dnstardk1 = -1/(B*(B+1)*B*N*ns^(-B-2) - B*(B-1)*B*ns^(-B-1))
  dnstardN = B*ns/(B*(B+1)*N - B*(B-1)*ns)

  return(cbind(
    dnstardN,
    dnstardk1,
    dnstarda,
    dnstardb,
    dnstardc
  ))
}



##' Finds best value of n to sample next
##'
##' @export
##' @name next_n
##' @description Recommends a value of `n` at which to next evaluate individual cost in order to most accurately estimate optimal holdout size. Currently only for use with a power-law parametrisation of k2.
##'
##' Approximately finds a set of n points which, given estimates of cost, minimise width of 95% confidence interval around OHS. Uses a greedy algorithm, so various parameters can be learned along the way.
##'
##' Given existing training set size/cost estimates `nset` and `d`, with `var_w[i]=variance(d[i])`, finds, for each candidate point `n[i]`, the median width of the 95% confidence interval for OHS if
##'
##' ``nset <- c(nset,n[i])``
##' ``var_w <- c(var_w,mean(var_w))``
##' ``d <- c(d,rnorm(powerlaw(n[i],theta),variance=mean(var_w)))``
##'
##' @param n Set of training set sizes to evaluate
##' @param nset Training set sizes for which a loss has been evaluated
##' @param d Loss at training set sizes `nset`
##' @param var_w Variance of error in loss estimate at each training set size.
##' @param N Total number of samples on which the model will be fitted/used
##' @param k1 Mean loss per sample with no predictive score in place
##' @param nmed number of times to re-evaluate d and confidence interval width.
##' @param ... Passed to `powersolve` and `powersolve_se`
##' @return Vector `out` of same length as `n`, where `out[i]` is the expected width of the 95% confidence interval for OHS should `n` be added to `nset`.
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
##' nstart=10
##' vwmin=0.001; vwmax=0.005
##' nset0=round(runif(nstart,1000,N/2))
##' var_w0=runif(nstart,vwmin,vwmax)
##' d0=rnorm(nstart,mean=powerlaw(nset0,theta_true),sd=sqrt(var_w0))
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
##' ## Add line corresponding to recommended new point. This is slow.
##' exp_imp <- next_n(n,nset=nset,d=d,var_w = var_w, N=N,k1=k1,nmed=10,
##'                      lower=theta_lower,upper=theta_upper,
##'                      init=theta_true)
##' abline(v=n[which.min(exp_imp)])
##'
##'
next_n=function(n,nset,d,var_w,N,k1,nmed=100,...) {
  out=rep(0,length(n))
  theta=powersolve(nset,d,y_var=var_w,...)$par
  for (i in 1:length(n)) {
    out_i=rep(0,nmed)
    for (j in 1:nmed) {
      # Data appended with new candidate point
      nsetx=c(nset,n[i]);
      var_wx=c(var_w,mean(var_w));
      dx=c(d,rnorm(1,mean=powerlaw(n[i],theta),sd=sqrt(mean(var_w))))

      # Parameter estimates with new candidate point
      thetax=powersolve(nsetx,dx,y_var=var_wx,...)$par
      covx=powersolve_se(nsetx,dx,y_var=var_wx,method="fisher",...)
      cov_all=matrix(0,5,5); cov_all[3:5,3:5]=covx

      if (!is.na(cov_all[1,1])) {
        cx=ci_ohs(N,k1,thetax,sigma=cov_all,mode="asymptotic",grad_nstar=grad_nstar_powerlaw)
        out_i[j]=max(cx)-min(cx)
      } else out[j]=Inf
    }
    out[i]=median(out_i)
  }
  return(out)
}



##' Measure of error for emulation-based OHS emulation
##'
##' @export
##' @name error_ohs_emulation
##' @description Measure of error for semiparametric (emulation) based estimation of optimal holdout set sizes.
##'
##' Returns a set of values of n for which a `1-alpha` credible interval for cost at includes a lower value than the cost at the estimated optimal holdout size.
##'
##' This is not a confidence interval, credible interval or credible set for the OHS, and is prone to misinterpretation.
##'
##' @keywords estimation,emulation
##' @param nset Training set sizes for which a loss has been evaluated
##' @param d Loss at training set sizes `nset`
##' @param var_w Variance of error in loss estimate at each training set size.
##' @param N Total number of samples on which the model will be fitted/used
##' @param k1 Mean loss per sample with no predictive score in place
##' @param alpha Use 1-alpha credible interval. Defaults to 0.1.
##' @param var_u Marginal variance for Gaussian process kernel. Defaults to 1e7
##' @param k_width Kernel width for Gaussian process kernel. Defaults to 5000
##' @param mean_fn Functional form governing expected loss per sample given sample size. Should take two parameters: n (sample size) and theta (parameters). Defaults to function `powerlaw`.
##' @param theta Current estimates of parameter values for mean_fn. Defaults to the MLE power-law solution corresponding to n,d, and var_w.
##' @param npoll Check npoll equally spaced values between 1 and N for minimum. If NULL, check all values (this can be slow). Defaults to 1000
##' @return Vector of values `n` for which 1-alpha credible interval for cost `l(n)` at n contains mean posterior loss at estimated optimal holdout size.
##' @examples
##'
##'  # Set seed
##' set.seed(57365)
##'
##' # Parameters
##' N=100000;
##' k1=0.3
##' A=8000; B=1.5; C=0.15; theta=c(A,B,C)
##'
##' # True mean function
##' k2_true=function(n) powerlaw(n,theta)
##'
##' # True OHS
##' nx=1:N
##' ohs_true=nx[which.min(k1*nx + k2_true(nx)*(N-nx))]
##'
##' # Values of n for which cost has been estimated
##' np=50 # this many points
##' nset=round(runif(np,1,N))
##' var_w=runif(np,0.001,0.0015)
##' d=rnorm(np,mean=k2_true(nset),sd=sqrt(var_w))
##'
##' # Compute OHS
##' res1=optimal_holdout_size_emulation(nset,d,var_w,N,k1)
##'
##' # Error estimates
##' ex=error_ohs_emulation(nset,d,var_w,N,k1)
##'
##' # Plot
##' plot(res1)
##'
##' # Add error
##' abline(v=ohs_true)
##' abline(v=ex,col=rgb(1,0,0,alpha=0.2))
##'
##' # Show justification for error
##' n=seq(1,N,length=1000)
##' mu=mu_fn(n,nset,d,var_w,N,k1); psi=pmax(0,psi_fn(n, nset, var_w, N)); Z=-qnorm(0.1/2)
##' lines(n,mu - Z*sqrt(psi),lty=2,lwd=2)
##' legend("topright",
##'     c("Err. region",expression(paste(mu(n)- "z"[alpha/2]*sqrt(psi(n))))),
##'     pch=c(16,NA),lty=c(NA,2),lwd=c(NA,2),col=c("pink","black"),bty="n")
error_ohs_emulation=function(nset,d,var_w,N,k1,
  alpha=0.1,
  var_u=1e7,
  k_width=5000,
  mean_fn=powerlaw,
  theta=powersolve(nset,d,y_var=var_w)$par,
  npoll=1000,
  ...){

  # Candidate values n
  if (!is.null(npoll)) n=seq(1,N,length=npoll) else n=1:N

  # mu and psi
  xmu=mu_fn(n,nset,d,var_w,N,k1,var_u,k_width,mean_fn,theta)
  xpsi=pmax(0,psi_fn(n, nset, var_w, N, var_u, k_width))

  # Compute minimum
  w=which.min(xmu)
  ohs=n[w]
  est_min=xmu[w]

  z=-qnorm(alpha/2) # Need to be this many standard deviations below the mean

  # Values of n for which 1-alpha credible interval for l(n) includes est_min
  n_cont=n[which(xmu - z*sqrt(xpsi) <= est_min)]

  return(n_cont)
}
