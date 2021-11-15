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

##' Estimate optimal holdout size
##'
##'
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
##' @param k2 Function governing expected cost to an individual sample given a predictive score fitted to n samples. Must take two arguments: n (number of samples) and theta (parameters). Defaults to a power-law form ``k2(n,c(a,b,c))=a n^(-b) + c``.
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
  k2 = function(n,par) par[1]*(n^(-par[2])) + par[3],
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

  ## Solve
  par_mat=cbind(N,k1,theta)
  np=dim(par_mat)[2]
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

  return(out)
}





##' Confidence interval for optimal holdout size
##'
##' @name ci_ohs
##' @description Compute confidence interval for optimal holdout size given a set of n_e estimates of parameters.
##'
##' This can be done either asymptotically, using a method analogous to the Fisher information matrix, or empirically (using bootstrap resampling)
##' @keywords estimation
##' @param N Vector of estimates of total number of samples on which the predictive score will be used/fitted. Can be a vector.
##' @param k1 Vector of estimates of cost value in the absence of a predictive score. Can be a vector.
##' @param theta Matrix of estimates of parameters for function k2(n) governing expected cost to an individual sample given a predictive score fitted to n samples. Can be a matrix of dimension n x n_par, where n_par is the number of parameters of k2.
##' @param alpha Construct 1-alpha confidence interval. Defaults to 0.05
##' @param k2 Function governing expected cost to an individual sample given a predictive score fitted to n samples. Must take two arguments: n (number of samples) and theta (parameters). Defaults to a power-law form ``k2(n,c(a,b,c))=a n^(-b) + c``.
##' @param grad_nstar Function giving partial derivatives of optimal holdout set, taking three arguments: N, k1, and theta. Only used for asymptotic confidence intervals. F NULL, estimated empirically
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
  k2 = function(n,par) par[1]*(n^(-par[2])) + par[3],
  grad_nstar=NULL,
  n_boot=10000,
  seed=NULL,
  mode="empirical",
  ...
) {

  ## Error handlers
  if (length(N)<5) stop("At least five samples necessary for estimation")
  if (!(is.numeric(N) & is.numeric(k1) & is.numeric(theta))) stop("Parameters N, k1 and theta must be numeric")
  if (!((length(N)==dim(theta)[1]) & (length(N)==length(k1)))) stop("Parameters N and k1 must have same length, which must match the first dimension of parameter theta")
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
    par_mat=cbind(N,k1,theta)
    mu=colMeans(par_mat)

    ## Determine function for grad_nstar if not specified
    if (is.null(grad_nstar)) {
      grad_nstar=function(N,k1,theta) {
        dx=1e-5 # estimate gradient by using this difference
        par2=outer(rep(1,2+length(theta)),c(N,k1,theta))
        par2=par2 + dx*diag(dim(par2)[1]) # parameters, shifted by dx one-at-a-time
        ohs2=optimal_holdout_size(N=par2[,1],k1=par2[,2],theta=par2[,3:dim(par2)[2]],k2=k2)$size
        ohs1=optimal_holdout_size(N,k1,theta,k2=k2)$size
        return((ohs2-ohs1)/dx)
      }
    }

    # Parameters for asymptotic confidence interval
    nstar=optimal_holdout_size(mu[1],mu[2],mu[3:length(mu)],k2=k2)$size
    beta_est=grad_nstar(mu[1],mu[2],mu[3:length(mu)])
    sigma_hat=var(par_mat)
    z_a=-qnorm(alpha/2)
    n_e=length(N)
    # Variation from estimated value in asymptotic confidence interval
    di=z_a*sqrt((t(beta_est) %*% sigma_hat %*% beta_est)/n_e)

    cx=nstar+c(-di,di)
    names(cx)=c("lower","upper")
    return(cx)
  }

  if (mode=="empirical") {
    par_mat=cbind(N,k1,theta)
    n_e=dim(par_mat)[1]

    # Random seed
    if (!is.null(seed)) set.seed(seed)

    # Populate with mean parameters from bootstrap resamples
    ci_mat=matrix(0,n_boot,2+dim(theta)[2])
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




