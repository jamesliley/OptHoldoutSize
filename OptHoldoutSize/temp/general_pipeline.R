################################################################################
## Setup                                                                      ##
################################################################################

# Set random seed
set.seed(21423)

# Load package
library(OptHoldoutSize)

# Suppose we have population size and cost-per-sample without a risk score as follows
N=100000
k1=0.4

# Suppose that true values of a,b,c are given by
theta_true=c(10000,1.2,0.2)
theta_lower=c(1,0.5,0.1) # lower bounds for estimating theta
theta_upper=c(20000,2,0.5) # upper bounds for estimating theta
theta_init=(theta_lower+theta_upper)/2 # We will start from this value when finding theta

# Kernel width and variance for Gaussian process
kw0=5000
vu0=1e7

# We will presume that these are the values of n for which cost can potentially be evaluated.
n=seq(1000,N,length=300)

# Include legend on plots or not; inclusion can obscure plot elements on small figures
inc_legend=FALSE


## True form of k2, option 1: parametric assumptions satisfied (true)
true_k2_pTRUE=function(n) powerlaw(n,theta_true)

## True form of k2, option 2: parametric assumptions NOT satisfied (false)
true_k2_pFALSE=function(n) powerlaw(n,theta_true) + (1e4)*dnorm(n,mean=4e4,sd=8e3)

## Plot
plot(0,type="n",xlim=range(n),ylim=range(true_k2_pTRUE(n)),
  xlab="n",ylab=expression(paste("k"[2],"(n)")))

lines(n,true_k2_pTRUE(n),type="l",col="black")
lines(n,true_k2_pFALSE(n),type="l",col="red")
legend("topright",c("Power law", "Not power law"),col=c("black","red"),lty=1)


nsamp=200 # Presume we have this many estimates of k_2(n), between 1000 and N
vwmin=0.001; vwmax=0.02 # Sample variances var_w uniformly between these values

nset=round(runif(nsamp,1000,N))
var_w=runif(nsamp,vwmin,vwmax)
d_pTRUE=rnorm(nsamp,mean=true_k2_pTRUE(nset),sd=sqrt(var_w))
d_pFALSE=rnorm(nsamp,mean=true_k2_pFALSE(nset),sd=sqrt(var_w))


nc=1000:N

true_ohs_pTRUE=nc[which.min(k1*nc + true_k2_pTRUE(nc)*(N-nc))]
true_ohs_pFALSE=nc[which.min(k1*nc + true_k2_pFALSE(nc)*(N-nc))]

print(true_ohs_pTRUE)

print(true_ohs_pFALSE)

print("Line 62")


# Estimate a,b, and c from values nset and d
est_abc_pTRUE=powersolve(nset,d_pTRUE,
  lower=theta_lower,upper=theta_upper,init=theta_init)$par
est_abc_pFALSE=powersolve(nset,d_pFALSE,
  lower=theta_lower,upper=theta_upper,init=theta_init)$par

# Estimate optimal holdout sizes using parametric method
param_ohs_pTRUE=optimal_holdout_size(N,k1,theta=est_abc_pTRUE)$size
param_ohs_pFALSE=optimal_holdout_size(N,k1,theta=est_abc_pFALSE)$size

# Estimate optimal holdout sizes using semi-parametric (emulation) method
emul_ohs_pTRUE=optimal_holdout_size_emulation(nset,d_pTRUE,var_w,theta=est_abc_pTRUE,N=N,k1=k1)$size
emul_ohs_pFALSE=optimal_holdout_size_emulation(nset,d_pFALSE,var_w,theta=est_abc_pFALSE,N=N,k1=k1)$size

# In the parametrised case, the parametric model is better
print(true_ohs_pTRUE)
print(param_ohs_pTRUE)
print(emul_ohs_pTRUE)

# In the misparametrised case, the semi-parametric model is better
print(true_ohs_pFALSE)
print(param_ohs_pFALSE)
print(emul_ohs_pFALSE)





################################################################################
## Resample values of d and recalculate OHS                                   ##
################################################################################

n_var=1000 # Resample d this many times to estimate OHS variance

ohs_resample=matrix(NA,n_var,4) # This will be populated with OHS estimates

for (i in 1:n_var) {
  set.seed(36253 + i)

  # Resample values d
  d_pTRUE_r=rnorm(nsamp,mean=true_k2_pTRUE(nset),sd=sqrt(var_w))
  d_pFALSE_r=rnorm(nsamp,mean=true_k2_pFALSE(nset),sd=sqrt(var_w))


  # Estimate a,b, and c from values nset and d
  est_abc_pTRUE_r=powersolve(nset,d_pTRUE_r,y_var=var_w,
    lower=theta_lower,upper=theta_upper,init=theta_init)$par
  est_abc_pFALSE_r=powersolve(nset,d_pFALSE_r,y_var=var_w,
    lower=theta_lower,upper=theta_upper,init=theta_init)$par

  # Estimate optimal holdout sizes using parametric method
  param_ohs_pTRUE_r=optimal_holdout_size(N,k1,theta=est_abc_pTRUE_r)$size
  param_ohs_pFALSE_r=optimal_holdout_size(N,k1,theta=est_abc_pFALSE_r)$size

  # Estimate optimal holdout sizes using semi-parametric (emulation) method
  emul_ohs_pTRUE_r=optimal_holdout_size_emulation(nset,d_pTRUE_r,theta=est_abc_pTRUE_r,var_w,N,k1)$size
  emul_ohs_pFALSE_r=optimal_holdout_size_emulation(nset,d_pFALSE_r,theta=est_abc_pTRUE_r,var_w,N,k1)$size

  ohs_resample[i,]=c(param_ohs_pTRUE_r, param_ohs_pFALSE_r, emul_ohs_pTRUE_r, emul_ohs_pFALSE_r)

  print(i)
}

colnames(ohs_resample)=c("param_pTRUE","param_pFALSE","emul_pTRUE", "emul_pFALSE")
save(ohs_resample,file="data/ohs_resample.RData")



################################################################################
## Next point using parametrisation                                           ##
################################################################################

## Choose an initial five training sizes at which to evaluate k2
set.seed(32424)
nstart=5
nset0=round(runif(nstart,1000,N/2))
var_w0=runif(nstart,vwmin,vwmax)
d0_pTRUE=rnorm(nstart,mean=true_k2_pTRUE(nset0),sd=sqrt(var_w0))
d0_pFALSE=rnorm(nstart,mean=true_k2_pFALSE(nset0),sd=sqrt(var_w0))


# These are our sets of training sizes and k2 estimates, which will be built up.
nset_pTRUE=nset0
d_pTRUE=d0_pTRUE
var_w_pTRUE=var_w0

nset_pFALSE=nset0
d_pFALSE=d0_pFALSE
var_w_pFALSE=var_w0


# Go up to this many points
max_points=200

while(length(nset_pTRUE)<= max_points) {
  set.seed(37261 + length(nset_pTRUE))

  # Estimate parameters
  theta_pTRUE=powersolve(nset_pTRUE,d_pTRUE,y_var=var_w_pTRUE,lower=theta_lower,upper=theta_upper,init=theta_init)$par
  theta_pFALSE=powersolve(nset_pFALSE,d_pFALSE,y_var=var_w_pFALSE,lower=theta_lower,upper=theta_upper,init=theta_init)$par


  # Two panels
  par(mfrow=c(1,2))
  yrange=c(0,100000)

  ## First panel
  plot(0,xlim=range(n),ylim=yrange,type="n",
    xlab="Training/holdout set size",
    ylab="Total cost (= num. cases)")
  points(nset_pTRUE,k1*nset_pTRUE + d_pTRUE*(N-nset_pTRUE),pch=16,col="purple")
  lines(n,k1*n + powerlaw(n,theta_pTRUE)*(N-n),lty=2)
  lines(n,k1*n + true_k2_pTRUE(n)*(N-n),lty=3,lwd=3)
  legend("topright",
    c("Par. est. cost",
      "True",
      "d",
      "Next pt."),
    lty=c(2,3,NA,NA),lwd=c(1,3,NA,NA),pch=c(NA,NA,16,124),pt.cex=c(NA,NA,1,1),
    col=c("black","black","purple","black"),bg="white",border=NA)

  # Add vertical line at next suggested point
  ci_pTRUE = next_n(n,nset_pTRUE,d=d_pTRUE,var_w=var_w_pTRUE,N=N,k1=k1,nmed=15)
  if (!all(is.na(ci_pTRUE))) nextn_pTRUE=n[which.min(ci_pTRUE)] else nextn_pTRUE=round(runif(1,1000,N))
  abline(v=nextn_pTRUE)


  ## Second panel
  plot(0,xlim=range(n),ylim=yrange,type="n",
    xlab="Training/holdout set size",
    ylab="Total cost (= num. cases)")
  points(nset_pFALSE,k1*nset_pFALSE + d_pFALSE*(N-nset_pFALSE),pch=16,col="purple")
  lines(n,k1*n + powerlaw(n,theta_pFALSE)*(N-n),lty=2)
  lines(n,k1*n + true_k2_pFALSE(n)*(N-n),lty=3,lwd=3)
  legend("topright",
    c("Par. est. cost",
      "True",
      "d",
      "Next pt."),
    lty=c(2,3,NA,NA),lwd=c(1,3,NA,NA),pch=c(NA,NA,16,124),pt.cex=c(NA,NA,1,1),
    col=c("black","black","purple","black"),bg="white",border=NA)

  # Add vertical line at next suggested point
  ci_pFALSE = next_n(n,nset_pFALSE,d=d_pFALSE,var_w=var_w_pFALSE,N=N,k1=k1,nmed=15)
  if (!all(is.na(ci_pFALSE))) nextn_pFALSE=n[which.min(ci_pFALSE)] else nextn_pFALSE=round(runif(1,1000,N))
  abline(v=nextn_pFALSE)



  # New estimates of k2
  var_w_new_pTRUE=runif(1,vwmin,vwmax)
  d_new_pTRUE=rnorm(1,mean=true_k2_pTRUE(nextn_pTRUE),sd=sqrt(var_w_new_pTRUE))

  var_w_new_pFALSE=runif(1,vwmin,vwmax)
  d_new_pFALSE=rnorm(1,mean=true_k2_pFALSE(nextn_pFALSE),sd=sqrt(var_w_new_pFALSE))


  # Update data
  nset_pTRUE=c(nset_pTRUE,nextn_pTRUE)
  d_pTRUE=c(d_pTRUE,d_new_pTRUE)
  var_w_pTRUE=c(var_w_pTRUE,var_w_new_pTRUE)

  nset_pFALSE=c(nset_pFALSE,nextn_pFALSE)
  d_pFALSE=c(d_pFALSE,d_new_pFALSE)
  var_w_pFALSE=c(var_w_pFALSE,var_w_new_pFALSE)

  print(length(nset_pFALSE))

  data_nextpoint_par=list(
    nset_pTRUE=nset_pTRUE,nset_pFALSE=nset_pFALSE,
    d_pTRUE=d_pTRUE,d_pFALSE=d_pFALSE,
    var_w_pTRUE=var_w_pTRUE,var_w_pFALSE=var_w_pFALSE)

  save(data_nextpoint_par,file="data/data_nextpoint_par.RData")

  #    Sys.sleep(10)

}

data_nextpoint_par=list(
  nset_pTRUE=nset_pTRUE,nset_pFALSE=nset_pFALSE,
  d_pTRUE=d_pTRUE,d_pFALSE=d_pFALSE,
  var_w_pTRUE=var_w_pTRUE,var_w_pFALSE=var_w_pFALSE)

save(data_nextpoint_par,file="data/data_nextpoint_par.RData")






################################################################################
## Next point using Bayesian emulation                                        ##
################################################################################

## Choose an initial five training sizes at which to evaluate k2
set.seed(32424) # start from same seed as before
nstart=5
nset0=round(runif(nstart,1000,N/2))
var_w0=runif(nstart,vwmin,vwmax)
d0_pTRUE=rnorm(nstart,mean=true_k2_pTRUE(nset0),sd=sqrt(var_w0))
d0_pFALSE=rnorm(nstart,mean=true_k2_pFALSE(nset0),sd=sqrt(var_w0))


# These are our sets of training sizes and k2 estimates, which will be built up.
nset_pTRUE=nset0
d_pTRUE=d0_pTRUE
var_w_pTRUE=var_w0

nset_pFALSE=nset0
d_pFALSE=d0_pFALSE
var_w_pFALSE=var_w0


# Go up to this many points
max_points=200

while(length(nset_pTRUE)<= max_points) {
  set.seed(46352 + length(nset_pTRUE))

  # Estimate parameters for parametric part of semi-parametric method
  theta_pTRUE=powersolve(nset_pTRUE,d_pTRUE,y_var=var_w_pTRUE,lower=theta_lower,upper=theta_upper,init=theta_init)$par
  theta_pFALSE=powersolve(nset_pTRUE,d_pFALSE,y_var=var_w_pTRUE,lower=theta_lower,upper=theta_upper,init=theta_init)$par

  # Mean and variance of emulator for cost function, parametric assumptions satisfied
  p_mu_pTRUE=mu_fn(n,nset=nset_pTRUE,d=d_pTRUE,var_w = var_w_pTRUE,theta=theta_pTRUE,N=N,k1=k1)
  p_var_pTRUE=psi_fn(n,nset=nset_pTRUE,N=N,var_w = var_w_pTRUE)

  # Mean and variance of emulator for cost function, parametric assumptions not satisfied
  p_mu_pFALSE=mu_fn(n,nset=nset_pFALSE,d=d_pFALSE,var_w = var_w_pFALSE,theta=theta_pFALSE,N=N,k1=k1)
  p_var_pFALSE=psi_fn(n,nset=nset_pFALSE,N=N,var_w = var_w_pFALSE)



  # Two panels
  par(mfrow=c(1,2))
  yrange=c(0,100000)


  ## First panel
  plot(0,xlim=range(n),ylim=yrange,type="n",
    xlab="Training/holdout set size",
    ylab="Total cost (= num. cases)")
  lines(n,p_mu_pTRUE,col="blue")
  lines(n,p_mu_pTRUE - 3*sqrt(pmax(0,p_var_pTRUE)),col="red")
  lines(n,p_mu_pTRUE + 3*sqrt(pmax(0,p_var_pTRUE)),col="red")
  points(nset_pTRUE,k1*nset_pTRUE + d_pTRUE*(N-nset_pTRUE),pch=16,col="purple")
  lines(n,k1*n + powerlaw(n,theta_pTRUE)*(N-n),lty=2)
  lines(n,k1*n + true_k2_pTRUE(n)*(N-n),lty=3,lwd=3)
  legend("topright",
    c(expression(mu(n)),
      expression(mu(n) %+-% 3*sqrt(psi(n))),
      "prior(n)",
      "True",
      "d",
      "Next pt."),
    lty=c(1,1,2,3,NA,NA),lwd=c(1,1,1,3,NA,NA),pch=c(NA,NA,NA,NA,16,124),pt.cex=c(NA,NA,NA,NA,1,1),
    col=c("blue","red","black","black","purple","black"),bg="white",border=NA)

  # Add vertical line at next suggested point
  exp_imp_em_pTRUE = exp_imp_fn(n,nset=nset_pTRUE,d=d_pTRUE,var_w = var_w_pTRUE, N=N,k1=k1)
  nextn_pTRUE = n[which.max(exp_imp_em_pTRUE)]
  abline(v=nextn_pTRUE)


  ## Second panel
  plot(0,xlim=range(n),ylim=yrange,type="n",
    xlab="Training/holdout set size",
    ylab="Total cost (= num. cases)")
  lines(n,p_mu_pFALSE,col="blue")
  lines(n,p_mu_pFALSE - 3*sqrt(pmax(0,p_var_pFALSE)),col="red")
  lines(n,p_mu_pFALSE + 3*sqrt(pmax(0,p_var_pFALSE)),col="red")
  points(nset_pFALSE,k1*nset_pFALSE + d_pFALSE*(N-nset_pFALSE),pch=16,col="purple")
  lines(n,k1*n + powerlaw(n,theta_pFALSE)*(N-n),lty=2)
  lines(n,k1*n + true_k2_pFALSE(n)*(N-n),lty=3,lwd=3)
  legend("topright",
    c(expression(mu(n)),
      expression(mu(n) %+-% 3*sqrt(psi(n))),
      "prior(n)",
      "True",
      "d",
      "Next pt."),
    lty=c(1,1,2,3,NA,NA),lwd=c(1,1,1,3,NA,NA),pch=c(NA,NA,NA,NA,16,124),pt.cex=c(NA,NA,NA,NA,1,1),
    col=c("blue","red","black","black","purple","black"),bg="white",border=NA)

  # Add vertical line at next suggested point
  exp_imp_em_pFALSE = exp_imp_fn(n,nset=nset_pFALSE,d=d_pFALSE,var_w = var_w_pFALSE, N=N,k1=k1)
  nextn_pFALSE = n[which.max(exp_imp_em_pFALSE)]
  abline(v=nextn_pFALSE)


  # New estimates of k2
  var_w_new_pTRUE=runif(1,vwmin,vwmax)
  d_new_pTRUE=rnorm(1,mean=true_k2_pTRUE(nextn_pTRUE),sd=sqrt(var_w_new_pTRUE))

  var_w_new_pFALSE=runif(1,vwmin,vwmax)
  d_new_pFALSE=rnorm(1,mean=true_k2_pFALSE(nextn_pFALSE),sd=sqrt(var_w_new_pFALSE))


  # Update data
  nset_pTRUE=c(nset_pTRUE,nextn_pTRUE)
  d_pTRUE=c(d_pTRUE,d_new_pTRUE)
  var_w_pTRUE=c(var_w_pTRUE,var_w_new_pTRUE)

  nset_pFALSE=c(nset_pFALSE,nextn_pFALSE)
  d_pFALSE=c(d_pFALSE,d_new_pFALSE)
  var_w_pFALSE=c(var_w_pFALSE,var_w_new_pFALSE)

  print(length(nset_pFALSE))


}

data_nextpoint_em=list(
  nset_pTRUE=nset_pTRUE,nset_pFALSE=nset_pFALSE,
  d_pTRUE=d_pTRUE,d_pFALSE=d_pFALSE,
  var_w_pTRUE=var_w_pTRUE,var_w_pFALSE=var_w_pFALSE)

save(data_nextpoint_em,file="data/data_nextpoint_em.RData")


################################################################################
## Analyse convergence rate of each algorithm                                 ##
################################################################################

# Function to resample values of d and regenerate OHS given nset and var_w
ntri=function(nset,var_w,k2,nx=100,method="MLE") {
  out=rep(0,nx)
  for (i in 1:nx) {
    d1=rnorm(length(nset),mean=k2(nset),sd=sqrt(var_w))
    theta1=powersolve(nset,d1,y_var=var_w,lower=theta_lower,upper=theta_upper,init=theta_true)$par
    if (method=="MLE") {
      out[i]=optimal_holdout_size(N,k1,theta1)$size
    } else {
      nn=seq(1000,N,length=1000)
      p_mu=mu_fn(nn,nset=nset,d=d1,var_w = var_w, N=N,k1=k1,theta=theta1)
      out[i]=nn[which.min(p_mu)]
    }
  }
  return(out)
}


# Load datasets of 'next point'
load("data/data_nextpoint_em.RData")
load("data/data_nextpoint_par.RData")

# Maximum number of training set sizes we will consider
n_iter=200

# Generate random 'next points'
set.seed(36279)
data_nextpoint_rand=list(
  nset_pTRUE=round(runif(nmax,1000,N)),
  nset_pFALSE=round(runif(nmax,1000,N)),
  var_w_pTRUE=runif(nmax,vwmin,vwmax),
  var_w_pFALSE=runif(nmax,vwmin,vwmax)
)

# Initialise matrices of records
# ohs_array[n,i,j,k,l] is
#  using the first n training set sizes
#  the ith resample
#  using the jth version of k2 (j=1: pTRUE, j=2: pFALSE)
#  using the kth algorithm (k=1: parametric, k=2: semiparametric/emulation)
#  using the lth method of selecting next points (l=1: random, l=2: systematic)
nr=200 # recalculate OHS this many times
ohs_array=array(NA,dim=c(n_iter,nr,2,2,2))

for (i in 5:n_iter) {
  set.seed(363726 + i)
  # Resamplings for parametric algorithm, random next point
  ohs_array[i,,1,1,1]=ntri(
    nset=data_nextpoint_rand$nset_pTRUE[1:i],
    var_w=data_nextpoint_rand$var_w_pTRUE[1:i],
    k2=true_k2_pTRUE,nx=nr,method="MLE")
  ohs_array[i,,2,1,1]=ntri(
    nset=data_nextpoint_rand$nset_pFALSE[1:i],
    var_w=data_nextpoint_rand$var_w_pFALSE[1:i],
    k2=true_k2_pFALSE,nx=nr,method="MLE")

  # Resamplings for semiparametric/emulation algorithm, random next point
  ohs_array[i,,1,2,1]=ntri(
    nset=data_nextpoint_rand$nset_pTRUE[1:i],
    var_w=data_nextpoint_rand$var_w_pTRUE[1:i],
    k2=true_k2_pTRUE,nx=nr,method="EM")
  ohs_array[i,,2,2,1]=ntri(
    nset=data_nextpoint_rand$nset_pFALSE[1:i],
    var_w=data_nextpoint_rand$var_w_pFALSE[1:i],
    k2=true_k2_pFALSE,nx=nr,method="EM")

  # Resamplings for parametric algorithm, nonrandom (systematic) next point
  ohs_array[i,,1,1,2]=ntri(
    nset=data_nextpoint_par$nset_pTRUE[1:i],
    var_w=data_nextpoint_par$var_w_pTRUE[1:i],
    k2=true_k2_pTRUE,nx=nr,method="MLE")
  ohs_array[i,,2,1,2]=ntri(
    nset=data_nextpoint_par$nset_pFALSE[1:i],
    var_w=data_nextpoint_par$var_w_pFALSE[1:i],
    k2=true_k2_pFALSE,nx=nr,method="MLE")

  # Resamplings for semiparametric/emulation algorithm, nonrandom (systematic) next point
  ohs_array[i,,1,2,2]=ntri(
    nset=data_nextpoint_em$nset_pTRUE[1:i],
    var_w=data_nextpoint_em$var_w_pTRUE[1:i],
    k2=true_k2_pTRUE,nx=nr,method="EM")
  ohs_array[i,,2,2,2]=ntri(
    nset=data_nextpoint_em$nset_pFALSE[1:i],
    var_w=data_nextpoint_em$var_w_pFALSE[1:i],
    k2=true_k2_pFALSE,nx=nr,method="EM")

  print(i)

  save(ohs_array,file="data/ohs_array.RData")
}





## Draw figure


# Load data
data(ohs_array)
# ohs_array[n,i,j,k,l] is
#  using the first n training set sizes
#  the ith resample
#  using the jth version of k2 (j=1: pTRUE, j=2: pFALSE)
#  using the kth algorithm (k=1: parametric, k=2: semiparametric/emulation)
#  using the lth method of selecting next points (l=1: random, l=2: systematic)

# Settings
alpha=0.1; # Plot 1-alpha confidence intervals
dd=3 # horizontal line spacing
n_iter=dim(ohs_array)[1] # X axis range
ymax=80000 # Y axis range

# Plot drawing function
plot_ci_convergence=function(title,key,M1,M2,ohs_true) {

  # Set up plot parameters
  par(mar=c(1,4,4,0.1))
  layout(mat=rbind(matrix(1,4,4),matrix(2,2,4)))


  # Initialise
  plot(0,xlim=c(5,n_iter),ylim=c(0,ymax),type="n",
    ylab="OHS and error",xaxt="n",main=title)
  abline(h=ohs_true,col="blue",lty=2)

  # Plot medians
  points(1:n_iter,rowMedians(M1,na.rm=T),pch=16,cex=0.5,col="black")
  points(1:n_iter,rowMedians(M2,na.rm=T),pch=16,cex=0.5,col="red")

  # CIs
  ci1=rbind(apply(M1,1,function(x) pmax(0,quantile(x,alpha/2,na.rm=T))),
    apply(M1,1,function(x) pmin(ymax,quantile(x,1-alpha/2,na.rm=T))))
  ci2=rbind(apply(M2,1,function(x) pmax(0,quantile(x,alpha/2,na.rm=T))),
    apply(M2,1,function(x) pmin(ymax,quantile(x,1-alpha/2,na.rm=T))))

  # Plot CIs
  segments(
    1:n_iter,ci1[1,],
    1:n_iter,ci1[2,],
    col="black"
  )
  segments(
    1:n_iter + 1/dd,ci2[1,],
    1:n_iter + 1/dd,ci2[2,],
    col="red"
  )

  # Add legend
  legend("topright",
    legend=key,
    col=c("black","red"),lty=1)


  # Bottom panel setup
  par(mar=c(4,4,0.1,0.1))
  plot(0,xlim=c(5,n_iter),ylim=c(0,5e5),type="n",
    ylab="CI width",xlab="Training set size")

  # Draw lines
  lines(1:n_iter,ci1[2,]-ci1[1,],col="black")
  lines(1:n_iter,ci2[2,]-ci2[1,],col="red")

}



# Extract matrices from aray
M111=ohs_array[1:n_iter,,1,1,1] # pTRUE, param algorithm, random nextpoint
M112=ohs_array[1:n_iter,,1,1,2] # pTRUE, param algorithm, systematic nextpoint

M211=ohs_array[1:n_iter,,2,1,1] # pFALSE, param algorithm, random nextpoint
M212=ohs_array[1:n_iter,,2,1,2] # pFALSE, param algorithm, systematic nextpoint

M121=ohs_array[1:n_iter,,1,1,1] # pTRUE, emul algorithm, random nextpoint
M122=ohs_array[1:n_iter,,1,1,2] # pTRUE, emul algorithm, systematic nextpoint

M221=ohs_array[1:n_iter,,2,1,1] # pFALSE, emul algorithm, random nextpoint
M222=ohs_array[1:n_iter,,2,1,2] # pFALSE, emul algorithm, systematic nextpoint


# True OHS
nc=1000:N
true_ohs_pTRUE=nc[which.min(k1*nc + true_k2_pTRUE(nc)*(N-nc))]
true_ohs_pFALSE=nc[which.min(k1*nc + true_k2_pFALSE(nc)*(N-nc))]


par(mfrow=c(2,2))
plot_ci_convergence("Params. satis, param. alg.",c("Rand. next n","Syst. next n"),M111,M112,true_ohs_pTRUE)
plot_ci_convergence("Params. not satis, param. alg.",c("Rand. next n","Syst. next n"),M211,M212,true_ohs_pFALSE)

plot_ci_convergence("Params. satis, emul. alg.",c("Rand. next n","Syst. next n"),M121,M122,true_ohs_pTRUE)
plot_ci_convergence("Params. not satis, emul. alg.",c("Rand. next n","Syst. next n"),M221,M222,true_ohs_pFALSE)




