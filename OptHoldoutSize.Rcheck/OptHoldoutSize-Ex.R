pkgname <- "OptHoldoutSize"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('OptHoldoutSize')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("add_aspre_interactions")
### * add_aspre_interactions

flush(stderr()); flush(stdout())

### Name: add_aspre_interactions
### Title: Add interaction terms corresponding to ASPRE model
### Aliases: add_aspre_interactions
### Keywords: aspre

### ** Examples


# Load ASPRE related data
data(params_aspre)

X=sim_random_aspre(1000,params_aspre)
Xnew=add_aspre_interactions(X)

print(colnames(X))
print(colnames(Xnew))




cleanEx()
nameEx("aspre")
### * aspre

flush(stderr()); flush(stdout())

### Name: aspre
### Title: Computes ASPRE score
### Aliases: aspre
### Keywords: aspre

### ** Examples


# Load ASPRE related data
data(params_aspre)

X=sim_random_aspre(1000,params_aspre)
Xnew=add_aspre_interactions(X)

aspre_score=aspre(Xnew)

plot(density(aspre_score))




cleanEx()
nameEx("aspre_k2")
### * aspre_k2

flush(stderr()); flush(stdout())

### Name: aspre_k2
### Title: Cost estimating function in ASPRE simulation
### Aliases: aspre_k2
### Keywords: aspre

### ** Examples


# Simulate
set.seed(32142)

N=1000; p=15
X=matrix(rnorm(N*p),N,p); PRE=rbinom(N,1,prob=logit(X%*% rnorm(p)))
aspre_k2(1000,X,PRE)



cleanEx()
nameEx("ci_ohs")
### * ci_ohs

flush(stderr()); flush(stdout())

### Name: ci_ohs
### Title: Confidence interval for optimal holdout size, when estimated
###   using parametric method
### Aliases: ci_ohs
### Keywords: estimation

### ** Examples

## We will assume that our observations of N, k1, and theta=(a,b,c) are distributed with mean mu_par and variance sigma_par
mu_par=c(N=10000,k1=0.35,A=3,B=1.5,C=0.1)
sigma_par=cbind(
  c(100^2,       1,      0,       0,       0),
  c(    1,  0.07^2,      0,       0,       0),
  c(    0,       0,  0.5^2,    0.05,  -0.001),
  c(    0,       0,   0.05,   0.4^2,  -0.002),
  c(    0,       0, -0.001,  -0.002,  0.02^2)
)

# Firstly, we make 500 observations
par_obs=rmnorm(500,mean=mu_par,varcov=sigma_par)

# Optimal holdout size and asymptotic and empirical confidence intervals
ohs=optimal_holdout_size(N=mean(par_obs[,1]),k1=mean(par_obs[,2]),theta=colMeans(par_obs[,3:5]))$size
ci_a=ci_ohs(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],alpha=0.05,seed=12345,mode="asymptotic")
ci_e=ci_ohs(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],alpha=0.05,seed=12345,mode="empirical")


# Assess cover at various n_e
n_e_values=c(20,30,50,100,150,200,300,500,750,1000,1500)
ntrial=5000
alpha_trial=0.1 # use 90% confidence intervals
nstar_true=optimal_holdout_size(N=mu_par[1],k1=mu_par[2],theta=mu_par[3:5])$size

## The matrices indicating cover take are included in this package but take around 30 minutes to generate. They are generated using the code below (including random seeds).
data(ci_cover_a_yn)
data(ci_cover_e_yn)

if (!exists("ci_cover_a_yn")) {
  ci_cover_a_yn=matrix(NA,length(n_e_values),ntrial) # Entry [i,j] is 1 if ith asymptotic CI for jth value of n_e covers true nstar
  ci_cover_e_yn=matrix(NA,length(n_e_values),ntrial) # Entry [i,j] is 1 if ith empirical CI for jth value of n_e covers true nstar

  for (i in 1:length(n_e_values)) {
    n_e=n_e_values[i]
    for (j in 1:ntrial) {
      # Set seed
      set.seed(j*ntrial + i + 12345)

      # Make n_e observations
      par_obs=rmnorm(n_e,mean=mu_par,varcov=sigma_par)
      ci_a=ci_ohs(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],alpha=alpha_trial,mode="asymptotic")
      ci_e=ci_ohs(N=par_obs[,1],k1=par_obs[,2],theta=par_obs[,3:5],alpha=alpha_trial,mode="empirical",n_boot=500)

      if (nstar_true>ci_a[1] & nstar_true<ci_a[2]) ci_cover_a_yn[i,j]=1 else ci_cover_a_yn[i,j]=0
      if (nstar_true>ci_e[1] & nstar_true<ci_e[2]) ci_cover_e_yn[i,j]=1 else ci_cover_e_yn[i,j]=0
    }
    print(paste0("Completed for n_e = ",n_e))
  }

}

# Cover at each n_e value and standard error
cover_a=rowMeans(ci_cover_a_yn)
cover_e=rowMeans(ci_cover_e_yn)
zse_a=2*sqrt(cover_a*(1-cover_a)/ntrial)
zse_e=2*sqrt(cover_e*(1-cover_e)/ntrial)


# Draw plot. Convergence to 1-alpha cover is evident. Cover is not far from alpha even at small n_e.

plot(0,type="n",xlim=range(n_e_values),ylim=c(0.7,1),xlab=expression("n"[e]),ylab="Cover")

# Asymptotic cover and 2*SE pointwise envelope
polygon(c(n_e_values,rev(n_e_values)),c(cover_a+zse_a,rev(cover_a-zse_a)),
  col=rgb(1,1,1,alpha=0.3),border=NA)
lines(n_e_values,cover_a,col="black")

# Empirical cover and 2*SE pointwiseenvelope
polygon(c(n_e_values,rev(n_e_values)),c(cover_e+zse_e,rev(cover_e-zse_e)),
  col=rgb(0,0,1,alpha=0.3),border=NA)
lines(n_e_values,cover_e,col="blue")

abline(h=1-alpha_trial,col="red")
legend("bottomright",c("Asym.","Emp.",expression(paste("1-",alpha))),lty=1,col=c("black","blue","red"))




cleanEx()
nameEx("cov_fn")
### * cov_fn

flush(stderr()); flush(stdout())

### Name: cov_fn
### Title: Covariance function for Gaussian process
### Aliases: cov_fn
### Keywords: emulation

### ** Examples


##' # We will sample from Gaussian processes GP(0,k(.,.)=cov_fn(.,.;var_u,theta)) at these values of n
nvals=seq(1,300,length=100)

# We will consider two theta values
kw1=10; kw2=30

# We will consider two var_u values
var1=1; var2=10

# Covariance matrices
cov11=outer(nvals,nvals,function(n,ndash) cov_fn(n,ndash,var_u=var1,k_width=kw1))
cov12=outer(nvals,nvals,function(n,ndash) cov_fn(n,ndash,var_u=var1,k_width=kw2))
cov21=outer(nvals,nvals,function(n,ndash) cov_fn(n,ndash,var_u=var2,k_width=kw1))
cov22=outer(nvals,nvals,function(n,ndash) cov_fn(n,ndash,var_u=var2,k_width=kw2))

# Dampen slightly to ensure positive definiteness
damp=1e-5
cov11=(1-damp)*(1-diag(length(nvals)))*cov11 + diag(length(nvals))*cov11
cov12=(1-damp)*(1-diag(length(nvals)))*cov12 + diag(length(nvals))*cov12
cov21=(1-damp)*(1-diag(length(nvals)))*cov21 + diag(length(nvals))*cov21
cov22=(1-damp)*(1-diag(length(nvals)))*cov22 + diag(length(nvals))*cov22

# Sample
set.seed(35243)
y11=rmnorm(1,mean=0,varcov=cov11)
y12=rmnorm(1,mean=0,varcov=cov12)
y21=rmnorm(1,mean=0,varcov=cov21)
y22=rmnorm(1,mean=0,varcov=cov22)

# Plot
rr=max(abs(c(y11,y12,y21,y22)))
plot(0,xlim=range(nvals),ylim=c(-rr,rr+10),xlab="n",ylab=expression("GP(0,cov_fn(.,.;var_u,theta))"))
lines(nvals,y11,lty=1,col="black")
lines(nvals,y12,lty=2,col="black")
lines(nvals,y21,lty=1,col="red")
lines(nvals,y22,lty=2,col="red")
legend("topright",c("k_width=10, var_u=1", "k_width=30, var_u=1", "k_width=10, var_u=10","k_width=30, var_u=10"),
  lty=c(1,2,1,2),col=c("black","black","red","red"))




cleanEx()
nameEx("error_ohs_emulation")
### * error_ohs_emulation

flush(stderr()); flush(stdout())

### Name: error_ohs_emulation
### Title: Measure of error for emulation-based OHS emulation
### Aliases: error_ohs_emulation
### Keywords: estimation,emulation

### ** Examples


 # Set seed
set.seed(57365)

# Parameters
N=100000;
k1=0.3
A=8000; B=1.5; C=0.15; theta=c(A,B,C)

# True mean function
k2_true=function(n) powerlaw(n,theta)

# True OHS
nx=1:N
ohs_true=nx[which.min(k1*nx + k2_true(nx)*(N-nx))]

# Values of n for which cost has been estimated
np=50 # this many points
nset=round(runif(np,1,N))
var_k2=runif(np,0.001,0.0015)
k2=rnorm(np,mean=k2_true(nset),sd=sqrt(var_k2))

# Compute OHS
res1=optimal_holdout_size_emulation(nset,k2,var_k2,N,k1)

# Error estimates
ex=error_ohs_emulation(nset,k2,var_k2,N,k1)

# Plot
plot(res1)

# Add error
abline(v=ohs_true)
abline(v=ex,col=rgb(1,0,0,alpha=0.2))

# Show justification for error
n=seq(1,N,length=1000)
mu=mu_fn(n,nset,k2,var_k2,N,k1); psi=pmax(0,psi_fn(n, nset, var_k2, N)); Z=-qnorm(0.1/2)
lines(n,mu - Z*sqrt(psi),lty=2,lwd=2)
legend("topright",
    c("Err. region",expression(paste(mu(n)- "z"[alpha/2]*sqrt(psi(n))))),
    pch=c(16,NA),lty=c(NA,2),lwd=c(NA,2),col=c("pink","black"),bty="n")



cleanEx()
nameEx("exp_imp_fn")
### * exp_imp_fn

flush(stderr()); flush(stdout())

### Name: exp_imp_fn
### Title: Expected improvement
### Aliases: exp_imp_fn

### ** Examples


# Set seed.
set.seed(24015)

# Kernel width and Gaussian process variance
kw0=5000
vu0=1e7

# Include legend on plots or not; inclusion can obscure plot elements on small figures
inc_legend=FALSE

# Suppose we have population size and cost-per-sample without a risk score as follows
N=100000
k1=0.4

# Suppose that true values of a,b,c are given by
theta_true=c(10000,1.2,0.2)
theta_lower=c(1,0.5,0.1) # lower bounds for estimating theta
theta_upper=c(20000,2,0.5) # upper bounds for estimating theta



# We start with five random holdout set sizes (nset0),
#  with corresponding cost-per-individual estimates k2_0 derived
#  with various errors var_k2_0
nstart=4
vwmin=0.001; vwmax=0.005
nset0=round(runif(nstart,1000,N/2))
var_k2_0=runif(nstart,vwmin,vwmax)
k2_0=rnorm(nstart,mean=powerlaw(nset0,theta_true),sd=sqrt(var_k2_0))

# We estimate theta from these three points
theta0=powersolve(nset0,k2_0,y_var=var_k2_0,lower=theta_lower,upper=theta_upper,init=theta_true)$par

# We will estimate the posterior at these values of n
n=seq(1000,N,length=1000)

# Mean and variance
p_mu=mu_fn(n,nset=nset0,k2=k2_0,var_k2 = var_k2_0, N=N,k1=k1,theta=theta0,k_width=kw0,var_u=vu0)
p_var=psi_fn(n,nset=nset0,N=N,var_k2 = var_k2_0,k_width=kw0,var_u=vu0)

# Plot
yrange=c(-30000,100000)
plot(0,xlim=range(n),ylim=yrange,type="n",
  xlab="Training/holdout set size",
  ylab="Total cost (= num. cases)")
lines(n,p_mu,col="blue")
lines(n,p_mu - 3*sqrt(p_var),col="red")
lines(n,p_mu + 3*sqrt(p_var),col="red")
points(nset0,k1*nset0 + k2_0*(N-nset0),pch=16,col="purple")
lines(n,k1*n + powerlaw(n,theta0)*(N-n),lty=2)
lines(n,k1*n + powerlaw(n,theta_true)*(N-n),lty=3,lwd=3)
if (inc_legend) {
  legend("topright",
    c(expression(mu(n)),
      expression(mu(n) %+-% 3*sqrt(psi(n))),
      "prior(n)",
      "True",
      "d"),
    lty=c(1,1,2,3,NA),lwd=c(1,1,1,3,NA),pch=c(NA,NA,NA,NA,16),pt.cex=c(NA,NA,NA,NA,1),
    col=c("blue","red","black","purple"),bg="white")
}

## Add line corresponding to recommended new point
exp_imp_em <- exp_imp_fn(n,nset=nset0,k2=k2_0,var_k2 = var_k2_0, N=N,k1=k1,theta=theta0,k_width=kw0,var_u=vu0)
abline(v=n[which.max(exp_imp_em)])




cleanEx()
nameEx("gen_base_coefs")
### * gen_base_coefs

flush(stderr()); flush(stdout())

### Name: gen_base_coefs
### Title: Coefficients for imperfect risk score
### Aliases: gen_base_coefs
### Keywords: simulation

### ** Examples


# See examples for model_predict



cleanEx()
nameEx("gen_preds")
### * gen_preds

flush(stderr()); flush(stdout())

### Name: gen_preds
### Title: Generate matrix of random observations
### Aliases: gen_preds
### Keywords: simulation

### ** Examples


# See examples for model_predict



cleanEx()
nameEx("gen_resp")
### * gen_resp

flush(stderr()); flush(stdout())

### Name: gen_resp
### Title: Generate response
### Aliases: gen_resp
### Keywords: simulation

### ** Examples


# See examples for model_predict



cleanEx()
nameEx("grad_nstar_powerlaw")
### * grad_nstar_powerlaw

flush(stderr()); flush(stdout())

### Name: grad_nstar_powerlaw
### Title: Gradient of optimal holdout size (power law)
### Aliases: grad_nstar_powerlaw
### Keywords: estimation

### ** Examples


# Evaluate optimal holdout set size for a range of values of k1, and compute derivative
N=10000;
k1=seq(0.1,0.5,length=20)
A=3; B=1.5; C=0.15; theta=c(A,B,C)

nstar=optimal_holdout_size(N,k1,theta)
grad_nstar=grad_nstar_powerlaw(N,k1,theta)

plot(0,type="n",ylim=c(-2000,500),xlim=range(k1),xlab=expression("k"[1]),ylab="Optimal holdout set size")
lines(nstar$k1,nstar$size,col="black")
lines(nstar$k1,grad_nstar[,2],col="red")
legend("bottomright",c(expression("n"["*"]),expression(paste(partialdiff[k1],"n"["*"]))),
    col=c("black","red"),lty=1)




cleanEx()
nameEx("logistic")
### * logistic

flush(stderr()); flush(stdout())

### Name: logistic
### Title: Logistic
### Aliases: logistic

### ** Examples


# Plot
x=seq(0,1,length=100)
plot(x,logistic(x),type="l")

# Logit and logistic are inverses
x=seq(-5,5,length=1000)
plot(x,logistic(logit(x)),type="l")



cleanEx()
nameEx("logit")
### * logit

flush(stderr()); flush(stdout())

### Name: logit
### Title: Logit
### Aliases: logit

### ** Examples


# Plot
x=seq(-5,5,length=1000)
plot(x,logit(x),type="l")



cleanEx()
nameEx("model_predict")
### * model_predict

flush(stderr()); flush(stdout())

### Name: model_predict
### Title: Make predictions
### Aliases: model_predict
### Keywords: simulation

### ** Examples


## Set seed for reproducibility
seed=1234
set.seed(seed)

# Initialisation of patient data
n_iter <- 500           # Number of point estimates to be calculated
nobs <- 5000            # Number of observations, i.e patients
npreds <- 7             # Number of predictors

# Model family
family="log_reg"

# Baseline behaviour is an oracle Bayes-optimal predictor on only one variable
max_base_powers <- 1
base_vars=1

# Check the following holdout size fractions
frac_ho = 0.1


# Set ground truth coefficients, and the accuracy at baseline
coefs_general <- rnorm(npreds,sd=1/sqrt(npreds))
coefs_base <- gen_base_coefs(coefs_general, max_base_powers = max_base_powers)

# Generate dataset
X <- gen_preds(nobs, npreds)

# Generate labels
newdata <- gen_resp(X, coefs = coefs_general)
Y <- newdata$classes

# Combined dataset
pat_data <- cbind(X, Y)
pat_data$Y = factor(pat_data$Y)

# For each holdout size, split data into intervention and holdout set
mask <- split_data(pat_data, frac_ho)
data_interv <- pat_data[!mask,]
data_hold <- pat_data[mask,]

# Train model
trained_model <- model_train(data_hold, model_family = family)
thresh <- 0.5

# Make predictions
class_pred <- model_predict(data_interv, trained_model,
                            return_type = "class",
                            threshold = 0.5, model_family = family)


# Simulate baseline predictions
base_pred <- oracle_pred(data_hold,coefs_base[base_vars, ], num_vars = base_vars)


# Contingency table for model-based predictor (on intervention set)
print(table(class_pred,data_interv$Y))

# Contingency table for model-based predictor (on holdout set)
print(table(base_pred,data_hold$Y))




cleanEx()
nameEx("model_train")
### * model_train

flush(stderr()); flush(stdout())

### Name: model_train
### Title: Train model (wrapper)
### Aliases: model_train
### Keywords: simulation

### ** Examples


# See examples for model_predict



cleanEx()
nameEx("mu_fn")
### * mu_fn

flush(stderr()); flush(stdout())

### Name: mu_fn
### Title: Updating function for mean.
### Aliases: mu_fn

### ** Examples


# Suppose we have population size and cost-per-sample without a risk score as follows
N=100000
k1=0.4

# Kernel width and variance for GP
k_width=5000
var_u=8000000

# Suppose we begin with k2() estimates at n-values
nset=c(10000,20000,30000)

# with cost-per-individual estimates
k2=c(0.35,0.26,0.28)

# and associated error on those estimates
var_k2=c(0.02^2,0.01^2,0.03^2)

# We estimate theta from these three points
theta=powersolve(nset,k2,y_var=var_k2)$par

# We will estimate the posterior at these values of n
n=seq(1000,50000,length=1000)

# Mean and variance
p_mu=mu_fn(n,nset=nset,k2=k2,var_k2 = var_k2, N=N,k1=k1,theta=theta,
           k_width=k_width,var_u=var_u)
p_var=psi_fn(n,nset=nset,N=N,var_k2 = var_k2,k_width=k_width,var_u=var_u)

# Plot
plot(0,xlim=range(n),ylim=c(20000,60000),type="n",
     xlab="Training/holdout set size",
     ylab="Total cost (= num. cases)")
lines(n,p_mu,col="blue")
lines(n,p_mu - 3*sqrt(p_var),col="red")
lines(n,p_mu + 3*sqrt(p_var),col="red")
points(nset,k1*nset + k2*(N-nset),pch=16,col="purple")
lines(n,k1*n + powerlaw(n,theta)*(N-n),lty=2)
segments(nset,k1*nset + (k2 - 3*sqrt(var_k2))*(N-nset),
         nset,k1*nset + (k2 + 3*sqrt(var_k2))*(N-nset))
legend("topright",
       c(expression(mu(n)),
         expression(mu(n) %+-% 3*sqrt(psi(n))),
         "prior(n)",
         "d",
         "3SD(d|n)"),
       lty=c(1,1,2,NA,NA),lwd=c(1,1,1,NA,NA),pch=c(NA,NA,NA,16,124),
       pt.cex=c(NA,NA,NA,1,1),
       col=c("blue","red","black","purple","black"),bg="white")



cleanEx()
nameEx("next_n")
### * next_n

flush(stderr()); flush(stdout())

### Name: next_n
### Title: Finds best value of n to sample next
### Aliases: next_n

### ** Examples


# Set seed.
set.seed(24015)

# Kernel width and Gaussian process variance
kw0=5000
vu0=1e7

# Include legend on plots or not; inclusion can obscure plot elements on small figures
inc_legend=FALSE

# Suppose we have population size and cost-per-sample without a risk score as follows
N=100000
k1=0.4

# Suppose that true values of a,b,c are given by
theta_true=c(10000,1.2,0.2)
theta_lower=c(1,0.5,0.1) # lower bounds for estimating theta
theta_upper=c(20000,2,0.5) # upper bounds for estimating theta



# We start with five random holdout set sizes (nset0),
#  with corresponding cost-per-individual estimates k2_0 derived
#  with various errors var_k2_0
nstart=10
vwmin=0.001; vwmax=0.005
nset0=round(runif(nstart,1000,N/2))
var_k2_0=runif(nstart,vwmin,vwmax)
k2_0=rnorm(nstart,mean=powerlaw(nset0,theta_true),sd=sqrt(var_k2_0))

# We estimate theta from these three points
theta0=powersolve(nset0,k2_0,y_var=var_k2_0,lower=theta_lower,upper=theta_upper,init=theta_true)$par

# We will estimate the posterior at these values of n
n=seq(1000,N,length=1000)

# Mean and variance
p_mu=mu_fn(n,nset=nset0,k2=k2_0,var_k2 = var_k2_0, N=N,k1=k1,theta=theta0,k_width=kw0,var_u=vu0)
p_var=psi_fn(n,nset=nset0,N=N,var_k2 = var_k2_0,k_width=kw0,var_u=vu0)

# Plot
yrange=c(-30000,100000)
plot(0,xlim=range(n),ylim=yrange,type="n",
  xlab="Training/holdout set size",
  ylab="Total cost (= num. cases)")
lines(n,p_mu,col="blue")
lines(n,p_mu - 3*sqrt(p_var),col="red")
lines(n,p_mu + 3*sqrt(p_var),col="red")
points(nset0,k1*nset0 + k2_0*(N-nset0),pch=16,col="purple")
lines(n,k1*n + powerlaw(n,theta0)*(N-n),lty=2)
lines(n,k1*n + powerlaw(n,theta_true)*(N-n),lty=3,lwd=3)
if (inc_legend) {
  legend("topright",
    c(expression(mu(n)),
      expression(mu(n) %+-% 3*sqrt(psi(n))),
      "prior(n)",
      "True",
      "d"),
    lty=c(1,1,2,3,NA),lwd=c(1,1,1,3,NA),pch=c(NA,NA,NA,NA,16),pt.cex=c(NA,NA,NA,NA,1),
    col=c("blue","red","black","purple"),bg="white")
}

## Add line corresponding to recommended new point. This is slow.
nn=seq(1000,N,length=20)
exp_imp <- next_n(nn,nset=nset0,k2=k2_0,var_k2 = var_k2_0, N=N,k1=k1,nmed=10,
                     lower=theta_lower,upper=theta_upper)
abline(v=nn[which.min(exp_imp)])





cleanEx()
nameEx("optimal_holdout_size")
### * optimal_holdout_size

flush(stderr()); flush(stdout())

### Name: optimal_holdout_size
### Title: Estimate optimal holdout size under parametric assumptions
### Aliases: optimal_holdout_size
### Keywords: estimation

### ** Examples


# Evaluate optimal holdout set size for a range of values of k1 and two values of N, some of which lead to infinite values
N1=10000; N2=12000
k1=seq(0.1,0.5,length=20)
A=3; B=1.5; C=0.15; theta=c(A,B,C)

res1=optimal_holdout_size(N1,k1,theta)
res2=optimal_holdout_size(N2,k1,theta)

par(mfrow=c(1,2))
plot(0,type="n",ylim=c(0,500),xlim=range(res1$k1),xlab=expression("k"[1]),ylab="Optimal holdout set size")
  lines(res1$k1,res1$size,col="black")
  lines(res2$k1,res2$size,col="red")
  legend("topright",as.character(c(N1,N2)),title="N:",col=c("black","red"),lty=1)
plot(0,type="n",ylim=c(1500,1600),xlim=range(res1$k1),xlab=expression("k"[1]),ylab="Minimum cost")
  lines(res1$k1,res1$cost,col="black")
  lines(res2$k1,res2$cost,col="red")
  legend("topleft",as.character(c(N1,N2)),title="N:",col=c("black","red"),lty=1)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("optimal_holdout_size_emulation")
### * optimal_holdout_size_emulation

flush(stderr()); flush(stdout())

### Name: optimal_holdout_size_emulation
### Title: Estimate optimal holdout size under semi-parametric assumptions
### Aliases: optimal_holdout_size_emulation
### Keywords: estimation,emulation

### ** Examples


# See examples for mu_fn()



cleanEx()
nameEx("oracle_pred")
### * oracle_pred

flush(stderr()); flush(stdout())

### Name: oracle_pred
### Title: Generate responses
### Aliases: oracle_pred
### Keywords: simulation

### ** Examples


# See examples for model_predict



cleanEx()
nameEx("plot.optholdoutsize")
### * plot.optholdoutsize

flush(stderr()); flush(stdout())

### Name: plot.optholdoutsize
### Title: Plot estimated cost function
### Aliases: plot.optholdoutsize
### Keywords: estimation

### ** Examples


# Simple example

N=100000;
k1=0.3
A=8000; B=1.5; C=0.15; theta=c(A,B,C)

res1=optimal_holdout_size(N,k1,theta)

plot(res1)




cleanEx()
nameEx("plot.optholdoutsize_emul")
### * plot.optholdoutsize_emul

flush(stderr()); flush(stdout())

### Name: plot.optholdoutsize_emul
### Title: Plot estimated cost function using emulation (semiparametric)
### Aliases: plot.optholdoutsize_emul
### Keywords: estimation

### ** Examples


# Simple example

# Parameters
N=100000;
k1=0.3
A=8000; B=1.5; C=0.15; theta=c(A,B,C)

# True mean function
k2_true=function(n) powerlaw(n,theta)

# Values of n for which cost has been estimated
np=50 # this many points
nset=round(runif(np,1,N))
var_k2=runif(np,0.001,0.002)
k2=rnorm(np,mean=k2_true(nset),sd=sqrt(var_k2))

# Compute OHS
res1=optimal_holdout_size_emulation(nset,k2,var_k2,N,k1)

# Plot
plot(res1)



cleanEx()
nameEx("powerlaw")
### * powerlaw

flush(stderr()); flush(stdout())

### Name: powerlaw
### Title: Power law function
### Aliases: powerlaw

### ** Examples


ncheck=seq(1000,10000)
plot(ncheck, powerlaw(ncheck, c(5e3,1.2,0.3)),type="l",xlab="n",ylab="powerlaw(n)")




cleanEx()
nameEx("powersolve")
### * powersolve

flush(stderr()); flush(stdout())

### Name: powersolve
### Title: Fit power law curve
### Aliases: powersolve
### Keywords: estimation,aspre

### ** Examples


# Retrieval of original values
A_true=2000; B_true=1.5; C_true=0.3; sigma=0.002

X=1000*abs(rnorm(10000,mean=4))
Y=A_true*(X^(-B_true)) + C_true + rnorm(length(X),sd=sigma)

c(A_true,B_true,C_true)
powersolve(X[1:10],Y[1:10])$par
powersolve(X[1:100],Y[1:100])$par
powersolve(X[1:1000],Y[1:1000])$par
powersolve(X[1:10000],Y[1:10000])$par



cleanEx()
nameEx("powersolve_general")
### * powersolve_general

flush(stderr()); flush(stdout())

### Name: powersolve_general
### Title: General solver for power law curve
### Aliases: powersolve_general
### Keywords: estimation,aspre

### ** Examples


# Retrieval of original values
A_true=2000; B_true=1.5; C_true=0.3; sigma=0.002

X=1000*abs(rnorm(10000,mean=4))
Y=A_true*(X^(-B_true)) + C_true + rnorm(length(X),sd=sigma)

c(A_true,B_true,C_true)
powersolve_general(X[1:10],Y[1:10])$par
powersolve_general(X[1:100],Y[1:100])$par
powersolve_general(X[1:1000],Y[1:1000])$par
powersolve_general(X[1:10000],Y[1:10000])$par



cleanEx()
nameEx("powersolve_se")
### * powersolve_se

flush(stderr()); flush(stdout())

### Name: powersolve_se
### Title: Standard error matrix for learning curve parameters (power law)
### Aliases: powersolve_se
### Keywords: estimation,aspre

### ** Examples


A_true=10; B_true=1.5; C_true=0.3; sigma=0.1

set.seed(31525)

X=1+3*rchisq(10000,df=5)
Y=A_true*(X^(-B_true)) + C_true + rnorm(length(X),sd=sigma)

# 'Observations' - 100 samples
obs=sample(length(X),100,rep=FALSE)
Xobs=X[obs]; Yobs=Y[obs]

# True covariance matrix of MLE of a,b,c on these x values
ntest=1000
abc_mat_xfix=matrix(0,ntest,3)
abc_mat_xvar=matrix(0,ntest,3)
E1=A_true*(Xobs^(-B_true)) + C_true
for (i in 1:ntest) {
  Y1=E1 + rnorm(length(Xobs),sd=sigma)
  abc_mat_xfix[i,]=powersolve(Xobs,Y1)$par # Estimate (a,b,c) with same X

  X2=1+3*rchisq(length(Xobs),df=5)
  Y2=A_true*(X2^(-B_true)) + C_true + rnorm(length(Xobs),sd=sigma)
  abc_mat_xvar[i,]=powersolve(X2,Y2)$par # Estimate (a,b,c) with variable X
}

Ve1=var(abc_mat_xfix) # empirical variance of MLE(a,b,c)|X
Vf=powersolve_se(Xobs,Yobs,method='fisher') # estimated SE matrix, asymptotic

Ve2=var(abc_mat_xvar) # empirical variance of MLE(a,b,c)
Vb=powersolve_se(Xobs,Yobs,method='bootstrap') # estimated SE matrix, bootstrap

cat("Empirical variance of MLE(a,b,c)|X\n")
print(Ve1)
cat("\n")
cat("Asymptotic variance of MLE(a,b,c)|X\n")
print(Vf)
cat("\n\n")
cat("Empirical variance of MLE(a,b,c)\n")
print(Ve2)
cat("\n")
cat("Bootstrap-estimated variance of MLE(a,b,c)\n")
print(Vb)
cat("\n\n")




cleanEx()
nameEx("psi_fn")
### * psi_fn

flush(stderr()); flush(stdout())

### Name: psi_fn
### Title: Updating function for variance.
### Aliases: psi_fn

### ** Examples


# See examples for `mu_fn`




cleanEx()
nameEx("sens10")
### * sens10

flush(stderr()); flush(stdout())

### Name: sens10
### Title: Sensitivity at theshold quantile 10%
### Aliases: sens10
### Keywords: aspre

### ** Examples


# Simulate
set.seed(32142)

N=1000
X=rnorm(N); Y=rbinom(N,1,prob=logit(X/2))

pi_int=0.1
q10=quantile(X,1-pi_int) # 10% of X values are above this threshold

print(length(which(Y==1 & X>q10))/length(which(X>q10)))
print(sens10(Y,X,pi_int))




cleanEx()
nameEx("sim_random_aspre")
### * sim_random_aspre

flush(stderr()); flush(stdout())

### Name: sim_random_aspre
### Title: Simulate random dataset similar to ASPRE training data
### Aliases: sim_random_aspre
### Keywords: aspre

### ** Examples


# Load ASPRE related data
data(params_aspre)

X=sim_random_aspre(1000,params_aspre)

print(c(median(X$age),params_aspre$age$median))

print(rbind(table(X$parity)/1000,params_aspre$parity$freq))




cleanEx()
nameEx("split_data")
### * split_data

flush(stderr()); flush(stdout())

### Name: split_data
### Title: Split data
### Aliases: split_data
### Keywords: simulation

### ** Examples


# See examples for model_predict



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
