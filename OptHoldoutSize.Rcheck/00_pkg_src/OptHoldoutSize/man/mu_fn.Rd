% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OptHoldoutSize_emulation.R
\name{mu_fn}
\alias{mu_fn}
\title{Updating function for mean.}
\usage{
mu_fn(
  n,
  nset,
  k2,
  var_k2,
  N,
  k1,
  var_u = 1e+07,
  k_width = 5000,
  k2form = powerlaw,
  theta = powersolve(nset, k2, y_var = var_k2)$par
)
}
\arguments{
\item{n}{Set of training set sizes to evaluate}

\item{nset}{Training set sizes for which k2() has been evaluated}

\item{k2}{Estimated k2() values at training set sizes \code{nset}}

\item{var_k2}{Variance of error in k2() estimate at each training set size.}

\item{N}{Total number of samples on which the model will be fitted/used}

\item{k1}{Mean cost per sample with no predictive score in place}

\item{var_u}{Marginal variance for Gaussian process kernel. Defaults to 1e7}

\item{k_width}{Kernel width for Gaussian process kernel. Defaults to 5000}

\item{k2form}{Functional form governing expected cost per sample given sample size. Should take two parameters: n (sample size) and theta (parameters). Defaults to function \code{powerlaw}.}

\item{theta}{Current estimates of parameter values for k2form. Defaults to the MLE power-law solution corresponding to n,k2, and var_k2.}
}
\value{
Vector Mu of same length of n where Mu_i=mean(posterior(cost(n_i)))
}
\description{
Posterior mean for emulator given points \code{n}.
}
\examples{

# Suppose we have population size and cost-per-sample without a risk score as follows
N=100000
k1=0.4

# Kernel width and variance for GP
k_width=5000
var_u=8000000

# Suppose we begin with k2() estimates at n-values
nset=c(10000,20000,30000)

# with cost-per-individual estimates
# (note that since empirical k2(n) is non-monotonic, it cannot be perfectly
#  approximated with a power-law function)
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
         expression(mu(n) \%+-\% 3*sqrt(psi(n))),
         "prior(n)",
         "d",
         "3SD(d|n)"),
       lty=c(1,1,2,NA,NA),lwd=c(1,1,1,NA,NA),pch=c(NA,NA,NA,16,124),
       pt.cex=c(NA,NA,NA,1,1),
       col=c("blue","red","black","purple","black"),bg="white")
}
