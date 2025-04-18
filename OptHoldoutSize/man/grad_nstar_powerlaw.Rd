% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OptHoldoutSize_estimation.R
\name{grad_nstar_powerlaw}
\alias{grad_nstar_powerlaw}
\title{Gradient of optimal holdout size (power law)}
\usage{
grad_nstar_powerlaw(N, k1, theta)
}
\arguments{
\item{N}{Total number of samples on which the predictive score will be used/fitted. Can be a vector.}

\item{k1}{Cost value in the absence of a predictive score. Can be a vector.}

\item{theta}{Parameters for function k2(n) governing expected cost to an individual sample given a predictive score fitted to n samples. Can be a matrix of dimension n x n_par, where n_par is the number of parameters of k2.}
}
\value{
List/data frame of dimension (number of evaluations) x 5 containing partial derivatives of nstar (optimal holdout size) with respect to N, k1, a, b, c respectively.
}
\description{
Compute gradient of optimal holdout size assuming a power-law form of k2

Assumes cost function is \verb{l(n;k1,N,theta) = k1 n + k2(n;theta) (N-n)} with \verb{k2(n;theta)=k2(n;a,b,c)= a n^(-b) + c}
}
\examples{

# Evaluate optimal holdout set size for a range of values of k1, and compute
#  derivative
N=10000;
k1=seq(0.1,0.5,length=20)
A=3; B=1.5; C=0.15; theta=c(A,B,C)

nstar=optimal_holdout_size(N,k1,theta)
grad_nstar=grad_nstar_powerlaw(N,k1,theta)

plot(0,type="n",ylim=c(-2000,500),xlim=range(k1),xlab=expression("k"[1]),
  ylab="Optimal holdout set size")
lines(nstar$k1,nstar$size,col="black")
lines(nstar$k1,grad_nstar[,2],col="red")
legend("bottomright",c(expression("n"["*"]),
    expression(paste(partialdiff[k1],"n"["*"]))),
    col=c("black","red"),lty=1)

}
\keyword{estimation}
