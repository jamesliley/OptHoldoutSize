% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OptHoldoutSize_estimation.R
\name{optimal_holdout_size}
\alias{optimal_holdout_size}
\title{Estimate optimal holdout size under parametric assumptions}
\usage{
optimal_holdout_size(
  N,
  k1,
  theta,
  k2form = powerlaw,
  round_result = FALSE,
  ...
)
}
\arguments{
\item{N}{Total number of samples on which the predictive score will be used/fitted. Can be a vector.}

\item{k1}{Cost value in the absence of a predictive score. Can be a vector.}

\item{theta}{Parameters for function k2form(n) governing expected cost to an individual sample given a predictive score fitted to n samples. Can be a matrix of dimension n x n_par, where n_par is the number of parameters of k2.}

\item{k2form}{Function governing expected cost to an individual sample given a predictive score fitted to n samples. Must take two arguments: n (number of samples) and theta (parameters). Defaults to a power-law form \verb{powerlaw(n,c(a,b,c))=a n^(-b) + c}.}

\item{round_result}{Set to TRUE to solve over integral sizes}

\item{...}{Passed to function \code{optimize}}
}
\value{
List/data frame of dimension (number of evaluations) x (4 + n_par) containing input data and results. Columns size and cost are optimal holdout size and cost at this size respectively. Parameters N, k1, theta.1, theta.2,...,theta.n_par are input data.
}
\description{
Compute optimal holdout size for updating a predictive score given appropriate parameters of cost function

Evaluates empirical minimisation of cost function \verb{l(n;k1,N,theta) = k1 n + k2form(n;theta) (N-n)}.

The function will return \code{Inf} if no minimum exists. It does not check if the minimum is unique, but this can be guaranteed using the assumptions for theorem 1 in the manuscript.

This calls the function \code{optimize} from package \code{stats}.
}
\examples{

# Evaluate optimal holdout set size for a range of values of k1 and two values of
#   N, some of which lead to infinite values
N1=10000; N2=12000
k1=seq(0.1,0.5,length=20)
A=3; B=1.5; C=0.15; theta=c(A,B,C)

res1=optimal_holdout_size(N1,k1,theta)
res2=optimal_holdout_size(N2,k1,theta)

oldpar=par(mfrow=c(1,2))
plot(0,type="n",ylim=c(0,500),xlim=range(res1$k1),xlab=expression("k"[1]),
  ylab="Optimal holdout set size")
  lines(res1$k1,res1$size,col="black")
  lines(res2$k1,res2$size,col="red")
  legend("topright",as.character(c(N1,N2)),title="N:",col=c("black","red"),lty=1)
plot(0,type="n",ylim=c(1500,1600),xlim=range(res1$k1),xlab=expression("k"[1]),
  ylab="Minimum cost")
  lines(res1$k1,res1$cost,col="black")
  lines(res2$k1,res2$cost,col="red")
  legend("topleft",as.character(c(N1,N2)),title="N:",col=c("black","red"),lty=1)

par(oldpar)
}
\keyword{estimation}
