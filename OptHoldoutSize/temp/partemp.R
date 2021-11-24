
spec=FALSE
## True mean cost function
if (spec==TRUE) {
  # Mean cost follows a power-law form
  true_mean=function(n) powerlaw(n,theta_true)
} else {
  # Mean cost follows a double-descent form
  true_mean=function(n) powerlaw(n,theta_true) + (2e4)*dnorm(n,mean=4e4,sd=8e3)
}
#

# Kernel width and Gaussian process variance
kw0=5000
vu0=1e7


nstart=4
vwmin=0.001; vwmax=0.005
nset0=round(runif(nstart,1000,N/2))
var_w0=runif(nstart,vwmin,vwmax)
d0=rnorm(nstart,mean=true_mean(nset0),sd=sqrt(var_w0))

nmin=c()






## Repeatedly run next part

# Now we will do the same repeatedly, adding new points to n each time.
n_iter=100 # total number of points

nset=nset0; d=d0; var_w=var_w0; theta=theta0
## Find next point to add
exp_imp_em <- exp_imp_fn(n,nset=nset,d=d,
  var_w = var_w, N=N,k1=k1,theta=theta,k_width=kw0,var_u=vu0)
n_new=n[which.max(exp_imp_em)]

## Estimate cost at new point with error
var_w_new=runif(1,vwmin,vwmax)
d_new=rnorm(1,mean=true_mean(n_new),sd=sqrt(var_w_new))

nset=c(nset,n_new)
var_w=c(var_w,var_w_new)
d=c(d,d_new)

nset_em=nset; d_em=d; var_w_em=var_w; theta1=theta

# Mean and variance
p_mu=mu_fn(n,nset=nset_em,d=d_em,var_w = var_w_em, N=N,k1=k1,theta=theta1,k_width=kw0,var_u=vu0)
p_var=psi_fn(n,nset=nset_em,N=N,var_w = var_w_em,k_width=kw0,var_u=vu0)

nmin=c(nmin,n[which.min(p_mu)])

theta=powersolve(nset,d,y_var=var_w,lower=theta_lower,upper=theta_upper,init=theta_true,control=list(parscale=theta_true))$par
nset0=nset; theta0=theta; var_w0=var_w; d0=d
yrange=c(-30000,100000)
plot(0,xlim=range(n),ylim=yrange,type="n",
  xlab="Training/holdout set size",
  ylab="Total cost (= num. cases)")
lines(n,p_mu,col="blue")
lines(n,p_mu - 3*sqrt(p_var),col="red")
lines(n,p_mu + 3*sqrt(p_var),col="red")
points(nset0,k1*nset0 + d0*(N-nset0),pch=16,col="purple")
lines(n,k1*n + powerlaw_mean_fn(n,theta0)*(N-n),lty=2)
lines(n,k1*n + true_mean(n)*(N-n),lty=3,lwd=3)
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
exp_imp_em <- exp_imp_fn(n,nset=nset0,d=d0,var_w = var_w0, N=N,k1=k1,theta=theta0,k_width=kw0,var_u=vu0)
abline(v=n[which.max(exp_imp_em)])

