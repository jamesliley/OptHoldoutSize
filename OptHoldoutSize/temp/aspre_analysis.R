################################################################################
## R script to estimate parameters of learning curve for ASPRE model          ##
################################################################################
##
## James Liley
## 21 October 2021
##
## Note: we do not have access to the script used to fit the model. THIS IS A
##  SIMULATION ONLY, INTENDED ONLY TO DEMONSTRATE HOW THE ALGORITHM COULD BE
##  APPLIED. DO NOT INTERPRET ANY RESULTS AS CLINICALLY USABLE.
##
## This simulation works by taking the final ASPRE model (Akolekar et al, 2013)
##  with >55,000 samples and considering it to have optimal performance.

######################################################################
## Parameters, switches and libraries                               ##
######################################################################

# Set random seed
seed=463825
set.seed(seed)

# Libraries
library(mvtnorm)
library(matrixStats)
library(mle.tools)

#### Parameters to do with ASPRE

# Total individuals in trial; all data
n_aspre_total=58794

# Population untreated PRE prevalence
pi_PRE = 1426/58974

# Maximum score sensitivity amongst highest 10%: assumed to be that of ASPRE
sens_max = (138+194)/2707  # = 0.122645 , from abstract of Rolnik 2017 Ultrasound in O&G

# Intervene with aspirin on this proportion of individuals
pi_intervention=0.1

# Aspirin reduces PRE risk by approximately this much
alpha=0.37
SE_alpha=0.09

#### Control parameters for simulation

# Check these potential training set sizes (eg holdout set sizes)
training_set_sizes=round(seq(1000,30000,length=20))

# Run this many simulations at each training set size
n_sim=50

# Save plot to file, or not
save_plot=FALSE



######################################################################
## Functions and scripts                                            ##
######################################################################

logit=function(x) 1/(1+exp(-x))
logistic=function(x) -log((1/x)-1)

##' Computes sensitivity of a risk score at a threshold at which
##'  10% of samples (or some proportion pi_int) are above it.
##'
##'  @name sens10
##'  @param Y true labels
##'  @param Ypred predictions
##'  @param pi compute sensitivity when a proportion pi_int of samples exceed threshold, default 0.1
##'  @return sensitivity at this threshold
sens10=function(Y,Ypred,pi_int=pi_intervention) {
  qa=quantile(Ypred,1-pi_int)
  return(mean(Y*(Ypred >= qa))/pi_int)
}



##' Fit power law curve
##' Find least-squares solution: MLE of (a,b,c) under model
##'  y= ax^-b + c + e, e~N(0,s^2)
##'
##' @name powersolve
##' @param x x values
##' @param y y values
##' @param init initial values of (a,b,c) to start
##' @param ... further parameters passed to optim
##' @return MLE values of (a,b,c)
powersolve=function(x,y,init=c(1,-0.8,0.0165),...) {
  fabc=function(abc) sum( (y- (abc[1]*(x^(abc[2])) + abc[3]))^2 )
  out=optim(par=init,fn=fabc,...)
  return(out)
}


######################################################################
## ASPRE_specific definitions and functions                         ##
######################################################################

# Distribution of parameters in ASPRE model
# Population distributions from https://www.nejm.org/doi/suppl/10.1056/NEJMoa1704559/suppl_file/nejmoa1704559_appendix.pdf
params_aspre=list(
  age=list(
    median= 31.45,
    IQR=c(27.1,35.8)),
  weight=list(
    median=70.5,
    IQR=c(61.15,84.05)),
  height=list(
    median= 163,
    IQR= c(159,167.5)),
  race=list(
    name=c("Caucasian","Afro_Caribbean","S_Asian","E_Asian","Mixed"),
    freq=c(0.671,0.253,0.0455,0.0175,0.013)),
  conception=list(
    name=c("Natural","Ovulation_assisted","In_vitro"),
    freq=c(0.942,0.0085,0.0495)
  ),
  fam_hx=list(
    freq=0.0865),
  htn_hist=list(
    freq=0.0675),
  sle_apl_hist=list(
    freq=0.005),
  dm_hist=list(
    freq=0.0155),
  parity=list(
    name=c("Nulliparous","Par_no_prev_PE","Par_prev_PE"),
    freq=c(0.685,0.206,0.109)),
  years_since_last_pregnancy=list( # this is simulated for all but will only be used if multiparous. Gravidity not reported.
    median=4.4,
    IQR=c(2.7,7.25)),
  last_delivery_gestational_age=list( # this is simulated for all but will only be used if multiparous
    median=39,
    IQR=c(36.5,40)
  ),
  map=list(
    median=1.079,
    IQR=c(1.022,1.140)),
  uterine_artery_PI=list(
    median=1.238,
    IQR=c(1.0435,1.474)),
  serum_PPA=list(
    median=0.778,
    IQR=c(0.5,1.191)),
  serum_PPG=list(
    median=0.694,
    IQR=c(0.507,0.916))
)



##' Generate random population of parameters where properties match that of params
##' Assume independence of parameter variation for simplicity.
##'
##' @param n size of population
##' @param pars list of parameter
sim_random=function(n,params=params_aspre) {
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


##' Add various interaction terms to X. Interaction terms correspond to those in ASPRE,
##'
##' @name add_aspre_interactions
##' @param X data frame
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


## Full ASPRE model, from https://www.nejm.org/doi/suppl/10.1056/NEJMoa1704559/suppl_file/nejmoa1704559_appendix.pdf
# Taken as ground truth
# Model is to predict gestational age at PE, eg higher=lower PE risk, so coefficients are negated for model to predict PE risk

##' Computes ASPRE model prediction on a matrix X of covariates
##'
##' @param X matrix, assumed to be output of sim_random with parameter params=params_aspre and transformed using add_aspre_interactions
##' @return vector of scores.
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



######################################################################
## Mockup of real data                                              ##
######################################################################

X=sim_random(n_aspre_total)
X1=add_aspre_interactions(X)

# Risk will be monotonic to ASPRE risk, but we will transform to match
#  population prevalence of PE and sensitivity of ASPRE score.
risk0=aspre(X1)

# Find a linear transformation ax+b of lrisk such that population prevalence
#  and expected sensitivity match. Suppose P(Y_i=1)=score_i
# Expected sensitivity = E_{Y|scores}(sens)
#                      = (1/(pi_intervention*n_aspre_total))*E{sum_{i:score(i)>thresh} [Y_i]}
#                      = (1/5879)*sum_{i:score(i)>thresh} [(score(i)])
lrisk0=logistic(risk0)
f_ab=function(ab) {
  a=ab[1]; b=ab[2]
  risk_ab=a*lrisk0 + b
  pop_prev=mean(logit(risk_ab))
  q_pi=quantile(risk_ab,0.9)
  sens=(1/(pi_intervention*n_aspre_total))*sum(logit(risk_ab)*(risk_ab>q_pi))
  return((pop_prev-pi_PRE)^2 + (sens - sens_max)^2)
}
abmin=optim(c(1,0),f_ab)$par
lrisk=abmin[1]*lrisk0 +abmin[2]
risk=logit(lrisk)

# Set PRE
PRE=rbinom(n_aspre_total,1,prob=risk) # ASPRE=ground truth



######################################################################
## Approximate ASPRE model with logistic model                      ##
######################################################################

## ASPRE model performance is close enough to well-approximated by a logistic model
if (FALSE) {
  set.seed(seed)
  train=sample(n_aspre_total,40000); test=setdiff(1:n_aspre_total,train)
  aspre_simple=glm(PRE~.,data=cbind(X,PRE)[train,],family=binomial(link="logit"))
  ytest=predict(aspre_simple,X[test,],type="response")
  pre_test=PRE[test]
  sens_max_logistic=sens10(pre_test,ytest) # close enough to ASPRE sensitivity
}


######################################################################
## Establish learning curve: E(cost per sample) with training size  ##
######################################################################

# Number of training set sizes to check
ntrial=length(training_set_sizes)

# This will be set to performance at each trial
performance=matrix(0,ntrial,n_sim)

# Run simulation
for (i in 1:ntrial) {
  for (j in 1:n_sim) {
    set.seed(seed + i*n_sim + j)
    sub=sample(n_aspre_total,training_set_sizes[i],rep=T) # subsample training data
    csub=setdiff(1:n_aspre_total,unique(sub)) # complement; use for testing performance
    Xsub=X[sub,]; PREsub=PRE[sub]
    Xc=X[csub,]; PREc=PRE[csub]
    g1=glm(PREsub~.,data=data.frame(Xsub,PREsub),family=binomial(link="logit"))
    score_p=(predict(g1,data.frame(Xc),type="response"))
    dif=sens10(PREc,score_p)
    performance[i,j]=dif
  }
  print(paste0("Completed simulations for training set size: ",training_set_sizes[i]))
}

# Transform sensitivity to cost
cost=pi_PRE - pi_intervention*alpha*performance


# Means
mean_cost=rowMeans(cost)




######################################################################
## Estimate power law parameters and variance                       ##
######################################################################

# Estimate power law parameters
fpar=powersolve(training_set_sizes,mean_cost)$par
fx=function(x) fpar[1]*(x^fpar[2]) + fpar[3]



######################################################################
## Draw figure for learning curve                                   ##
######################################################################

if (save_plot) pdf("./learning_curve_estimate.pdf",width=4,height=4)

plot(0,xlim=range(training_set_sizes),ylim=range(c(cost)),type="n",
  xlab="Training set size",
  ylab=expression(paste("Indiv. cost ", "(", "",
    "", phantom() %prop% phantom(), " sens.", ")", "")))
for (i in 1:ntrial) points(rep(training_set_sizes[i],n_sim),cost[i,],pch=16,cex=0.5)
lines(training_set_sizes,mean_cost,col="red")
xx=seq(min(training_set_sizes),max(training_set_sizes))
lines(xx,fx(xx),col="blue")
legend("topright",c("Values","Means","Fitted","Minimum"),bty="n",
  lty=c(NA,1,1,1),pch=c(16,NA,NA,NA),pt.cex=c(0.5,NA,NA,NA),
  col=c("black","red","blue","black"))

if (save_plot) dev.off()



######################################################################
## Estimate optimum holdout size and variance                       ##
######################################################################

# Parameters
N=400000; SE_N=1500
A=fpar[1]; SE_A=0
B=-fpar[2]; SE_B=0
C=fpar[3]; SE_C=0

# Parameter calculation for k1
NICE_sensitivity=0.2
pi_1=NICE_sensitivity*(239/8875)/pi_intervention
pi_0=(1-NICE_sensitivity)*(239/8875)/(1-pi_intervention)
SE_pi_1=sqrt(pi_1*(1-pi_1)/(8875*0.1))
SE_pi_0=sqrt(pi_0*(1-pi_0)/(8875*0.9))
k1=pi_0*(1-pi_intervention) + pi_1*pi_intervention*alpha

# Standard error for k1
pi_1_s=rnorm(1000,mean=pi_1,sd=SE_pi_1)
pi_0_s=rnorm(1000,mean=pi_0,sd=SE_pi_0)
alpha_s=rnorm(1000,mean=alpha,sd=SE_alpha)
SE_k1=sd(pi_0_s*(1-pi_intervention) + pi_1_s*pi_intervention*alpha_s)




# Total loss function (per population)
fn=function(n,Nx=N,k1x=k1,Ax=A,Bx=B,Cx=C)  k1x*n + (Ax *n^(-Bx) + Cx)*(Nx - n)

# Optimal holdout set size and cost
optim_aspre=optimize(fn,c(1000,100000))
OHS_ASPRE=optim_aspre$minimum
Min_cost_ASPRE=optim_aspre$objective

# Errors: Monte Carlo
n_trial=1000
OHS_vec=rep(0,n_trial); cost_vec=rep(0,n_trial)
Nvec=rnorm(n_trial,mean=N,sd=SE_N)
k1vec=rnorm(n_trial,mean=k1,sd=SE_k1)
Avec=rnorm(n_trial,mean=A,sd=SE_A)
Bvec=rnorm(n_trial,mean=B,sd=SE_B)
Cvec=rnorm(n_trial,mean=C,sd=SE_C)
for (i in 1:n_trial) {
  optim_se=optimize(fn,c(1000,100000),Nx=Nvec[i],k1x=k1vec[i],Ax=Avec[i],Bx=Bvec[i],Cx=Cvec[i])
  OHS_vec[i]=optim_se$minimum; cost_vec[i]=optim_se$objective
}

# IQR for ASPRE
IQR_OHS_ASPRE=quantile(OHS_vec,c(0.25,0.75))
IQR_min_cost_ASPRE=quantile(cost_vec,c(0.25,0.75))




######################################################################
## Draw figure for total cost                                       ##
######################################################################

if (save_plot) pdf("./cost_function_aspre.pdf",width=4,height=4)

x=seq(1000,100000,length=100)
y=fn(x)
k1v=c(0.021,0.025)
Nv=c(N-5000,N+5000)

y_k1=matrix(0,length(x),length(k1v))
y_N=matrix(0,length(x),length(Nv))
min_k1_x=rep(0,length(k1v)); min_k1_y=rep(0,length(k1v))
min_N_x=rep(0,length(k1v)); min_N_y=rep(0,length(k1v))
for (i in 1:length(k1v)) {
  y_k1[,i]=fn(x,k1x=k1v[i])
  ox=optimize(fn,c(100,100000),k1x=k1v[i])
  min_k1_x[i]=ox$minimum; min_k1_y[i]=ox$objective
}
for (i in 1:length(Nv)) {
  y_N[,i]=fn(x,Nx=Nv[i])
  ox=optimize(fn,c(100,100000),Nx=Nv[i])
  min_N_x[i]=ox$minimum; min_N_y[i]=ox$objective
}

plot(0,xlim=range(x),ylim=range(cbind(y,y_k1,y_N)),type="n",
  xlab="Training/holdout set size",
  ylab="Total cost (= num. cases)")
for (i in 1:length(k1v)) {
  lines(x,y_k1[,i],col="red")
}
for (i in 1:length(Nv)) {
  lines(x,y_N[,i],col="blue")
}
lines(x,y,lwd=2)
points(x=OHS_ASPRE,y=Min_cost_ASPRE,pch=16,cex=1)
points(x=min_k1_x,y=min_k1_y,pch=16,cex=1,col="red")
points(x=min_N_x,y=min_N_y,pch=16,cex=1,col="blue")
legend("topright",c("Estimated",expression("Varying k"[1]),"Varying N","Minimum"),
  lty=c(1,1,1,NA),lwd=c(2,1,1,NA),pch=c(NA,NA,NA,16),pt.cex=c(NA,NA,NA,1),
  col=c("black","red","blue","black"),bg="white")

if (save_plot) dev.off()


