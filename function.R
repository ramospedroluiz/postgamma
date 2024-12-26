require(coda)
require(MCMCpack)
require(VGAM)     


###########################################################################
### Gibbs with Metropolis-Hasting algorithm ###
### R: Iteration Number; burn: Burn in; jump: Jump size; b= Control generation values ###
### log_posteriori: logarithm of posteriori density ###
###########################################################################

MCMC<-function(x,delta,R,burn,jump,b=1) {
igamma<-function(a,b){(pgamma(a,b, lower=FALSE)*gamma(b)) }
posterior <- function (phi,mu) {
posterior<-log(phi*trigamma(phi)-1)-0.5*log(phi)+(d*phi*log(mu))-n*lgamma(phi)+((phi-1)*sum(delta*log(t)))-mu*sum(delta*t)+sum((1-delta)*log(igamma(mu*t,phi)))
return(posterior) }
  mu<-length(R+1); vphi<-length(R+1); n<-length(x)
  d<-sum(delta)
  #Compute the Initial values
  v1<-sum(t*log(t))/sum(t)
  vphi[1]<-(((d-2))-sum((1-delta)*log(t))+(n-d)*v1)/(n*v1-sum(log(t)))  
  mu[1]<-((n*(vphi[1]-1)+d)/sum(t)) 
  c1<-rep(0,times=R)    
  c2<-rep(0,times=R)    
  ## Realizando o M-H hibrido
  a1<-0; a2<-0; i<-1; c10<-0
  try(
    while (i<=R) {
      if(i<1) i<-2
      prop1<-rgamma(1,shape=b*vphi[i],rate=b)
      ratio1<-posterior(prop1,mu[i])-posterior(vphi[i],mu[i])+dgamma(vphi[i],shape=b*prop1,rate=b,log=TRUE)-dgamma(prop1,shape=b*vphi[i],rate=b,log=TRUE)
      alpha1<-min(1,exp(ratio1)); u1<-runif(1)
      if (u1<alpha1 & alpha1!="NaN" & alpha1!="Inf" & alpha1!="-Inf") {vphi[i+1]<-prop1 ; c1[i]<-0 ; a1<-0} else {vphi[i+1]<-vphi[i] ; a1<-a1+1; c1[i]<-1}
      
      prop2<-rgamma(1,shape=b*mu[i],rate=b)
      ratio2<-posterior(vphi[i],prop2)-posterior(vphi[i],mu[i])+dgamma(mu[i],shape=b*prop2,rate=b,log=TRUE)-dgamma(prop2,shape=b*mu[i],rate=b,log=TRUE)
      alpha2<-min(1,exp(ratio2)); u2<-runif(1)
      if (u2<alpha2 & alpha2!="NaN" & alpha2!="Inf" & alpha2!="-Inf") {mu[i+1]<-prop2 ; c2[i]<-0 ; a2<-0} else {mu[i+1]<-mu[i] ; a2<-a2+1; c2[i]<-1}
      if(a1==40) {i<-i-50; a1=0}
      i<-i+1 ; c10<-c10+1
      if(c10==500000) {i=R+1; alfa<-rep(0,times=R); vphi<-rep(0,times=R)}
    })
  try(vvphi<-vphi[seq(burn,R,jump)])
  try(vmu<- mu[seq(burn,R,jump)])
  ace1<- (1-sum(c1)/length(c1))
  atc1<- mean(acf(vvphi,plot=F)$acf)
  atc2<- mean(acf(vmu,plot=F)$acf)
  ge1<-abs(geweke.diag(vvphi)$z[1])
  ge2<-abs(geweke.diag(vmu)$z[1])
  prai<-quantile(vvphi, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
  pras<-quantile(vvphi, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
  prli<-quantile(vmu, probs = 0.025, na.rm = FALSE,names = FALSE,type = 7)
  prls<-quantile(vmu, probs = 0.975, na.rm = FALSE,names = FALSE,type = 7)
  auxvphi<-mean(vvphi); auxalfa<-mean(vmu); 
  return(list(acep=(1-sum(c1)/length(c1)),vphi=auxvphi,alpha=auxalfa, CIL_vphi=prai,CIS_vphi=pras, CIL_alpha=prli,CIS_alpha=prls, Geweke.statistics=ge1))
}

################################################################
## Example ###
## t: data vector ###
pphi <- 3       # Parametro mu
pmu <- 2       # Parametro alfa 
n <-50      # Tamanho da amostra
maux<-5.12
################################################################
t<-c()
y<-rgamma(n, pphi, rate=pmu) 
delta<-rep(1,n) # Gera vetor de tamanho n com valores 1
cax<-runif(n,0,maux)
for (i in 1:n) {
  if (y[i]<=cax[i]) {t[i]<-y[i] ;delta[i]<-1 }
  if (y[i]>cax[i]) {t[i]<-cax[i] ;delta[i]<-0 }
}
d<-sum(delta)

MCMC(t,delta,R=5500,burn=500,jump=5)
