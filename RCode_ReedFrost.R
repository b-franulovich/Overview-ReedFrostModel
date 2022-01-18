########################
### Reed-Frost Model ###
########################
##    Deterministic   ##
########################
rf.deterministic<-function(p,I0,S0){
  q<-1-p
  I=rep(0,15)
  S=rep(0,15)
  S[1]=S0
  I[1]=I0
  for (t in 1:15){
    I[t+1]=S[t]*(1-q^(I[t]))
    S[t+1]=S[t]-I[t+1]
  }
  list(I0=I0,S0=S0,I=I,S=S)
}

determ.results<-rf.deterministic(p=0.01,I0=1,S0=200)
  
plot(determ.results$I,ylab="Number of Infected Individuals",xlab="Time",
     main="Deterministic Model of Infected Individuals")
lines(determ.results$I,lty=2,col="#6666FF")
plot(determ.results$S,ylab="Number of Susceptible Individuals",xlab="Time",
     main="Deterministic Model of Susceptible Individuals")
lines(determ.results$S,lty=2,col="#9900CC")

########################
##      Stochastic    ##
########################
##       One Rep      ##
########################
p=0.01
q=1-p
II=rep(0,15)
SS=rep(0,15)
SS[1]=100
II[1]=1

rf.deterministic<-func
for (i in 1:10){
  II[i+1]=rbinom(1,SS[i],(1-q^(II[i])))
  SS[i+1]=SS[i]-II[i+1]
}

plot(II,type="l")

###########################
##      Stochastic       ##
###########################
##    Many Simulations   ##
###########################
# S0 initial number of susceptibles
# I0 initial number of infectious
# n number of trials

reedfrost <- function(p, I0, S0, n) {
  q<-1-p
  s<-St<-rep(S0,n)
  i<-It<-rep(I0,n)

  time<-0

  while (sum(It)>0){
    It<-rbinom(n,St,1-(q^It))
    St<-St-It
    i<-rbind(i,It)
    s<-rbind(s,St)
    time<-time+1
  }
  i<-as.matrix(i)
  s<-as.matrix(s)
  list(I0=I0,S0=S0,p=p,n=n,i=i,s=s)
}


## Number of infected plot ## 
infect.plot<-function(rfsim){
  t<-nrow(rfsim$i)-1 # time it takes to reach zero infectious
  plot(x=(0:t),y=rfsim$i[,1],type="l",ylim=c(0,max(rfsim$i)),
       xlab="Time",ylab="Number of Infectious Individuals")
  n<-ncol(rfsim$i)
  if (ncol(rfsim$i)>1){
    for (i in 2:n){
      lines(x=0:t,y=rfsim$i[,i],col=i)
    }
  }
}

## Number of susceptible plot ## 
sus.plot<-function(rfsim){
  t<-nrow(rfsim$s)-1 # time it takes to reach zero infectious
  plot(x=(0:t),y=rfsim$s[,1],type="l",ylim=c(0,max(rfsim$s)),
       xlab="Time",ylab="Number of Susceptible Individuals")
  n<-ncol(rfsim$s)
  if (ncol(rfsim$s)>1){
    for (i in 2:n){
      lines(x=0:t,y=rfsim$s[,i],col=i)
    }
  }
}

## Example ## 
rf.result <- reedfrost(p=0.01, I0=1, S0=200, n=10)
infect.plot(rf.result)
sus.plot(rf.result)

###########################
##         MLE           ##
##     estimating p      ##
###########################
##    simulated data     ##
###########################
## p=0.05, I0=3, S0=500 ##
###########################
rf.result.determ<-rf.deterministic(p=0.05,I0=3,S0=500)
rf.result.stoch <- reedfrost(p=0.05, I0=3, S0=500, n=10)
I1<-c(3,65,425,10,0) 
S1<-c(500,435,10,0,0)

rf.llike<-function(q,I,S){
  LL<-rep(0,length(I)-1)
  for (t in 1:(length(I)-1)){
  
    LL[t]<-lfactorial(S[t])-lfactorial(I[t+1])-lfactorial(S[t]-I[t+1])+
      (I[t+1])*log(1-q^I[t])+
        I[t]*(S[t]-I[t+1])*log(q)
  }
  LL<-list(LL)
  LL<-rapply( LL, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
  LL<-unlist(LL)
  sum(LL)
}


MLE=optimize(f=rf.llike,I=I1,S=S1,lower=0,upper=1,
            maximum=TRUE)
qhat=MLE$maximum                    #  ML estimate of q
qhat
phat=1-qhat
phat
llikeval=MLE$objective


qs<-seq(0,1,by=0.001)
rf.llikevals<-sapply(qs,rf.llike,I=I1,S=S1)

plot(qs,rf.llikevals,type="l",main="Log-Likelihood of Reed-Frost Model",
     xlab="q (1-p)",
     ylab="Log-Likelihood")
abline(v=qhat,lty=2)
abline(h=llikeval,lty=2)







