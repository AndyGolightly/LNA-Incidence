## Inference

#Evaluate log of prior density
logpriorOPM=function(param,a,b)
{
 # a,b = prior hyper-parameters
 # param = (\log\gamma,...,\log\lambda)  
 p=length(a)
 return(sum(dnorm(param[1:(p-1)],a[1:(p-1)],b[1:(p-1)],log=TRUE))+dunif(exp(param[p]),a[p],b[p],log=TRUE)+param[p])
}

# Forward filter (binomial obs model)
filter = function(x0,ODEdt,param,a,C,F,data,afun,bfun,Ffun,S,inter)
{
 # a = mean of n0
 # C = variance of n0
 # F = picks out observed component of N_t
 # inter = inter-observation time
 p = length(param)
 d = length(a)
 ODEend = inter/ODEdt +1
 n = length(data) #no. obs
 mt = rep(0,d); Vt = matrix(0, ncol=d, nrow=d); Gt = matrix(0, ncol=d, nrow=d)
 ll = 0
 #loop
 for(i in 1:n)
 {
  #update m and V
  lnastep = tstep(inter,ODEdt,x0,a,C,param,afun,bfun,Ffun,S)    
  mt = lnastep[[1]][,ODEend]
  Vt = lnastep[[2]][,,ODEend]
  Gt = lnastep[[4]][,,ODEend]
  mtd = mt-a
  Vtd = Vt+C-Gt%*%C-C%*%t(Gt)
  #update marginal likelihood
  ll = ll+dnorm(data[i],param[p]*t(F)%*%mtd,sqrt((param[p]^2)*t(F)%*%Vtd%*%F+param[p]*(1-param[p])*t(F)%*%mtd),log=T)
  #update posterior
  covRY = Vt-Gt%*%C; covYR = t(covRY)
  a = mt+param[p]*covRY%*%F%*%solve((param[p]^2)*t(F)%*%Vtd%*%F+param[p]*(1-param[p])*t(F)%*%mtd)%*%(data[i]-param[p]*t(F)%*%mtd)
  C = Vt-(param[p]^2)*covRY%*%F%*%solve((param[p]^2)*t(F)%*%Vtd%*%F+param[p]*(1-param[p])*t(F)%*%mtd)%*%t(F)%*%covYR
  }
 return(ll)
}

# FFBS (binomial obs model)
filterFFBS = function(x0,ODEdt,param,a,C,F,data,afun,bfun,Ffun,S,inter)
{
 d = length(a)
 dx = length(x0)
 p = length(param)
 ODEend = inter/ODEdt +1
 n = length(data) #no. obs
 mt = rep(0,d); Vt=matrix(0, ncol=d, nrow=d); Gt=matrix(0, ncol=d, nrow=d)
 rsamp = matrix(0,ncol=n,nrow=d); xsamp=matrix(0,ncol=n,nrow=d)
 mtA = array(0,dim=c(d,1,n)); VtA = array(0,dim=c(d,d,n)); #Store m and V here
 atA = array(0,dim=c(d,1,n)); CtA = array(0,dim=c(d,d,n)); #Store a and C here
 GtA = array(0,dim=c(d,d,n));  #Store G here
 ll = 0
 atA[,,1] = a; CtA[,,1] = C

 #loop
 for(i in 1:n)
 {
  #update m and V
  lnastep=tstep(inter,ODEdt,x0,a,C,param,afun,bfun,Ffun,S)    
  mt = lnastep[[1]][,ODEend]
  Vt = lnastep[[2]][,,ODEend]
  Gt = lnastep[[4]][,,ODEend]
  mtd = mt-a
  Vtd = Vt+C-Gt%*%C-C%*%t(Gt)
  mtA[,,(i)] = mt; VtA[,,(i)] = Vt; GtA[,,(i)] = Gt
  #update marginal likelihood
  ll = ll+dnorm(data[i],param[p]*t(F)%*%mtd,sqrt((param[p]^2)*t(F)%*%Vtd%*%F+param[p]*(1-param[p])*t(F)%*%mtd),log=T)
  #update posterior
  covRY = Vt-Gt%*%C; covYR = t(covRY)
  a = mt+param[p]*covRY%*%F%*%solve((param[p]^2)*t(F)%*%Vtd%*%F+param[p]*(1-param[p])*t(F)%*%mtd)%*%(data[i]-param[p]*t(F)%*%mtd)
  C = Vt-(param[p]^2)*covRY%*%F%*%solve((param[p]^2)*t(F)%*%Vtd%*%F+param[p]*(1-param[p])*t(F)%*%mtd)%*%t(F)%*%covYR
  atA[,,(i)] = a; CtA[,,(i)] = C
  }
  #Backward sweep
  rsamp[,n] = rmvn(a,C)
  xsamp[1:(dx-1),n] = x0[1:(dx-1)]+S%*%rsamp[1:(d-1),n]; xsamp[dx,n] = rsamp[d,n]
  #Backward sampler
  for(i in (n-1):1)
  {
   #calculate ahat and Chat
   ahat = atA[,,i] + CtA[,,(i)]%*%t(GtA[,,(i+1)])%*%solve(VtA[,,(i+1)])%*%(rsamp[,i+1]-mtA[,,(i+1)])
   Chat = CtA[,,i] - CtA[,,i]%*%t(GtA[,,(i+1)])%*%solve(VtA[,,(i+1)])%*%GtA[,,(i+1)]%*%CtA[,,(i)]
   rsamp[,i]=rmvn(ahat,Chat)
   xsamp[1:(dx-1),i] = x0[1:(dx-1)]+S%*%rsamp[1:d,i]; xsamp[dx,i] = rsamp[d,i]
  }   
 xsamp = cbind(x0,xsamp)
 return(list(rsamp,xsamp,ll))
}

# MH scheme

LNAinfOPM = function(iters,sigma,x0,x,n0,V0,a,b,init,afun,bfun,Ffun,S,inter,ODEdt,Fobs)
{
  # iters = no. of iterations
  # sigma = tuning matrix in random walk proposal
  # n0,V0 = mean and variance of initial incidence value (with \log\beta0 appended)
  # a,b = prior hyper-parameters
  # init = initial parameter values
  n = length(x) #rows = no. obs
  p = length(init)  #no. params
  mat = matrix(0,ncol=p+1,nrow=iters)
  curr = log(init) #current param value (log-scale)
  loglikecurr = -1e10 #current log-likelihood value 
  mat[1,] = c(curr,0) #params on log-scale and candidate log-likelihood	
  count = 0
  for (i in 2:iters) 
   {
    #update params   
    can = rmvn(curr,sigma)
    if(exp(can[p])>1)
    {
     laprob = -1e10 #reject if lambda>1
    }else
    {
     loglikecan = filter(x0,ODEdt,exp(can),n0,V0,Fobs,x,afun,bfun,Ffun,S,inter)
     laprob = loglikecan+logpriorOPM(can,a,b)-loglikecurr-logpriorOPM(curr,a,b)
     if((abs(laprob)>1e10)||(is.nan(laprob)==TRUE)){laprob=-1e10} #hack
    } 
    u1 = runif(1)
    if (log(u1) < laprob)
    { 
      curr = can
      loglikecurr = loglikecan 
      count = count+1  #track no. acceptances
    }   
    mat[i,] = c(curr,loglikecurr)
  }
  print(count/(iters-1))
  return(mat[2:iters,])
}






