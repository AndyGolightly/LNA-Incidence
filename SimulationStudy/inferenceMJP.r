## MJP for SIR incidence processes

# source("helper.r")

#Simulate MJP at discrete times

Epidisc=function(param,x,tincr,maxtime)
{    
 # param = parameters
 # x = initial state / current prevalence (St,It)
 # tincr = time increment
 # maxtime = end time
    
 len = maxtime/tincr+1; i=1
 Scurr = numeric(len); Icurr=numeric(len); discretetime=numeric(len); 
 r1 = numeric(len); r2=numeric(len);
 y = rep(0,2) # current cumulative incidence
 time = 0; nextreactiontime=0
 tol=1e-7
 while (time <= maxtime+tol)
  {
   #Calculate the overall rate
   rate = param[1]*x[1]*x[2]+param[2]*x[2]
   if (rate == 0)
   {
    #No infectives left
    nextreactiontime=maxtime+1
    while (time < nextreactiontime & time <= (maxtime+0.01))
    {
     Scurr[i] = x[1]; Icurr[i] = x[2]; discretetime[i] = time
     r1[i] = y[1]; r2[i] = y[2];
     time = time + tincr
     i = i+1
    }	  
    break
   }
   else
   {
    nextreactiontime = nextreactiontime + rexp(1, rate)
   }
   #Store values at discrete intervals
   while (time < nextreactiontime & time <= (maxtime+0.01))
   {
    Scurr[i] = x[1]; Icurr[i] = x[2]; discretetime[i] = time
    r1[i] = y[1]; r2[i] = y[2]; 
    time = time + tincr
    i = i+1
   }
   #Choose which reaction happens
   u = runif(1)
   if (u<(param[1]*x[1]*x[2])/rate)
   {
    #infection
    x[1] = x[1]-1
    x[2] = x[2]+1
    y[1] = y[1]+1
   }
   else
   {
    #Removal
    x[2] = x[2]-1
    y[2] = y[2]+1
   }
  }
 return(list(time=discretetime,Svec=Scurr,Ivec=Icurr,r1=r1,r2=r2))
}


#Weighted resampling on a single inter-obs interval
wrMJP=function(N,x0,r0,yT,T,theta,S,uni)
{
  # N = no. of particles
  # x0 = initial state (S0,I0)
  # r0 = sample of cumulative incidences up to current time t
  # yT = observed incidence (infections) over (t,t+T]
  # T = inter observation time
  # theta = parameters
  # S = stoichiometry matrix
  # uni = indep. U(0,1) draw
  flag = 0
  d = dim(r0)[2] #cols = no. components	
  mat = matrix(0,nrow=N,ncol=d) #store values at end here
  wts = rep(0,N)
  #calculate x state using cummulative incidences
  xmat = x0+S%*%t(r0)
  #generate path and calculate weight (Binomial obs model)
  if(N==1){
    end = Epidisc(theta[1:2],xmat,T,T)
    wts = exp(dbinom(yT,end$r1[length(end$r1)],theta[3],log=TRUE))
    mat[1,] = c(end$r1[length(end$r1)]+r0[1,1],end$r2[length(end$r2)]+r0[1,2])
    indices = 1
  }else{
   for(i in 1:N) 
   {
     end = Epidisc(theta[1:2],xmat[,i],T,T)
     wts[i] = exp(dbinom(yT,end$r1[length(end$r1)],theta[3],log=TRUE))
     if(is.nan(wts[i])==TRUE){wts[i] = 0.0}
     mat[i,] = c(end$r1[length(end$r1)]+r0[i,1],end$r2[length(end$r2)]+r0[i,2])
   }
   if(sum(wts)==0.0){
    wts = rep(0.0,N)
    indices = 1:N
    flag = 1
   }else{
    sorted=esort(mat)
    mat=mat[sorted,]; wts=wts[sorted]
    #systematic resampling
    indices = sysresamp(wts,N,uni)
   }
 } 
 return(list(mat[indices,],mean(wts),mat,flag))
}

#MH Scheme

pmmhMJP=function(iters,sigma,x0,r0,x,a,b,init,inter,N,S)
{
  # iters = no. of iterations
  # sigma = 3 x 3 tuning matrix in random walk proposal
  # x0 = initial condition (S0,I0)
  # r0 = sample of initial incidences (N x 2 matrix)
  # x = data (noisy cumulative incidence in (t,t+inter])
  # a,b = prior hyper-parameters (length-3 vectors)
  # init = initial parameter values (length-3 vector)
  # inter = inter-obs time
  # N = no. of particles
  n = length(x) #no. obs
  d = dim(r0)[2] #cols = no. components
  p = length(init)  #no. params
  mat = matrix(0,ncol=p+1,nrow=iters) #store output here
  curr = log(init) #current param value (log-scale)
  loglikecurr = -1000000 #log-likelihood of current value (to accept first candidate)
  mat[1,] = c(curr,0) #params on log-scale and candidate log-like	
  count = 0 #count acceptances
  for (i in 2:iters) {
    xi = r0
    can = rmvn(curr,sigma)
    if(exp(can[3])>1)
    {
     laprob = -1e10
     loglikecan = -1e10
    }else{
     loglikecan = 0
     uprop=runif(n)
     for(j in 1:n)
     {
      wrstuff = wrMJP(N,x0,xi,x[j],inter,exp(can),S,uprop[j])
      if(wrstuff[[4]]==1) #Reject move if all weights are 0
      {
       loglikecan = -1e10
       j = n
      }else{ 
       loglikecan = loglikecan+log(wrstuff[[2]])
       xi = wrstuff[[1]]
      }
     }
     laprob = loglikecan+logprior(can,a,b)-loglikecurr-logprior(curr,a,b)
    }
    u = runif(1)
    if (log(u) < laprob)
    { 
      curr = can
      loglikecurr = loglikecan 
      count=count+1  #track no. acceptances
    }   
    mat[i,] = c(curr,loglikecan)
  } 
  print(count/iters)
  return(mat[2:iters,])
}





