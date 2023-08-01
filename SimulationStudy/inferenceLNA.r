## LNA for SIR incidence processes
## Analytic integration of the latent process

# source("helper.r")

# Forward filter (binomial obs model)
filter = function(x0,ODEdt,param,a,C,F,data,afun,bfun,Ffun,S,inter)
{
 # x0 = initial state (S0,I0)
 # param = parameters (\beta,\gamma,\lambda)
 # ODEdt = Euler time step
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
  ll = ll+dnorm(data[i],param[p]*t(F)%*%mtd,sqrt((param[p]^2)*t(F)%*%Vtd%*%F+param[p]*(1-param[p])*t(F)%*%mtd),log=TRUE)
  #update posterior
  covRY = Vt-Gt%*%C; covYR = t(covRY)
  a = mt+param[p]*covRY%*%F%*%solve((param[p]^2)*t(F)%*%Vtd%*%F+param[p]*(1-param[p])*t(F)%*%mtd)%*%(data[i]-param[p]*t(F)%*%mtd)
  C = Vt-(param[p]^2)*covRY%*%F%*%solve((param[p]^2)*t(F)%*%Vtd%*%F+param[p]*(1-param[p])*t(F)%*%mtd)%*%t(F)%*%covYR
  }
 return(ll)
}

# MH scheme (FFMH)
LNAinf = function(iters,sigma,x0,x,n0,V0,a,b,init,afun,bfun,Ffun,S,inter,ODEdt,Fobs)
{
  # iters = no. of iterations
  # sigma = tuning matrix in random walk proposal (3 x 3 matrix)
  # x0 = initial condition (S0,I0)
  # x = data vector (noisy cumulative incidence in (t,t+inter])
  # n0,V0 = mean and variance matrix of initial incidence value
  # a,b = prior hyper-parameters (length-3 vectors)
  # init = initial parameter values (length-3 vector)
  # afun, bfun = drift and diffusion coefficients
  # Ffun = Jacobian matrix
  # inter = inter-obs time
  # ODEdt = ODE time step
  # Fobs = 2-vector to pick out observed state
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
     laprob = loglikecan+logprior(can,a,b)-loglikecurr-logprior(curr,a,b)
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











