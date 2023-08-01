## LNA for SIR incidence processes
## Integration of latent process via (correlated) pseudo-marginal MH 

# source("helper.r")

# Forward simulation of the LNA
forwardInc = function(x0,r0,deltat,inter,afun,bfun,Ffun,S,theta,ivec)
{
  # x0 = initial condition (S0,I0)
  # r0 = initial incidence
  # deltat = ODE time step
  # inter = length of time interval over which LNA is applied
  # afun, bfun, Ffun = drift, diffusion, Jacobian functions
  # S = stoichiometry matrix 
  # theta = parameters (\beta,\gamma,\lambda)
  # ivec = d-vector of N(0,1) quantities / innovations, d is no. of incidence components
  d = length(r0)
  tol = 1e-4
  ODEend = inter/deltat+1
  lnastep = tstep(inter,deltat,x0,r0,diag(c(0,0)),theta,afun,bfun,Ffun,S)
  z = lnastep[[1]][,ODEend]; V = lnastep[[2]][,,ODEend]
  # Generate forward simulation 
  rT = z+sqrtmat(V)%*%ivec
  for(j in 1:d){
     if(rT[j]<tol){
      rT[j] = tol #hack 
     }
    }
  return(rT)
}

#Weighted resampling on a single inter-obs interval
wrInc=function(N,x0,r0,yT,deltat,T,afun,bfun,Ffun,S,theta,imat,uni)
{
  # N = no. of particles
  # x0 = initial state (S0,I0)
  # r0 = sample of cumulative incidences up to current time t
  # yT = observed incidence (infections) over (t,t+T]
  # deltat = LNA ODE time step
  # T = inter observation time
  # afun, bfun, Ffun = drift, diffusion and Jacobian functions
  # theta = parameters
  # S = stoichiometry matrix
  # imat = N times d matrix of N(0,1) innovations
  # uni = indep. U(0,1) draw
  flag = 0
  d = dim(r0)[2] #cols = no. components	
  mat = matrix(0,nrow=N,ncol=d) #store values at end here
  wts = rep(0,N)
  #generate path and calculate weight (Binomial obs model)
  if(N==1){
    end = forwardInc(x0,r0,deltat,T,afun,bfun,Ffun,S,theta,imat)
    wts = exp(dbinom(yT,round(end[1]-r0[1,]),theta[3],log=TRUE))
    mat[1,] = end
    indices = 1
  }else{
   for(i in 1:N) 
   {
     end = forwardInc(x0,r0[i,],deltat,T,afun,bfun,Ffun,S,theta,imat[i,])
     wts[i] = exp(dbinom(yT,round(end[1]-r0[i,1]),theta[3],log=TRUE))
     if(is.nan(wts[i])==TRUE){wts[i] = 0.0}
     mat[i,] = end
   }
   if(sum(wts)==0.0){ 
    wts = rep(0.0,N)
    indices = 1:N
    flag = 1
   }else{
    sorted = esort(mat) #Euclidean sorting
    mat = mat[sorted,]; wts = wts[sorted]
    #systematic resampling
    indices = sysresamp(wts,N,uni)
   }
 } 
 return(list(mat[indices,],mean(wts),mat,flag))
}

#(C)PMMH Scheme

cpmmhLNA=function(iters,sigma,x0,r0,x,a,b,init,afun,bfun,Ffun,S,inter,ODEdt,rho,N)
{
  # iters = no. of iterations
  # sigma = 3 x 3 tuning matrix in random walk proposal
  # x0 = initial condition (S0,I0)
  # r0 = sample of initial incidences (N x 2 matrix)
  # x = data (noisy cumulative incidence in (t,t+inter])
  # a,b = prior hyper-parameters (length-3 vectors)
  # afun, bfun, Ffun = drift, diffusion and Jacobian functions
  # S = stoichiometry matrix
  # init = initial parameter values (length-3 vector)
  # inter = inter-obs time
  # ODEdt = LNA ODE time step
  # rho \in [0,1] = correlation between innovations at successive iterations (0 for PMMH)  
  # N = no. of particles
  n = length(x) #no. obs
  d = dim(r0)[2] #cols = no. components
  iarray = array(rnorm(N*d*n,0,1),dim=c(N,d,n)) # no. particles x no. incidence components x no. obs intervals 
  uvec = runif(n) #vector of uniforms, one per resampling step
  p = length(init)  #no. params
  mat = matrix(0,ncol=p+1,nrow=iters) #store output here
  curr = log(init) #current param value (log-scale)
  loglikecurr = -1000000 #log-likelihood of current value (to accept first candidate)
  mat[1,] = c(curr,0) #params on log-scale and candidate log-like	
  count = 0 #count acceptances
  for (i in 2:iters) {
    xi=r0
    can = rmvn(curr,sigma)
    if(exp(can[3])>1)
    {
     laprob = -1e10
     loglikecan = -1e10
    }else{
     iarrayprop = rho*iarray+sqrt(1-rho^2)*rnorm(N*d*n,0,1) #update innovations
     uprop = pnorm(rho*qnorm(uvec)+sqrt(1-rho^2)*rnorm(n)) #update uniforms for resampling step
     loglikecan = 0
     for(j in 1:n)
     {
      wrstuff = wrInc(N,x0,xi,x[j],ODEdt,inter,afun,bfun,Ffun,S,exp(can),iarrayprop[,,j],uprop[j])
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
      iarray = iarrayprop
      uvec = uprop
      count = count+1  #track no. acceptances
    }   
    mat[i,] = c(curr,loglikecan)
  } 
  print(count/iters)
  return(mat[2:iters,])
}






