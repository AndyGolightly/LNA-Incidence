##Helper functions

#Square root representation of a matrix (uses svd)
sqrtmat = function(V)
{
 spec = svd(V)
 return(spec$u%*%diag(sqrt(spec$d))%*%t(spec$u))
}

# Generates a single draw from a multivariate N(m,V) distn
rmvn = function(m,V)
{
  p = length(m)
  z = rnorm(p)
  return(m+sqrtmat(V)%*%z)
}

#Evaluate log of prior density
logprior=function(param,a,b)
{
 # a,b = prior hyper-parameters
 # param = (\log\beta,\log\gamma,\log\lambda)   
 p=length(a)
 return(sum(dgamma(exp(param[1:(p-1)]),a[1:(p-1)],b[1:(p-1)],log=TRUE))+dunif(exp(param[p]),a[p],b[p],log=TRUE)+sum(param[1:p]))
}

# Stoichiometry matrix
Sepi = matrix(c(-1,0,1,-1),nrow=2,byrow=TRUE)

## LNA helper functions

# Drift function
alpha=function(x,theta)
{
 # x = state (S,I)
 # theta = parameters (\beta,\gamma,\lambda)
 return(c(theta[1]*x[1]*x[2],theta[2]*x[2]))
}

# Diffusion matrix
beta=function(x,theta)
{
 # x = state (S,I)
 # theta = parameters (\beta,\gamma,\lambda)
 mat = matrix(0,ncol=2,nrow=2,byrow=T)
 mat[1,1] = theta[1]*x[1]*x[2]
 mat[2,2] = theta[2]*x[2]
 return(mat)
}

# Jacobian matrix
Fmat=function(x,z,theta)
{
 # x = initial state (S0,I0)
 # z = deterministic incidence process (n_{SI},n_{IR})
 # theta = parameters (\beta,\gamma,\lambda)
 F = matrix(0,ncol=2,nrow=2)
 F[1,] = c(theta[1]*(x[1]-x[2]-2*z[1]+z[2]),theta[1]*(-x[1]+z[1]))
 F[2,] = c(theta[2], -theta[2])
 return(F)
}

# Euler solve of LNA ODE system
tstep=function(T,dt,x0,initz,initV,theta,afun,bfun,Ffun,S)
{
 # T = final time
 # dt = time step
 # x0 = initial state (S0,I0)
 # initz = initial incidence state z0
 # initV = V0
 # theta = parameters (\beta,\gamma,\lambda)
 # afun, bfun, Ffun = drift, diffusion, Jacobian functions
 # S = stoichiometry matrix 
 n = T/dt
 d = length(initz)
 z = matrix(0,nrow=d,ncol=n+1) #store z solution here
 x = matrix(0,nrow=d,ncol=n+1) #state process
 V = array(0,dim=c(d,d,n+1)) #store V solution here
 G = array(0,dim=c(d,d,n+1)) #store G solution here
 #Initialise
 z[,1] = initz
 V[,,1] = initV
 G[,,1] = diag(rep(1,d))
 x[,1] = x0
 for(i in 2:(n+1))
 {
  # Time step LNA ODEs
  Fm = Ffun(x0,z[,i-1],theta)
  G[,,i] = G[,,i-1]+Fm%*%G[,,i-1]*dt
  z[,i] = z[,i-1]+afun(x[,i-1],theta)*dt
  V[,,i] = V[,,i-1]+(V[,,i-1]%*%t(Fm)+bfun(x[,i-1],theta)+Fm%*%V[,,i-1])*dt
  # Update prevalence x using incidence z
  x[,i] = x0+S%*%z[,i]
 }
 return(list(z,V,x,G))
}

## (C)PMMH helper functions

#Euclidean distance
euclid=function(vec)
{
return(sqrt(sum(vec^2)))
}

#returns indices of Euclidean sorted (smallest to largest) vectors stored in the rows of xmat
esort=function(xmat)
{
 N = dim(xmat)[1]
 indices = rep(0,N)
 iset = 1:N
 indices[1] = which.min(xmat[,1])
 for(j in 2:N)
 {
  xstar = xmat[indices[j-1],]
  iset = iset[iset!=indices[(j-1)]]
  dist = rep(0,length(iset))
  for(i in 1:length(iset))
  {
   dist[i] = euclid(xstar-xmat[iset[i],])
  }
  ind = which.min(dist)
  indices[j] = iset[ind]
 }
 return(indices)
}

#Systematic resampling - returns indices of resamples
sysresamp=function(wts,N,uni)
{
 # wts = unormalised weight vector
 # N = number of resamples required
 # uni = uniform random draw
 vec = rep(0,N)
 wsum = sum(wts)
 k = 1
 u = uni/N
 wsumtarg = u
 wsumcurr = wts[k]/wsum
 delta = 1/N
 for(i in 1:N)
 {
  while (wsumcurr<wsumtarg)
  {
   k = k+1
   wsumcurr = wsumcurr+wts[k]/wsum
  }   
  vec[i] = k 
  wsumtarg = wsumtarg+delta
 }
 return(vec)
}

