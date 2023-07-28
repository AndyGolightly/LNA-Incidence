## LNA for SIR and SIRS incidence processes
## Drift and diffusion functions, RHS of ODE solution

## SIR model

# Stoichiometry matrix
Sepi = matrix(c(-1,0,1,-1),nrow=2,byrow=TRUE)

# Drift function
alpha = function(x,theta)
{
 # x = state (S,I,\log\beta)
 # theta = parameters (\gamma,\sigma_{\beta},\lambda)
 c(exp(x[3])*x[1]*x[2],theta[1]*x[2],0)
}

# Diffusion matrix
beta = function(x,theta)
{
 # x = state (S,I,\log\beta)
 # theta = parameters (\gamma,\sigma_{\beta},\lambda)
 mat=matrix(0,ncol=3,nrow=3,byrow=T)
 mat[1,1]=exp(x[3])*x[1]*x[2]
 mat[2,2]=theta[1]*x[2]
 mat[3,3]=theta[2]^2
 return(mat)
}

# Jacobian matrix
Fmat = function(x,z,theta)
{
 # x = initial state (S0,I0,\log\beta0)
 # z = deterministic incidence process appended with \log\beta (n_{SI},n_{IR},\log\beta)
 # theta = parameters (\gamma,\sigma_{\beta},\lambda)
 F = matrix(0,ncol=3,nrow=3)
 F[1,] = c(exp(z[3])*(x[1]-x[2]-2*z[1]+z[2]),exp(z[3])*(-x[1]+z[1]),exp(z[3])*(x[1]-z[1])*(x[2]+z[1]-z[2]))
 F[2,] = c(theta[1], -theta[1],0)
 return(F)
}

## SIRS 

# Stoichiometry matrix
Sepi2 = matrix(c(-1,0,1,1,-1,0),nrow=2,byrow=TRUE)

# Drift function (includes total population size N as argument)
alpha2 = function(x,theta,N=40000)
{
 # x = state (S,I,\log\beta)
 # theta = parameters (\gamma,\kappa,\sigma_{\beta},\lambda)
 c(exp(x[3])*x[1]*x[2],theta[1]*x[2],theta[2]*(N-(x[1]+x[2])),0)
}

# Diffusion function (includes total population size N as argument)
beta2=function(x,theta,N=40000)
{
 # x = state (S,I,\log\beta)
 # theta = parameters (\gamma,\kappa,\sigma_{\beta},\lambda)
 mat = matrix(0,ncol=4,nrow=4,byrow=T)
 mat[1,1] = exp(x[3])*x[1]*x[2]
 mat[2,2] = theta[1]*x[2]
 mat[3,3] = theta[2]*(N-(x[1]+x[2]))
 mat[4,4] = theta[3]^2
 return(mat)
}

# Jacobian matrix
Fmat2 = function(x,z,theta)
{
 # x = imitial state (S0,I0,\log\beta0)
 # z = deterministic incidence process appended with \log\beta (n_{SI},n_{IR},n_{RS},\log\beta)
 # theta = parameters (\gamma,\kappa,\sigma_{\beta},\lambda)
 F = matrix(0,ncol=4,nrow=4)
 F[1,] = c(exp(z[4])*(x[1]-x[2]-2*z[1]+z[2]+z[3]),exp(z[4])*(-x[1]+z[1]-z[3]),exp(z[4])*(x[2]+z[1]-z[2]),exp(z[4])*(x[1]-z[1]+z[3])*(x[2]+z[1]-z[2]))
 F[2,] = c(theta[1], -theta[1],0,0)
 F[3,] = c(0,theta[2],-theta[2],0) 
 return(F)
}

# Euler solve of LNA ODE system (SIR or SIRS)
tstep = function(T,dt,x0,initz,initV,theta,afun,bfun,Ffun,S)
{
 # T = final time
 # dt = time step
 # x0 = initial state (S0,I0,\log\beta0)
 # initz = initial incidence state z0 (appended with \log\beta0)
 # initV = V0
 # theta = parameters
 # afun, bfun, Ffun = drift, diffusion, Jacobian functions
 # S = stoichiometry matrix 
 n = T/dt
 d = length(initz)
 dx = length(x0)
 z = matrix(0,nrow=d,ncol=n+1) #store z solution here
 x = matrix(0,nrow=dx,ncol=n+1) #state process
 V = array(0,dim=c(d,d,n+1)) #store V solution here
 G = array(0,dim=c(d,d,n+1)) #store G solution here
 #Initialise
 z[,1] = initz; V[,,1] = initV; G[,,1] = diag(rep(1,d)); x[,1] = x0
 for(i in 2:(n+1))
 {
  # Time step LNA ODEs
  Fm = Ffun(x0,z[,i-1],theta)
  G[,,i] = G[,,i-1]+Fm%*%G[,,i-1]*dt
  z[,i] = z[,i-1]+afun(x[,i-1],theta)*dt
  V[,,i] = V[,,i-1]+(V[,,i-1]%*%t(Fm)+bfun(x[,i-1],theta)+Fm%*%V[,,i-1])*dt
  # Update prevalence x using incidence z
  x[1:(dx-1),i]= x0[1:(dx-1)]+S%*%z[1:(d-1),i]
  x[dx,i] = z[d,i]
 }
 return(list(z,V,x,G))
}






