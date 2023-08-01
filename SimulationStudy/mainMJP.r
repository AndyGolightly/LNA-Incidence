## Run MJP inference code

source("helper.r")
source("inferenceMJP.r")

# Example: data set D_2

# Read in data
obs2 = scan("obs2.dat") 

# Tuning matrix for random walk proposal
sigtuneMJP = matrix(c(0.038128127, 0.03665185, -0.003598614,
0.036651851, 0.06166513, 0.018442814,
-0.003598614, 0.01844281,  0.028652206), ncol=3)

iters = 5000 # no. iterations
x0 = c(359,1) # initial state (S0,I0)
a = c(10,10,0) # prior hyper-parameters
b = c(10^4,30,1)
init1 = c(exp(-7),3*exp(-2.5),0.8) #ground truth
inter = 10 # inter-observation time
N = 125 # no. particles
r0 = matrix(0,ncol=2,nrow=N) #intial cumulative incidence sample 

# Run MH scheme (posterior output on log scale)
set.seed(1)
system.time(outMJP<-pmmhMJP(iters,sigtuneMJP,x0,r0,obs2,a,b,init1,inter,N,Sepi))





