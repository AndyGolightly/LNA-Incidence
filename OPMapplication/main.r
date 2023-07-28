## Run LNA inference code

source("helper.r")
source("SIRS.r")
source("inference.r")

# Read in data
# OPM removals (trees), Richmond Park
OPMdata = c(1024,1414,958,540,594,557,587,1029)

## SIR model

# Initial parameter values (log scale) and tuning matrix for SIR model 
init1 = c(-0.1117046, -0.6068124, -0.4847999)
sigtune1 = matrix(c(0.13097550, -0.08678061, -0.08860525, -0.08678061, 
0.21274618,  0.06128517, -0.08860525,  0.06128517,  0.07105732), ncol=3)

iters = 1000 # no. iterations
x0 = c(40000-1400,1400,-10) # initial state (S0,I0,\log\beta0)
n0 = c(0,0,-10) # initial incidence mean (with \log\beta0 appended)
V0 = diag(c(0,0,0)) # initial incidence variance
a = c(0,1,0) # prior hyper-parameters
b = c(0.5,1,1)
inter = 1 # inter-observation time
ODEdt = 0.01 # Euler time step
Fobs = c(0,1,0) # pick out observed compnent -- removals in this case

# Run MH scheme (posterior output on log scale)
set.seed(1)
system.time(outOPM<-LNAinfOPM(iters,0.9*sigtune1,x0,OPMdata,n0,V0,a,b,exp(init1),alpha,beta,Fmat,Sepi,inter,ODEdt,Fobs))

# Calculate DIC
DIC1 = -4*mean(outOPM[,4])+2*filter(x0,ODEdt,apply(exp(outOPM[,1:3]),2,mean),n0,V0,Fobs,OPMdata,alpha,beta,Fmat,Sepi,inter)

## SIRS model 

# Initial parameter values (log scale) and tuning matrix for SIRS model 
init2 = c(0.01872092, -0.96796443, -0.82862965, -0.61552812)
sigtune2 = matrix(c(0.09880592, -0.06196876, -0.07273215, -0.07194142, 
-0.06196876, 1.02800709, 0.06982288, 0.04398925, 
-0.07273215, 0.06982288, 0.14309788, 0.05415371, 
-0.07194142, 0.04398925, 0.05415371, 0.05813329),ncol=4)

iters = 1000
x0 = c(40000-1400,1400,-10)
n0 = c(0,0,0,-10)
V0 = diag(c(0,0,0,0))
a = c(0,0,1,0)
b = c(0.5,1,1,1)
inter = 1
ODEdt = 0.01
Fobs = c(0,1,0,0)

# Run MH scheme (posterior output on log scale)
set.seed(1)
system.time(outOPM2<-LNAinfOPM(iters,0.8*sigtune2,x0,OPMdata,n0,V0,a,b,exp(init2),alpha2,beta2,Fmat2,Sepi2,inter,ODEdt,Fobs))

# Calculate DIC
DIC2 = -4*mean(outOPM2[,5])+2*filter(x0,ODEdt,apply(exp(outOPM2[,1:4]),2,mean),n0,V0,Fobs,OPMdata,alpha2,beta2,Fmat2,Sepi2,inter)



