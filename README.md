# LNA-Incidence

Computer code for implementing the ensemble MCMC method of [https://arxiv.org/abs/2303.15371](https://arxiv.org/abs/2303.15371)

The article considers parameter inference for compartment models of epidemics constrained using incidence data. A linear noise approximation of the most natural Markov jump process representation of cumulative incidence is developed. This is combined with a Gaussian approximation of various observation models to give a computationally efficient inference scheme. Specifically, an analytic approximation of the observed data likelihood is used to drive an MCMC scheme. 

The code can be used to reproduce the results for the real data application given in Section 4. 
