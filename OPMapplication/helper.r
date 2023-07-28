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


