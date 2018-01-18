# Basic Helper Functions
# Version 1.1
# Last Updated on Jan 10, 2017

# dependencies
suppressMessages(require(pracma))
suppressMessages(require(compiler))
suppressMessages(require(pbivnorm))

compressCount <- function(y){
  # compress multivariate count vectors into c(count,value) pairs
  # Args:
  #   y: vector or matrix of (non-negative integer) quantiles. 
  #      If x is a matrix, each row is taken to be a quantile.
  
  y = as.matrix(y) # force 1-d vector into a column vector
  
  pool = unique(y, MARGIN = 1) # grab unique count vector values
  
  ret = matrix(0, ncol(pool) + 1, nrow(pool))
  for (i in 1:nrow(pool)){
    ret[, i] = c(sum(apply(y, 1, function(x) {
      all(x == pool[i,])
    })), pool[i,])
  }
  
  return(ret)
}
compressCount <- cmpfun(compressCount)

lowerOffDiagonal <- function(S) {
  # Transform the lower off-diagonal elements of a symmetric matrix into 
  # a vector.
  # Args:
  #   S: a symmetric matrix.
  
  p = ncol(S)
  ret = NULL
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      ret = c(ret, S[i,j])
    }
  }
  return(ret)
}
lowerOffDiagonal <- cmpfun(lowerOffDiagonal)

vec2Corr <- function(rhos, p) {
  # Convert a vector into a correlation matrix.
  # Args:
  #   rhos: a vector containing the off-diagonal elements of the correlation matrix.
  #   p: dimension of the correlation matrix.
  
  if (length(rhos) != p*(p-1)/2) {
    stop("number of off-diagonal elements does not match the target dimension!")
  }
  
  ret = matrix(1,p,p)
  counter = 0
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      counter = counter+1
      ret[i,j] = rhos[counter]
      ret[j,i] = ret[i,j]
    }
  }
  return(ret)
}
vec2Corr <- cmpfun(vec2Corr)

genLLik <- function(compressed_y, genIndividualLik, ...){
  # Compute the data log-likelihood.
  # Args:
  #   compressed_y: matrix of compressed multivariate counts. Each column is
  #                 taken to be a quantile.
  #   genIndividualLik: function that computes the likelihood for a quantile.
  
  ret = 0
  for (i in 1:ncol(compressed_y)){
    ret = ret+compressed_y[1, i]*
      log(genIndividualLik(y=compressed_y[-1, i], ...))
  }
  return(ret)
}
genLLik <- cmpfun(genLLik)

genLLikGrad <- function(compressed_y, genIndividualLik, genIndividualLikGrad, ...){
  # Compute the gradient of the data log-likelihood.
  # Args:
  #   compressed_y: matrix of compressed multivariate counts. Each column is
  #                 taken to be a quantile.
  #   genIndividualLik: function that computes the likelihood for a quantile.
  #   genIndividualLikGrad: function that computes the gradient of the likelihood
  #                         for a quantile.
  
  ret = 0
  for (i in 1:ncol(compressed_y)){
    ret = ret+compressed_y[1,i]*genIndividualLikGrad(y=compressed_y[-1, i], ...)/
      genIndividualLik(y=compressed_y[-1, i], ...)
  }
  return(ret)
}
genLLikGrad <- cmpfun(genLLikGrad)

forcePSD <- function(S) {
  # Force a matrix to be positive semi-definite by truncating the negative 
  # singular values to zero. This is the best positive semi-definite in a 
  # spectrum norm sense.
  # Args:
  #   S: matrix
  
  r <- eigen(S)
  D <- diag(as.numeric(r$values>0)*r$values)
  return(list(is_psd = all(r$values>0),
              C = r$vectors%*%D%*%solve(r$vectors)))
}

Cov2CorrMat <- function(S) {
  # Convert a covariance matrix to a correlation matrix.
  # Args:
  #   S: p.s.d. matrix.
  
  d = diag(S)
  return(diag(1/sqrt(d))%*%S%*%diag(1/sqrt(d)))
}

pbivnormBM <- function(X,rho) {
  # wrapper function for pbivnorm.
  # Args:
  #   X: matrix of upper and lower integration limits for the CDF.
  #   rho: correlation parameter.
  
  ret = NULL
  for (i in 1:nrow(X)) {
    if (X[i,1]==-Inf | X[i,2]==-Inf) {
      ret = c(ret, 0)
    } else if (X[i,1]==Inf & X[i,2]==Inf) {
      ret = c(ret, 1)
    } else {
      ret = c(ret, pbivnorm(x=X[i,1],y=X[i,2],rho=rho))
    }
  }
  return(ret)
}
pbivnormBM <- cmpfun(pbivnormBM)

rTruncatedNormal <- function(n, l, u, mu, s) {
  # sample from truncated normal distribution.
  # Args:
  #   n: number of samples.
  #   l: lower bound of the truncation.
  #   u: upper bound of the truncation.
  #   mu: mean.
  #   s: variance.
  
  pl = pnorm(l, mu, sqrt(s))
  pu = pnorm(u, mu, sqrt(s))
  
  # handle edge cases
  if (pu==0) {
    return(u)
  }
  if (pl==1) {
    return(l)
  }
  
  psample = runif(n, pl, pu)
  return(qnorm(psample, mu, sqrt(s)))
}

quad2Log <- function(logIntegrand, n, xa, xb, ya, yb, ...) {
  # Calculate the log integral via two-dimensional Gaussian quadrature.
  # Args:
  #   logIntegrand: the log of the target integrand function.
  #   n: number of nodes used per direction.
  #   xa, ya: lower limits of integration; must be finite.
  #   xb, yb: upper limits of integration; must be finite.
  #   ...: additional arguments to be passed to logIntegrand.
  
  cx <- gaussLegendre(n, xa, xb)
  cy <- gaussLegendre(n, ya, yb)
  xygrid = expand.grid(x=cx$x, y=cy$x)
  wxygrid = expand.grid(x=cx$w, y=cy$w)
  
  tmp = log(wxygrid$x)+log(wxygrid$y)+logIntegrand(xygrid$x, xygrid$y, ...)
  return(min(tmp)+log(sum(exp(tmp-min(tmp)))))
}