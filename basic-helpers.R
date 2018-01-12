# Basic Helper Functions
# Version 1.1
# Last Updated on Jan 10, 2017

# dependencies
suppressMessages(require(compiler))

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
