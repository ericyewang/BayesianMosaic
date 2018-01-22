# Helper Functions for Samplers
# Version 1.1
# Last Updated on Jan 13, 2018

# dependencies
source("~/Documents/yw_git/bayesian_mosaic/basic-helpers.R")
# source("/home/collabor/yw104/BayesianMosaic/basic-helpers.R")
suppressMessages(require(mvtnorm))
suppressMessages(require(MCMCpack))
suppressMessages(require(parallel))
ncores = detectCores() - 1

# likelihood functions via numerical integration, 1-d or 2-d
genLLikPoissonLogNormal <- function(y, mu, v, ...){
  # compute the log-likelihood function value for a individual observation from a
  # Poisson log-normal distribution via numerical integration.
  # args:
  #   y: observation.
  #   mu: mean.
  #   v: variance.
  
  integrand <- function(x){
    return(dpois(y, exp(x))*dnorm(x,mean=mu,sd=sqrt(v)))
  }
  log(integrate(integrand, lower = min(log(y+0.1), mu)-4*sqrt(v), 
            upper = max(log(y+0.1), mu)+4*sqrt(v), stop.on.error=FALSE)$value)
}

genLLikPoissonLogGaussian2D <- function(y, mu, upper_tri, ...){
  # compute the log-likelihood function value for a individual observation from a multivariate 
  # Poisson log-Gaussian distribution via numerical integration.
  # args:
  #   y: observation.
  #   mu: mean.
  #   upper_tri: the upper triangular Cholesky factorization of the covariance matrix.
  
  integrand <- function(z1, z2) {
    # change of variables
    x1 = mu[1]+upper_tri[1,1]*z1+upper_tri[2,1]*z2
    x2 = mu[2]+upper_tri[1,2]*z1+upper_tri[2,2]*z2
    dpois(y[1], exp(x1), log=TRUE) + dpois(y[2], exp(x2), log=TRUE) + 
      dnorm(z1, log=TRUE) + dnorm(z2, log=TRUE)
    # note that dmvnorm(cbind(c(x1),c(x2)), mu, S)*det(upper_tri) equals dnorm(z1)*dnorm(z2)
  }
  # finite lower and upper bounds have to be provided for quad2d
  # in our case we set it to be 10 sd away from the center to make sure that
  # the function value is negligible
  quad2dLog(integrand, 64, -10, 10, -10, 10)
}

genLLikBinomialLogNormal <- function(y, N, mu, v, ...){
  # compute the log-likelihood function value for a individual observation from a
  # Poisson log-normal distribution via numerical integration
  # args:
  #   y: observation
  #   N: number of trials
  #   mu: mean
  #   v: variance
  
  integrand <- function(x){
    return(dbinom(y, N, 1/(1+exp(-x)))*dnorm(x,mean=mu,sd=sqrt(v)))
  }
  log(integrate(integrand, lower = min(log((y+0.1)/N), mu)-4*sqrt(v), 
            upper = max(log((y+0.1)/N), mu)+4*sqrt(v), stop.on.error=FALSE)$value)
}

genLLikBinomialLogGaussian2D <- function(y, Ns, mu, upper_tri, ...){
  # compute the log-likelihood function value for a individual observation from a multivariate 
  # binomial log-Gaussian distribution via numerical integration.
  # args:
  #   y: (vector) observations.
  #   Ns: (vector) number of trials.
  #   mu: mean.
  #   upper_tri: the upper triangular Cholesky factorization of the covariance matrix.
  
  integrand <- function(z1, z2) {
    x1 = mu[1]+upper_tri[1,1]*z1+upper_tri[2,1]*z2
    x2 = mu[2]+upper_tri[1,2]*z1+upper_tri[2,2]*z2
    dbinom(y[1], Ns[1], 1/(1+exp(-x1)), log=TRUE)+
      dbinom(y[2], Ns[2], 1/(1+exp(-x2)), log=TRUE)+
      dnorm(z1, log=TRUE)+dnorm(z2, log=TRUE)
  }
  
  quad2dLog(integrand, 64, -10, 10, -10, 10)
}

genLLikRoundedNormal <- function(y, mu, v, ...){
  # compute the log-likelihood function value for a individual observation from a
  # rounded normal via numerical integration.
  # args:
  #   y: observation.
  #   mu: mean.
  #   v: variance.
  
  if (y == 0) {
    return(pnorm(0, mean=mu, sd=sqrt(v), log=TRUE))
  } else {
    py = pnorm(y, mean=mu, sd=sqrt(v), log=TRUE)
    pym1 = pnorm(y-1, mean=mu, sd=sqrt(v), log=TRUE)
    return(pym1+log(exp(py-pym1)-1))
  }
}

genLLikRoundedGaussian2D <- function(y, mu, v, rho, ...){
  # compute the log-likelihood function value for a individual observation from a 
  # rounded bivariate Gaussian via numerical integration.
  # args:
  #   y: observation.
  #   mu: mean.
  #   upper_tri_inv: inverse of the upper triangular Cholesky factorization of 
  #                  the covariance matrix.
  
  ur = y # upper right corner of the support
  ll = y-1 # lower left corner of the support
  if (y[1]==0) {
    ll[1] = -Inf
  }
  if (y[2]==0) {
    ll[2] = -Inf
  }
  ul = c(ll[1],ur[2]) # upper left corner of the support
  lr = c(ur[1],ll[2]) # lower right corner of the support
  
  # change of variables
  X = rbind((ur-mu)/sqrt(v), (ul-mu)/sqrt(v), 
            (lr-mu)/sqrt(v), (ll-mu)/sqrt(v))
  pbis = pbivnormBM(X,rho)
  return(log(pbis[1]-pbis[2]-pbis[3]+pbis[4]))
}

genGroupLLik <- function(group_compressed_y, genIndividualLLik, group_mus, ...){
  # Compute the data likelihood for each group.
  # Args:
  #   group_compressed_y: a list of matrices. Each matrix contains compressed 
  #                       multivariate counts for a group of data. Within each
  #                       matrix, each column is taken to be a quantile.
  #   genIndividualLLik: function that computes the log-likelihood for a quantile.
  
  group_mus = as.matrix(group_mus)
  ret = NULL
  for (i in 1:nrow(group_mus)){
    ret = c(ret, genLLik(compressed_y=group_compressed_y[[i]], 
                         genIndividualLLik=genIndividualLLik, mu=group_mus[i,], 
                         ...))
  }
  return(ret)
}
genGroupLLik <- cmpfun(genGroupLLik)

# gradient w.r.t. rho of the likelihood functions via numerical integration, 2-d
genLikPoissonLogGaussian2DGradient <- function(y, mu, v, upper_tri, rho, ...){
  # likelihood function value for a individual observation via numerical integration
  # args:
  #   y: observation.
  #   mu: (vector) mean.
  #   v: (vector) diagonal elements of the covariance matrix.
  #   upper_tri: the upper triangular Cholesky factorization of the covariance matrix.
  #   rho: correlation.
  
  integrand <- function(z1, z2) {
    x1 = mu[1]+upper_tri[1,1]*z1+upper_tri[2,1]*z2
    x2 = mu[2]+upper_tri[1,2]*z1+upper_tri[2,2]*z2
    x1_tilde = (x1 - mu[1]) / sqrt(v[1])
    x2_tilde = (x2 - mu[2]) / sqrt(v[2])
    gradient_term = (rho * (1 - rho^2) + x1_tilde * x2_tilde * (1 + rho^2) - rho * (x1_tilde^2 + x2_tilde^2)) / (1 - rho^2)^2
    return(dpois(y[1], exp(x1))*dpois(y[2], exp(x2))*dnorm(z1)*dnorm(z2)*gradient_term)
  }
  
  quad2d(integrand, -10, 10, -10, 10, 64)
}

genLikBinomialLogGaussian2DGradient <- function(y, Ns, mu, v, upper_tri, rho, ...){
  # likelihood function value for a individual observation via numerical integration
  # args:
  #   y: observation.
  #   Ns: (vector) number of trials.
  #   mu: (vector) mean.
  #   v: (vector) diagonal elements of the covariance matrix.
  #   upper_tri: the upper triangular Cholesky factorization of the covariance matrix.
  #   rho: correlation.
  
  integrand <- function(z1, z2) {
    x1 = mu[1]+upper_tri[1,1]*z1+upper_tri[2,1]*z2
    x2 = mu[2]+upper_tri[1,2]*z1+upper_tri[2,2]*z2
    x1_tilde = (x1 - mu[1]) / sqrt(v[1])
    x2_tilde = (x2 - mu[2]) / sqrt(v[2])
    gradient_term = (rho * (1 - rho^2) + x1_tilde * x2_tilde * (1 + rho^2) - rho * (x1_tilde^2 + x2_tilde^2)) / (1 - rho^2)^2
    return(dbinom(y[1], Ns[1], 1/(1+exp(-x1)))*dbinom(y[2], Ns[2], 1/(1+exp(-x2)))*dnorm(z1)*dnorm(z2)*gradient_term)
  }
  quad2d(integrand, -10, 10, -10, 10, 64)
}

genLikRoundedGaussian2DGradient <- function(y, mu, v, rho, ...){
  # likelihood function value for a individual observation via numerical integration
  # args:
  #   y: observation.
  #   Ns: (vector) number of trials.
  #   mu: (vector) mean.
  #   v: (vector) diagonal elements of the covariance matrix.
  #   upper_tri: the upper triangular Cholesky factorization of the covariance matrix.
  #   rho: correlation.
  
  small_value = 1E-8
  (exp(genLLikRoundedGaussian2D(y=y, mu=mu, v=v, rho=rho+small_value, ...))-
      exp(genLLikRoundedGaussian2D(y=y, mu=mu, v=v, rho=rho-small_value, ...)))/2/small_value
}

getFcn <- function(model, fcn_name) {
  # Get the correct function for specified model
  # Args:
  #   model: (string) model name.
  #   fcn_name: (string) function name.
  
  if (model=="mvtPoisson") {
    genIndividualLLik1D = genLLikPoissonLogNormal
    genIndividualLLik2D = genLLikPoissonLogGaussian2D
    genIndividualLikGrad2D = genLikPoissonLogGaussian2DGradient
  } else if (model=="mvtBinomial") {
    Ns = args[[1]]$Ns
    genIndividualLLik1D = genLLikBinomialLogNormal
    genIndividualLLik2D = genLLikBinomialLogGaussian2D
    genIndividualLikGrad2D = genLikBinomialLogGaussian2DGradient
  } else if (model=="roundedGaussian") {
    genIndividualLLik1D = genLLikRoundedNormal
    genIndividualLLik2D = genLLikRoundedGaussian2D
    genIndividualLikGrad2D = genLikRoundedGaussian2DGradient
  }
  
  if (fcn_name=="IndividualLLik1D") {
    return(genIndividualLLik1D)
  } else if (fcn_name=="IndividualLLik2D") {
    return(genIndividualLLik2D)
  } else if (fcn_name=="IndividualLikGrad2D") {
    return(genIndividualLikGrad2D)
  }
}

# Sample rho via Laplacian Approximation
findZeroGradient <- function(fx) {
  # find the zero location of a monotone function
  # args:
  #   fx: some monotone function that crosses zero
  
  left = -.99
  right = .99
  counter = 0
  while (right - left >= 1E-2) {
    counter = counter + 1
    if (counter > 100) {
      stop("findZeroGradient not converging!")
    }
    mid = (right + left) / 2
    if (fx(mid) > 0) {
      right = mid
    } else {
      left = mid
    }
  }
  return(c(left, (fx(right) - fx(left)) / (right - left)))
}

sampleLaplaceApprox <- function(compressed_y, model, mu, vvars, ...) {
  # sample from the Laplace approximation
  # args:
  #   group_summary: aggregated summaries for the dataset
  #   mu: mean vector
  #   vvars: vector of variances
  
  genIndividualLLik2D = getFcn(model, "IndividualLLik2D")
  genIndividualLikGrad2D = getFcn(model, "IndividualLikGrad2D")
  x_sample = -Inf
  
  llikGradient <- function(rho) {
    Sigma = diag(sqrt(vvars)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sqrt(vvars))
    upper_tri = chol(Sigma)
    return(-genLLikGrad(compressed_y=compressed_y, genIndividualLik=genIndividualLLik2D, 
                        genIndividualLikGrad=genIndividualLikGrad2D, mu=mu, 
                        v=vvars, upper_tri=upper_tri, rho=rho)-1/(1+rho)+1/(1-rho))
  }
  
  # maximize the log-likelihood
  optim_res = findZeroGradient(llikGradient)
  
  # laplace approximation
  laplace_mean = optim_res[1]
  laplace_var = 1 / c(optim_res[2])
  
  if (laplace_var <= 0) {
    return(runif(1, min=-1, max=1))
  }
  
  while (abs(x_sample) >= 1) {
    # sample once
    x_sample = laplace_mean + sqrt(laplace_var) * rnorm(1)
  }
  
  return(x_sample)
}

# Data-augmented MCMC to directly sample from the actual posterior
auxMatrix <- function(Sigma) {
  # TODO
  
  p = nrow(Sigma)
  aux_mat = diag(diag(Sigma))
  for (j in 1:p) {
    aux_mat[j, -j] = Sigma[j, -j] %*% solve(Sigma[-j, -j])
    aux_mat[j, j] = Sigma[j, j] - Sigma[j, -j] %*% solve(Sigma[-j, -j]) %*% Sigma[-j, j]
  }
  return (aux_mat)
}

# Impute missing latent x
imputeLatentX <- function(imcomplete_x, sample_mu, sample_diag, sample_corr) {
  # TODO
  
  ns = nrow(sample_diag)
  p = ncol(sample_diag)
  pred_id = which(is.na(imcomplete_x))
  preds = NULL
  for (i in 1:ns) {
    Sigma = diag(sqrt(sample_diag[i,]))%*%vec2Corr(sample_corr[i,],p)%*%
      diag(sqrt(sample_diag[i,]))
    R = Sigma[pred_id,-pred_id] %*% ginv(Sigma[-pred_id,-pred_id])
    V = Sigma[pred_id,pred_id] - Sigma[pred_id,-pred_id]%*%
      ginv(Sigma[-pred_id,-pred_id])%*%Sigma[-pred_id,pred_id]
    if (V<0 & V>-1E-8){
      V = 0 # handle numerical error when Sigma is rank deficit
    } else if (V<=-1E-8) {
      stop("conditional variance negative!")
    }
    preds = c(preds, sqrt(V)*rnorm(1)+
                sample_mu[i,pred_id]+sum(R*(imcomplete_x[-pred_id]-sample_mu[i,-pred_id])))
  }
  
  return(list(mean=mean(preds),
              sd=sd(preds),
              interv=quantile(preds,c(.025,.975)),
              preds=preds))
}

evalPredAccuracy <- function(nsim, mu, S, sample_mu, sample_diag, sample_corr) {
  # TODO
  
  p = ncol(S)
  errs = NULL
  coverages = NULL
  for (i in 1:nsim) {
    imcomplete_x = rmvnorm(1,mu,S)
    pred_id = sample(1:p,1)
    targ = imcomplete_x[pred_id]
    imcomplete_x[pred_id] = NA
    pred = imputeLatentX(imcomplete_x, sample_mu, sample_diag, sample_corr)
    errs = c(errs, pred$mean-targ)
    coverages = c(coverages, targ>pred$interv[1]&targ<pred$interv[2])
  }
  
  return(c(sqrt(mean(errs^2)), mean(coverages)))
}
