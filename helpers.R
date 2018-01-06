# Bayesian Mosaic
# Version 1.1
# Last Updated on Dec 12, 2017
suppressMessages(require(rstan))
suppressMessages(require(mvtnorm))
suppressMessages(require(pracma))
suppressMessages(require(MCMCpack))
suppressMessages(require(compiler))
suppressMessages(require(ars))
suppressMessages(require(parallel))
ncores = detectCores() - 1

# Helper Functions----
compressCount <- function(y){
  # compress multivariate count vectors into c(count,value) pairs
  # args:
  #   y: count observations
  
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

Corr2Vec <- function(S) {
  # TODO
  
  p = ncol(S)
  ret = NULL
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      ret = c(ret, S[i,j])
    }
  }
  return(ret)
}

Vec2Corr <- function(rhos, p) {
  # TODO
  
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

# likelihood functions via numerical integration, 1-d or 2-d
genLikPoissonLogNormal <- function(y, mu, v){
  # compute the likelihood function value for a individual observation from a
  # Poisson log-normal distribution via numerical integration
  # args:
  #   y: observation
  #   mu: mean
  #   v: variance
  
  integrand <- function(x){
    return( exp(y * x - lfactorial(y) - (x - mu)^2 / 2 / v - exp(x)) * (2 * pi * v)^-.5 )
  }
  integrate(integrand, lower = min(log(y + 0.1), mu) - 4 * sqrt(v), 
            upper = max(log(y + 0.1), mu) + 4 * sqrt(v), stop.on.error=FALSE)$value
}

genLikPoissonLogGaussian2D <- function(y, mu, Sigma, Sinv){
  # compute the likelihood function value for a individual observation from a multivariate 
  # Poisson log-Gaussian distribution via numerical integration
  # args:
  #   y: observation
  #   mu: mean
  #   Sigma: covariance matrix
  #   Sinv: precision matrix (avoid redundant matrix inversion)
  
  upper_tri = chol(Sigma)
  integrand <- function(z1, z2) {
    # change of variables
    x1 = mu[1]+upper_tri[1,1]*z1+upper_tri[2,1]*z2
    x2 = mu[2]+upper_tri[1,2]*z1+upper_tri[2,2]*z2
    dpois(y[1], exp(x1))*dpois(y[2], exp(x2))*dnorm(z1)*dnorm(z2)
    # note that dmvnorm(cbind(c(x1),c(x2)), mu, S)*det(upper_tri) equals dnorm(z1)*dnorm(z2)
  }
  # finite lower and upper bounds have to be provided for quad2d
  # in our case we set it to be 10 sd away from the center to make sure that
  # the function value is negligible
  quad2d(integrand, -10, 10, -10, 10, 64)
}

genLikBinomialLogNormal <- function(y, N, mu, v){
  # compute the likelihood function value for a individual observation from a
  # Poisson log-normal distribution via numerical integration
  # args:
  #   y: observation
  #   N: number of trials
  #   mu: mean
  #   v: variance
  
  integrand <- function(x){
    return( exp(y * x - (x - mu)^2 / 2 / v) * (2 * pi * v)^-.5 / (1 + exp(x))^N )
  }
  integrate(integrand, lower = min(log((y + 0.1) / N), mu) - 4 * sqrt(v), 
            upper = max(log((y + 0.1) / N), mu) + 4 * sqrt(v), stop.on.error=FALSE)$value
}

genLikBinomialLogGaussian2D <- function(y, Ns, mu, Sigma, Sinv){
  # TODO
  
  upper_tri = chol(Sigma)
  integrand <- function(z1, z2) {
    x1 = mu[1]+upper_tri[1,1]*z1+upper_tri[2,1]*z2
    x2 = mu[2]+upper_tri[1,2]*z1+upper_tri[2,2]*z2
    dbinom(y[1], Ns[1], 1/(1+exp(-x1)))*dbinom(y[2], Ns[2], 1/(1+exp(-x2)))*dnorm(z1)*dnorm(z2)
  }
  
  quad2d(integrand, -10, 10, -10, 10, 64)
}

genLogLikelihood <- function(compressed_y, likelihood, mu, log = TRUE, ...){
  # TODO
  
  ret = 0
  
  if (likelihood == "PoissonLogNormal") {
    for (i in 1:ncol(compressed_y)){
      ret = ret + compressed_y[1, i] * log(genLikPoissonLogNormal(y = compressed_y[-1, i], mu = mu, ...))
    }
  } else if (likelihood == "PoissonLogGaussian2D") {
    for (i in 1:ncol(compressed_y)){
      ret = ret + compressed_y[1, i] * log(genLikPoissonLogGaussian2D(y = compressed_y[-1, i], mu = mu, ...))
    }
  } else if (likelihood == "BinomialLogNormal") {
    for (i in 1:ncol(compressed_y)){
      ret = ret + compressed_y[1, i] * log(genLikBinomialLogNormal(y = compressed_y[-1, i], mu = mu, ...))
    }
  } else if (likelihood == "BinomialLogGaussian2D") {
    for (i in 1:ncol(compressed_y)){
      ret = ret + compressed_y[1, i] * log(genLikBinomialLogGaussian2D(y = compressed_y[-1, i], mu = mu, ...))
    }
  }
  
  if (log == FALSE){
    ret = exp(ret)
  }
  return(ret)
}

genGroupLogLikhood <- function(group_compressed_y, likelihood, group_mus, log = TRUE, ...){
  # TODO
  
  ret = NULL
  for (i in 1:ncol(group_mus)){
    ret = c(ret, genLogLikelihood(compressed_y = group_compressed_y[[i]], 
                                  likelihood = likelihood, mu = group_mus[,i], ...))
  }
  if (log == FALSE){
    ret = exp(ret)
  }
  return(ret)
}

# gradient w.r.t. rho of the likelihood functions via numerical integration, 2-d
genLikPoissonLogGaussian2DGradient <- function(y, mu, Sigma, Sinv){
  # likelihood function value for a individual observation via numerical integration
  # args:
  #   y: observation
  #   mu: mean
  #   Sigma: covariance matrix
  #   Sinv: precision matrix
  
  upper_tri = chol(Sigma)
  v = diag(Sigma)
  rho = Sigma[1, 2] / sqrt(prod(v))
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

genLikBinomialLogGaussian2DGradient <- function(y, Ns, mu, Sigma, Sinv){
  # TODO
  
  upper_tri = chol(Sigma)
  v = diag(Sigma)
  rho = Sigma[1, 2] / sqrt(prod(v))
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

genLogLikelihoodGradient <- function(compressed_y, likelihood, ...){
  # TODO
  
  ret = 0
  
  if (likelihood == "PoissonLogGaussian2D") {
    for (i in 1:ncol(compressed_y)){
      ret = ret + compressed_y[1, i] * genLikPoissonLogGaussian2DGradient(y = compressed_y[-1, i], ...) / 
        genLikPoissonLogGaussian2D(y = compressed_y[-1, i], ...)
    }
  } else if (likelihood == "BinomialLogGaussian2D") {
    for (i in 1:ncol(compressed_y)){
      ret = ret + compressed_y[1, i] * genLikBinomialLogGaussian2DGradient(y = compressed_y[-1, i], ...) / 
        genLikBinomialLogGaussian2D(y = compressed_y[-1, i], ...)
    }
  }
  
  return(ret)
}

# Complie the Functions for Efficiency----
compressCount <- cmpfun(compressCount)
genLogLikelihood <- cmpfun(genLogLikelihood)
genGroupLogLikhood <- cmpfun(genGroupLogLikhood)
genLogLikelihoodGradient <- cmpfun(genLogLikelihoodGradient)

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
  
  x_sample = -Inf
  
  if (model == "mvtPoisson") {
    llikGradient <- function(rho) {
      Sigma = diag(sqrt(vvars)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sqrt(vvars))
      Sinv = solve(Sigma)
      return(-genLogLikelihoodGradient(compressed_y = compressed_y, likelihood = "PoissonLogGaussian2D", 
                                        mu = mu, Sigma = Sigma, Sinv = Sinv) - 1 / (1 + rho) + 1 / (1 - rho))
    }
  } else if (model == "mvtBinomial") {
    llikGradient <- function(rho) {
      Sigma = diag(sqrt(vvars)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sqrt(vvars))
      Sinv = solve(Sigma)
      return(-genLogLikelihoodGradient(compressed_y = compressed_y, likelihood = "BinomialLogGaussian2D", 
                                       mu = mu, Sigma = Sigma, Sinv = Sinv, ...) - 1 / (1 + rho) + 1 / (1 - rho))
    }
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

sampleViaDAMCMC <- function(Y, ns, inits, stop_time=NULL, verbose=FALSE, parallel=FALSE){
  # TODO
  
  start_time = proc.time()
  
  n = nrow(Y)
  p = ncol(Y)
  # initialization
  latent_x = matrix(0, n, p)
  if (is.null(inits)) {
    mu = rep(0, p)
    Sigma = diag(1, p)
  } else {
    mu = inits$mu
    Sigma = inits$Sigma
  }
  Sinv = solve(Sigma)
  
  # outputs
  sample_mu = NULL
  sample_diag_Sigma = NULL
  sample_Correlation = NULL
  
  # helpers
  f <- function(x, y, mu, s) {
    y * x - exp(x) - (x - mu)^2 / 2 / s
  }
  fprima <- function(x, y, mu, s) {
    y - exp(x) - (x - mu) / s
  }
  
  # sampling
  for (iter in 1:ns) {
    # early stop
    if (!is.null(stop_time) & (proc.time()-start_time)[3]>stop_time) {
      return(list(sample_mu=sample_mu, sample_diag_Sigma=sample_diag_Sigma,
                  sample_Correlation=sample_Correlation, ealry_stop=TRUE))
    }
    
    # print intermediate sampling info
    if (verbose && (iter %% max(floor(ns/100), 1) == 0)) {
      cat("iteration: ", iter, "\n")
    }
    
    # sample latent x
    aux_mat = auxMatrix(Sigma)
    sample_latent_x <- function(i) {
      x = latent_x[i,]
      for (j in 1:p) {
        mux = mu[j] + sum(aux_mat[j, -j]*(x[-j]-mu[-j]))
        # subtract the global mean from latent x's
        x[j] = ars(n = 1, f = f, fprima = fprima, x = c(-50, 1, 50), ns = 100,
                   m = 3, emax = 64, lb = FALSE, ub = FALSE,
                   y = Y[i,j], mu = mux, s = aux_mat[j, j])
      }
      return( x )
    }
    if (parallel) {
      latent_x = matrix(unlist(mclapply(1:n, FUN=sample_latent_x, mc.cores=ncores)), n, p, byrow = TRUE)
    } else {
      for (i in 1:n) {
        latent_x[i,] = sample_latent_x(i)
      }
    }
    stopifnot(is.matrix(latent_x) && nrow(latent_x) == n)
    
    # sample mu & Sigma
    mu = apply(latent_x, 2, mean)
    latent_x_demeaned = latent_x-matrix(rep(mu,n),n,p,byrow=TRUE)
    Sigma = riwish(p + 2 + n, diag(1, p) + t(latent_x_demeaned)%*%latent_x_demeaned)
    
    # prepare outputs
    sample_mu = rbind(sample_mu, c(mu))
    sample_diag_Sigma = rbind(sample_diag_Sigma, diag(Sigma))
    Corr = diag(sqrt(1/diag(Sigma)))%*%Sigma%*%diag(sqrt(1/diag(Sigma)))
    sample_Correlation = rbind(sample_Correlation, Corr2Vec(Corr))
  }
  
  return(list(sample_mu=sample_mu, sample_diag_Sigma=sample_diag_Sigma,
              sample_Correlation=sample_Correlation, ealry_stop=FALSE))
}

# MCMC based on 2-d integration
sampleVia2DInt <- function(y, ns, proposal_var, likelihood, stop_time=NULL,
                           ada_prop=0.5, verbose=FALSE, ...) {
  #   ada_prop: # the proportion of burn-ins for which adapting the proposal variance is based on
  # TODO
  
  acr_target = 0.25 # target acceptance rate
  acr_tol = 0.05 # tolerance for acceptance rate
  dcf_pv = 0.9 # discount factor for adapting the proposal variance
  dcf_multiplier = 0.9 # discount factor for adapting the multiplier
  dcf_acr = 0.99 # discount factor for averaging the historical acceptance rate
  
  start_time = proc.time()
  
  # compress the counts
  compressed_y = compressCount(y)
  
  # initialization
  rho = 0
  mu1 = 0
  s1 = 1
  mu2 = 0
  s2 = 1
  prv_llik = -Inf
  historical_acr_mu = 1
  historical_acr_s = 1
  historical_acr_rho = 1
  
  # storing the results
  sample_rho = NULL
  sample_mu1 = NULL
  sample_s1 = NULL
  sample_mu2 = NULL
  sample_s2 = NULL
  
  # preparation for adapting
  proposal_var_rho = proposal_var
  window_memory_rho = rep(0, 100)
  window_memory_rho[100] = rho
  proposal_var_mu1 = proposal_var
  window_memory_mu1 = rep(0, 100)
  window_memory_mu1[100] = mu1
  proposal_var_s1 = proposal_var
  window_memory_s1 = rep(0, 100)
  window_memory_s1[100] = s1
  proposal_var_mu2 = proposal_var
  window_memory_mu2 = rep(0, 100)
  window_memory_mu2[100] = mu2
  proposal_var_s2 = proposal_var
  window_memory_s2 = rep(0, 100)
  window_memory_s2[100] = s2
  adaptive_multiplier_mu = 1
  adaptive_multiplier_s = 1
  adaptive_multiplier_rho = 1
  
  # sampling
  for (iter in 1:ns) {
    # early stop
    if (!is.null(stop_time) & (proc.time()-start_time)[3]>stop_time) {
      return(list(sample_rho=sample_rho, sample_mu1=sample_mu1,
                  sample_s1=sample_s1, sample_mu2=sample_mu2,
                  sample_s2=sample_s2, ealry_stop=TRUE))
    }
    
    # print intermediate sampling info
    if (verbose && (iter %% max(floor(ns/100), 1) == 0)) {
      cat("iteration: ", iter, "\n")
    }
    
    # update mu's
    new_mu1 = mu1 + sqrt(proposal_var_mu1)*rnorm(1)
    new_mu2 = mu2 + sqrt(proposal_var_mu2)*rnorm(1)
    Sigma = diag(c(sqrt(s1), sqrt(s2)))%*%matrix(c(1,rho,rho,1),2,2)%*%
      diag(c(sqrt(s1), sqrt(s2)))
    Sinv = solve(Sigma)
    new_llik = genLogLikelihood(compressed_y=compressed_y, likelihood=likelihood, 
                                mu=c(new_mu1, new_mu2), log=TRUE, Sigma=Sigma,
                                Sinv=Sinv, ...)
    acr = min(exp(new_llik-prv_llik), 1) # symmetric proposal
    historical_acr_mu = dcf_acr*historical_acr_mu + (1-dcf_acr)*acr
    if (runif(1) <= acr) {
      mu1 = new_mu1
      mu2 = new_mu2
      prv_llik = new_llik
    }
    
    # update s
    new_s1 = s1 + sqrt(proposal_var_s1)*rnorm(1)
    new_s2 = s2 + sqrt(proposal_var_s2)*rnorm(1)
    if (new_s1>0 & new_s2>0) {
      new_Sigma = diag(c(sqrt(new_s1), sqrt(new_s2)))%*%matrix(c(1,rho,rho,1),2,2)%*%
        diag(c(sqrt(new_s1), sqrt(new_s2)))
      new_Sinv = solve(new_Sigma)
      new_llik = genLogLikelihood(compressed_y=compressed_y, likelihood=likelihood, 
                                  mu=c(new_mu1, new_mu2), log=TRUE, Sigma=new_Sigma,
                                  Sinv=new_Sinv, ...)
      new_lpost = new_llik + dgamma(1/new_s1, shape=2, rate=2, log=TRUE) + 
        dgamma(1/new_s2, shape=2, rate=2, log=TRUE)
      prv_lpost = prv_llik + dgamma(1/s1, shape=2, rate=2, log=TRUE) + 
        dgamma(1/s2, shape=2, rate=2, log=TRUE)
      acr = min(exp(new_lpost-prv_lpost), 1) # symmetric proposal
      historical_acr_s = dcf_acr*historical_acr_s + (1-dcf_acr)*acr
      if (runif(1) <= acr) {
        s1 = new_s1
        s2 = new_s2
        prv_llik = new_llik
      }
    } else {
      historical_acr_s = dcf_acr*historical_acr_s
    }
    
    # update rho
    new_rho = rho + sqrt(proposal_var_rho)*rnorm(1)
    if (abs(new_rho)<1) {
      new_Sigma = diag(c(sqrt(s1), sqrt(s2)))%*%matrix(c(1,new_rho,new_rho,1),2,2)%*%
        diag(c(sqrt(s1), sqrt(s2)))
      new_Sinv = solve(new_Sigma)
      new_llik = genLogLikelihood(compressed_y=compressed_y, likelihood=likelihood, 
                                  mu=c(new_mu1, new_mu2), log=TRUE, Sigma=new_Sigma,
                                  Sinv=new_Sinv, ...)
      acr = min(exp(new_llik-prv_llik), 1) # symmetric proposal
      historical_acr_rho = dcf_acr*historical_acr_rho + (1-dcf_acr)*acr
      if (runif(1) <= acr) {
        rho = new_rho
        prv_llik = new_llik
      }
    } else {
      historical_acr_rho = dcf_acr*historical_acr_rho
    }
    
    # store samples
    if (iter > nb & (iter - nb) %% njump == 0) {
      sample_rho = c(sample_rho, rho)
      sample_mu1 = c(sample_mu1, mu1)
      sample_s1 = c(sample_s1, s1)
      sample_mu2 = c(sample_mu2, mu2)
      sample_s2 = c(sample_s2, s2)
    }
    
    # adapting proposal sd
    if (iter<(nb*(1+ada_prop)/2) & iter>=(nb*(1-ada_prop)/2)) {
      # mu
      window_memory_mu1[1:99] = window_memory_mu1[2:100]
      window_memory_mu1[100] = mu1
      window_memory_mu2[1:99] = window_memory_mu2[2:100]
      window_memory_mu2[100] = mu2
      if (historical_acr_mu > (acr_target + acr_tol)) {
        # to shrink acr, make the multipler larger
        adaptive_multiplier_mu = dcf_multiplier*adaptive_multiplier_mu + 
          (1-dcf_multiplier) * (adaptive_multiplier_mu+1)
      }
      if (historical_acr_mu < (acr_target - acr_tol)) {
        # to increase acr, make the multipler smaller
        adaptive_multiplier_mu = dcf_multiplier*adaptive_multiplier_mu + 
          (1-dcf_multiplier) * max(adaptive_multiplier_mu-1, 1)
      }
      proposal_var_mu1 = dcf_pv*proposal_var_mu1 + 
        (1-dcf_pv)*adaptive_multiplier_mu*var(window_memory_mu1, na.rm = TRUE)
      proposal_var_mu2 = dcf_pv*proposal_var_mu2 + 
        (1-dcf_pv)*adaptive_multiplier_mu*var(window_memory_mu2, na.rm = TRUE)
      
      # s
      window_memory_s1[1:99] = window_memory_s1[2:100]
      window_memory_s1[100] = s1
      window_memory_s2[1:99] = window_memory_s2[2:100]
      window_memory_s2[100] = s2
      if (historical_acr_s > (acr_target + acr_tol)) {
        # to shrink acr, make the multipler larger
        adaptive_multiplier_s = dcf_multiplier*adaptive_multiplier_s + 
          (1-dcf_multiplier) * (adaptive_multiplier_s+1)
      }
      if (historical_acr_s < (acr_target - acr_tol)) {
        # to increase acr, make the multipler smaller
        adaptive_multiplier_s = dcf_multiplier*adaptive_multiplier_s + 
          (1-dcf_multiplier) * max(adaptive_multiplier_s-1, 1)
      }
      proposal_var_s1 = dcf_pv*proposal_var_s1 + 
        (1-dcf_pv)*adaptive_multiplier_s*var(window_memory_s1, na.rm = TRUE)
      proposal_var_s2 = dcf_pv*proposal_var_s2 + 
        (1-dcf_pv)*adaptive_multiplier_s*var(window_memory_s2, na.rm = TRUE)
      
      # rho
      window_memory_rho[1:99] = window_memory_rho[2:100]
      window_memory_rho[100] = rho
      if (historical_acr_rho > (acr_target + acr_tol)) {
        # to shrink acr, make the multipler larger
        adaptive_multiplier_rho = dcf_multiplier*adaptive_multiplier_rho + 
          (1-dcf_multiplier) * (adaptive_multiplier_rho+1)
      }
      if (historical_acr_rho < (acr_target - acr_tol)) {
        # to increase acr, make the multipler smaller
        adaptive_multiplier_rho = dcf_multiplier*adaptive_multiplier_rho + 
          (1-dcf_multiplier) * max(adaptive_multiplier_rho-1, 1)
      }
      proposal_var_rho = dcf_pv*proposal_var_rho + 
        (1-dcf_pv)*adaptive_multiplier_rho*var(window_memory_rho, na.rm = TRUE)
    }
  }
  
  return(list(sample_rho=sample_rho, sample_mu1=sample_mu1,
              sample_s1=sample_s1, sample_mu2=sample_mu2,
              sample_s2=sample_s2, ealry_stop=FALSE))
}

