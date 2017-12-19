# Bayesian Mosaic
# Version 1.1
# Last Updated on Dec 12, 2017
suppressMessages(require(rstan))
suppressMessages(require(mvtnorm))
suppressMessages(require(pracma))
suppressMessages(require(MCMCpack))
suppressMessages(require(compiler))
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

genLikPoissonLogGuassian2D <- function(y, mu, Sigma, Sinv){
  # compute the likelihood function value for a individual observation from a multivariate 
  # Poisson log-Gaussian distribution via numerical integration
  # args:
  #   y: observation
  #   mu: mean
  #   Sigma: covariance matrix
  #   Sinv: precision matrix (avoid redundant matrix inversion)
  
  integrand <- function(x1, x2) {
    return(exp(y[1] * x1 + y[2] * x2 - lfactorial(y[1]) - lfactorial(y[2]) -
                 ((x1 - mu[1])^2 * Sinv[1, 1] + 2 * (x1 - mu[1]) * (x2 - mu[2]) * Sinv[1, 2] + (x2 - mu[2])^2 * Sinv[2, 2]) / 2 -
                 exp(x1) - exp(x2)) * (2 * pi)^(-1) * det(Sigma)^-.5)
  }
  
  return(quad2d(integrand, min(log(y[1] + 0.1), mu[1]) - 4 * sqrt(Sigma[1, 1]), 
                max(log(y[1] + 0.1), mu[1]) + 4 * sqrt(Sigma[1, 1]), 
                min(log(y[2] + 0.1), mu[2]) - 4 * sqrt(Sigma[2, 2]), 
                max(log(y[2] + 0.1) , mu[2]) + 4 * sqrt(Sigma[2, 2]), 64))
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

genLikBinomialLogGuassian2D <- function(y, Ns, mu, Sigma, Sinv){
  # TODO
  
  integrand <- function(x1, x2) {
    return(exp(y[1] * x1 + y[2] * x2 - ((x1 - mu[1])^2 * Sinv[1, 1] + 2 * (x1 - mu[1]) * (x2 - mu[2]) * Sinv[1, 2] + (x2 - mu[2])^2 * Sinv[2, 2]) / 2)
           * (2 * pi)^(-1) * det(Sigma)^-.5 / ((1 + exp(x1))^Ns[1]) / ((1 + exp(x2))^Ns[2]))
  }
  
  return(quad2d(integrand, min(log((y[1] + 0.1) / Ns[1]), mu[1]) - 4 * sqrt(Sigma[1, 1]),
                max(log((y[1] + 0.1) / Ns[1]), mu[1]) + 4 * sqrt(Sigma[1, 1]), 
                min(log((y[2] + 0.1) / Ns[2]), mu[2]) - 4 * sqrt(Sigma[2, 2]),
                max(log((y[2] + 0.1) / Ns[2]) , mu[2]) + 4 * sqrt(Sigma[2, 2]), 64))
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
      ret = ret + compressed_y[1, i] * log(genLikBinomialLogGuassian2D(y = compressed_y[-1, i], mu = mu, ...))
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
genLikPoissonLogGuassian2DGradient <- function(y, mu, Sigma, Sinv){
  # likelihood function value for a individual observation via numerical integration
  # args:
  #   y: observation
  #   mu: mean
  #   Sigma: covariance matrix
  #   Sinv: precision matrix
  
  v = diag(Sigma)
  rho = Sigma[1, 2] / sqrt(prod(v))
  integrand <- function(x1, x2) {
    x1_tilde = (x1 - mu[1]) / sqrt(v[1])
    x2_tilde = (x2 - mu[2]) / sqrt(v[2])
    gradient_term = (rho * (1 - rho^2) + x1_tilde * x2_tilde * (1 + rho^2) - rho * (x1_tilde^2 + x2_tilde^2)) / (1 - rho^2)^2
    return(exp(y[1] * x1 + y[2] * x2 - lfactorial(y[1]) - lfactorial(y[2]) -
                 ((x1 - mu[1])^2 * Sinv[1, 1] + 2 * (x1 - mu[1]) * (x2 - mu[2]) * Sinv[1, 2] + (x2 - mu[2])^2 * Sinv[2, 2]) / 2 -
                 exp(x1) - exp(x2)) * (2 * pi)^(-1) * det(Sigma)^-.5 * gradient_term)
  }
  return(quad2d(integrand, min(log(y[1] + 0.1), mu[1])-4 * sqrt(v[1]), 
                max(log(y[1] + 0.1), mu[1]) + 4 * sqrt(v[1]), min(log(y[2] + 0.1), mu[2]) - 4 * sqrt(v[2]), 
                max(log(y[2] + 0.1) , mu[2]) + 4 * sqrt(v[2]), 64))
}

genLikBinomialLogGuassian2DGradient <- function(y, Ns, mu, Sigma, Sinv){
  # TODO
  
  v = diag(Sigma)
  rho = Sigma[1, 2] / sqrt(prod(v))
  integrand <- function(x1, x2) {
    x1_tilde = (x1 - mu[1]) / sqrt(v[1])
    x2_tilde = (x2 - mu[2]) / sqrt(v[2])
    gradient_term = (rho * (1 - rho^2) + x1_tilde * x2_tilde * (1 + rho^2) - rho * (x1_tilde^2 + x2_tilde^2)) / (1 - rho^2)^2
    return(exp(y[1] * x1 + y[2] * x2 - ((x1 - mu[1])^2 * Sinv[1, 1] + 2 * (x1 - mu[1]) * (x2 - mu[2]) * Sinv[1, 2] + (x2 - mu[2])^2 * Sinv[2, 2]) / 2)
           * (2 * pi)^(-1) * det(Sigma)^-.5 / ((1 + exp(x1))^Ns[1]) / ((1 + exp(x2))^Ns[2]) * gradient_term)
  }
  return(quad2d(integrand, min(log((y[1] + 0.1) / Ns[1]), mu[1]) - 4 * sqrt(v[1]), 
                max(log((y[1] + 0.1) / Ns[1]), mu[1]) + 4 * sqrt(v[1]), min(log((y[2] + 0.1) / Ns[2]), mu[2]) - 4 * sqrt(v[2]), 
                max(log((y[2] + 0.1) / Ns[2]) , mu[2]) + 4 * sqrt(v[2]), 64))
}

genLogLikelihoodGradient <- function(compressed_y, likelihood, ...){
  # TODO
  
  ret = 0
  
  if (likelihood == "PoissonLogGaussian2D") {
    for (i in 1:ncol(compressed_y)){
      ret = ret + compressed_y[1, i] * genLikPoissonLogGuassian2DGradient(y = compressed_y[-1, i], ...) / 
        genLikPoissonLogGuassian2D(y = compressed_y[-1, i], ...)
    }
  } else if (likelihood == "BinomialLogGaussian2D") {
    for (i in 1:ncol(compressed_y)){
      ret = ret + compressed_y[1, i] * genLikBinomialLogGuassian2DGradient(y = compressed_y[-1, i], ...) / 
        genLikBinomialLogGuassian2D(y = compressed_y[-1, i], ...)
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
    return(NA)
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

sampleViaDAMCMC <- function(Y, ns, verbose = FALSE){
  # TODO
  
  n = nrow(Y)
  p = ncol(Y)
  # initialization
  latent_x = matrix(0, n, p)
  mu = rep(0, p)
  Sigma = diag(1, p)
  Sinv = solve(Sigma)
  
  # outputs
  sample_mu = NULL
  sample_Sigma = NULL
  
  # helpers
  f <- function(x, y, mu, s) {
    y * x - exp(x) - (x - mu)^2 / 2 / s
  }
  fprima <- function(x, y, mu, s) {
    y - exp(x) - (x - mu) / s
  }
  
  aux_mat = auxMatrix(Sigma)
  
  # sampling
  for (iter in 1:ns) {
    # print intermediate sampling info
    if (verbose && (iter %% max(floor(ns / 100), 1) == 0)) {
      cat("iteration: ", iter, "\n")
    }
    
    # sample latent x
    latent_x = matrix(unlist(mclapply(1:n, FUN = function(i){
      x = latent_x[i,]
      for (j in 1:p) {
        mux = mu[j] + sum(aux_mat[j, -j] * (x[-j] - mu[-j]))
        # subtract the global mean from latent x's
        x[j] = ars(n = 1, f = f, fprima = fprima, x = c(-50, 1, 50), ns = 100,
                   m = 3, emax = 64, lb = FALSE, ub = FALSE,
                   y = Y[i,j], mu = mux, s = aux_mat[j, j])
      }
      return( x )
    },mc.cores = ncores)), n, p, byrow = TRUE)
    stopifnot(is.matrix(latent_x) && nrow(latent_x) == n)
    
    # sample mu & Sigma
    mu = apply(latent_x, 2, mean)
    Sigma = riwish(p + 2 + n, diag(1, p) + t(mu) %*% mu)
    
    # prepare outputs
    sample_mu = cbind(sample_mu, c(mu))
    sample_Sigma = cbind(sample_Sigma, c(Sigma))
  }
  
  return(list(sample_mu = sample_mu, sample_Sigma = sample_Sigma))
}

# Bayesian Mosaic
sampleKnot <- function(y, nb, ns, njump, proposal_var, model, verbose = FALSE, ...) {
  # TODO
  
  acr_target = 0.25 # target acceptance rate
  acr_tol = 0.05 # tolerance for acceptance rate
  dcf_pv = 0.9 # discount factor for adapting the proposal variance
  dcf_multiplier = 0.9 # discount factor for adapting the multiplier
  dcf_acr = 0.99 # discount factor for averaging the historical acceptance rate
  ada_prop = 0.8 # the proportion of burn-ins for which adapting the proposal variance is based on
  
  # compress the counts
  compressed_y = compressCount(y)
  
  # initialization
  s = 1
  mu = 0
  prv_llik = -Inf
  historical_acr = 1
  
  # storing the results
  sample_s = NULL
  sample_mu = NULL
  
  # preparation for adapting
  proposal_var_mu = proposal_var
  window_memory_mu = rep(0, 100)
  window_memory_mu[100] = mu
  proposal_var_s = proposal_var
  window_memory_s = rep(0, 100)
  window_memory_s[100] = s
  adaptive_multiplier = 1
  
  # sampling
  for (iter in 1:(nb + ns * njump)) {
    # print intermediate sampling info
    if (verbose && (iter %% max(floor((nb + ns * njump) / 100), 1)) == 0) {
      cat("iteration: ", iter, "\n")
    }
    
    # update mu & s
    new_mu = mu + sqrt(proposal_var_mu) * rnorm(1)
    new_s = s + sqrt(proposal_var_s) * rnorm(1)
    if (new_s > 0) {
      if (model == "mvtPoisson") {
        new_llik = genLogLikelihood(compressed_y = compressed_y, likelihood = "PoissonLogNormal", 
                                    mu = new_mu, log = TRUE, v = new_s)
      } else if (model == "mvtBinomial") {
        new_llik = genLogLikelihood(compressed_y = compressed_y, likelihood = "BinomialLogNormal", 
                                    mu = new_mu, log = TRUE, v = new_s, ...)
      }
      new_lpost = new_llik + dgamma(1 / new_s, shape = 2, rate = 2, log = TRUE)
      prv_lpost = prv_llik + dgamma(1 / s, shape = 2, rate = 2, log = TRUE)
      acr = min(exp(new_lpost - prv_lpost), 1) # symmetric proposal
      historical_acr = dcf_acr * historical_acr + (1 - dcf_acr) * acr
      if (runif(1) <= acr) {
        mu = new_mu
        s = new_s
        prv_llik = new_llik
      }
    } else {
      historical_acr = dcf_acr * historical_acr
    }
    
    # store samples
    if (iter > nb & (iter - nb) %% njump == 0) {
      sample_s = c(sample_s, s)
      sample_mu = c(sample_mu, mu)
    }
    
    # adapting proposal sd
    if (iter <= (nb * ada_prop)) {
      # mu and s
      window_memory_mu[1:99] = window_memory_mu[2:100]
      window_memory_mu[100] = mu
      window_memory_s[1:99] = window_memory_s[2:100]
      window_memory_s[100] = s
      if (historical_acr > (acr_target + acr_tol)) {
        # to shrink acr, make the multipler larger
        adaptive_multiplier = dcf_multiplier * adaptive_multiplier + 
          (1 - dcf_multiplier) * (adaptive_multiplier + 1)
      }
      if (historical_acr < (acr_target - acr_tol)) {
        # to increase acr, make the multipler smaller
        adaptive_multiplier = dcf_multiplier * adaptive_multiplier + 
          (1 - dcf_multiplier) * max(adaptive_multiplier - 1, 1)
      }
      proposal_var_mu = dcf_pv * proposal_var_mu + 
        (1 - dcf_pv) * adaptive_multiplier * var(window_memory_mu, na.rm = TRUE)
      proposal_var_s = dcf_pv * proposal_var_s + 
        (1 - dcf_pv) * adaptive_multiplier * var(window_memory_s, na.rm = TRUE)
    }
  }
  
  return(list(sample_s = sample_s,
              sample_mu = sample_mu,
              historical_acr = historical_acr))
}

sampleTile <- function(y, sample_mu, sample_s, model, verbose = FALSE, ...) {
  # TODO
  
  ns = ncol(sample_mu)
  # compress the counts
  compressed_y = compressCount(y)
  
  # storing the results
  sample_rho = NULL
  
  # sampling
  for (iter in 1:ns) {
    # print intermediate sampling info
    if (verbose) {
      cat("iteration: ", iter, "\n")
    }
    
    sample_rho = c(sample_rho, sampleLaplaceApprox(compressed_y, model, sample_mu[,iter], sample_s[,iter], ...))
  }
  
  return(sample_rho)
}

bayesianMosaic <- function(Y, nb, ns, njump, proposal_var, model, verbose = FALSE, Ns = NULL) {
  # TODO
  
  if (verbose) {
    cat("Bayesian Mosaic\n")
    cat("learning posterior knots...\n")
  }
  
  p = ncol(Y)
  knots = mclapply(1:p, FUN = function(k){
    return( sampleKnot(y = Y[,k], nb = nb, ns = ns, njump = njump, 
                       proposal_var = proposal_var, model = model, verbose = verbose, N = Ns[k]) )
  }, mc.cores = ncores)
  
  if (verbose) {
    cat("learning posterior tiles...\n")
  }
  
  pairs = NULL
  for (x in 1:(p - 1)) {
    for (y in (x + 1):p) {
      pairs = rbind(pairs, c(x, y))
    }
  }
  n_pair = nrow(pairs)
  
  tiles = mclapply(1:n_pair, FUN = function(k){
    s = pairs[k, 1]
    t = pairs[k, 2]
    y = cbind(Y[,s], Y[,t])
    sample_mu = rbind(knots[[k]]$sample_mu, knots[[t]]$sample_mu)
    sample_s = rbind(knots[[k]]$sample_s, knots[[t]]$sample_s)
    return( sampleTile(y = y, sample_mu = sample_mu, sample_s = sample_s, 
                       model = model, verbose = verbose, Ns = Ns[c(s, t)]) )
  },mc.cores = ncores)
  
  # write outputs
  if (verbose) {
    cat("writing outputs...\n")
  }
  
  sample_mu = matrix(0, p, ns)
  sample_diag_Sigma = matrix(0, p, ns)
  sample_Correlation = array(1, c(p, p, ns))
  
  for (j in 1:p) {
    sample_mu[j,] = knots[[j]]$sample_mu
    sample_diag_Sigma[j,] = knots[[j]]$sample_s
  }
  
  for (k in 1:n_pair) {
    s = pairs[k, 1]
    t = pairs[k, 2]
    sample_Correlation[s, t,] = tiles[[k]]
    sample_Correlation[t, s,] = sample_Correlation[s, t,]
  }
  
  return(list(sample_mu = sample_mu,
              sample_diag_Sigma = sample_diag_Sigma,
              sample_Correlation = sample_Correlation))
}
