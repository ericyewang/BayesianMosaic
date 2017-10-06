# Bayesian Mosaic Simulation Study
# hierarchical Poisson log Gaussian

rm(list=ls())

# dependencies
suppressMessages(require(MCMCpack))
suppressMessages(require(mvtnorm))
suppressMessages(require(pracma))
suppressMessages(require(LaplacesDemon))
suppressMessages(require(Rmisc))
suppressMessages(require(ggplot2))

# global helper functions
aggregateData <- function(y){
  # aggregate data
  # since the count data only takes a small number of values
  # args:
  #   y: count observations
  
  y = as.matrix(y)
  
  pool = unique(y, MARGIN = 1)
  ret = matrix(0, ncol(pool) + 1, nrow(pool))
  for (i in 1:nrow(pool)){
    ret[, i] = c(pool[i,], sum(apply(y, 1, function(x) {
      all(x == pool[i,])
    })))
  }
  
  return(ret)
}

# helper functions for sampling tile posteriors
getLik <- function(y, mu, v){
  # likelihood function value for a individual observation via numerical integration
  # transform the variable so that the range is 0 to 1, have tried -Inf to Inf 
  # and there is some bizzare numerical issue there.
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

getGroupLik <- function(group_summary, mu, v, logarithm = TRUE){
  # get (the logarithm of) of the likelihood of a group
  # args:
  #   group_summary: summary of the observations in the group
  #   mu: mean
  #   v: variance
  #   logarithm: whether or not to return the log likelihood
  
  ret = 0
  for (i in 1:ncol(group_summary)){
    ret = ret + group_summary[2, i] * log(getLik(group_summary[1, i], mu, v))
  }
  if (logarithm == FALSE){
    ret = exp(ret)
  }
  return(ret)
}

getGroupLiks <- function(group_summaries, group_mus, v, logarithm = TRUE){
  # get (the logarithm of) of the likelihoods of all groups
  # args:
  #   group_summaries: a list of group summaries
  #   group_mus: a vector of group means
  #   v: variance
  #   logarithm: whether or not to return the log likelihood
  
  ret = NULL
  for (i in 1:length(group_mus)){
    ret = c(ret, getGroupLik(group_summaries[[i]], group_mus[i], v))
  }
  if (logarithm == FALSE){
    ret = exp(ret)
  }
  return(ret)
}

# sample from the knot posteriors
sampleKnot <- function(y, shoot, nb, ns, njump, lam, a, b, 
                       proposal_var, verbose = FALSE) {
  # sample from the knot posterior
  # Args:
  #   y: observation
  #   shoot: group index of each observation
  #   nb: number of burn-ins
  #   ns: number of samples to collect
  #   njump: thinning parameter
  #   lam: parameter in NIW that controls the prior concentration of mu
  #   a: shape parameter of inverse gamma
  #   b: scale parameter of inverse gamma
  #   proposal_var: initial variance of the proposal normal distribution
  #   verbose: whether or not to print intermediate sampling info
  
  # summary statistics
  K = max(shoot)
  group_summaries = list()
  for( k in 1:K) {
    group_summaries[[k]] = aggregateData(y[shoot == k])
  }
  
  # initialization
  group_mus = rep(0, K)
  s1 = 1
  mu0 = 0
  s2 = 1
  prv_lliks = rep(-Inf, K)
  ars = rep(0, K + 1)
    
  # storing the results
  mgroup_mus = NULL
  vs1 = NULL
  vmu0 = NULL
  vs2 = NULL
  mar = NULL
  
  # preparation for adapting
  proposal_vars_group_mus = rep(proposal_var, K)
  window_memory_group_mus = matrix(NA, 100, K)
  window_memory_group_mus[100,] = group_mus
  proposal_var_s1 = proposal_var
  window_memory_s1 = c(rep(NA, 99), s1)
  
  # sampling
  for (iter in 1:(nb + ns * njump)) {
    # print intermediate sampling info
    if (verbose && (iter %% floor((nb + ns * njump) / 10)) == 0) {
      cat("iteration: ", iter, "\n")
    }
    
    # update group_mus
    new_group_mus = group_mus + sqrt(proposal_vars_group_mus) * rnorm(K)
    new_lliks = getGroupLiks(group_summaries, new_group_mus, s1, logarithm = TRUE)
    new_lposts = new_lliks + dnorm(new_group_mus, mu0, sqrt(s2), log = TRUE)
    prv_lposts = prv_lliks + dnorm(group_mus, mu0, sqrt(s2), log = TRUE)
    for (i in 1:K) {
      ars[i] = min(exp(new_lposts[i] - prv_lposts[i]), 1) # symmetric proposal
      if (runif(1) <= ars[i]) {
        group_mus[i] = new_group_mus[i]
        prv_lliks[i] = new_lliks[i]
      }
    }
    
    # update s1
    new_s1 = s1 + sqrt(proposal_var_s1) * rnorm(1)
    if (new_s1 > 0) {
      new_lliks = getGroupLiks(group_summaries, group_mus, new_s1, logarithm = TRUE)
      new_lpost = sum(new_lliks) + dgamma(1 / new_s1, shape = a, rate = b, log = TRUE)
      prv_lpost = sum(prv_lliks) + dgamma(1 / s1, shape = a, rate = b, log = TRUE)
      ars[K + 1] = min(exp(new_lpost - prv_lpost), 1)
      if (runif(1) <= ars[K + 1]) {
        s1 = new_s1
        prv_lliks = new_lliks
      }
    } else {
      ars[K + 1] = 0
    }
    
    # update mu0 and s2
    s2 = 1 / rgamma(1, shape = a + K / 2, rate = b + ((K - 1) * var(group_mus) + lam * K * mean(group_mus)^2 / (lam + K)) / 2)
    mu0 = sum(group_mus) / (lam + K) + sqrt(s2 / (lam + K)) * rnorm(1)
    
    # store samples
    if (iter > nb & (iter - nb) %% njump == 0) {
      mgroup_mus = cbind(mgroup_mus, group_mus)
      vs1 = c(vs1, s1)
      vmu0 = c(vmu0, mu0)
      vs2 = c(vs2, s2)
      mar = cbind(mar, ars)
    }
    
    # adapting proposal sd
    if (iter <= nb) {
      # only adapt the proposal variance during burn-in to ensure ergodicity
      window_memory_group_mus[1:99,] = window_memory_group_mus[2:100,]
      window_memory_group_mus[100,] = group_mus
      window_memory_s1[1:99] = window_memory_s1[2:100]
      window_memory_s1[100] = s1
      historical_var_s1 = var(window_memory_s1, na.rm = TRUE)
      historical_var_group_mus = apply(window_memory_group_mus, 2, function(x) {
        var(x, na.rm = TRUE)
      })
      proposal_vars_group_mus[historical_var_group_mus > 0] = historical_var_group_mus[historical_var_group_mus > 0]
      if (historical_var_s1 > 0) {
        proposal_var_s1 = historical_var_s1
      }
    }
  }
  
  return(list(mgroup_mus = mgroup_mus, vs1 = vs1, vmu0 = vmu0, vs2 = vs2, mar = mar))
}

# helper functions for sampling tile posteriors
getLik2D <- function(y, mu, Sigma, Sinv){
  # likelihood function value for a individual observation via numerical integration
  # args:
  # args:
  #   y: observation
  #   mu: mean
  #   Sigma: covariance matrix
  #   Sinv: precision matrix
  
  v = diag(Sigma)
  integrand <- function(x1, x2) {
    return(exp(y[1] * x1 + y[2] * x2 - lfactorial(y[1]) - lfactorial(y[2]) -
                 ((x1 - mu[1])^2 * Sinv[1, 1] + 2 * (x1 - mu[1]) * (x2 - mu[2]) * Sinv[1, 2] + (x2 - mu[2])^2 * Sinv[2, 2]) / 2 -
                 exp(x1) - exp(x2)) * (2 * pi)^(-1) * det(Sigma)^-.5)
  }
  return(quad2d(integrand, min(log(y[1] + 0.1), mu[1])-4 * sqrt(v[1]), 
                max(log(y[1] + 0.1), mu[1]) + 4 * sqrt(v[1]), min(log(y[2] + 0.1), mu[2]) - 4 * sqrt(v[2]), 
                max(log(y[2] + 0.1) , mu[2]) + 4 * sqrt(v[2]), 64))
}

getGroupLik2D <- function(group_summary, mu, v, rho, logarithm = TRUE){
  # get (the logarithm of) of the likelihood of a group
  # args:
  #   group_summary: summary of the observations in the group
  #   mu: mean vector
  #   v: marginal variance vector
  #   rho: correlation parameter
  #   logarithm: whether or not to return the log likelihood
  
  Sigma = diag(sqrt(v))%*%matrix(c(1,rho,rho,1),2,2)%*%diag(sqrt(v))
  Sinv = solve(Sigma)
  ret = 0
  for (i in 1:ncol(group_summary)){
    ret = ret + group_summary[3, i] * log(getLik2D(group_summary[1:2, i], mu, Sigma, Sinv))
  }
  if (logarithm == FALSE){
    ret = exp(ret)
  }
  return(ret)
}

getGroupLiks2D <- function(group_summaries, group_mus, v, rho, logarithm = TRUE){
  # get (the logarithm of) of the likelihoods of all groups
  # args:
  #   group_summaries: a list of group summaries
  #   group_mus: a matrix of group mean vectors
  #   v: variance
  #   rho: correlation parameter
  #   logarithm: whether or not to return the log likelihood
  
  ret = NULL
  for (i in 1:ncol(group_mus)){
    ret = c(ret, getGroupLik2D(group_summaries[[i]], group_mus[, i], v, rho))
  }
  if (logarithm == FALSE){
    ret = exp(ret)
  }
  return(ret)
}

# sample from the tile posteriors given knots fixed at their posterior means
sampleTile <- function(ys, yt, shoot, samples_knot_s, samples_knot_t, nb, ns, njump,
                       proposal_var, verbose = FALSE){
  # sample from the tile posterior
  # Args:
  #   ys: dimension s of y
  #   yt: dimension t of y
  #   shoot: group index of each observation
  #   samples_knot_s: posterior samples from the knot posterior of dimension s
  #   samples_knot_t: posterior samples from the knot posterior of dimension t
  #   nb: number of burn-ins
  #   ns: number of samples to collect
  #   njump: thinning parameter
  #   nu: d.f. of inverse-wishart
  #   Psi: scale matrix of inverse-wishart
  #   proposal_var: initial variance of the proposal normal distribution
  #   verbose: whether or not to print intermediate sampling info
  
  small_var = 0.001
  
  # extract knot info
  vs1 = c(mean(samples_knot_s$vs1), mean(samples_knot_t$vs1))
  vs2 = c(mean(samples_knot_s$vs2), mean(samples_knot_t$vs2))
  vmu0 = c(mean(samples_knot_s$vmu0), mean(samples_knot_t$vmu0))
  proposal_vars_group_mus = array(0, c(2, 2, K))
  group_mus = matrix(0, 2, K)
  for (k in 1:K) {
    tmp_mus = cbind(samples_knot_s$mgroup_mus[k,], samples_knot_t$mgroup_mus[k,])
    # use the intermediate samples from knot posteriors to calculate the proposal covariance matrix
    proposal_vars_group_mus[,, k] = cov(tmp_mus)
    # and the mean as the initial values of these augmented variables
    group_mus[, k] = apply(tmp_mus, 2, mean)
  }
  
  # summary statistics
  K = max(shoot)
  group_summaries = list()
  for( k in 1:K) {
    group_summaries[[k]] = aggregateData(y[shoot == k, ])
  }
  
  # initialization
  rho1 = 0
  rho2 = 0
  prv_lliks = rep(-Inf, K)
  ars = rep(0, K + 2)
  new_group_mus = matrix(0, 2, K)
  
  # storing the results
  mgroup_mus = NULL
  vrho1 = NULL
  vrho2 = NULL
  mar = NULL
  
  # preparation for adapting
  window_memory_group_mus = array(NA, c(100, 2, K))
  window_memory_group_mus[100,,] = group_mus
  proposal_var_rho1 = proposal_var
  window_memory_rho1 = c(rep(NA, 99), rho1)
  proposal_var_rho2 = proposal_var
  window_memory_rho2 = c(rep(NA, 99), rho2)
  
  # sampling
  for (iter in 1:(nb + ns * njump)) {
    # print intermediate sampling info
    if (verbose && (iter %% floor((nb + ns * njump) / 50)) == 0) {
      cat("iteration: ", iter, "\n")
    }
    
    # covariance matrices of different levels
    Sigma1 = diag(sqrt(vs1)) %*% matrix(c(1, rho1, rho1, 1), 2, 2) %*% diag(sqrt(vs1))
    Sigma2 = diag(sqrt(vs2)) %*% matrix(c(1, rho2, rho2, 1), 2, 2) %*% diag(sqrt(vs2))
    
    # update rho1
    new_rho1 = rho1 + sqrt(proposal_var_rho1) * rnorm(1)
    new_Sigma1 = diag(sqrt(vs1)) %*% matrix(c(1, new_rho1, new_rho1, 1), 2, 2) %*% diag(sqrt(vs1))
    if (abs(new_rho1) < 1) {
      new_lliks = getGroupLiks2D(group_summaries, group_mus, vs1, new_rho1)
      new_lpost = sum(new_lliks) + dbeta((new_rho1 + 1) / 2, 2, 2, log = TRUE)
      prv_lpost = sum(prv_lliks) + dbeta((rho1 + 1) / 2, 2, 2, log = TRUE)
      ars[K + 1] = min(exp(new_lpost - prv_lpost), 1)
      if (runif(1) <= ars[K + 1]) {
        rho1 = new_rho1
        prv_lliks = new_lliks
      }
    } else {
      ars[K + 1] = 0
    }
    
    # update group_mus
    for (k in 1:K) {
      new_group_mus[, k] = rmvnorm(1, group_mus[, k], proposal_vars_group_mus[,, k])
    }
    new_lliks = getGroupLiks2D(group_summaries, new_group_mus, vs1, rho1)
    new_lposts = new_lliks + apply(new_group_mus, 2, function(x) {
      dmvnorm(x, vmu0, Sigma2, log = TRUE)
    })
    prv_lposts = prv_lliks + apply(group_mus, 2, function(x) {
      dmvnorm(x, vmu0, Sigma2, log = TRUE)
    })
    for (k in 1:K) {
      ars[k] = min(exp(new_lposts[k] - prv_lposts[k]), 1) # symmetric proposal
      if (runif(1) <= ars[k]) {
        group_mus[, k] = new_group_mus[, k]
        prv_lliks[k] = new_lliks[k]
      }
    }
    
    # update rho2
    new_rho2 = rho2 + sqrt(proposal_var_rho2) * rnorm(1)
    new_Sigma2 = diag(sqrt(vs2)) %*% matrix(c(1, new_rho2, new_rho2, 1), 2, 2) %*% diag(sqrt(vs2))
    if (abs(new_rho2) < 1) {
      new_lliks_lvl2 = sum(dmvnorm(t(group_mus), vmu0, new_Sigma2, log = TRUE))
      new_lpost_lvl2 = new_lliks_lvl2 + dbeta((new_rho2 + 1) / 2, 2, 2, log = TRUE)
      prv_lpost_lvl2 = sum(dmvnorm(t(group_mus), vmu0, Sigma2, log = TRUE)) + dbeta((rho2 + 1) / 2, 2, 2, log = TRUE)
      ars[K + 2] = min(exp(new_lpost_lvl2 - prv_lpost_lvl2), 1)
      if (runif(1) <= ars[K + 2]) {
        rho2 = new_rho2
      }
    } else {
      ars[K + 2] = 0
    }
    
    # store samples
    if (iter > nb & (iter - nb) %% njump == 0) {
      mgroup_mus = cbind(mgroup_mus, c(group_mus))
      vrho1 = c(vrho1, rho1)
      vrho2 = c(vrho2, rho2)
      mar = cbind(mar, ars)
    }
    
    # adapting proposal sd
    if (iter <= nb) {
      # only adapt the proposal variance during burn-in to ensure ergodicity
      window_memory_group_mus[1:99,,] = window_memory_group_mus[2:100,,]
      window_memory_group_mus[100,,] = group_mus
      window_memory_rho1[1:99] = window_memory_rho1[2:100]
      window_memory_rho1[100] = rho1
      window_memory_rho2[1:99] = window_memory_rho2[2:100]
      window_memory_rho2[100] = rho2
      historical_var_rho1 = var(window_memory_rho1, na.rm = TRUE)
      historical_var_rho2 = var(window_memory_rho2, na.rm = TRUE)
      n_not_na = sum(!is.na(window_memory_rho1))
      smoothing_weight = n_not_na / (n_not_na + 1)
      for (k in 1:K) {
        tmp = na.omit(window_memory_group_mus[,, k])
        cov_tmp = cov(tmp)
        if (all(!is.na(cov_tmp)) && all(eigen(cov_tmp)$values > 0)) {
          # if cov_tmp is well defined and is positive definite
          proposal_vars_group_mus[,, k] = smoothing_weight * cov_tmp + (1 - smoothing_weight) * proposal_vars_group_mus[,, k]
        }
      }
      if (historical_var_rho1 > 0) {
        proposal_var_rho1 = smoothing_weight * historical_var_rho1 + (1 - smoothing_weight) * proposal_var_rho1
        # to prevent proposal_var_rho1 to be too small
        proposal_var_rho1 = proposal_var_rho1 + (proposal_var_rho1 < small_var) * (small_var - proposal_var_rho1)
      }
      if (historical_var_rho2 > 0) {
        proposal_var_rho2 = smoothing_weight * historical_var_rho2 + (1 - smoothing_weight) * proposal_var_rho2
        # to prevent proposal_var_rho2 to be too small
        proposal_var_rho2 = proposal_var_rho2 + (proposal_var_rho2 < small_var) * (small_var - proposal_var_rho2)
      }
    }
  }
  
  return(list(mgroup_mus = mgroup_mus, vrho1 = vrho1, vrho2 = vrho2, mar = mar))
}

# run one experiment
experimentOnce <- function(n, rho, lam, nu, Psi){
  # Args:
  #   n: sample size
  #   rho: correlation parameter
  #   lam: parameter in NIW that controls the prior concentration of mu
  #   nu: d.f. of inverse-wishart
  #   Psi: scale matrix of inverse-wishart
  
  # simulate dataset
  p = 2 # without loss of generosity, we focus on 2D case
  mu = rnorm(p)
  diag_Sigma = 0.5+rnorm(2)^2
  Sigma = diag(sqrt(diag_Sigma)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sqrt(diag_Sigma))
  y = rmvnorm(n, mu, Sigma)
  
  # compute summary statistics
  ybar = apply(y, 2, mean)
  S = cov(y)*(n-1)
  
  # posterior distribution of the full model
  lam_f = lam + n
  nu_f = nu + n
  mu_f = (n * ybar) / lam_f
  Psi_f = Psi + S + lam * n / lam_f * ybar %*% t(ybar)
  
  # knot posteriors
  # knot posteriors are the same as the true marginal posteriors, hence we only need to sample
  # marginal distribution of IW, refer to https://en.wikipedia.org/wiki/Inverse-Wishart_distribution#Theorems
  ns_knot = 1000 # number of posterior samples
  samples_knots = array(0, c(2, ns_knot, p))
  for (s in 1:p) {
    tmp_sigma = rinvgamma(ns_knot, nu_f / 2, Psi_f[s, s] / 2)
    samples_knots[1, , s] = tmp_sigma
    samples_knots[2, , s] = mu_f[s] + sqrt(tmp_sigma / lam_f) * rnorm(ns_knot)
  }
  
  # prepare for sampling tiles
  njump_tile = 20 # thinning parameter
  nb_tile = 1000 # number of burn-ins
  ns_tile = 1000 # number of posterior samples to collect
  knots_to_feed_tiles = array(0, c(2, nb_tile + njump_tile * ns_tile, p))
  for (s in 1:p) {
    knots_to_feed_tiles[,,s] = samples_knots[, sample(1:ns_knot, nb_tile + njump_tile * ns_tile, replace = TRUE), s]
  }
  
  # tile posteriors
  tile12 = sampleTile(y[,1], y[,2], knots_to_feed_tiles[,,1], knots_to_feed_tiles[,,2], 
                      nb_tile, ns_tile, njump_tile, lam, nu, Psi, 0.1)
  
  return(list(
    rhos = tile12$rhos, 
    mus = tile12$mus, 
    vars = tile12$vars, 
    ars = tile12$ars,
    lam_f = lam_f,
    nu_f = nu_f,
    mu_f = mu_f,
    Psi_f = Psi_f
  ))
}