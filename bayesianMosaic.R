# Bayesian Mosaic
# Version 1.0

suppressMessages(require(testthat))
suppressMessages(require(parallel))
ncores = detectCores() - 1

# Helper Functions----
ilogit <- function(x){
  # inverse logit function
  # args:
  #   x: input value
  # outputs:
  #   inverse logit of x
  
  return( 1/(1+exp(-x)) )
}

aggregateCount <- function(y) {
  # aggregate sparse count data
  # args:
  #   y: vector of count observations
  # outputs:
  #   2-by-m matrix with 1st row being values and 2nd row being counts
  
  pool = unique(y)
  summaries = NULL
  for (i in 1:length(pool)){
    summaries = cbind(summaries, c(pool[i], sum(y==pool[i])))
  }
  
  return(summaries)
}

aggregateCount2d <- function(y){
  # aggregate 2-d sparse count data
  # args:
  #   y: n-by-2 matrix of count observations
  # outputs:
  #   3-by-m matrix with 1st & 2nd row being values and 3rd row being counts
  
  k = nrow(y)
  n = ncol(y)
  summaries = NULL
  while(n > 0){
    tmp = y[, 1]
    idx = apply(y, 2, FUN = function(x) {
      all(x == tmp)
    })
    summaries = cbind(summaries, c(tmp, sum(idx)))
    y = y[, !idx, drop=FALSE]
    n = n - sum(idx)
  }
  
  return(summaries)
}

getlik_mvtpois <- function(y, mu, v){
  # Compute the log-likelihood of a Poisson log-normal
  # likelihood function value for a individual observation via numerical integration
  # transform the variable so that the range is 0 to 1, have tried -Inf to Inf 
  # and there is some bizzare numerical issue there.
  # args:
  #   y: vector of oberservations
  #   mu: gaussian mean
  #   v: gaussian variance
  # outputs:
  #   likelihood
  
  integrand <- function(x){
    return(exp(y * x - lfactorial(y) - (x - mu)^2 / 2 / v - exp(x)) * (2 * pi * v)^-.5)
  }
  
  integrate(integrand, lower = min(log(y + 0.1), mu) - 4 * sqrt(v), 
            upper = max(log(y + 0.1), mu) + 4 * sqrt(v), stop.on.error = FALSE)$value
}

getLikMvtpois <- function(dataSummary, mu, v, logarithm = TRUE){
  # get (the logarithm of) the data likelihood
  # args:
  #   y: vector of oberservations
  #   mu: gaussian mean
  #   v: gaussian variance
  #   logarithm: bool indicate whether or not to return the logarithm
  # outputs:
  #   (the logarithm of) the data likelihood
  
  ret = 0
  for (i in 1:length(dataSummary)){
    ret = ret + dataSummary[2, i] * log(getlik_mvtpois(dataSummary[1, i], mu, v))
  }
  if (logarithm == FALSE){
    ret = exp(ret)
  }
  return(ret)
}

getlik_mvtpois_2d <- function(y, mu, Sigma, Sinv, v){
  # get the log-likelihood of a 2-d Poisson log-normal
  # likelihood function value for a individual observation via numerical integration
  # TODO
  
  integrand <- function(x1, x2){
    return(exp(y[1] * x1 + y[2] * x2 - lfactorial(y[1]) - lfactorial(y[2]) -
                  ((x1 - mu[1])^2 * Sinv[1, 1] + 2 * (x1 - mu[1]) * (x2 - mu[2]) * Sinv[1, 2] + (x2 - mu[2])^2 * Sinv[2, 2]) / 2 -
                  exp(x1) - exp(x2)) * (2 * pi)^(-1) * det(Sigma)^-.5)
  }
  return(quad2d(integrand, min(log(y[1] + 0.1), mu[1]) - 4 * sqrt(v[1]), 
                 max(log(y[1] + 0.1), mu[1]) + 4 * sqrt(v[1]), min(log(y[2] + 0.1), mu[2]) - 4 * sqrt(v[2]), 
                 max(log(y[2] + 0.1), mu[2]) + 4 * sqrt(v[2]), 64))
}

getLikMvtpois2d <- function(dataSummary, mu, v, rho){
  # TODO
  
  p = length(mu)
  Sigma = diag(sqrt(v)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sqrt(v))
  Sinv = solve(Sigma)
  sum(sapply(1:ncol(dataSummary), FUN = function(k) {
    dataSummary[p + 1, k] * log(getlik_mvtpois_2d(dataSummary[1:p, k], mu, Sigma, Sinv, v))
  }))
}

# Model 1: Multivariate Poisson (mvtpois)----
genKnotMvtpois <- function(y, verbose = FALSE){
  # M-H sampler based on numerical integration
  # args:
  #   y: vector of oberservations
  #   verbose: TODO
  # outputs:
  #   TODO
  
  # hyper-parameters
  nb = 1000 # number of burn-ins
  nc = 1000 # number of MCMC samples
  sig=1e3 # prior variance for gaussian mean
  a=2 # gamma prior shape parameter
  b = 2 # gamma prior scale parameter
  ps_m = 1 # TODO
  ps_s = 1 # TODO
  njump = 1 # TODO
  
  dataSummaries = aggregateCount(y)
  # Initialize mu and s
  m = 0
  s = 1
  # storing the results
  vs = NULL
  vm = NULL
  ars = rep(0,nb+nc)
  
  pllik = getLikMvtpois(dataSummaries, m, s) + dnorm(m, 0, sqrt(sig), TRUE) +
    dgamma(s, a, b, log=TRUE)
  # iterative
  for (iter in 1:(nb + nc)){
    if (verbose && iter%%ceiling((nb+nc) / 100) == 0){
      cat("iteration", iter, "\n")
    }
    
    # propose new m, s and v
    newm = m + sqrt(ps_m) * rnorm(1)
    news = rlnorm(1, log(s) - ps_s / 2, sqrt(ps_s))
    # evaluate new log-likelihood
    nllik = getLikMvtpois(dataSummaries, newm, news) + dnorm(newm, 0, sqrt(sig), TRUE) +
      dgamma(news, a, b, log = TRUE)
    # acceptance rate
    ars[iter] = min(1,exp(nllik - pllik +
                            dlnorm(s, log(news) - ps_s / 2, sqrt(ps_s), log = TRUE) -
                            dlnorm(news, log(s) - ps_s / 2, sqrt(ps_s), log = TRUE)))
    if (runif(1) <= ars[iter]){
      m = newm
      s = news
      pllik = nllik
    }
    
    vm = c(vm, m)
    vs = c(vs, s)
  }
  
  if (njump > 1){
    vm = vm[seq(1, nc, njump)]
    vs = vs[seq(1, nc, njump)]
  }
  
  return(list(vmu = vm, vs = vs, ars = ars))
}

genTileMvtpois <- function(y, knot1, knot2, verbose = FALSE) {
  # Adaptive Metropolis (AM) sampler based on numerical integration
  # TODO
  
  # hyper-parameters
  nb = 1000 # number of burn-ins
  nc = 1000 # number of MCMC samples
  sd_mh = 1 # standard deviation of the mh proposal
  sd_rho = 5 # TODO
  njump = 1 # TODO
  
  dataSummaries = aggregateCount2d(y)
  p = nrow(vy)
  t0 = floor(nb / 4)
  t1 = floor(nb / 2)
  nm = nb + nc
  nMargin = ncol(vmu)
  # storing the log acceptance ratios
  ars = rep(0, nb + nc)
  
  ldens <- function(x, mu, s){
    # log density function of the target distribution
    # up to a normalizing constant
    return(getLikMvtpois2d(dataSummaries, mu, s, 1.9 * ilogit(x) - 0.95) +
              dnorm(x, 0, sd_rho, log = TRUE))
  }
  
  # preparation
  sp = 0
  sps = NULL
  mean_history = 0
  var_history = sd_mh^2
  count_history = 5
  for (iter in 1:nm){
    if (verbose && iter %% ceiling((nm) / 100) == 0){
      cat("iteration", iter, "\n")
    }
    
    idx_empirical = sample(1:nMargin, 1)
    cur_mu = c(knot1$vmu[idx_empirical], knot2$vmu[idx_empirical])
    cur_s = c(knot1$vs[idx_empirical], knot2$vs[idx_empirical])
    
    # propose new sample
    if (iter <= t1){
      nsp = sp + sd_mh * rnorm(1)
    }else{
      nsp = sp + sqrt(var_history) * rnorm(1)
    }
    ars[iter] = ldens(nsp, cur_mu, cur_s) - ldens(sp, cur_mu, cur_s)
    if (log(runif(1)) < ars[iter]){
      sp = nsp
    }
    # update the proposal variance
    if (iter >= t0 && iter <= nb){
      tmp = var_history * (count_history - 1) + count_history * mean_history^2
      mean_history = (count_history * mean_history + sp) / (count_history + 1)
      count_history = count_history + 1
      if (count_history == 1){
        var_history = 0
      }else{
        var_history = (tmp + sp^2 - count_history * mean_history^2) / (count_history - 1)
      }
    }
    
    sps = c(sps, sp)
  }
  
  vrho = 1.9 * ilogit(sps) - 0.95
  if (njump > 1){
    vrho = vrho[seq(1, nc, njump)]
  }
  return(list(vrho = vrho, ars = ars))
}

# model caller----
callModel <- function(y, X, stage, model, knot1 = NULL, knot2 = NULL) {
  # wrapper function that calles model function at a certain stage
  # args:
  #   y: TODO
  #   X: TODO
  #   knot1: TODO
  #   knot2: TODO
  #   stage: TODO
  #   model: TODO
  # outputs:
  #   TODO
  
  if (model == 'mvtpois') {
    if (stage = 'knot') {
      genKnotMvtpois(y, FALSE)
    } else {
      genTileMvtpois(y, knot1, knot2, FALSE)
    }
  }
}

# main functions----
genKnots <- function(y, X = NULL, model) {
  # generate the knots that connect titles
  # args:
  #   y: a data frame (or object coercible by as.data.frame to a data frame) containing the data for all dimensions
  #   X: a data frame (or object coercible by as.data.frame to a data frame) containing extra features
  #   model: a string indicating the model to apply Bayesian Mosaic on
  # output:
  # TODO
  
  knots = list()
  dim_y = dim(y)
  # validation
  if (dim_y[1] < dim_y[2]) {
    stop("sample size smaller than number of dimensions")
  }
  nd = dim_y[2]
  
  # genKnot for each dimension
  for (j in 1:nd) {
    knots[[1]] = callModel(y[, j], X, 'knot', model)
  }
  
  # return knots
  return(knots)
}

genTiles <- function(y, X = NULL, knots, model, idx_pairs, knots, parallel = FALSE) {
  # generate tiles that later will be collaged into a mosaic
  # args:
  #   y: a data frame (or object coercible by as.data.frame to a data frame) containing the data for all dimensions
  #   X: a data frame (or object coercible by as.data.frame to a data frame) containing extra features
  #   knots: TODO
  #   model: a string indicating the model to apply Bayesian Mosaic on
  #   indices: a vector of indices, tiles of all pairs of these indices will be generated
  #   knots: TODO
  #   parallel: bool indicating whether or not to parallelize generating these tiles
  # output:
  # TODO
  
  tiles = list()
  num_pairs = nrow(idx_pairs)
  # validation
  if (ncol(idx_pairs) != 2) {
    stop("idx_pairs should have two columns")
  }
  # TODO: validate knots is consistent with y
  
  # generate tile for all pairs of given indices
  for (p in 1:num_pairs) {
    idx_pair = idx_pairs[p, ]
    tiles[[p]] = callModel(y[, idx_pair], X, 'tile', model, knots[[idx_pair[1]]], knots[[idx_pair[2]]])
  }
  
  # return tiles
  return(tiles)
}