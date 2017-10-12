# Bayesian Mosaic Simulation Study 2
# Bivariate log-normal mixture of Poisson

rm(list=ls())

# dependencies
suppressMessages(require(MCMCpack))
suppressMessages(require(mvtnorm))
suppressMessages(require(pracma))
suppressMessages(require(ggplot2))
suppressMessages(require(compiler))
suppressMessages(require(ars))

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

# sample from the knot marginal
sampleKnot <- function(y, nb, ns, njump, proposal_var, verbose = FALSE) {
  # sample from the knot marginal
  # Args:
  #   y: observation
  #   shoot: group index of each observation
  #   nb: number of burn-ins
  #   ns: number of samples to collect
  #   njump: thinning parameter
  #   proposal_var: initial variance of the proposal normal distribution
  #   verbose: whether or not to print intermediate sampling info
  
  # summary statistics
  group_summary = aggregateData(y)
  
  # initialization
  s1 = 1
  mu0 = 0
  prv_llik = -Inf
  ars = rep(0, 2)
    
  # storing the results
  mgroup_mus = NULL
  vs1 = NULL
  vmu0 = NULL
  mar = NULL
  
  # preparation for adapting
  proposal_var_mu0 = proposal_var
  window_memory_mu0 = c(rep(NA, 99), mu0)
  proposal_var_s1 = proposal_var
  window_memory_s1 = c(rep(NA, 99), s1)
  
  # sampling
  for (iter in 1:(nb + ns * njump)) {
    # print intermediate sampling info
    if (verbose && (iter %% floor((nb + ns * njump) / 10)) == 0) {
      cat("iteration: ", iter, "\n")
    }
    
    # update mu0
    new_mu0 = mu0 + sqrt(proposal_var_mu0) * rnorm(1)
    new_llik = getGroupLik(group_summary, new_mu0, s1)
    new_lpost = new_llik + dnorm(new_mu0, 0, 1E4, log = TRUE)
    prv_lpost = prv_llik + dnorm(mu0, 0, 1E4, log = TRUE)
    ars[1] = min(exp(new_lpost - prv_lpost), 1) # symmetric proposal
    if (runif(1) <= ars[1]) {
      mu0 = new_mu0
      prv_llik = new_llik
    }
    
    # update s1
    new_s1 = s1 + sqrt(proposal_var_s1) * rnorm(1)
    if (new_s1 > 0) {
      new_llik = getGroupLik(group_summary, mu0, new_s1)
      new_lpost = new_llik + dgamma(1 / new_s1, shape = 2, rate = 2, log = TRUE)
      prv_lpost = prv_llik + dgamma(1 / s1, shape = 2, rate = 2, log = TRUE)
      ars[2] = min(exp(new_lpost - prv_lpost), 1)
      if (runif(1) <= ars[2]) {
        s1 = new_s1
        prv_llik = new_llik
      }
    } else {
      ars[2] = 0
    }
    
    # store samples
    if (iter > nb & (iter - nb) %% njump == 0) {
      vs1 = c(vs1, s1)
      vmu0 = c(vmu0, mu0)
      mar = cbind(mar, ars)
    }
    
    # adapting proposal sd
    if (iter <= nb) {
      # only adapt the proposal variance during burn-in to ensure ergodicity
      window_memory_mu0[1:99] = window_memory_mu0[2:100]
      window_memory_mu0[100] = mu0
      window_memory_s1[1:99] = window_memory_s1[2:100]
      window_memory_s1[100] = s1
      historical_var_mu0 = var(window_memory_mu0, na.rm = TRUE)
      historical_var_s1 = var(window_memory_s1, na.rm = TRUE)
    }
  }
  
  return(list(vs1 = vs1, vmu0 = vmu0, mar = mar))
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

getLik2DGradient <- function(y, mu, Sigma, Sinv){
  # likelihood function value for a individual observation via numerical integration
  # args:
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

getGroupLLik2DGradient <- function(group_summary, mu, v, rho){
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
    ret = ret + group_summary[3, i] * getLik2DGradient(group_summary[1:2, i], mu, Sigma, Sinv) / getLik2D(group_summary[1:2, i], mu, Sigma, Sinv)
  }
  return(ret)
}

findZeroGradient <- function(fx) {
  # find the zero location of a monotone function
  # args:
  #   fx: some monotone function that crosses zero
  
  left = -.99
  right = .99
  counter = 0
  while (right - left >= 1E-4) {
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

sampleLaplaceApprox <- cmpfun(function(data_summary, mu, vvars) {
  # sample from the Laplace approximation
  # args:
  #   group_summary: aggregated summaries for the dataset
  #   mu: mean vector
  #   vvars: vector of variances
  
  x_sample = -Inf
  
  while (abs(x_sample) >= 1) {
    llikGradient <- function(x) {
      -getGroupLLik2DGradient(data_summary, mu, vvars, x) - 1 / (1 + x) + 1 / (1 - x)
    }
    
    # maximize the log-likelihood
    optim_res = findZeroGradient(llikGradient)
    
    # laplace approximation
    laplace_mean = optim_res[1]
    laplace_var = 1 / c(optim_res[2])
    
    if (laplace_var <= 0) {
      return(NA)
    }
    
    # sample once
    x_sample = laplace_mean + sqrt(laplace_var) * rnorm(1)
  }

  return(x_sample)
})

# sample from the tile conditionals
sampleTile <- function(Y, ms1, mmu0, verbose = FALSE){
  # sample from the tile posterior
  # Args:
  #   Y: dimension s of y
  #   ms1: posterior draws of the variances from the knot marginal
  #   mmu0: posterior draws of the means from the knot marginal
  #   verbose: whether or not to print intermediate sampling info
  
  # summary statistics
  data_summary = aggregateData(Y)
  vrho = NULL
  
  # sampling
  for (iter in 1:nrow(ms1)) {
    # print intermediate sampling info
    if (verbose && (iter %% floor(nrow(ms1) / 100)) == 0) {
      cat("iteration: ", iter, "\n")
    }
    
    vrho = c(vrho, sampleLaplaceApprox(data_summary, mmu0[iter,], ms1[iter,]))
  }
  
  return(vrho)
}

# DA MCMC to directly sample from the actual posterior
sampleViaDAMCMC <- function(Y, ns, verbose = FALSE){
  n = nrow(Y)
  # initialization
  Sigma = diag(1, 2)
  mu = rep(0, 2)
  latent_x = matrix(0, n, 2)
  
  # outputs
  mv = NULL
  vrho = NULL
  mmu = NULL
  
  # helpers
  f <- function(x, y, mu, s) {
    y * x - exp(x) - (x - mu)^2 / 2 / s
  }
  fprima <- function(x, y, mu, s) {
    y - exp(x) - (x - mu) / s
  }
  
  # sampling
  for (iter in 1:ns) {
    # print intermediate sampling info
    if (verbose && (iter %% max(floor(ns / 100), 1) == 0)) {
      cat("iteration: ", iter, "\n")
    }
    
    # sample latent x
    for (i in 1:n) {
      mux1 = mu[1]+Sigma[1,2]/Sigma[2,2]*(latent_x[i,2]-mu[2])
      sx1 = Sigma[1,1] - Sigma[1,2]^2/Sigma[2,2]
      latent_x[i,1] = ars(n=1,f=f,fprima=fprima,x=c(-10,1,10),ns=100,m=3,emax=64,lb=FALSE,ub=FALSE,
                          y=Y[i,1],mu=mux1,s=sx1)
      mux2 = mu[2]+Sigma[1,2]/Sigma[1,1]*(latent_x[i,1]-mu[1])
      sx2 = Sigma[2,2] - Sigma[1,2]^2/Sigma[1,1]
      latent_x[i,2] = ars(n=1,f=f,fprima=fprima,x=c(-10,1,10),ns=100,m=3,emax=64,lb=FALSE,ub=FALSE,
                          y=Y[i,2],mu=mux2,s=sx2)
    }
    
    # sample mu
    Sigma_mu = solve(n * solve(Sigma) + .01*diag(1,2))
    mu = rmvnorm(1, Sigma_mu%*%apply(latent_x%*%solve(Sigma), 2, sum), Sigma_mu)
    
    # sample Sigma
    latent_x_demeaned = apply(latent_x, 1, function(x) {
      x - mu
    })
    Sigma = riwish(6 + n, diag(1,2) + latent_x_demeaned%*%t(latent_x_demeaned))
    
    mv = cbind(mv, diag(Sigma))
    vrho = c(vrho, Sigma[1,2]/sqrt(Sigma[1,1]*Sigma[2,2]))
    mmu = cbind(mmu, c(mu))
  }
  
  return(list(mv=mv,vrho=vrho,mmu=mmu))
}

# simulation
set.seed(2018)
n = 1000
p = 2 # without loss of generosity, we focus on 2D case
mu = c(-2, -3)
diag_Sigma = 0.5+0.5 * rnorm(2)^2
rho = .8
Sigma = diag(sqrt(diag_Sigma)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sqrt(diag_Sigma))
x = rmvnorm(n, mu, Sigma)
y = matrix(rpois(2 * n, exp(c(x))), n, 2)

# fit Bayesian mosaic
nb = 2000
ns = 500
njump = 10
proposal_var = .1
ct_knot = proc.time()
knot1 = sampleKnot(y[, 1], nb, ns, njump, proposal_var, verbose = FALSE)
knot2 = sampleKnot(y[, 2], nb, ns, njump, proposal_var, verbose = FALSE)
ct_knot = proc.time() - ct_knot
ct_tile = proc.time()
tile12 = sampleTile(y, cbind(knot1$vs1, knot2$vs1), cbind(knot1$vmu0, knot2$vmu0))
ct_tile = proc.time() - ct_tile

ct_damcmc = proc.time()
ns_damcmc = 60000
res_damcmc = sampleViaDAMCMC(y, ns_damcmc, verbose = TRUE)
ct_damcmc = proc.time() - ct_damcmc

# compare the traceplot of log likelihood
nb_damcmc = 10000
v_da = apply(res_damcmc$mv[, (nb_damcmc+1):ns_damcmc], 1, mean)
v_bm = c(mean(knot1$vs1), mean(knot2$vs1))
mu_da = apply(res_damcmc$mmu[, (nb_damcmc+1):ns_damcmc], 1, mean)
mu_bm = c(mean(knot1$vmu0), mean(knot2$vmu0))
rho_da = mean(res_damcmc$vrho[(nb_damcmc+1):ns_damcmc])
rho_bm = mean(tile12)
group_summary = aggregateData(y)
getGroupLLik2D(group_summary, mu_da, v_da, rho_da)
getGroupLLik2D(group_summary, mu_bm, v_bm, rho_bm)

# traceplots
par(mar = c(2.5, 3, 3, 3), mfrow=c(2,1))
ratio_ct = (ct_knot[3] / 2 + ct_tile[3]) / ct_damcmc[3]
plot(1:500, tile12, 'l', main = "Bayesian mosaic", xlab = "", ylab = "rho")
ns_damcmc = length(res_damcmc$vrho)
plot(1:floor(ns_damcmc * ratio_ct), res_damcmc$vrho[1:floor(ns_damcmc * ratio_ct)], 'l', main = "DAMCMC", xlab = "", ylab = "rho")
par(mfrow=c(1,1))

# compare the effective sample size
effectiveSize(tile12)
effectiveSize(res_damcmc$vrho[1:floor(ns_damcmc * ratio_ct)])
effectiveSize(res_damcmc$vrho)
