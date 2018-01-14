# Other Samplers
# Version 1.1
# Last Updated on Jan 7, 2018

# dependencies
source("~/Documents/yw_git/bayesian_mosaic/sampler-helpers.R")
# source("/home/collabor/yw104/BayesianMosaic/sampler-helpers.R")
suppressMessages(require(ars))

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
  sample_diag = NULL
  sample_corr = NULL
  
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
      return(list(sample_mu=sample_mu, sample_diag=sample_diag,
                  sample_corr=sample_corr, ealry_stop=TRUE))
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
    sample_diag = rbind(sample_diag, diag(Sigma))
    Corr = diag(sqrt(1/diag(Sigma)))%*%Sigma%*%diag(sqrt(1/diag(Sigma)))
    sample_corr = rbind(sample_corr, lowerOffDiagonal(Corr))
  }
  
  return(list(sample_mu=sample_mu, sample_diag=sample_diag,
              sample_corr=sample_corr, ealry_stop=FALSE))
}

# MCMC based on 2-d integration
sampleVia2DInt <- function(y, ns, proposal_var, genIndividualLik, stop_time=NULL,
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
    new_llik = genLLik(compressed_y=compressed_y, genIndividualLik=genIndividualLik, 
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
      new_llik = genLLik(compressed_y=compressed_y, genIndividualLik=genIndividualLik, 
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
      new_llik = genLLik(compressed_y=compressed_y, genIndividualLik=genIndividualLik, 
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

# DA-MCMC on rounded Gaussian
damcmcRoundedGaussian <- function(Y, ns, inits, stop_time=NULL, verbose=FALSE, parallel=FALSE){
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
  sample_diag = NULL
  sample_corr = NULL
  
  # sampling
  for (iter in 1:ns) {
    # early stop
    if (!is.null(stop_time) & (proc.time()-start_time)[3]>stop_time) {
      return(list(sample_mu=sample_mu, sample_diag=sample_diag,
                  sample_corr=sample_corr, ealry_stop=TRUE))
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
    sample_diag = rbind(sample_diag, diag(Sigma))
    Corr = diag(sqrt(1/diag(Sigma)))%*%Sigma%*%diag(sqrt(1/diag(Sigma)))
    sample_corr = rbind(sample_corr, lowerOffDiagonal(Corr))
  }
  
  return(list(sample_mu=sample_mu, sample_diag=sample_diag,
              sample_corr=sample_corr, ealry_stop=FALSE))
}
