# Bayesian Mosaic
# Version 1.1
# Last Updated on Dec 31, 2017

suppressMessages(require(rstan))
# source("~/Documents/yw_git/bayesian_mosaic/sampler-helpers.R")
source("/home/collabor/yw104/BayesianMosaic/sampler-helpers.R")

# # Covariance Matrix Enforcer----
# code_corr_correct = '
# data {
# int<lower=0> N; // number of draws from Bayesian mosaic
# int<lower=2> K; // number of species
# int<lower=2> K_choose_2;
# vector[N] corr_mats[K, K]; // predictor matrix
# }
# parameters {
# corr_matrix[K] Omega; // correlation matrix
# vector<lower=0>[K_choose_2] L_sigma;
# }
# model {
# for (x in 1:(K - 1)) {
# for (y in (x + 1):K) {
# corr_mats[x, y] ~ normal(Omega[x, y], sqrt(L_sigma[(K - x + 1 + K - 1) * (x - 1) / 2 + y - x]));
# }
# }
# Omega ~ lkj_corr(2.0); // regularize to unit correlation
# L_sigma ~ gamma(2, 2);
# }
# '
# 
# corr_correct = stan_model(model_name = "correlation_correction", model_code = code_corr_correct)
# 
# correctCorrelation <- function(sample_corr, ns) {
#   # TODO
#   
#   p = dim(sample_corr)[1]
#   n = dim(sample_corr)[3]
#   
#   complete_idx = complete.cases(matrix(c(sample_corr), n, p*p, byrow=TRUE))
#   n_complete = length(complete_idx)
#   
#   data = list(N=n_complete, K=p, K_choose_2=p*(p-1)/2, corr_mats=sample_corr[,,complete_idx])
#   fit <- sampling(corr_correct, data=data, chains=1, iter=ns*2, cores=1)
#   return(apply(extract(fit, pars = 'Omega')$'Omega', c(2, 3), mean))
# }
# 
# fij <- function(nu, pm_ij, pm_ii, pm_jj, var_ij) {
#   #TODO
#   
#   (sqrt(((nu+4)*pm_ij^2+(nu+2)*pm_ii*pm_jj)/nu/(nu+3))-sqrt(var_ij))^2
# }
# 
# findBestIWishart <- function(sample_diag, sample_corr, corrected_Correlation) {
#   # TODO
#   
#   p = ncol(sample_diag)
#   pm_diag_Sigma = apply(sample_diag, 2, mean)
#   corrected_Covariance = diag(sqrt(pm_diag_Sigma))%*%corrected_Correlation%*%diag(sqrt(pm_diag_Sigma))
#   
#   
#   fn <- function(nu) {
#     ret = 0
#     for (x in 1:p) {
#       for (y in x:p) {
#         ret = ret + fij(nu, corrected_Covariance[x, y], corrected_Covariance[x, x], 
#                         corrected_Covariance[y, y], 
#                         var(sample_corr[x,y,]*sqrt(sample_diag[,x]*sample_diag[,y]), na.rm=TRUE))
#       }
#     }
#     return(ret)
#   }
#   
#   nu_optim = optim(par = 1, fn = fn, method = "L-BFGS-B", lower = 0)$par
#   
#   return(list(v=nu_optim+p+1, 
#               S=corrected_Covariance*nu_optim))
# }

# LKJ Sampler----
code_lkj = '
data {
int<lower=2> K;
}
parameters {
corr_matrix[K] Omega; // correlation matrix
}
model {
Omega ~ lkj_corr(2.0);
}
'

lkj = stan_model(model_name = "lkj_sampler", model_code = code_lkj)

sampleLKJ <- function(n, p) {
  data = list(K = p)
  fit <- sampling(lkj, data = data, chains = 1, iter = n * 2)
  return(extract(fit, pars = 'Omega')$'Omega')
}

# Bayesian Mosaic
sampleKnot <- function(y, nb, ns, njump, proposal_var, inits, genIndividualLik, verbose = FALSE, ...) {
  # TODO
  
  acr_target = 0.25 # target acceptance rate
  acr_tol = 0.05 # tolerance for acceptance rate
  dcf_pv = 0.9 # discount factor for adapting the proposal variance
  dcf_multiplier = 0.9 # discount factor for adapting the multiplier
  dcf_acr = 0.99 # discount factor for averaging the historical acceptance rate
  ada_prop = 0.5 # the proportion of burn-ins for which adapting the proposal variance is based on
  
  # compress the counts
  compressed_y = compressCount(y)
  
  # initialization
  mu = inits$mu
  s = inits$s
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
      new_llik = genLLik(compressed_y=compressed_y, genIndividualLik=genIndividualLik,
                         mu = new_mu, v = new_s, ...)
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
    if (iter<(nb*(1+ada_prop)/2) & iter>=(nb*(1-ada_prop)/2)) {
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

bayesianMosaic <- function(Y, nb, ns, njump, proposal_var, model, 
                           verbose=0, parallel=FALSE, ...) {
  # TODO
  
  args = as.list(sys.call())
  p = ncol(Y)
  genIndividualLik1D = getFcn(model, "IndividualLik1D")
  genIndividualLik2D = getFcn(model, "IndividualLik2D")
  initial_mu = apply(Y, 2, FUN=function(y){
    if (model=="mvtPoisson") {
      mu = log(mean(y)+1E-4)
    } else if (model=="mvtBinomial") {
      pr = mean(y)/args$N + 1E-4
      mu = log(pr/(1-pr))
    }
  })
  
  # posterior knots
  verb = FALSE
  if (verbose > 0) {
    cat("learning posterior knots...\n")
    if (verbose > 1) {
      verb = TRUE
    }
  }
  
  genSampleKnot <- function(k) {
    if (is.null(args$Ns)) {N=NULL} else {N=args$Ns[k]}
    sampleKnot(y=Y[,k], nb=nb, ns=ns, njump=njump, proposal_var=proposal_var, 
               inits = list(mu=initial_mu[k], s=1),
               genIndividualLik=genIndividualLik1D, verbose = verb, N=N)
  }
  if (parallel) {
    knots = mclapply(1:p, FUN=genSampleKnot, mc.cores=ncores)
  } else {
    knots = list()
    for (k in 1:p) {
      knots[[k]] = genSampleKnot(k)
    }
  }
  
  # posterior tiles
  if (verbose > 0) {
    cat("learning posterior tiles...\n")
  }
  
  pairs = NULL
  for (x in 1:(p-1)) {
    for (y in (x+1):p) {
      pairs = rbind(pairs, c(x, y))
    }
  }
  n_pair = nrow(pairs)
  
  genSampleTile <- function(k){
    s = pairs[k, 1]
    t = pairs[k, 2]
    y = cbind(Y[,s], Y[,t])
    sample_mu = rbind(knots[[s]]$sample_mu, knots[[t]]$sample_mu)
    sample_s = rbind(knots[[s]]$sample_s, knots[[t]]$sample_s)
    if (model=="mvtPoisson") {
      return( sampleTile(y = y, sample_mu = sample_mu, sample_s = sample_s, 
                         model="mvtPoisson", verbose = verb) )
    } else if (model=="mvtBinomial") {
      return( sampleTile(y = y, sample_mu = sample_mu, sample_s = sample_s, 
                         model="mvtBinomial", verbose = verb, Ns = Ns[c(s, t)]) )
    }
  }
  
  if (parallel) {
    tiles = mclapply(1:n_pair, FUN=genSampleTile, mc.cores=ncores)
  } else {
    tiles = list()
    for (k in 1:n_pair) {
      tiles[[k]] = genSampleTile(k)
    }
  }
  
  # write outputs
  if (verbose > 0) {
    cat("writing outputs...\n")
  }
  
  sample_corr = matrix(unlist(tiles), ns, n_pair)
  sample_mu = matrix(0, ns, p)
  sample_diag = matrix(0, ns, p)
  for (j in 1:p) {
    sample_mu[,j] = knots[[j]]$sample_mu
    sample_diag[,j] = knots[[j]]$sample_s
  }
  
  # handle covariance matrix constraint
  if (verbose > 0) {
    cat("enforcing constraints...\n")
  }
  
  # Corr_mats = array(0,c(p,p,ns))
  # for (i in 1:ns) {
  #   Corr_mats[,,i] = vec2Corr(sample_corr[i,], p)
  # }
  # corrected_Correlation = correctCorrelation(Corr_mats, 1000)
  # best_iw_param = findBestIWishart(sample_diag, Corr_mats, corrected_Correlation)
  
  counter_c = 0
  sample_diag_c = NULL
  sample_corr_c = NULL
  for (i in 1:ns) {
    tmp_Sigma = diag(sqrt(sample_diag[i,]))%*%vec2Corr(sample_corr[i,], p)%*%
      diag(sqrt(sample_diag[i,]))
    correct_res = forcePSD(tmp_Sigma)
    counter_c = counter_c+correct_res$is_psd
    sample_diag_c = rbind(sample_diag_c, diag(correct_res$C))
    sample_corr_c = rbind(sample_corr_c, 
                          lowerOffDiagonal(Cov2CorrMat(correct_res$C)))
  }
  
  # return(list(model=model,
  #             sample_mu = sample_mu,
  #             sample_diag = sample_diag,
  #             sample_corr = sample_corr,
  #             best_iw_v=best_iw_param$v,
  #             best_iw_S=best_iw_param$S))
  
  return(list(model=model,
              sample_mu = sample_mu,
              sample_diag = sample_diag,
              sample_corr = sample_corr,
              sample_diag_c=sample_diag_c,
              sample_corr_c=sample_corr_c,
              perc_psd = counter_c/ns))
}
