# Bayesian Mosaic
# Version 1.2
# Last Updated on April 6, 2018

suppressMessages(require(rstan))
source("~/Documents/yw_git/bayesian_mosaic/sampler-helpers.R")
# source("/home/collabor/yw104/BayesianMosaic/sampler-helpers.R")

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
sampleKnot <- function(y, nb, ns, njump, inits, genIndividualLLik, verbose = FALSE, ...) {
  # TODO
  
  # compress the counts
  compressed_y = compressCount(y)
  
  # construct a proposal distribution for independence MH
  hull = constructHull(compressed_y, genIndividualLLik, mu_range=c(-10,0), 
                s_range=c(0,5), n_grid=200, ...)
  
  # initialization
  mu = inits$mu
  s = inits$s
  prv_llik = -Inf
  ld_hull = 0
  
  # storing the results
  sample_s = NULL
  sample_mu = NULL
  
  # sampling
  for (iter in 1:(nb + ns * njump)) {
    # print intermediate sampling info
    if (verbose && (iter %% max(floor((nb + ns * njump) / 100), 1)) == 0) {
      cat("iteration: ", iter, "\n")
    }
    
    # update mu & s
    # independence MH
    tmp = sampleHull(hull)
    new_mu = tmp[1]
    new_s = tmp[2]
    new_ld_hull = tmp[3]
    if (new_s > 0 & new_s < 10 & abs(new_mu) < 100) {
      new_llik = genLLik(compressed_y=compressed_y, 
                         genIndividualLLik=genIndividualLLik,
                         mu = new_mu, v = new_s, ...)
      new_lpost = new_llik + ld_hull
      prv_lpost = prv_llik + new_ld_hull
      acr = min(exp(new_lpost - prv_lpost), 1) # symmetric proposal
      if (runif(1) <= acr) {
        mu = new_mu
        s = new_s
        prv_llik = new_llik
        ld_hull = new_ld_hull
      }
    }
    
    # store samples
    if (iter > nb & (iter - nb) %% njump == 0) {
      sample_s = c(sample_s, s)
      sample_mu = c(sample_mu, mu)
    }
  }
  
  return(list(sample_s = sample_s,
              sample_mu = sample_mu))
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

bayesianMosaic <- function(Y, nb, ns, njump, model, 
                           verbose=0, parallel=FALSE, ...) {
  # TODO
  
  args = as.list(sys.call())
  p = ncol(Y)
  genIndividualLLik1D = getFcn(model, "IndividualLLik1D")
  genIndividualLLik2D = getFcn(model, "IndividualLLik2D")
  initial_mu = apply(Y, 2, FUN=function(y){
    if (model=="mvtPoisson") {
      mu = log(mean(y)+1E-4)
    } else if (model=="mvtBinomial") {
      pr = mean(y)/args$N + 1E-4
      mu = log(pr/(1-pr))
    } else if (model=="roundedGaussian") {
      mu = mean(y)
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
    sampleKnot(y=Y[,k], nb=nb, ns=ns, njump=njump,
               inits = list(mu=initial_mu[k], s=1),
               genIndividualLLik=genIndividualLLik1D, verbose = verb, N=N)
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
    if (is.null(args$Ns)) {Ns=NULL} else {Ns=args$Ns[c(s, t)]}
    return(sampleTile(y=y, sample_mu=sample_mu, sample_s=sample_s, 
                       model=model, verbose=verb, Ns=Ns))
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
  
  return(list(model=model,
              sample_mu = sample_mu,
              sample_diag = sample_diag,
              sample_corr = sample_corr,
              sample_diag_c=sample_diag_c,
              sample_corr_c=sample_corr_c,
              perc_psd = counter_c/ns))
}
