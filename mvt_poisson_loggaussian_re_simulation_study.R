# Bayesian Mosaic Simulation Study
# hierarchical Poisson log Gaussian

rm(list=ls())

# dependencies
suppressMessages(require(MCMCpack))
suppressMessages(require(mvtnorm))
suppressMessages(require(LaplacesDemon))
suppressMessages(require(Rmisc))
suppressMessages(require(ggplot2))

# helpers
ilogit <- function(x){
  # inverse logit function
  return( 1/(1+exp(-x)) )
}

aggregateData <- function(y){
  # aggregate data into list([value,count])
  # since the count data only takes a small number of values
  # args:
  #   y: count observations
  
  pool = unique(y)
  ret = matrix(0, 2, length(pool))
  for (i in 1:length(pool)){
    ret[, i] = c(pool[i], sum(y == pool[i]))
  }
  return(ret)
}

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
    ret = c(ret, getGroupLik(group_summaries[[i]], group_mus[i], v, logarithm = TRUE))
  }
  if (logarithm == FALSE){
    ret = exp(ret)
  }
  return(ret)
}

# sample from the knot posteriors
sampleKnot <- function(y, group_idx, nb, ns, njump, lam, a, b, proposal_var, verbose = FALSE) {
  # sample from the knot posterior
  # Args:
  #   y: observation
  #   group_idx: group index of each observation
  #   nb: number of burn-ins
  #   ns: number of samples to collect
  #   njump: thinning parameter
  #   lam: parameter in NIW that controls the prior concentration of mu
  #   a: shape parameter of inverse gamma
  #   b: scale parameter of inverse gamma
  #   proposal_var: initial variance of the proposal normal distribution
  #   verbose: whether or not to print intermediate sampling info
  
  # summary statistics
  K = max(group_idx)
  group_summaries = list()
  for( k in 1:K) {
    group_summaries[[k]] = aggregateData(y[group_idx == k])
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
  proposal_vars_s1 = proposal_var
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
    new_s1 = s1 + sqrt(proposal_vars_s1) * rnorm(1)
    if (new_s1 > 0) {
      new_lliks = getGroupLiks(group_summaries, group_mus, new_s1, logarithm = TRUE)
      new_lpost = sum(new_lliks) + dgamma(new_s1, 2, .5, log = TRUE)
      prv_lpost = sum(prv_lliks) + dgamma(s1, 2, .5, log = TRUE)
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
    window_memory_group_mus[1:99,] = window_memory_group_mus[2:100,]
    window_memory_group_mus[100,] = group_mus
    window_memory_s1[1:99] = window_memory_s1[2:100]
    window_memory_s1[100] = s1
    historical_var_s1 = sd(window_memory_s1, na.rm = TRUE)
    historical_var_group_mus = apply(window_memory_group_mus, 2, function(x) {
      sd(x, na.rm = TRUE)
    })
    proposal_vars_group_mus[historical_var_group_mus > 0] = historical_var_group_mus[historical_var_group_mus > 0]
    if (historical_var_s1 > 0) {
      proposal_var_s1 = historical_var_s1
    }
  }
  
  return(list(mgroup_mus = mgroup_mus, vs1 = vs1, vmu0 = vmu0, vs2 = vs2, mar = mar))
}

# sample from the tile posteriors
sampleTile <- function(ys, yt, samples_knot_s, samples_knot_t, nb, ns, njump, lam, nu, Psi, proposal_var){
  # sample from the tile posterior
  # Args:
  #   ys: dimension s of y
  #   yt: dimension t of y
  #   samples_knot_s: posterior samples from the knot posterior of dimension s
  #   samples_knot_t: posterior samples from the knot posterior of dimension t
  #   nb: number of burn-ins
  #   ns: number of samples to collect
  #   njump: thinning parameter
  #   lam: parameter in NIW that controls the prior concentration of mu
  #   nu: d.f. of inverse-wishart
  #   Psi: scale matrix of inverse-wishart
  #   proposal_var: initial variance of the proposal normal distribution
  
  # initialization
  Y = cbind(ys, yt)
  rho = 0
  prv_lpost = -Inf
  
  # summary statistics
  YY = t(Y) %*% Y
  YS = apply(Y, 2, sum)
  
  # preparation for adapting
  window_memory = c(rep(NA, 99), rho)
  
  # outputs
  ars = NULL
  rhos = NULL
  mus = NULL
  vars = NULL
  
  for (iter in 1:(nb + ns * njump)) {
    mu_tmp = c(samples_knot_s[2, iter], samples_knot_t[2, iter])
    diag_Sigma_tmp = c(samples_knot_s[1, iter], samples_knot_t[1, iter])
    
    new_rho = rho + sqrt(proposal_var) * rnorm(1)
    if (abs(new_rho) < 1){ # reject if proposal is out of (-1,1)
      new_Sigma_tmp = diag(sqrt(diag_Sigma_tmp)) %*% matrix(c(1, new_rho, new_rho, 1), 2, 2) %*% diag(sqrt(diag_Sigma_tmp))
      new_Sigma_tmp_inv = solve(new_Sigma_tmp)
      new_llik = -sum(diag(YY %*% new_Sigma_tmp_inv)) / 2 + YS %*% new_Sigma_tmp_inv %*% mu_tmp - n * t(mu_tmp) %*% new_Sigma_tmp_inv %*% mu_tmp / 2 - 
        n * as.numeric(determinant(2 * pi * new_Sigma_tmp)$modulus) / 2
      new_lprior = dmvnorm(mu_tmp, mean = rep(0, 2), sigma = new_Sigma_tmp / lam, log = TRUE) + log(diwish(new_Sigma_tmp, nu, Psi))
      new_lpost = new_llik + new_lprior
      # new_lpost = new_llik
      
      ar = min(exp(new_lpost - prv_lpost), 1) # symmetric proposal
      if (runif(1) <= ar) {
        rho = new_rho
        prv_lpost = new_lpost
      }
    } else {
      ar = 0
    }
    
    
    if (iter > nb & (iter - nb) %% njump == 0) {
      rhos = c(rhos, rho)
      mus = cbind(mus, mu_tmp)
      vars = cbind(vars, diag_Sigma_tmp)
      ars = c(ars, ar)
    }
    
    # adapting proposal sd
    window_memory[1:99] = window_memory[2:100]
    window_memory[100] = rho
    historical_var = sd(window_memory, na.rm = TRUE)
    if (historical_var > 0) {
      proposal_var = historical_var
    }
  }
  
  return(list(rhos = rhos, mus = mus, vars = vars, ars = ars))
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

# run experiments
lam = 0.1
nu = 4
Psi = diag(1,2)
experiment_results = list()
counter = 0
set.seed(2017)
for (n in 100 * 2^(0:7)) {
  for (rho in c(-0.9, -0.5, 0, 0.5, 0.9)) {
    for (j in 1:50) {
      counter = counter + 1
      cat(counter, '\n')
      
      res = experimentOnce(n, rho, lam, nu, Psi)
      
      # draw samples from the true posteriors
      ns = 1000
      Sigma_samples_true = matrix(0, 3, ns)
      mu_samples_true = matrix(0, 2, ns)
      for (i in 1:ns) {
        tmp = riwish(res$nu_f, res$Psi_f)
        mu_samples_true[, i] = rmvnorm(1, res$mu_f, tmp / res$lam_f)
        Sigma_samples_true[1:2, i] = diag(tmp)
        Sigma_samples_true[3, i] = tmp[1, 2] / sqrt(tmp[1, 1] * tmp[2, 2])
      }
      
      mu_mean_diffs = apply(res$mus, 1, mean) - apply(mu_samples_true, 1, mean)
      var_mean_diffs = c(apply(res$vars, 1, mean), mean(res$rhos)) - apply(Sigma_samples_true, 1, mean)
      mu_sd_diffs = apply(res$mus, 1, sd) - apply(mu_samples_true, 1, sd)
      var_sd_diffs = c(apply(res$vars, 1, sd), sd(res$rhos)) - apply(Sigma_samples_true, 1, sd)
      
      suppressWarnings(ks_stat_rho <- ks.test(res$rhos, Sigma_samples_true[3,])$statistic)
      
      experiment_results[[counter]] = list(
        n = n,
        rho = rho,
        mu_mean_diffs = abs(mu_mean_diffs),
        var_mean_diffs = abs(var_mean_diffs),
        mu_sd_diffs = abs(mu_sd_diffs),
        var_sd_diffs = abs(var_sd_diffs),
        ks_stat_rho = ks_stat_rho
      )
    }
  }
}

# compair actual marginal posterior density of rho and the Bayesian mosaic posterior density of rho
df_rho = data.frame(rhos = c(res$rhos, Sigma_samples_true[3,]), 
                    group = as.factor(c(rep("Bayesian Mosaic", length(res$rhos)), rep("True Posterior", ncol(Sigma_samples_true)))))
ggplot(df_rho) + geom_density(aes(rhos, fill = group), alpha = 0.3) + theme_bw() + xlab('rho')

# reshape the results
x_n = NULL
x_rho = NULL
y_rho_mean_diff = NULL
y_rho_sd_diff = NULL
y_var1_mean_diff = NULL
y_var1_sd_diff = NULL
y_var2_mean_diff = NULL
y_var2_sd_diff = NULL
for (i in 1:counter) {
  x_n = c(x_n, experiment_results[[i]]$n)
  x_rho = c(x_rho, experiment_results[[i]]$rho)
  y_rho_mean_diff = c(y_rho_mean_diff, experiment_results[[i]]$var_mean_diffs[3])
  y_rho_sd_diff = c(y_rho_sd_diff, experiment_results[[i]]$var_sd_diffs[3])
  y_var1_mean_diff = c(y_var1_mean_diff, experiment_results[[i]]$var_mean_diffs[1])
  y_var1_sd_diff = c(y_var1_sd_diff, experiment_results[[i]]$var_sd_diffs[1])
  y_var2_mean_diff = c(y_var2_mean_diff, experiment_results[[i]]$var_mean_diffs[2])
  y_var2_sd_diff = c(y_var2_sd_diff, experiment_results[[i]]$var_sd_diffs[2])
}
df2vis = data.frame(
  n = rep(x_n, 3),
  rho = as.factor(rep(x_rho, 3)),
  mean_diff = c(y_rho_mean_diff, y_var1_mean_diff, y_var2_mean_diff),
  sd_diff = c(y_rho_sd_diff, y_var1_sd_diff, y_var2_sd_diff),
  estimand = as.factor(c(rep("rho", counter), rep("var1", counter), rep("var2", counter)))
)

# extract the sd of the differences
mean_diff_summary = summarySE(df2vis, measurevar="mean_diff", groupvars=c("n", "estimand"))
sd_diff_summary = summarySE(df2vis, measurevar="sd_diff", groupvars=c("n", "estimand"))

# visualize the trend of differences versus the sample size
pd <- position_dodge(0.2) # move them .05 to the left and right
ggplot(mean_diff_summary, aes(x = log(n), y = mean_diff, color = estimand)) + 
  geom_errorbar(aes(ymin = mean_diff - sd, ymax = mean_diff + sd), width = .1, position=pd) +
  geom_line(position=pd) + geom_point(position=pd) + theme_bw() + xlab("log sample size") + ylab("difference in posterior mean")
ggplot(sd_diff_summary, aes(x = log(n), y = sd_diff, color = estimand)) + 
  geom_errorbar(aes(ymin = sd_diff - sd, ymax = sd_diff + sd), width = .1, position=pd) +
  geom_line(position=pd) + geom_point(position=pd) + theme_bw() + xlab("log sample size") + ylab("difference in posterior sd")


