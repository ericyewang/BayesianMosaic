# Bayesian Mosaic Simulation Study
# 2D Gaussian with Normal Inverse Wishart (NIW) Prior
# Definition of NIW, see https://en.wikipedia.org/wiki/Normal-inverse-Wishart_distribution#Posterior_distribution_of_the_parameters

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


