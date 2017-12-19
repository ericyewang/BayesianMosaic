# Bayesian Mosaic Simulation Study 1
# Bivariate Gaussian

rm(list=ls())

# dependencies
source("~/Documents/yw_git/bayesian_mosaic/helpers.R")

# helper
sampleLaplaceApprox <- function(n, ys, yt, mu, vvars) {
  # sample from the Laplace approximation
  # args:
  #   n: number of samles
  #   ys: dimension s of y
  #   yt: dimension t of y
  #   mu: mean vector
  #   vvars: vector of variances
  
  # simulation study
  ys = (ys - mu[1]) / sqrt(vvars[1])
  yt = (yt - mu[2]) / sqrt(vvars[2])
  a = length(ys) / 2
  b = sum(ys^2 + yt^2) / 2
  c = -sum(ys * yt)
  
  llikfcn <- function(x) {
    rho = 2 / (exp(-x) + 1) - 1
    aux = 1 - rho^2
    -a * log(aux) - b / aux - c * rho / aux
  }
  
  llikGradient <- function(x) {
    rho = 2 / (exp(-x) + 1) - 1
    aux = 1 - rho^2
    llik_gradient_rho = (2 * a * rho * aux - 2 * b * rho - c - c * rho^2) / aux^2
    return(llik_gradient_rho * 2 * exp(x) / (1 + exp(x))^2)
  }
  
  # maximize the log-likelihood
  optim_res = optim(par = 0, fn = function(x){-llikfcn(x)}, gr = function(x){-llikGradient(x)}, 
                    method = "BFGS", lower = -Inf, upper = Inf, hessian = TRUE)
  
  # laplace approximation
  laplace_mean = optim_res$par
  laplace_var = 1 / c(optim_res$hessian)
  
  # sample once
  x_sample = laplace_mean + sqrt(laplace_var) * rnorm(n)
  
  return(2 / (exp(-x_sample) + 1) - 1)
}

# sample from the tile conditional
sampleTile <- function(ys, yt, samples_knot_s, samples_knot_t,
                       ns, verbose = FALSE){
  # sample from the tile conditional
  # Args:
  #   ys: y_s
  #   yt: y_t
  #   samples_knot_s: posterior samples from the knot marginal of dimension s
  #   samples_knot_t: posterior samples from the knot marginal of dimension t
  #   ns: number of samples to collect
  #   verbose: whether or not to print intermediate sampling info
  
  # initialization
  rho = 0
  
  # outputs
  rhos = NULL
  
  for (iter in 1:ns) {
    # print intermediate sampling info
    if (verbose && (iter %% floor((nb + ns * njump) / 100)) == 0) {
      cat("iteration: ", iter, "\n")
    }
    
    mu_tmp = c(samples_knot_s[2, iter], samples_knot_t[2, iter])
    diag_Sigma_tmp = c(samples_knot_s[1, iter], samples_knot_t[1, iter])
    
    rho = sampleLaplaceApprox(1, ys, yt, mu_tmp, diag_Sigma_tmp)
    
    rhos = c(rhos, rho)
  }
  
  return(rhos)
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
  ns = 1000 # number of posterior samples
  samples_knots = array(0, c(2, ns, p))
  for (s in 1:p) {
    tmp_sigma = rinvgamma(ns, nu_f / 2, Psi_f[s, s] / 2)
    samples_knots[1, , s] = tmp_sigma
    samples_knots[2, , s] = mu_f[s] + sqrt(tmp_sigma / lam_f) * rnorm(ns)
  }
  
  # tile posteriors
  tile12 = sampleTile(y[,1], y[,2], samples_knots[,,1], 
                      samples_knots[,,2], ns, FALSE)
  
  return(list(
    rhos = tile12, 
    mus = t(samples_knots[2, , s]), 
    vars = t(samples_knots[1, , s]),
    lam_f = lam_f,
    nu_f = nu_f,
    mu_f = mu_f,
    Psi_f = Psi_f
  ))
}

# run experiments
lam = 0.1
nu = 6
Psi = diag(1,2)
experiment_results = list()
counter = 0
set.seed(2017)
for (n in 100 * 2^(0:8)) {
  for (j in 1:100) {
    counter = counter + 1
    cat(counter, '\n')
    rho = runif(1, min = -1, max = 1)
    res = experimentOnce(n, rho, lam, nu, Psi)
    
    # draw samples from the true posteriors
    rho_samples_true = NULL
    for (jj in 1:1000) {
      tmp = riwish(res$nu_f, res$Psi_f)
      rho_samples_true = c(rho_samples_true, tmp[1, 2]/ sqrt(tmp[1, 1] * tmp[2, 2]))
    }
    
    experiment_results[[counter]] = list(
      n = n,
      rho_mean_diff = mean(res$rhos) - mean(rho_samples_true),
      rho_sd_diff = sd(res$rhos) - sd(rho_samples_true)
    )
  }
}

# reshape the results
x_n = NULL
y_rho_mean_diff = NULL
y_rho_sd_diff = NULL
for (i in 1:length(experiment_results)) {
  x_n = c(x_n, experiment_results[[i]]$n)
  y_rho_mean_diff = c(y_rho_mean_diff, experiment_results[[i]]$rho_mean_diff)
  y_rho_sd_diff = c(y_rho_sd_diff, experiment_results[[i]]$rho_sd_diff)
}
df2vis = data.frame(sample_size = as.factor(round(log(x_n), 1)),
                    y_rho_mean_diff = y_rho_mean_diff,
                    y_rho_sd_diff = sqrt(x_n) * y_rho_sd_diff)

# visualize the trend of differences versus the sample size
ggplot(df2vis, aes(x = sample_size, y = y_rho_mean_diff)) + geom_boxplot() + 
  stat_summary(fun.y=mean, geom="line", aes(group=1))  + 
  stat_summary(fun.y=mean, geom="point") + theme_bw() + xlab("log sample size") + ylab("difference in posterior mean")

ggplot(df2vis, aes(x = sample_size, y = y_rho_sd_diff)) + geom_boxplot() + 
  stat_summary(fun.y=mean, geom="line", aes(group=1))  + 
  stat_summary(fun.y=mean, geom="point") + theme_bw() + xlab("log sample size") + ylab("standardized difference in posterior sd")


# compare the posterior density from Bayesian mosaic to the truth
lam = 0.1
nu = 6
Psi = diag(1,2)
facets = list()
set.seed(2018)
facet_counter = 0
for (n in c(100, 1000, 10000)) {
  for (rho in c(0, 0.5, 0.9)) {
    facet_counter = facet_counter + 1
    
    res = experimentOnce(n, rho, lam, nu, Psi)
    
    n_samples = length(res$rhos)
    # draw samples from the true posteriors
    rho_samples_true = NULL
    for (jj in 1:n_samples) {
      tmp = riwish(res$nu_f, res$Psi_f)
      rho_samples_true = c(rho_samples_true, tmp[1, 2]/ sqrt(tmp[1, 1] * tmp[2, 2]))
    }
    
    group_n = c(group_n, rep(n, n_samples))
    group_rho = c(group_rho, rep(rho, n_samples))
    
    df_tmp = data.frame(rho_samples = c(res$rhos, rho_samples_true),
                        method = c(rep("Bayesian Mosaic", n_samples), rep("True Posterior", n_samples)))
    facets[[facet_counter]] = ggplot(df_tmp) + geom_density(aes(rho_samples, fill = method), alpha = 0.3) + 
      theme_bw() + xlab(paste("n=", n, ", rho=", rho, sep = "")) + theme(legend.position="none") + ylab("")
  }
}

suppressMessages(require(gridExtra))
grid.arrange(facets[[1]], facets[[2]], facets[[3]],
             facets[[4]], facets[[5]], facets[[6]],
             facets[[7]], facets[[8]], facets[[9]], nrow = 3)
