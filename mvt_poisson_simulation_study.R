# Bayesian Mosaic Simulation Study
# Multivariate log-Gaussian mixture of Poisson
# Version 1.2
# Last Updated on April 6, 2018
rm(list = ls())

# dependencies
suppressMessages(require(coda))
source("~/Documents/yw_git/bayesian_mosaic/bayesian-mosiac.R")
source("~/Documents/yw_git/bayesian_mosaic/other-samplers.R")
# source("/home/collabor/yw104/BayesianMosaic/bayesian-mosiac.R")
# source("/home/collabor/yw104/BayesianMosaic/other-samplers.R")

# helper
experimentOnce <- function(n, p, mu, diag, corr_mat, nb, ns, 
                           njump, seednum, parallel=FALSE) {
  # TODO
  
  set.seed(seednum)
  # 1. simulate dataset
  cat("Simulating data...\n")
  Sigma = diag(sqrt(diag)) %*% corr_mat %*% diag(sqrt(diag))
  x = rmvnorm(n, mu, Sigma)
  Y = matrix(rpois(p * n, exp(c(x))), n, p)
  
  # 2. fit Bayesian mosaic
  cat("Bayesian Mosaic...\n")
  ct_bm = proc.time()
  res = bayesianMosaic(Y, nb, ns, njump, model="mvtPoisson", 
                       verbose=0, parallel=parallel)
  ct_bm = proc.time() - ct_bm
  
  # 3. fit a DAMCMC using the same amount of time
  cat("Data-augmented MCMC...\n")
  ns_damcmc = 10000
  res_damcmc = sampleViaDAMCMC(Y, ns_damcmc, list(mu=mu,Sigma=Sigma), 
                               ct_bm[3], verbose=FALSE, parallel=parallel)
  
  # 4. gather performance info
  cat("Evaluating Performance...\n")
  n_pair = ncol(res$sample_corr)
  true_corr = lowerOffDiagonal(corr_mat)
  
  # bayesian mosaic
  perf_bm = list(perc_psd = res$perc_psd, # percentage of posterior correlation matrix that satisfies psd
                 pred_accuracy = evalPredAccuracy(500, mu, Sigma, res$sample_mu, 
                                                  res$sample_diag_c, 
                                                  res$sample_corr_c),
                 err_mu = apply(res$sample_mu,2,mean)-mu,
                 err_diag = apply(res$sample_diag,2,mean)-diag,
                 err_diag_c = apply(res$sample_diag_c,2,mean)-diag,
                 err_corr = apply(res$sample_corr,2,mean)-true_corr,
                 err_corr_c = apply(res$sample_corr_c,2,mean)-true_corr,
                 ess_mu = effectiveSize(res$sample_mu),
                 ess_diag = effectiveSize(res$sample_diag),
                 ess_corr = effectiveSize(res$sample_corr),
                 coverage_mu = sapply(1:p,FUN=function(j){
                   return(mu[j]>quantile(res$sample_mu[,j],.025)&
                            mu[j]<quantile(res$sample_mu[,j],.975))
                 }),
                 coverage_diag = sapply(1:p,FUN=function(j){
                   return(diag[j]>quantile(res$sample_diag[,j],.025)&
                            diag[j]<quantile(res$sample_diag[,j],.975))
                 }),
                 coverage_diag_c = sapply(1:p,FUN=function(j){
                   return(diag[j]>quantile(res$sample_diag_c[,j],.025)&
                            diag[j]<quantile(res$sample_diag_c[,j],.975))
                 }),
                 coverage_corr = sapply(1:n_pair,FUN=function(j){
                   return(true_corr[j]>quantile(res$sample_corr[,j],.025)&
                            true_corr[j]<quantile(res$sample_corr[,j],.975))
                 }),
                 coverage_corr_c = sapply(1:n_pair,FUN=function(j){
                   return(true_corr[j]>quantile(res$sample_corr_c[,j],.025)&
                            true_corr[j]<quantile(res$sample_corr_c[,j],.975))
                 }),
                 intlen_mu = sapply(1:p,FUN=function(j){
                   return(quantile(res$sample_mu[,j],.975)-
                            quantile(res$sample_mu[,j],.025))
                 }),
                 intlen_diag = sapply(1:p,FUN=function(j){
                   return(quantile(res$sample_diag[,j],.975)-
                            quantile(res$sample_diag[,j],.025))
                 }),
                 intlen_diag_c = sapply(1:p,FUN=function(j){
                   return(quantile(res$sample_diag_c[,j],.975)-
                            quantile(res$sample_diag_c[,j],.025))
                 }),
                 intlen_corr = sapply(1:n_pair,FUN=function(j){
                   return(quantile(res$sample_corr[,j],.975)-
                            quantile(res$sample_corr[,j],.025))
                 }),
                 intlen_corr_c = sapply(1:n_pair,FUN=function(j){
                   return(quantile(res$sample_corr_c[,j],.975)-
                            quantile(res$sample_corr_c[,j],.025))
                 }))
  
  # damcmc
  perf_damcmc = list(pred_accuracy = evalPredAccuracy(500, mu, Sigma, res_damcmc$sample_mu, 
                                                      res_damcmc$sample_diag, res_damcmc$sample_corr),
                     err_mu = apply(res_damcmc$sample_mu,2,mean)-mu,
                     err_diag = apply(res_damcmc$sample_diag,2,mean)-diag,
                     err_corr = apply(res_damcmc$sample_corr,2,mean)-true_corr,
                     ess_mu = effectiveSize(res_damcmc$sample_mu),
                     ess_diag = effectiveSize(res_damcmc$sample_diag),
                     ess_corr = effectiveSize(res_damcmc$sample_corr),
                     coverage_mu = sapply(1:p,FUN=function(j){
                       return(mu[j]>quantile(res_damcmc$sample_mu[,j],.025)&
                                mu[j]<quantile(res_damcmc$sample_mu[,j],.975))
                     }),
                     coverage_diag = sapply(1:p,FUN=function(j){
                       return(diag[j]>quantile(res_damcmc$sample_diag[,j],.025)&
                                diag[j]<quantile(res_damcmc$sample_diag[,j],.975))
                     }),
                     coverage_corr = sapply(1:n_pair,FUN=function(j){
                       return(true_corr[j]>quantile(res_damcmc$sample_corr[,j],.025)&
                                true_corr[j]<quantile(res_damcmc$sample_corr[,j],.975))
                     }),
                     intlen_mu = sapply(1:p,FUN=function(j){
                       return(quantile(res_damcmc$sample_mu[,j],.975)-
                                quantile(res_damcmc$sample_mu[,j],.025))
                     }),
                     intlen_diag = sapply(1:p,FUN=function(j){
                       return(quantile(res_damcmc$sample_diag[,j],.975)-
                                quantile(res_damcmc$sample_diag[,j],.025))
                     }),
                     intlen_corr = sapply(1:n_pair,FUN=function(j){
                       return(quantile(res_damcmc$sample_corr[,j],.975)-
                                quantile(res_damcmc$sample_corr[,j],.025))
                     }))
  
  return(list(eigenvalues_corr=eigen(corr_mat,only.values=TRUE)$values,
              perf_bm=sapply(perf_bm,FUN=function(x){as.vector(x)}), 
              perf_damcmc=sapply(perf_damcmc,FUN=function(x){as.vector(x)})))
}

# simulation study
mu0 = -4
diag0 = 0.5
p = 3
n_experiment = 2

n = 10000
nb = 10
ns = 10
njump = 1
parallel = TRUE

perfs = list()
corr_mats = sampleLKJ(1000, p) # simulate a bag of correlation matrix
corr_ids = sample(1:1000, n_experiment, replace=FALSE)
set.seed(2018)
seeds = sample(1:10000, n_experiment, replace=FALSE)
for (iter in 1:n_experiment) {
  cat("experiment ", iter, "\n")
  corr_mat = corr_mats[corr_ids[iter],,]
  mu = mu0+runif(p)
  diag = diag0+0.5*runif(p)
  seednum = seeds[iter]
  perfs[[iter]] = experimentOnce(n, p, mu, diag, corr_mat, nb, ns, 
                                 njump, seednum, parallel=parallel)
}

# # sparse Poisson counts
# n_experiment = 100
# for (iter in 1:n_experiment) {
#   
# }
# 
# # compare the traceplot of log likelihood
# nb_damcmc = 10000
# v_da = apply(res_damcmc$mv[, (nb_damcmc+1):ns_damcmc], 1, mean)
# v_bm = c(mean(knot1$vs1), mean(knot2$vs1))
# mu_da = apply(res_damcmc$mmu[, (nb_damcmc+1):ns_damcmc], 1, mean)
# mu_bm = c(mean(knot1$vmu0), mean(knot2$vmu0))
# rho_da = mean(res_damcmc$vrho[(nb_damcmc+1):ns_damcmc])
# rho_bm = mean(tile12)
# group_summary = aggregateData(y)
# getGroupLLik2D(group_summary, mu_da, v_da, rho_da)
# getGroupLLik2D(group_summary, mu_bm, v_bm, rho_bm)
# 
# # traceplots
# par(mar = c(2.5, 3, 3, 3), mfrow=c(2,1))
# ratio_ct = (ct_knot[3] / 2 + ct_tile[3]) / ct_damcmc[3]
# plot(1:500, tile12, 'l', main = "Bayesian mosaic", xlab = "", ylab = "rho")
# ns_damcmc = length(res_damcmc$vrho)
# plot(1:floor(ns_damcmc * ratio_ct), res_damcmc$vrho[1:floor(ns_damcmc * ratio_ct)], 'l', main = "DA-MCMC", xlab = "", ylab = "rho")
# par(mfrow=c(1,1))
# 
# # compare the effective sample size
# effectiveSize(tile12)
# effectiveSize(res_damcmc$vrho[1:floor(ns_damcmc * ratio_ct)])
# effectiveSize(res_damcmc$vrho)
