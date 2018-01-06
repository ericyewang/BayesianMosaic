# Bayesian Mosaic Simulation Study
# Multivariate log-Gaussian mixture of Poisson
# Version 1.1
# Last Updated on Dec 19, 2017
rm(list = ls())

# dependencies
suppressMessages(require(coda))
source("~/Documents/yw_git/bayesian_mosaic/basian_mosaic_functions.R")

# helper
experimentOnce <- function(n, p, mu, diag_Sigma, corr_mat, nb, ns, 
                           njump, parallel=FALSE) {
  # TODO
  
  # 1. simulate dataset
  cat("Simulating data...\n")
  Sigma = diag(sqrt(diag_Sigma)) %*% corr_mat %*% diag(sqrt(diag_Sigma))
  x = rmvnorm(n, mu, Sigma)
  Y = matrix(rpois(p * n, exp(c(x))), n, p)
  
  # 2. fit Bayesian mosaic
  cat("Bayesian Mosaic...\n")
  ct_bm = proc.time()
  res = bayesianMosaic(Y, nb, ns, njump, proposal_var=0.01, model="mvtPoisson", 
                       verbose=2, parallel=parallel)
  ct_bm = proc.time() - ct_bm
  
  # 3. fit a DAMCMC using the same amount of time
  cat("Data-augmented MCMC...\n")
  ns_damcmc = 10000
  res_damcmc = sampleViaDAMCMC(Y, ns_damcmc, list(mu=mu,Sigma=Sigma), 
                               ct_bm[3], verbose=TRUE, parallel=parallel)
  
  # 4. gather performance info
  cat("Evaluating Performance...\n")
  n_pair = ncol(res$sample_Correlation)
  true_corr = Corr2Vec(corr_mat)
  
  # bayesian mosaic
  biw_diag_Sigma = NULL
  biw_Correlation = NULL
  for (i in 1:1000) {
    tmp = riwish(res$best_iw_v, res$best_iw_S)
    tmp_corr = diag(sqrt(1/diag(tmp)))%*%tmp%*%diag(sqrt(1/diag(tmp)))
    biw_diag_Sigma = rbind(biw_diag_Sigma, diag(tmp))
    biw_Correlation = rbind(biw_Correlation, Corr2Vec(tmp_corr))
  }
  perf_bm = list(err_mu = apply(res$sample_mu,2,mean)-mu,
                 err_diag_Sigma = apply(res$sample_diag_Sigma,2,mean)-diag_Sigma,
                 err_biw_diag_Sigma = apply(biw_diag_Sigma,2,mean)-diag_Sigma,
                 err_Corr = apply(res$sample_Correlation,2,mean)-true_corr,
                 err_biw_Corr = apply(biw_Correlation,2,mean)-true_corr,
                 ess_mu = effectiveSize(res$sample_mu),
                 ess_diag_Sigma = effectiveSize(res$sample_diag_Sigma),
                 ess_corr = effectiveSize(res$sample_Correlation),
                 coverage_mu = sapply(1:p,FUN=function(j){
                   return(mu[j]>quantile(res$sample_mu[,j],.025)&
                            mu[j]<quantile(res$sample_mu[,j],.975))
                 }),
                 coverage_diag_Sigma = sapply(1:p,FUN=function(j){
                   return(diag_Sigma[j]>quantile(res$sample_diag_Sigma[,j],.025)&
                            diag_Sigma[j]<quantile(res$sample_diag_Sigma[,j],.975))
                 }),
                 coverage_biw_diag_Sigma = sapply(1:p,FUN=function(j){
                   return(diag_Sigma[j]>quantile(biw_diag_Sigma[,j],.025)&
                            diag_Sigma[j]<quantile(biw_diag_Sigma[,j],.975))
                 }),
                 coverage_Corr = sapply(1:n_pair,FUN=function(j){
                   return(true_corr[j]>quantile(res$sample_Correlation[,j],.025)&
                            true_corr[j]<quantile(res$sample_Correlation[,j],.975))
                 }),
                 coverage_biw_Corr = sapply(1:n_pair,FUN=function(j){
                   return(true_corr[j]>quantile(biw_Correlation[,j],.025)&
                            true_corr[j]<quantile(biw_Correlation[,j],.975))
                 }),
                 intlen_mu = sapply(1:p,FUN=function(j){
                   return(quantile(res$sample_mu[,j],.975)-
                            quantile(res$sample_mu[,j],.025))
                 }),
                 intlen_diag_Sigma = sapply(1:p,FUN=function(j){
                   return(quantile(res$sample_diag_Sigma[,j],.975)-
                            quantile(res$sample_diag_Sigma[,j],.025))
                 }),
                 intlen_biw_diag_Sigma = sapply(1:p,FUN=function(j){
                   return(quantile(biw_diag_Sigma[,j],.975)-
                            quantile(biw_diag_Sigma[,j],.025))
                 }),
                 intlen_Corr = sapply(1:n_pair,FUN=function(j){
                   return(quantile(res$sample_Correlation[,j],.975)-
                            quantile(res$sample_Correlation[,j],.025))
                 }),
                 intlen_biw_Corr = sapply(1:n_pair,FUN=function(j){
                   return(quantile(biw_Correlation[,j],.975)-
                            quantile(biw_Correlation[,j],.025))
                 }))
  
  # damcmc
  perf_damcmc = list(err_mu = apply(res_damcmc$sample_mu,2,mean)-mu,
                     err_diag_Sigma = apply(res_damcmc$sample_diag_Sigma,2,mean)-diag_Sigma,
                     err_cor = apply(res_damcmc$sample_Correlation,2,mean)-true_corr,
                     ess_mu = effectiveSize(res_damcmc$sample_mu),
                     ess_diag_Sigma = effectiveSize(res_damcmc$sample_diag_Sigma),
                     ess_corr = effectiveSize(res_damcmc$sample_Correlation),
                     coverage_mu = sapply(1:p,FUN=function(j){
                       return(mu[j]>quantile(res_damcmc$sample_mu[,j],.025)&
                                mu[j]<quantile(res_damcmc$sample_mu[,j],.975))
                     }),
                     coverage_diag_Sigma = sapply(1:p,FUN=function(j){
                       return(diag_Sigma[j]>quantile(res_damcmc$sample_diag_Sigma[,j],.025)&
                                diag_Sigma[j]<quantile(res_damcmc$sample_diag_Sigma[,j],.975))
                     }),
                     coverage_Corr = sapply(1:n_pair,FUN=function(j){
                       return(true_corr[j]>quantile(res_damcmc$sample_Correlation[,j],.025)&
                                true_corr[j]<quantile(res_damcmc$sample_Correlation[,j],.975))
                     }),
                     intlen_mu = sapply(1:p,FUN=function(j){
                       return(quantile(res_damcmc$sample_mu[,j],.975)-
                                quantile(res_damcmc$sample_mu[,j],.025))
                     }),
                     intlen_diag_Sigma = sapply(1:p,FUN=function(j){
                       return(quantile(res_damcmc$sample_diag_Sigma[,j],.975)-
                                quantile(res_damcmc$sample_diag_Sigma[,j],.025))
                     }),
                     intlen_Corr = sapply(1:n_pair,FUN=function(j){
                       return(quantile(res_damcmc$sample_Correlation[,j],.975)-
                                quantile(res_damcmc$sample_Correlation[,j],.025))
                     }))
  
  return(list(perf_bm=sapply(perf_bm,FUN=function(x){as.vector(x)}), 
              perf_damcmc=sapply(perf_damcmc,FUN=function(x){as.vector(x)})))
}

# simulate a bag of correlation matrix
n=10000
p = 3
corr_mats = sampleLKJ(1000, p)
nb=10
ns=10
njump=1
parallel=TRUE
corr_mat = corr_mats[1,,]
mu = -2 + rnorm(p)
diag_Sigma = 0.5 + 0.5 * rnorm(p)^2
res = experimentOnce(n, p, mu, diag_Sigma, corr_mat, nb, ns, njump, parallel=parallel)

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
