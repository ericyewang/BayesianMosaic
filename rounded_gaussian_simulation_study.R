# Bayesian Mosaic Simulation Study
# Rounded Gaussian
# Version 1.2
# Last Updated on May 20, 2018

rm(list = ls())

# change the dir to yours
source("~/Documents/yw_git/bayesian_mosaic/bayesian_mosiac.R")
source("~/Documents/yw_git/bayesian_mosaic/other_samplers.R")
suppressMessages(require(coda))

# helper
experimentOnce <- function(n, p, mu, diag, corr_mat, nb, ns, 
                           njump, parallel=FALSE) {
  # TODO
  
  # 1. simulate dataset
  cat("Simulating data...\n")
  Sigma = diag(sqrt(diag)) %*% corr_mat %*% diag(sqrt(diag))
  x = rmvnorm(n, mu, Sigma)
  Y = apply(x,2,FUN=function(x) {
    (x<0)*0+(x>0)*ceiling(x)
  })
  
  # 2. fit Bayesian mosaic
  cat("Bayesian Mosaic...\n")
  ct_bm = proc.time()
  res = bayesianMosaic(Y, nb, ns, njump, proposal_var=0.01, 
                       model="roundedGaussian", verbose=0, parallel=parallel)
  ct_bm = proc.time() - ct_bm
  
  # 3. fit a DAMCMC using the same amount of time
  cat("Data-augmented MCMC...\n")
  ns_damcmc = 10000
  res_damcmc = damcmcRoundedGaussian(Y, ns_damcmc, list(mu=mu,Sigma=Sigma), 
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
mu0 = 4
diag0 = 1
p = 4
n_experiment = 100

n = 10000
nb = 50000
ns = 500
njump = 30
parallel = TRUE

perfs = list()
corr_mats = sampleLKJ(1000, p) # simulate a bag of correlation matrix
corr_ids = sample(1:1000, n_experiment, replace=FALSE)
for (iter in 1:n_experiment) {
  cat("experiment ", iter, "\n")
  corr_mat = corr_mats[corr_ids[iter],,]
  mu = mu0+runif(p)
  diag = diag0+0.5*runif(p)
  perfs[[iter]] = experimentOnce(n, p, mu, diag, corr_mat, nb, ns, 
                                 njump, parallel=parallel)
}
