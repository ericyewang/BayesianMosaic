# Test Helpers for Samplers
# Version 1.2
# Last Updated on May 20, 2018

rm(list = ls())

# change the dir to yours
source("~/Documents/yw_git/bayesian_mosaic/sampler_helpers.R")

# test Poisson----
# test genLLikPoissonLogNormal
res_int = genLLikPoissonLogNormal(1, -1.2, .75)
res_mc = log(mean(dpois(1, exp(-1.2 + sqrt(.75) * rnorm(100000)))))
if (abs(res_int - res_mc) / res_int > .01) {
  cat("genLLikPoissonLogNormal failed unit test!\n")
}

# test genLLikPoissonLogGaussian2D
Sigma = matrix(c(1.2, .96, .96, 1.4), 2, 2)
upper_tri = chol(Sigma)
res_int = genLLikPoissonLogGaussian2D(c(0, 0), c(-4, -4), upper_tri)
latent_x = rmvnorm(10000, c(-4, -4), Sigma)
res_mc = log(mean(dpois(0, exp(latent_x[,1]))*dpois(0, exp(latent_x[,2]))))
if (abs(res_int - res_mc) / res_int > .01) {
  cat("genLLikPoissonLogGaussian2D failed unit test!\n")
}

# test genLikPoissonLogGaussian2DGradient
v = c(1.2, 1.4)
rho = 0.8
Sigma = diag(sqrt(v)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sqrt(v))
upper_tri = chol(Sigma)
res_int = genLikPoissonLogGaussian2DGradient(c(0, 0), c(-4,-4), diag(Sigma), 
                                             upper_tri, rho)

Sigma_u = diag(sqrt(v))%*%matrix(c(1, rho+1E-4, rho+1E-4, 1), 2, 2)%*%diag(sqrt(v))
upper_tri_u = chol(Sigma_u)
Sigma_l = diag(sqrt(v))%*%matrix(c(1, rho-1E-4, rho-1E-4, 1), 2, 2)%*%diag(sqrt(v))
upper_tri_l = chol(Sigma_l)
res_test = (exp(genLLikPoissonLogGaussian2D(c(0, 0), c(-4, -4), upper_tri_u))-
              exp(genLLikPoissonLogGaussian2D(c(0, 0), c(-4, -4), upper_tri_l)))/(2E-4)

if (abs(res_int - res_test) / res_int > .01) {
  cat("genLikPoissonLogGaussian2DGradient failed unit test!\n")
}

# test Binomial----
# test genLLikBinomialLogNormal
res_int = genLLikBinomialLogNormal(1, 8, -1.2, .75)
res_mc = log(mean(dbinom(1, 8, 1/(1+exp(-(-1.2+sqrt(.75)*rnorm(100000)))))))
if (abs(res_int - res_mc) / res_int > .01) {
  cat("genLLikBinomialLogNormal failed unit test!\n")
}

# test genLLikBinomialLogGaussian2D
Sigma = matrix(c(1.2, .96, .96, 1.4), 2, 2)
upper_tri = chol(Sigma)
res_int = genLLikBinomialLogGaussian2D(c(0, 0), c(7, 8), c(-4, -4), upper_tri)
latent_x = rmvnorm(10000, c(-4, -4), Sigma)
res_mc = log(mean(dbinom(0, 7, 1 / (1 + exp(-latent_x[,1])))*dbinom(0, 8, 1 / (1 + exp(-latent_x[,2])))))
if (abs(res_int - res_mc) / res_int > .01) {
  cat("genLLikBinomialLogGaussian2D failed unit test!\n")
}

# test genLikBinomialLogGaussian2DGradient
v = c(1.2, 1.4)
rho = 0.8
Sigma = diag(sqrt(v)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sqrt(v))
upper_tri = chol(Sigma)
res_int = genLikBinomialLogGaussian2DGradient(c(0, 0), c(7, 8), c(-4, -4), 
                                              diag(Sigma), upper_tri, rho)

Sigma_u = diag(sqrt(v))%*%matrix(c(1, rho+1E-4, rho+1E-4, 1), 2, 2)%*%diag(sqrt(v))
upper_tri_u = chol(Sigma_u)
Sigma_l = diag(sqrt(v))%*%matrix(c(1, rho-1E-4, rho-1E-4, 1), 2, 2)%*%diag(sqrt(v))
upper_tri_l = chol(Sigma_l)
res_test = (exp(genLLikBinomialLogGaussian2D(c(0, 0), c(7, 8), c(-4, -4), upper_tri_u))-
              exp(genLLikBinomialLogGaussian2D(c(0, 0), c(7, 8), c(-4, -4), upper_tri_l)))/(2E-4)

if (abs(res_int - res_test) / res_int > .01) {
  cat("genLikBinomialLogGaussian2DGradient failed unit test!\n")
}

# test Rounded Gaussian----
# test genLLikRoundedNormal
res_int = genLLikRoundedNormal(0, -1.2, .75)
res_true = pnorm(0, -1.2, sqrt(0.75), log=TRUE)
if (res_int!=res_true) {
  cat("genLLikPoissonLogNormal failed unit test!\n")
}
res_int = genLLikRoundedNormal(4, -1.2, .75)
res_true = log(pnorm(4, -1.2, sqrt(0.75))-pnorm(3, -1.2, sqrt(0.75)))
if (abs(res_int-res_true)>1E-8) {
  cat("genLLikPoissonLogNormal failed unit test!\n")
}

# test genLLikRoundedGaussian2D
Sigma = matrix(c(1.2, .96, .96, 1.4), 2, 2)
res_int = genLLikRoundedGaussian2D(c(0, 0), c(-4, -4), diag(Sigma), 
                                  Sigma[1,2]/sqrt(Sigma[1,1])/sqrt(Sigma[2,2]))
latent_x = rmvnorm(100000, c(-4, -4), Sigma)
res_mc = log(mean(latent_x[,1]<0&latent_x[,2]<0))
if (abs(res_int - res_mc) / res_int > .01) {
  cat("genLLikRoundedGaussian2D failed unit test!\n")
}

# test genLikRoundedGaussian2DGradient
v = c(1.2, 1.4)
rho = 0.8
Sigma = diag(sqrt(v)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sqrt(v))
res_int = genLikRoundedGaussian2DGradient(c(0, 0), c(-4,-4), diag(Sigma), rho)

Sigma_u = diag(sqrt(v))%*%matrix(c(1, rho+1E-8, rho+1E-8, 1), 2, 2)%*%diag(sqrt(v))
Sigma_l = diag(sqrt(v))%*%matrix(c(1, rho-1E-8, rho-1E-8, 1), 2, 2)%*%diag(sqrt(v))
res_test = (exp(genLLikRoundedGaussian2D(c(0, 0), c(-4, -4), diag(Sigma), rho+1E-8))-
              exp(genLLikRoundedGaussian2D(c(0, 0), c(-4, -4), diag(Sigma), rho-1E-8)))/(2E-8)

if (res_int!=res_test) {
  cat("genLikRoundedGaussian2DGradient failed unit test!\n")
}

# test imputeLatentX----
p = 5
mu = 5*rnorm(p)
S = 5*matrix(rnorm(p,1),p,1)
S = S%*%t(S)+diag(.1,p)
Corr = diag(sqrt(1/diag(S)))%*%S%*%diag(sqrt(1/diag(S)))
sample_x1 = NULL
for (i in 1:10000) {
  imcomplete_x = rmvnorm(1,mu,S)
  imcomplete_x[1] = NA
  sample_x1 = c(sample_x1, imputeLatentX(imcomplete_x, matrix(mu,1,p), matrix(diag(S),1,p), 
                matrix(lowerOffDiagonal(Corr),1,p*(p-1)/2))$preds)
}
if (abs((mean(sample_x1)-mu[1])/mu[1])>.05 | abs((var(sample_x1)-S[1,1])/S[1,1])>.05) {
  cat("imputeLatentX failed unit test!\n")
}

# test evalPredAccuracy----
p = 5
mu = 5*rnorm(p)
S = 5*matrix(rnorm(p,1),p,1)
S = S%*%t(S)+diag(.1,p)
Corr = diag(sqrt(1/diag(S)))%*%S%*%diag(sqrt(1/diag(S)))

sample_mu = matrix(rep(mu,500),500,p,byrow=TRUE)
sample_diag = matrix(rep(diag(S),500),500,p,byrow=TRUE)
sample_corr = matrix(rep(lowerOffDiagonal(Corr),500),500,p*(p-1)/2,byrow=TRUE)
res = evalPredAccuracy(500, mu, S, sample_mu, sample_diag, sample_corr)
if (abs(res[2]-.95)>.05) {
  cat("evalPredAccuracy failed unit test!\n")
}
