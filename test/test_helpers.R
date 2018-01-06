# Test Helpers

rm(list = ls())

# load helpers
source("~/Documents/yw_git/bayesian_mosaic/helpers.R")

# test Poisson----
# test genLikPoissonLogNormal
res_int = genLikPoissonLogNormal(1, -1.2, .75)
res_mc = mean(dpois(1, exp(-1.2 + sqrt(.75) * rnorm(100000))))
if (abs(res_int - res_mc) / res_int > .01) {
  cat("genLikPoissonLogNormal failed unit test!\n")
}

# test genLikPoissonLogGaussian2D
Sigma = matrix(c(1.2, .96, .96, 1.4), 2, 2)
Sinv = solve(Sigma)
res_int = genLikPoissonLogGaussian2D(c(0, 0), c(-4, -4), Sigma, Sinv)
latent_x = rmvnorm(10000, c(-4, -4), Sigma)
res_mc = mean(dpois(0, exp(latent_x[,1]))*dpois(0, exp(latent_x[,2])))
if (abs(res_int - res_mc) / res_int > .01) {
  cat("genLikPoissonLogGaussian2D failed unit test!\n")
}

# test genLikPoissonLogGaussian2DGradient
v = c(1.2, 1.4)
rho = 0.8
Sigma = diag(sqrt(v)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sqrt(v))
Sinv = solve(Sigma)
res_int = genLikPoissonLogGaussian2DGradient(c(0, 0), c(-4,-4), Sigma, Sinv)

Sigma_u = diag(sqrt(v)) %*% matrix(c(1, rho + 1E-4, rho + 1E-4, 1), 2, 2) %*% diag(sqrt(v))
Sinv_u = solve(Sigma_u)
Sigma_l = diag(sqrt(v)) %*% matrix(c(1, rho - 1E-4, rho - 1E-4, 1), 2, 2) %*% diag(sqrt(v))
Sinv_l = solve(Sigma_l)
res_test = (genLikPoissonLogGaussian2D(c(0, 0), c(-4, -4), Sigma_u, Sinv_u) - 
              genLikPoissonLogGaussian2D(c(0, 0), c(-4, -4), Sigma_l, Sinv_l)) / (2E-4)

if (abs(res_int - res_test) / res_int > .01) {
  cat("genLikPoissonLogGaussian2DGradient failed unit test!\n")
}

# test Binomial----
# test genLikBinomialLogNormal
res_int = genLikBinomialLogNormal(1, 8, -1.2, .75)
res_mc = mean(dbinom(1, 8, 1 / (1 + exp(-(-1.2 + sqrt(.75) * rnorm(100000))))) / choose(8, 1))
if (abs(res_int - res_mc) / res_int > .01) {
  cat("genLikBinomialLogNormal failed unit test!\n")
}

# test genLikBinomialLogGaussian2D
Sigma = matrix(c(1.2, .96, .96, 1.4), 2, 2)
Sinv = solve(Sigma)
res_int = genLikBinomialLogGaussian2D(c(0, 0), c(7, 8), c(-4, -4), Sigma, Sinv)
latent_x = rmvnorm(10000, c(-4, -4), Sigma)
res_mc = mean(dbinom(0, 7, 1 / (1 + exp(-latent_x[,1])))*dbinom(0, 8, 1 / (1 + exp(-latent_x[,2]))))
if (abs(res_int - res_mc) / res_int > .01) {
  cat("genLikBinomialLogGaussian2D failed unit test!\n")
}

# test genLikBinomialLogGaussian2DGradient
v = c(1.2, 1.4)
rho = 0.8
Sigma = diag(sqrt(v)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sqrt(v))
Sinv = solve(Sigma)
res_int = genLikBinomialLogGaussian2DGradient(c(0, 0), c(7, 8), c(-4, -4), Sigma, Sinv)

Sigma_u = diag(sqrt(v)) %*% matrix(c(1, rho + 1E-4, rho + 1E-4, 1), 2, 2) %*% diag(sqrt(v))
Sinv_u = solve(Sigma_u)
Sigma_l = diag(sqrt(v)) %*% matrix(c(1, rho - 1E-4, rho - 1E-4, 1), 2, 2) %*% diag(sqrt(v))
Sinv_l = solve(Sigma_l)
res_test = (genLikBinomialLogGaussian2D(c(0, 0), c(7, 8), c(-4, -4), Sigma_u, Sinv_u) - 
              genLikBinomialLogGaussian2D(c(0, 0), c(7, 8), c(-4, -4), Sigma_l, Sinv_l)) / (2E-4)

if (abs(res_int - res_test) / res_int > .01) {
  cat("genLikBinomialLogGaussian2DGradient failed unit test!\n")
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
                matrix(Corr2Vec(Corr),1,p*(p-1)/2))$preds)
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
sample_corr = matrix(rep(Corr2Vec(Corr),500),500,p*(p-1)/2,byrow=TRUE)
res = evalPredAccuracy(500, mu, S, sample_mu, sample_diag, sample_corr)
if (abs(res[2]-.95)>.05) {
  cat("evalPredAccuracy failed unit test!\n")
}