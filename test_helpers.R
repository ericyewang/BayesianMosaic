# Test Helpers

rm(list = ls())

# load helpers
source("~/Documents/yw_git/bayesian_mosaic/helpers.R")

# test binomial
# test genLikBinomialLogNormal
res_int = genLikBinomialLogNormal(1, 8, -1.2, .75)
res_mc = mean(dbinom(1, 8, 1 / (1 + exp(-(-1.2 + sqrt(.75) * rnorm(100000))))) / choose(8, 1))
if (abs(res_int - res_mc) / res_int > .01) {
  cat("genLikBinomialLogNormal failed unit test!\n")
}

# test genLikBinomialLogGuassian2D
Sigma = matrix(c(1.2, .96, .96, 1.4), 2, 2)
Sinv = solve(Sigma)
res_int = genLikBinomialLogGuassian2D(c(0, 0), c(7, 8), c(-4, -4), Sigma, Sinv)
latent_x = rmvnorm(10000, c(-4, -4), Sigma)
res_mc = mean(dbinom(0, 7, 1 / (1 + exp(-latent_x[,1]))) / choose(7, 0) * 
                dbinom(0, 8, 1 / (1 + exp(-latent_x[,2]))) / choose(8, 0))
if (abs(res_int - res_mc) / res_int > .01) {
  cat("genLikBinomialLogGuassian2D failed unit test!\n")
}

# test genLikBinomialLogGuassian2DGradient
v = c(1.2, 1.4)
rho = 0.8
Sigma = diag(sqrt(v)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sqrt(v))
Sinv = solve(Sigma)
res_int = genLikBinomialLogGuassian2DGradient(c(0, 0), c(7, 8), c(-4, -4), Sigma, Sinv)

Sigma_u = diag(sqrt(v)) %*% matrix(c(1, rho + 1E-4, rho + 1E-4, 1), 2, 2) %*% diag(sqrt(v))
Sinv_u = solve(Sigma_u)
Sigma_l = diag(sqrt(v)) %*% matrix(c(1, rho - 1E-4, rho - 1E-4, 1), 2, 2) %*% diag(sqrt(v))
Sinv_l = solve(Sigma_l)
res_test = (genLikBinomialLogGuassian2D(c(0, 0), c(7, 8), c(-4, -4), Sigma_u, Sinv_u) - 
              genLikBinomialLogGuassian2D(c(0, 0), c(7, 8), c(-4, -4), Sigma_l, Sinv_l)) / (2E-4)

if (abs(res_int - res_test) / res_int > .01) {
  cat("genLikBinomialLogGuassian2DGradient failed unit test!\n")
}
