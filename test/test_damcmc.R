# test sampleViaDAMCMC
rm(list = ls())

# dependencies
source("~/Documents/yw_git/bayesian_mosaic/helpers.R")

# simulation test
n = 10000
p = 3
mu = 1 + rnorm(p)
Sigma = riwish(6,diag(1,p))
x = rmvnorm(n, mu, Sigma)
Y = matrix(rpois(p * n, exp(c(x))), n, p)

ns_damcmc = 10000
res_damcmc = sampleViaDAMCMC(Y, ns_damcmc, verbose=TRUE, parallel=TRUE)

# how good the mixing is?
par(mfrow=c(2,1))
plot(1:ns_damcmc, res_damcmc$sample_mu[1,],'l')
plot(1:ns_damcmc, res_damcmc$sample_Sigma[1,],'l')
dev.off()

# how accurate is posterior certering at?
apply(res_damcmc$sample_mu, 1, mean)
mu

matrix(apply(res_damcmc$sample_Sigma, 1, mean), p, p)
Sigma

# conclusion: the posterior is centered at the correct values, hence the code
# should be bug free!
