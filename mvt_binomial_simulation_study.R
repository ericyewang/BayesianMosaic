# Bayesian Mosaic Simulation Study 2
# Bivariate log-normal mixture of Binomial
# Version 1.2
# Last Updated on May 20, 2018

rm(list = ls())

# change the dir to yours
source("~/Documents/yw_git/bayesian_mosaic/bayesian_mosiac.R")
source("~/Documents/yw_git/bayesian_mosaic/other_samplers.R")
suppressMessages(require(coda))

# simulation
# set.seed(2018)
n = 10000
p = 2 # without loss of generosity, we focus on 2D case
mu = c(-4, -4)
# diag_Sigma = 0.5+0.5 * rnorm(2)^2
diag_Sigma = c(1, 1)
rho = 0.9
Sigma = diag(sqrt(diag_Sigma)) %*% matrix(c(1, rho, rho, 1), 2, 2) %*% diag(sqrt(diag_Sigma))
x = rmvnorm(n, mu, Sigma)
Y = matrix(rbinom(2 * n, 10, 1 / (1 + exp(-c(x)))), n, 2)

# fit Bayesian mosaic
nb = 1000
ns = 100
njump = 2
proposal_var = 0.01
model = "mvtBinomial"
ct_bm = proc.time()
res = bayesianMosaic(Y, nb, ns, njump, proposal_var, model, verbose = TRUE, Ns = c(10, 10))
ct_bm = proc.time() - ct_bm

ct_damcmc = proc.time()
ns_damcmc = 60000
res_damcmc = sampleViaDAMCMC(y, ns_damcmc, verbose = TRUE)
ct_damcmc = proc.time() - ct_damcmc

# compare the traceplot of log likelihood
nb_damcmc = 10000
v_da = apply(res_damcmc$mv[, (nb_damcmc+1):ns_damcmc], 1, mean)
v_bm = c(mean(knot1$vs1), mean(knot2$vs1))
mu_da = apply(res_damcmc$mmu[, (nb_damcmc+1):ns_damcmc], 1, mean)
mu_bm = c(mean(knot1$vmu0), mean(knot2$vmu0))
rho_da = mean(res_damcmc$vrho[(nb_damcmc+1):ns_damcmc])
rho_bm = mean(tile12)
group_summary = aggregateData(y)
getGroupLLik2D(group_summary, mu_da, v_da, rho_da)
getGroupLLik2D(group_summary, mu_bm, v_bm, rho_bm)

# traceplots
par(mar = c(2.5, 3, 3, 3), mfrow=c(2,1))
ratio_ct = (ct_knot[3] / 2 + ct_tile[3]) / ct_damcmc[3]
plot(1:500, tile12, 'l', main = "Bayesian mosaic", xlab = "", ylab = "rho")
ns_damcmc = length(res_damcmc$vrho)
plot(1:floor(ns_damcmc * ratio_ct), res_damcmc$vrho[1:floor(ns_damcmc * ratio_ct)], 'l', main = "DA-MCMC", xlab = "", ylab = "rho")
par(mfrow=c(1,1))

# compare the effective sample size
effectiveSize(tile12)
effectiveSize(res_damcmc$vrho[1:floor(ns_damcmc * ratio_ct)])
effectiveSize(res_damcmc$vrho)
