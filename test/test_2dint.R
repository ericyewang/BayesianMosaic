# test sampleVia2DInt
rm(list = ls())

# dependencies
source("~/Documents/yw_git/bayesian_mosaic/helpers.R")

# simulation test
n = 10000
p = 2
mu = c(-3,-3)
diag_Sigma = c(1, 1)
rho = 0.9
Sigma = diag(sqrt(diag_Sigma))%*%matrix(c(1,rho,rho,1),2,2)%*%diag(sqrt(diag_Sigma))
x = rmvnorm(n, mu, Sigma)
Y = matrix(rpois(p * n, exp(c(x))), n, p)

nb = 5000
ns = 5000
njump = 1
proposal_var = 0.05^2
res = sampleVia2DInt(y=Y, nb=nb, ns=ns, njump=njump, proposal_var=proposal_var,
               likelihood="PoissonLogGaussian2D", ada_prop=0.5, verbose=TRUE)

# compressed_y = compressCount(Y)
# 
# tmp = NULL
# for (rho in c(-.9999,-.999,-.99,-.98,-.9,-.8,.5,0,.5,.8,.9,.98,.99,.999,.9999)) {
#   new_Sigma = diag(sqrt(diag(Sigma)))%*%matrix(c(1,rho,rho,1),2,2)%*%
#     diag(sqrt(diag(Sigma)))
#   new_Sinv = solve(new_Sigma)
#   tmp = c(tmp, genLogLikelihood(compressed_y=compressed_y, 
#                                 likelihood="PoissonLogGaussian2D", 
#                                 mu=mu, log=TRUE, Sigma=new_Sigma, Sinv=new_Sinv))
# }
# data.frame(rho=c(-.9999,-.999,-.99,-.98,-.9,-.8,.5,0,.5,.8,.9,.98,.99,.999,.9999),
#            loglik=tmp)
# plot(c(-.9999,-.999,-.99,-.98,-.9,-.8,.5,0,.5,.8,.9,.98,.99,.999,.9999), tmp)
