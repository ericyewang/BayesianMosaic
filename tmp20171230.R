n = 10000
p = 3
corr_mat = corr_mats[1,,]
nb = 5000
ns = 1000
njump = 5
res = experimentOnce(n, p, corr_mat, nb, ns, njump, TRUE)

apply(res$bm_fit$sample_mu, 1, FUN=function(x){
  c(mean(x), quantile(x, c(.025,.975)))
})
res$mu_true

apply(res$bm_fit$sample_diag_Sigma, 1, FUN=function(x){
  c(mean(x), quantile(x, c(.025,.975)))
})
diag(res$Sigma_true)

apply(res$bm_fit$sample_Correlation, c(1,2), FUN=mean)
apply(res$bm_fit$sample_Correlation, c(1,2), FUN=function(x){
  quantile(x, .975)
})
apply(res$bm_fit$sample_Correlation, c(1,2), FUN=function(x){
  quantile(x, .025)
})
diag(sqrt(1/diag(res$Sigma_true)))%*%res$Sigma_true%*%diag(sqrt(1/diag(res$Sigma_true)))

plot(1:3000, res$damcmc_fit$sample_mu[1,],'l')
plot(1:3000, res$damcmc_fit$sample_Sigma[1,],'l')
apply(res$damcmc_fit$sample_mu[,1501:3000],1,mean)
matrix(apply(res$damcmc_fit$sample_Sigma[,1501:3000],1,mean),p,p)

