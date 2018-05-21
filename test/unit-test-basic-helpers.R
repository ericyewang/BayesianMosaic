# Test Basic Helpers
# Version 1.2
# Last Updated on May 20, 2018

rm(list = ls())

# change the dir to yours
source("~/Documents/yw_git/bayesian_mosaic/basic_helpers.R")
suppressMessages(require(rbenchmark))

# test compressCount
if (!all(c(compressCount(c(rep(0,10), rep(1,2))))==c(10,0,2,1))) {
  cat("compressCount failed unit test!\n")
}

# test lowerOffDiagonal
if (!all(lowerOffDiagonal(matrix(c(1,.9,.9,1),2,2))==0.9)) {
  cat("lowerOffDiagonal failed unit test!\n")
}

# test vec2Corr
if (!all(c(vec2Corr(0.9,2))==c(1,0.9,0.9,1))) {
  cat("vec2Corr failed unit test!\n")
}

# test forcePSD
res = forcePSD(matrix(c(1,.9,.9,1),2,2))
if (res$is_psd & !all(abs(c(res$C)-c(1,.9,.9,1))<1E-8)) {
  cat("forcePSD failed unit test!\n")
}
res = forcePSD(matrix(c(1,1.2,1.2,1),2,2))
if (!res$is_psd & !all(abs(c(res$C)-c(1.1,1.1,1.1,1.1))<1E-8)) {
  cat("forcePSD failed unit test!\n")
}

# test Cov2CorrMat
if (!all(abs(c(Cov2CorrMat(matrix(c(2,1.8,1.8,2),2,2)))-c(1,0.9,0.9,1))<1E-8)) {
  cat("Cov2CorrMat failed unit test!\n")
}

# test pbivnormBM
if (!all(pbivnormBM(matrix(c(Inf,Inf,-Inf,-Inf,Inf,0,0,0),4,2,byrow=TRUE),0)==
  c(1,0,0.5,0.25))) {
  cat("pbivnormBM failed unit test!\n")
}

# test rTruncatedNormal
if (abs(mean(rTruncatedNormal(100000, 0, 6, 3, 1))-3)>0.01) {
  cat("rTruncatedNormal failed unit test!\n")
}
if (rTruncatedNormal(1, -Inf, 0, 20, 0.1)!=0) {
  cat("rTruncatedNormal failed unit test!\n")
}
if (rTruncatedNormal(1, 0, Inf, -20, 0.1)!=0) {
  cat("rTruncatedNormal failed unit test!\n")
}

# test quad2dLog
integrand <- function(x,y,c) {
  c*x^2*y^4
}
logIntegrand <- function(x,y,c) {
  log(c)+2*log(abs(x))+4*log(abs(y))
}
if(!abs(quad2dLog(logIntegrand=logIntegrand, n=64, xa=-2, xb=2, ya=-3, yb=4, c=3)-
  log(quad2d(f=integrand, n=64, xa=-2, xb=2, ya=-3, yb=4, c=3)))<1E-8) {
  cat("quad2Log failed unit test!\n")
}

# benchmark("quad2dLog" = {
#   quad2dLog(logIntegrand=logIntegrand, n=64, xa=-2, xb=2, ya=-3, yb=4, c=3)
# },
# "log(quad2d)" = {
#   log(quad2d(f=integrand, n=64, xa=-2, xb=2, ya=-3, yb=4, c=3))
# },
# replications = 5000,
# columns = c("test", "replications", "elapsed",
#             "relative", "user.self", "sys.self"))
# test replications elapsed relative user.self sys.self
# 2 log(quad2d)         5000  14.481    1.000    14.137    0.320
# 1   quad2dLog         5000  18.978    1.311    18.540    0.392
