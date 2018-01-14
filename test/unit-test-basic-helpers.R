# Test Basic Helpers
# Version 1.1
# Last Updated on Jan 13, 2018

rm(list = ls())

# load helpers
source("~/Documents/yw_git/bayesian_mosaic/basic-helpers.R")

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
if (abs(mean(rTruncatedNormal(n=100000, mu=3))-3)>0.01) {
  cat("rTruncatedNormal failed unit test!\n")
}
