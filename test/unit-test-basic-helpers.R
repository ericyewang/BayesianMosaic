# Test Basic Helpers
# Version 1.1
# Last Updated on Jan 11, 2018

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
