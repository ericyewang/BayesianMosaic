genKnots <- function(y, X, model) {
  # generate the knots that connect titles
  # args:
  #   y: a data frame (or object coercible by as.data.frame to a data frame) containing the data for all dimensions
  #   X: a data frame (or object coercible by as.data.frame to a data frame) containing extra features
  #   model: a string indicating the model to apply Bayesian Mosaic on
  # output:
  # TODO
  
  # validation
  dim_y = dim(y)
  if (dim_y[1] < dim_y[2]) {
    stop("sample size smaller than number of dimensions")
  }
  nd = dim_y[2]
  
  # genKnot for each dimension
  for (j in 1:nd) {
    # TODO
  }
  
  # return knots
  return(list())
}

genMosaic <- function(y, X, model) {}