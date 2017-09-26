# Helper----
aggregateData <- function(vy){
  # aggregate data into list([value,count])
  # since the count data only takes a small number of values
  pool = unique(vy)
  ret = list()
  for (i in 1:length(pool)){
    ret[[i]] = c(pool[i],sum(vy==pool[i]))
  }
  return(ret)
}

# Compute the log-likelihood of a Poisson log-normal----
getLik <- function(y,mu,v){
  # likelihood function value for a individual observation via numerical integration
  # transform the variable so that the range is 0 to 1, have tried -Inf to Inf 
  # and there is some bizzare numerical issue there.
  integrand <- function(x){
    return( exp(y*x-lfactorial(y)-(x-mu)^2/2/v-exp(x))*(2*pi*v)^-.5 )
  }
  integrate(integrand, lower = min(log(y+0.1),mu)-4*sqrt(v), 
            upper = max(log(y+0.1),mu)+4*sqrt(v),stop.on.error=FALSE)$value
}

# integrate the 1st layer random effects out----
getLikelihood <- function(data_compressed, mu, s, logarithm = TRUE){
  # get (the logarithm of) of the data likelihood
  ret = 0
  for (i in 1:length(dataSummary)){
    ret = ret + data_compressed[i, 2] * log(getLik(data_compressed[i, 1], mu, s))
  }
  if (logarithm == FALSE){
    ret = exp(ret)
  }
  return(ret)
}

# full conditional samplers----
sampleMujOnce <- function(data_compressed, prev_muj, prev_llik = NULL,
                          s, muj_mu, muj_s, proposal_s) {
  # TODO
  
  if (is.null(prev_llik)) { # if the previous log-likelihood is not provided
    prev_llik = getLikelihood(data_compressed, prev_muj, s)
  }
  
  new_muj = prev_muj + sqrt(proposal_s) * rnorm(1)
  cur_llik = getLikelihood(data_compressed, new_muj, s)
  
  # log acceptance ratio (symmetric proposal)
  lar = cur_llik - prev_llik - dnorm(prev_muj, muj_mu, sqrt(muj_s), log = TRUE) + 
    dnorm(new_muj, muj_mu, sqrt(muj_s), log = TRUE)
  
  
}
