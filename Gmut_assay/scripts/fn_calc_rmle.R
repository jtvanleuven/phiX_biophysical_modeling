###################################################################
# Fisher scoring iterations to compute restricted ML estimate 
# per Rubin et al. 
# (this is basically R version of their python code:
# github.com/FowlerLab/Enrich2/blob/master/enrich2/random_effects.py)

# Input: y is a numeric vector of replicate scores,
#        sigma2i is a numeric vector of std error squared
#        iteration is a numeric
###################################################################

calc_rmle <- function(y, sigma2i, iterations){
  w <- 1/sigma2i
  sw <- sum(w)
  beta0 <- sum(y*w)/sw
  sigma2ML <- sum((y - mean(y))^2/(length(y)-1))
  eps <- 0
  betaML <- 0
  for (k in 1:iterations){
    w <- 1/(sigma2i + sigma2ML)
    sw <- sum(w)
    sw2 <- sum(w^2)
    betaML <- sum(y*w)/sw
    sigma2ML_new <- sigma2ML*sum(((y-betaML)^2)*(w^2))/(sw-(sw2/sw))
    eps <- abs(sigma2ML - sigma2ML_new)
    sigma2ML <- sigma2ML_new
  }
  var_betaML <- 1/sum(1/(sigma2i+sigma2ML))
  out <- list()
  out[["betaML"]] <- betaML
  out[["var_betaML"]] <- var_betaML
  out[["eps"]] <- eps
  return(out)
}
