library(MASS)
library(abind)
library(magic)
library(matrixcalc)
library(transport)
library(expm)

one_hot <- function(position, dimension)
{
  vector = matrix(0, nrow = dimension, ncol = 1)
  
  if(position<=dimension){
    vector[position] = 1 
  }
  
  return(vector)
}

Phi<- function(mu, Sigma)
{
  Sigma = (Sigma+t(Sigma))/2
  
  return(tlrmvnmvt::pmvn(upper = mu[,1], sigma = Sigma)[1])
}

# Phi_SMC <- function(mu, Sigma,
#                 algorithm = tlrmvnmvt::GenzBretz()) {
#   
#   mu <- as.matrix(mu)
#   Sigma <- as.matrix(Sigma)
#   Sigma <- (Sigma + t(Sigma)) / 2
#   
#   p <- ncol(mu)
#   
#   if (!all(dim(Sigma) == c(p, p))) {
#     stop("Sigma must be a p x p matrix, where p = ncol(mu).")
#   }
#   
#   vapply(
#     seq_len(nrow(mu)),
#     FUN.VALUE = numeric(1),
#     FUN = function(i) {
#       as.numeric(
#         tlrmvnmvt::pmvn(
#           upper = mu[i, ],       # row i, not mu[, 1]
#           mean = numeric(p),
#           sigma = Sigma,
#           algorithm = algorithm
#         )
#       )
#     }
#   )
# }

# tlrmvnmvt::pmvn uses 9,980 samples (20 batches x 499)
Phi_batch <- function(mu, Sigma, M = 9980, seed = NULL) {
  
  mu <- as.matrix(mu)
  Sigma <- as.matrix(Sigma)
  Sigma <- (Sigma + t(Sigma)) / 2
  
  N <- nrow(mu)
  p <- ncol(mu)
  
  if (!all(dim(Sigma) == c(p, p))) {
    stop("Sigma must be a p x p matrix, where p = ncol(mu).")
  }
  
  # chol() returns an upper-triangular factor R:
  # t(R) %*% R = Sigma
  # lpmvnorm wants the corresponding lower-triangular factor.
  Cmat <- t(chol(Sigma))
  
  C <- mvtnorm::ltMatrices(
    Cmat[lower.tri(Cmat, diag = TRUE)],
    diag = TRUE
  )
  
  lower <- matrix(-Inf, nrow = p, ncol = N)
  upper <- t(mu)   # lpmvnorm uses one observation per column
  
  log_probabilities <- mvtnorm::lpmvnorm(
    lower = lower,
    upper = upper,
    mean = 0,
    chol = C,
    logLik = FALSE,
    M = M,
    seed = seed
  )
  
  exp(log_probabilities)
}

V_minus_l <- function(l, L)
{
  if(L>1){
    if(l!=0)
    {
      return(t(cbind(- one_hot(l,L), sapply(seq(1, L)[-l], function(i) one_hot(i,L) - one_hot(l,L)))))
    }
    
    else
    {
      return(diag(L))
    }
  }
  else{
    if(l!=0)
    {
      return(-diag(L))
    }
    
    else
    {
      return(diag(L))
    }
  }
}

X_minus_l <- function(x, l, L)
{
  if(L>1){
    if(l!=0)
    {
      return(t(cbind(x[,l], sapply(seq(1, L)[-l], function(i) x[,l] - x[,i]))))
    }
    
    else
    {
      return(t(sapply(seq(1, L), function(i) -x[,i])))
    }
  }
  else{
    if(l!=0)
    {
      return((t(x)))
    }
    
    else
    {
      return((t(-x)))
    }
  }
  
}

is_positive_definite <- function(matrix) {
  is.positive.definite(matrix, tol=1e-8)
}

sim_name <- function(name_init, L, Time, p, expected_m_t)
{
  string = paste(name_init, p, sep = "")
  string = paste(string, "_mt_", sep = "")
  string = paste(string, expected_m_t, sep = "")
  string = paste(string, "_L_", sep = "")
  string = paste(string, L, sep = "")
  string = paste(string, "_Time_", sep = "")
  string = paste(string, Time, sep = "")
  
  dir_string = paste(string, "/", sep ="")
  
  return(list(name = string, name_for_dir = dir_string))
}

wasserstein_distance_gaussian <- function(mu1, mu2, Sigma1, Sigma2) {
  term1 <- sum((mu1 - mu2)^2)
  sqrt_Sigma1 <- sqrtm(Sigma1)
  term2 <- sum(diag(Sigma1)) + sum(diag(Sigma2)) - 2*sum(diag(sqrtm(sqrt_Sigma1 %*% Sigma2 %*% sqrt_Sigma1)))
  return(sqrt(term1 + term2))
}

mean_wasserstein_gaussian <- function(output_1, output_2)
{
  Time = dim(output_1$a_history)[1]
  return(mean(sapply(seq(1, Time), function(t) wasserstein_distance_gaussian(output_1$a_history[t,], output_2$a_history[t,], 
                                                                             output_1$P_history[t,,], output_2$P_history[t,,]))))
}

mean_wasserstein_sampling <- function(sampling_1, sampling_2)
{
  Time = dim(sampling_1)[1]
  p = dim(sampling_1)[2]
  
  return(mean(sapply(seq(1, Time), function(t) mean(sapply(seq(1, p), 
                                                           function(i) wasserstein1d(sampling_1[t,i,], sampling_2[t,i,], p = 2))))))
}
