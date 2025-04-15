sim_name <- function(name_init="filt_smoot_lik_p_", L, Time, p, expected_m_t)
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