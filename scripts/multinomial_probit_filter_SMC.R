SMC_sampling_filter <- function(data, parameters, MC_samples)
{
  start_time <- Sys.time()
  L = dim(parameters$Sigma)[1]
  Time = length(data$Y)
  
  a_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1]))
  P_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1], dim(parameters$a_0)[1]))
  
  sampling_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1], MC_samples))
  
  a_history[1,]  = parameters$a_0
  P_history[1,,] = parameters$P_0
  log_likelihood = 0
  
  sampling_history[1,,] = t(mvrnorm(n = MC_samples, parameters$a_0, parameters$P_0))
  
  for(t in 1:Time)
  {
    sampling_prediction = parameters$G %*% sampling_history[t,,] + t(mvrnorm(n = MC_samples, 0*parameters$a_0, parameters$W))
    
    y_t = data$Y[[t]]
    x_t = data$X[[t]]
    
    m_t = dim(x_t)[1]

    transformed_X = lapply(seq_len(m_t), function(i) {
      X_minus_l(x_t[i,,], y_t[i], L)
    })
    transformed_Sigma = lapply(seq_len(m_t), function(i) {
      V_i = V_minus_l(y_t[i], L)
      V_i %*% parameters$Sigma %*% t(V_i)
    })

    particle_log_weight <- function(s) {
      sum(vapply(seq_len(m_t), function(i) {
        probability = Phi(
          transformed_X[[i]] %*% sampling_prediction[, s, drop = FALSE],
          transformed_Sigma[[i]]
        )
        log(pmax(probability, .Machine$double.xmin))
      }, numeric(1)))
    }

    log_weights = vapply(seq_len(MC_samples), particle_log_weight, numeric(1))
    max_log_weight = max(log_weights)
    if (!is.finite(max_log_weight)) {
      stop("All SMC particle weights are numerically zero at time ", t, ".")
    }
    scaled_weights = exp(log_weights - max_log_weight)
    weights_norm = scaled_weights / sum(scaled_weights)

    log_likelihood = log_likelihood +
      max_log_weight + log(mean(scaled_weights))

    ancestors = sample.int(MC_samples, MC_samples, replace = TRUE, prob = weights_norm)
    
    sampling_history[t+1,,] = sampling_prediction[,ancestors]
    
    a_t = t(t(apply(sampling_history[t+1,,], 1, mean)))
    P_t = cov(t(sampling_history[t+1,,]))
    
    a_history[t+1,]  = a_t
    P_history[t+1,,] = P_t
    
  }
  
  return(list(log_likelihood = log_likelihood, 
              a_history = a_history, P_history = P_history,
              sampling_history = sampling_history[2:(Time+1),,, drop = FALSE],
              running_time = as.numeric(Sys.time() - start_time, units = "secs")))
}
