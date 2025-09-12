gaussian_prior <- function(data, parameters)
{
  L = dim(parameters$Sigma)[1]
  Time = length(data$Y)
  
  a_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1]))
  P_history = array(NA, dim = c(Time+1, dim(parameters$a_0)[1], dim(parameters$a_0)[1]))
  
  a_history[1,]  = parameters$a_0
  P_history[1,,] = parameters$P_0
  
  for(t in 1:Time)
  {
    #print(t)
    ADF_prediction_parameters = ADF_prediction(a_history[t,], P_history[t,,], parameters)
    a_t_tm1 = ADF_prediction_parameters$a_t_tm1
    P_t_tm1 = ADF_prediction_parameters$P_t_tm1
    
    
    a_history[t+1,]  = a_t_tm1
    P_history[t+1,,] = P_t_tm1
    
  }
  
  return(list(a_history = a_history, P_history = P_history))
}


multinomial_probit_EM_invSigma <- function(data, 
                                           filter, type, 
                                           initial_parameters, target_trace, 
                                           EM_iterations, eps_sample_size,
                                           log = list(FALSE, c(NA, NA)))
{
  Time = length(data$Y)
  sum_m_t = sum(sapply(data$Y, length)) #Reduce("+", lapply(seq(1, Time), function(t) length(data$Y[[t]])))
  
  parameters_current = initial_parameters
  
  log_likelihood_history = array(NA, dim = c(EM_iterations+1))

  Sigma_history = array(NA, dim = c(EM_iterations+1, dim(parameters_current$Sigma)[1], dim(parameters_current$Sigma)[1]))
  
  for(iteration in 1:EM_iterations)
  {
    print("##########################################################")
    print(paste("Iteration nr ", iteration))
    print("Matrix Sigma")
    print(parameters_current$Sigma)
    Sigma_history[iteration,,] = parameters_current$Sigma
    
    filter_output = filter(data, parameters_current)
    smoother_output = Gaussian_smoother(data, parameters_current, filter_output)
    log_likelihood_history[iteration] = filter_output$log_likelihood
    print(paste("The current log-likelihood is ", filter_output$log_likelihood))

    eps_history = epsilon_smoother(data, parameters_current, filter_output, smoother_output, type)
      
    eps_SUN_sample = lapply(seq(1, Time), function(t) lapply(seq(1, length(data$Y[[t]])), function(i) rSUN(eps_history[[t]][[i]], eps_sample_size)))
      
    eps_second_moment = lapply(seq(1, Time),
                                 function(t) lapply(seq(1, length(data$Y[[t]])),
                                                    function(i) cov(t(eps_SUN_sample[[t]][[i]])) + (t(t(apply(eps_SUN_sample[[t]][[i]], 1, mean))) %*% t(apply(eps_SUN_sample[[t]][[i]], 1, mean)))))

    V = Reduce("+", lapply(seq(1, Time), function(t) Reduce("+", eps_second_moment[[t]])))

    # Function that returns A as a function of c
    A_func <- function(c) {
      return((1/sum_m_t) * (V - 2 * c * diag(dim(V)[1])))
    }
    
    # Function to compute trace of inverse of A(c)
    trace_inv_A <- function(c) {
      A <- A_func(c)
      
      if (det(A) <= 1e-10) {
        return(Inf)  # Avoid singular or near-singular matrices
      }
      
      A_inv <- solve(A)
      return(sum(diag(A_inv)))
    }
    
    # Root-finding function: solve trace_inv_A(c) = 5
    root_function <- function(c) {
      return(trace_inv_A(c) - target_trace)
    }
    
    upper_lim = min(eigen(V, only.values = T)$values)/2
    
    result <- uniroot(root_function, lower = -10e+10, upper = upper_lim)
    
    lambda_opt <- result$root
    
    print(paste("Found root with lambda ", lambda_opt))

    new_Sigma = (1/sum_m_t) * (V - 2 * lambda_opt * diag(1,L,L))
    
    parameters_current$Sigma = (new_Sigma + t(new_Sigma))/2
    
    if(log[[1]])
    {
      log_likelihood_history_subset = log_likelihood_history[1:iteration]
      Sigma_history_subset = Sigma_history[1:iteration, , ]
      
      EM_subset_output = list( log_likelihood_history = log_likelihood_history_subset, 
                               Sigma_history = Sigma_history_subset)
      
      save(EM_subset_output, file = paste(paste(log[[2]][1], sep = "", log[[2]][2]), sep = "", "_subset_EM.Rdata"))
      
    }
    
  }
  
  print("##########################################################")
  print("END")
  print("Matrix Sigma")
  print(parameters_current$Sigma)
  Sigma_history[iteration+1,,] = parameters_current$Sigma
  
  filter_output = filter(data, parameters_current)
  smoother_output = Gaussian_smoother(data, parameters_current, filter_output)
  log_likelihood_history[iteration+1] = filter_output$log_likelihood
  print(filter_output$log_likelihood)
  
  return(list(log_likelihood_history = log_likelihood_history, 
              Sigma_history = Sigma_history))
}
