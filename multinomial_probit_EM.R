# library(einsum)

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

multinomial_probit_EM <- function(data, filter, initial_parameters, EM_iterations = 100)
{
  Time = length(data$Y)
  sum_m_t = sum(sapply(data$Y, length)) #Reduce("+", lapply(seq(1, Time), function(t) length(data$Y[[t]])))
  
  parameters_current = initial_parameters
  
  for(iteration in 1:EM_iterations)
  {
    print(paste("Iteration nr ", iteration))
    filter_output = filter(data, parameters_current)
    prior_output  = gaussian_prior(data, parameters)
    smoother_output = Gaussian_smoother(data, parameters_current, filter_output, prior_output, "smoother")
    print(filter_output$log_likelihood)
    # print("Matrix G")
    # print(parameters_current$G)
    # print("Matrix W")
    # print(parameters_current$W)
    print("Matrix Sigma")
    print(parameters_current$Sigma)
    
    a_0_Tm1_T = smoother_output$a_history[1:Time,]
    a_1_T_T   = smoother_output$a_history[2:(Time+1),]
    
    P_0_Tm1_T = smoother_output$P_history[1:Time,,]
    P_1_T_T   = smoother_output$P_history[2:(Time+1),,]
    
    J_0_T_T   = smoother_output$J_history[1:Time,,]
    
    sum_barGP_a_t_tm1_T = Reduce("+", lapply(seq(1, Time), 
                                             function(t) J_0_T_T[t,,] %*% P_1_T_T[t,,] + a_1_T_T[t,] %*% t(a_0_Tm1_T[t,])))
    sum_P_a_tm1_tm1_T   = Reduce("+", lapply(seq(1, Time), 
                                             function(t) P_0_Tm1_T[t,,] + a_0_Tm1_T[t,] %*% t(a_0_Tm1_T[t,])))
    
    new_G = sum_barGP_a_t_tm1_T %*% solve(sum_P_a_tm1_tm1_T)

    new_W = (1/Time)*Reduce("+", lapply(seq(1, Time),
                                        function(t) P_1_T_T[t,,] + a_1_T_T[t,] %*% t(a_1_T_T[t,])
                                                    - J_0_T_T[t,,] %*% P_1_T_T[t,,] %*% t(new_G) - a_1_T_T[t,] %*% t(a_0_Tm1_T[t,]) %*% t(new_G)
                                                    - t(J_0_T_T[t,,] %*% P_1_T_T[t,,] %*% t(new_G) + a_1_T_T[t,] %*% t(a_0_Tm1_T[t,]) %*% t(new_G))
                                                    + new_G %*% (P_0_Tm1_T[t,,] + a_0_Tm1_T[t,] %*% t(a_0_Tm1_T[t,])) %*% t(new_G)))

    eps_SUN_sample = lapply(seq(1, Time), function(t) lapply(seq(1, length(data$Y[[t]])), function(i) rSUN(smoother_output$eps_history[[t]][[i]], 5000)))
    eps_second_moment = lapply(seq(1, Time), function(t) lapply(seq(1, length(data$Y[[t]])), function(i) cov(t(eps_SUN_sample[[t]][[i]])) + (t(t(apply(eps_SUN_sample[[t]][[i]], 1, mean))) %*% t(apply(eps_SUN_sample[[t]][[i]], 1, mean)))))
    
    # eps_second_moment = lapply(seq(1, Time), function(t) lapply(seq(1, length(data$Y[[t]])), function(i) einsum("in,jn->ij", eps_SUN_sample[[t]][[i]], eps_SUN_sample[[t]][[i]])/1000))
    
    new_Sigma = Reduce("+", lapply(seq(1, Time), function(t) Reduce("+", eps_second_moment[[t]])))/sum_m_t

    new_W = parameters_current$W
    new_G = parameters_current$G
    # new_Sigma = parameters_current$Sigma
    
    # parameters_current$G = new_G
    # parameters_current$W = (new_W + t(new_W))/2
    parameters_current$Sigma = (new_Sigma + t(new_Sigma))/2
    
  }
  
  return(parameters_current)
}

