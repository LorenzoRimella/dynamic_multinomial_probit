
multinomial_probit_EM <- function(data, filter, initial_parameters, EM_iterations = 100)
{
  Time = length(data$Y)
  
  parameters_current = initial_parameters
  
  for(iteration in 1:EM_iterations)
  {
    print(iteration)
    filter_output = filter(data, parameters_current)
    smoother_output = Gaussian_smoother(parameters_current, filter_output)
    print(filter_output$log_likelihood)
    print(parameters_current$G)
    print(parameters_current$W)
    
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

    # new_W = parameters_current$W
    
    parameters_current$G = new_G
    parameters_current$W = (new_W + t(new_W))/2
  }
  
  return(parameters_current)
}

