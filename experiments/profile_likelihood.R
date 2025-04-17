rm(list=ls())
# Set the directory to the path of the current script
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# check the values of P_0 etc. to make sure that there are not deterministic cases
# repeat the EM algorithm and report the covergence traceplots 

source("../scripts/utils.R")
source("../scripts/model.R")
source("../scripts/multinomial_probit_filter_EP.R")
source("../scripts/multinomial_probit_filter_SMC.R")

source("../scripts/multinomial_probit_filter_ADF.R")
source("../scripts/multinomial_probit_filter_exact.R")
source("../scripts/multinomial_probit_gaussian_smoothing.R")

library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(transport)
library(expm) # for sqrtm (matrix square root)
library(future.apply)

library(dplyr)
library(tidyr)


# SET UP THE SIMULATION
L = 3 # number of choices
Time= 10 # number of time steps

p_X = 2 # number of subject dependent covariates
p_Y = 3 # number of choice dependent covariates
p = p_X*L + p_Y

expected_m_t = 15

name_init="filt_smoot_lik_p_"

output_name = sim_name(name_init, L, Time, p, expected_m_t)
name_simulation = output_name$name
directory_name  = paste("../data/", output_name$name_for_dir, sep ="")

# Create the directory if needed
if (!dir.exists(directory_name)) {
  dir.create(directory_name)
}

set.seed(678)

# GENERATE PARAMETERS

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = 0.1*diag(beta_size)

W = 0.05*diag(beta_size)

G = 0.5*diag(beta_size)
for(i in 2:beta_size)
{
  sign = sample(c(-1, 1), 1)
  G[i,i-1] = sign*0.25
  G[i-1,i] = sign*0.25
}

Sigma = 0.2*diag(L)

# GENERATE DATA

parameters = list(a_0 = a_0, P_0 = L*P_0/sum(diag(Sigma)), G = G, W = L*W/sum(diag(Sigma)), Sigma = L*Sigma/sum(diag(Sigma)))

hyperparameters = list(lambda_M = expected_m_t, N_max = 4*expected_m_t, Y_max = 20, prediction_bool= FALSE)

X_cov = matrix(runif(hyperparameters$N_max*p_X, -0.5, 0.5), nrow = hyperparameters$N_max, ncol = p_X)
Y_cov = matrix(runif(hyperparameters$Y_max*p_Y, -0.5, 0.5), nrow = hyperparameters$Y_max, ncol = p_Y)

data_multinomial_probit = simulator(X_cov, Y_cov, parameters, hyperparameters, Time)

# Check mean variance of the fixed effect
Sigma_cumul = parameters$P_0
sum_var     = array(0, dim= c(L,L))
elem = 0 
for(t in 1:Time)
{
  Sigma_cumul = parameters$G%*%Sigma_cumul%*%t(parameters$G) + parameters$W
  
  elem = elem + dim(data_multinomial_probit$X[[t]])[1]
  
  list_variances = lapply(seq(1, dim(data_multinomial_probit$X[[t]])[1]), function(i) t(data_multinomial_probit$X[[t]][i,,])%*%Sigma_cumul%*%data_multinomial_probit$X[[t]][i,,])

  sum_var = sum_var + Reduce("+", list_variances)
}

df <- expand.grid(x = 1:nrow(sum_var/elem), y = 1:ncol(sum_var/elem))
df$z <- as.vector(sum_var/elem)  # Assign values from the matrix

ggplot(df, aes(x, y, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c() + 
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal()

# PLOT DATA
# Convert Beta matrix to long format
df_beta      = as.data.frame(data_multinomial_probit$Beta)
df_beta$Time = seq_len(Time+1)  

df_beta_long = df_beta %>%
  gather(Column, Value, -Time) %>%
  mutate(Type = "Beta")  

time_steps <- seq_along(data_multinomial_probit$Y)  
data_list <- lapply(seq_along(data_multinomial_probit$Y), 
                    function(i) {data.frame(table(data_multinomial_probit$Y[[i]])) %>%
                        mutate(time = i)})

df_y <- do.call(rbind, data_list)
names(df_y) <- c("Category", "Count", "Time")  
df_y$Time <- as.factor(df_y$Time)
df_y$Category <- factor(df_y$Category, levels = as.character(sort(as.numeric(levels(df_y$Category)))))

plot_1 = ggplot(df_y, aes(x = Time, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +  
  scale_x_discrete(breaks = c(1, seq(0, 100, by = 10)[2:11])) +
  labs(x = "Time", y = "Count", fill = "Category") +
  theme_minimal() 

plot_2 =ggplot(df_beta_long, aes(x = Time, y = Value, color = Type)) +
  geom_line(alpha = 0.8) +
  scale_color_manual(values = c("blue", "red")) +
  labs(x = "Time", y = "") +
  scale_x_discrete(breaks = c(1, seq(0, 100, by = 10)[2:11])) +
  theme_minimal() 

grid.arrange(plot_1, plot_2, ncol = 1)


MC_samples = 5000

Sigma_comp <- function(a,b)
{
  Sigma = diag(L)
  for(i in 2:L)
  {
    Sigma[i,i-1] = b
    Sigma[i-1,i] = b
  }
  Sigma[1,L] = a
  Sigma[L,1] = a
  
  return(Sigma)
}

Sigma_comp("a", "b")

grid_size = 20

a_grid = seq(-1, 1, length = grid_size) 
b_grid = seq(-1, 1, length = grid_size) 

parameters_temp = parameters

log_likelihood_surface = matrix(NA, nrow=grid_size, ncol=grid_size)

for(i in 1:grid_size)
{
  for(j in 1:grid_size)
  {
    
    parameters_temp$Sigma = Sigma_comp(a_grid[i], b_grid[j])
    print(c(i,j))
    print(parameters_temp$Sigma)
    if(is_positive_definite(parameters_temp$Sigma))
    {
      print("Run ADF")
      ADF_filter_sampler_output   = ADF_sampling_filter(data_multinomial_probit, parameters_temp, MC_samples)
      log_likelihood_surface[i,j] = ADF_filter_sampler_output$log_likelihood
    }
    
  }
}

save(log_likelihood_surface, file = "log_likelihood_surface_test.Rdata")

df <- expand.grid(a = a_grid, b = b_grid)
df$log_likelihood <- as.vector(log_likelihood_surface)

ggplot(df, aes(x = a, y = b, fill = log_likelihood)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Log-Likelihood Surface", x = "a", y = "b", fill = "Log-Likelihood")








# COMPUTE FILTERS

MC_samples = 1000

sliding_window = 10

# SAME DATA but different replicates

n_replicates = 5
data_list <- list()
i = 0

for(repl in 1:n_replicates)
{
  print(repl)
  print("EP")
  EP_filter_sampler_output   = EP_sampling_filter(data_multinomial_probit, parameters, EP_max_iterations, MC_samples, tol)
  EP_smoother_sampler_output = Gaussian_smoother(data_multinomial_probit, parameters, EP_filter_sampler_output)
  
  print("ADF")
  ADF_filter_sampler_output   = ADF_sampling_filter(data_multinomial_probit, parameters, MC_samples)
  ADF_smoother_sampler_output = Gaussian_smoother(data_multinomial_probit, parameters, ADF_filter_sampler_output)
  
  print("exact")
  exact_filter_sampler_output   = exact_sampling_filter(data_multinomial_probit, parameters, sliding_window, MC_samples)
  exact_smoother_sampler_output = Gaussian_smoother(data_multinomial_probit, parameters, exact_filter_sampler_output)
  
  print("SMC")
  SMC_filter_sampler_output   = SMC_sampling_filter(data_multinomial_probit, parameters, MC_samples)
  SMC_smoother_sampler_output = Gaussian_smoother(data_multinomial_probit, parameters, SMC_filter_sampler_output)
  
  df_EP = data.frame(type ="EP", 
                     filter_running_time   = EP_filter_sampler_output$running_time,
                     smoother_running_time = EP_smoother_sampler_output$running_time,
                     log_likelihood        = EP_filter_sampler_output$log_likelihood,
                     filter_distance_from_exact = mean_wasserstein_sampling(EP_filter_sampler_output$sampling_history,
                                                                            exact_filter_sampler_output$sampling_history),
                     smoothing_distance_from_exact = mean_wasserstein_gaussian(EP_smoother_sampler_output, 
                                                                               exact_smoother_sampler_output))
  
  df_ADF = data.frame(type ="ADF", 
                      filter_running_time   = ADF_filter_sampler_output$running_time,
                      smoother_running_time = ADF_smoother_sampler_output$running_time,
                      log_likelihood        = ADF_filter_sampler_output$log_likelihood,
                      filter_distance_from_exact = mean_wasserstein_sampling(ADF_filter_sampler_output$sampling_history,
                                                                             exact_filter_sampler_output$sampling_history),
                      smoothing_distance_from_exact = mean_wasserstein_gaussian(ADF_smoother_sampler_output, 
                                                                                exact_smoother_sampler_output))
  
  df_SMC = data.frame(type ="SMC", 
                      filter_running_time   = SMC_filter_sampler_output$running_time,
                      smoother_running_time = SMC_smoother_sampler_output$running_time,
                      log_likelihood        = SMC_filter_sampler_output$log_likelihood,
                      filter_distance_from_exact = mean_wasserstein_sampling(SMC_filter_sampler_output$sampling_history,
                                                                             exact_filter_sampler_output$sampling_history),
                      smoothing_distance_from_exact = mean_wasserstein_gaussian(SMC_smoother_sampler_output, 
                                                                                exact_smoother_sampler_output))
  
  df_exact = data.frame(type ="exact", 
                        filter_running_time   = exact_filter_sampler_output$running_time,
                        smoother_running_time = exact_smoother_sampler_output$running_time,
                        log_likelihood        = exact_filter_sampler_output$log_likelihood,
                        filter_distance_from_exact = 0,
                        smoothing_distance_from_exact = 0)
  
  data_list[[i+1]] = (df_EP)
  data_list[[i+2]] = (df_ADF)
  data_list[[i+3]] = (df_SMC)
  data_list[[i+4]] = (df_exact)
  i = i+4 
  
  data_frame_output = bind_rows(data_list)
  
  save(data_frame_output, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_data.Rdata"))
  
}

# Compute mean and standard deviation for each type
summary_table <- data_frame_output %>%
  group_by(type) %>%
  summarise(across(
    where(is.numeric),
    list(mean = ~mean(.), sd = ~sd(.)),
    .names = "{.col}_{.fn}"
  )) %>%
  mutate(across(
    ends_with("_mean"),
    ~paste0(round(.x, 3), " ± ", round(get(sub("_mean$", "_sd", cur_column())), 3)),
    .names = "{.col}_formatted"
  )) %>%
  select(type, ends_with("_formatted")) %>%
  rename_with(~sub("_mean_formatted$", "", .x))

# Print table
print(summary_table, row.names = FALSE)

save(data_multinomial_probit, file    = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_data_multinomial_probit.Rdata"))
save(EP_filter_sampler_output, file    = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_EP_filter.Rdata"))
save(ADF_filter_sampler_output, file   = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_ADF_filter.Rdata"))
save(exact_filter_sampler_output, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_exact_filter.Rdata"))

color_list <- c("blue", "black", "red", "green", "brown")

# PLOT FILTERING
# Create a data frame for plotting
data_list <- list()
for (j in 1:4) {
  i =1
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = EP_filter_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = EP_filter_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(EP_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = EP_filter_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(EP_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i],
    facet_row = j,
    facet_col = "EP"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = ADF_filter_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = ADF_filter_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_filter_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+1],
    facet_row = j,
    facet_col = "ADF"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = exact_filter_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = exact_filter_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(exact_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = exact_filter_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(exact_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+1],
    facet_row = j,
    facet_col = "Exact"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
}

plot_data <- bind_rows(data_list)

# Create the ggplot
ggplot(plot_data, aes(x = x)) +
  geom_line(aes(y = y1, color = color), alpha = 0.5) +
  geom_line(aes(y = y2, color = color), alpha = 1.0) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = color), alpha = 0.5) +
  facet_grid(facet_row ~ facet_col) +
  scale_color_identity() +
  scale_fill_identity() +
  theme_minimal() +
  theme( panel.spacing = unit(1, "lines"))+
  labs(x = "Time", y = "Beta")    

# PLOT SMOOTHING
# Create a data frame for plotting
data_list <- list()
for (j in 1:4) {
  i =1
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = EP_smoother_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = EP_smoother_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(EP_smoother_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = EP_smoother_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(EP_smoother_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i],
    facet_row = j,
    facet_col = "EP"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = ADF_smoother_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = ADF_smoother_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_smoother_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_smoother_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_smoother_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+1],
    facet_row = j,
    facet_col = "ADF"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = exact_smoother_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = exact_smoother_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(exact_smoother_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = exact_smoother_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(exact_smoother_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+1],
    facet_row = j,
    facet_col = "Exact"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
}

plot_data <- bind_rows(data_list)

# Create the ggplot
ggplot(plot_data, aes(x = x)) +
  geom_line(aes(y = y1, color = color), alpha = 0.5) +
  geom_line(aes(y = y2, color = color), alpha = 1.0) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = color), alpha = 0.5) +
  facet_grid(facet_row ~ facet_col) +
  scale_color_identity() +
  scale_fill_identity() +
  theme_minimal() +
  theme( panel.spacing = unit(1, "lines"))+
  labs(x = "Time", y = "Beta")    
