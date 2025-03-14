source("R_directory/dynamic_multinomial_probit/scripts/model.R")
source("R_directory/dynamic_multinomial_probit/scripts/multinomial_probit_filter_EP.R")
source("R_directory/dynamic_multinomial_probit/scripts/multinomial_probit_filter_SMC.R")
source("R_directory/dynamic_multinomial_probit/scripts/multinomial_probit_filter_ADF.R")
source("R_directory/dynamic_multinomial_probit/scripts/multinomial_probit_filter_exact.R")
source("R_directory/dynamic_multinomial_probit/scripts/multinomial_probit_gaussian_smoothing.R")

library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(transport)
library(expm) # for sqrtm (matrix square root)
library(future.apply)

name_simulation = "small_dim_small_data_small_T"
directory_name  = "R_directory/dynamic_multinomial_probit/data/small_dim_small_data_small_T/"

# Create the directory if needed
if (!dir.exists(directory_name)) {
  dir.create(directory_name)
}

set.seed(123)

# GENERATE PARAMETERS
# check the loglikelihood for exact

L = 3 # number of choices
Time= 10 # number of time steps

p = 5 # number of subject dependent covariate

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = 0.1*diag(beta_size)

G = 0.5*diag(beta_size)
for(i in 2:beta_size)
{
  sign = sample(c(-1, 1), 1)
  G[i,i-1] = sign*0.25
  G[i-1,i] = sign*0.25
}

# Plot the asymptotic G
G_cumul = G
for(t in 1:Time){G_cumul = G%*%G_cumul}

df <- expand.grid(x = 1:nrow(G_cumul), y = 1:ncol(G_cumul))
df$z <- as.vector(G_cumul)  # Assign values from the matrix

ggplot(df, aes(x, y, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c() + 
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal()

W = 0.1*diag(beta_size)

Sigma = 0.5*diag(L)
for(i in 2:L)
{
  sign = sample(c(-1, 1), 1)
  Sigma[i,i-1] = sign*0.2
  Sigma[i-1,i] = sign*0.2
}

# GENERATE DATA

parameters = list(a_0 = a_0, P_0 = P_0, G = G, W = W, Sigma = Sigma)

hyperparameters = list(lambda_M = 5, N_max = 20, prediction_bool= FALSE)

X_cov = matrix(runif(hyperparameters$N_max*p, -0.5, 0.5), nrow = hyperparameters$N_max, ncol = p)

data_multinomial_probit = simulator(X_cov, parameters, hyperparameters, Time)

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

# COMPUTE FILTERS
EP_iterations = 5
MC_samples = 1000

sliding_window = 11

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

n_replicates = 100
data_list <- list()


EP_filter_sampler_output   = EP_sampling_filter(data_multinomial_probit, parameters)
EP_smoother_sampler_output = Gaussian_smoother(data_multinomial_probit, parameters, EP_filter_sampler_output)

ADF_filter_sampler_output   = ADF_sampling_filter(data_multinomial_probit, parameters, MC_samples)
ADF_smoother_sampler_output = Gaussian_smoother(data_multinomial_probit, parameters, ADF_filter_sampler_output)

exact_filter_sampler_output   = exact_sampling_filter(data_multinomial_probit, parameters, sliding_window, MC_samples)
exact_smoother_sampler_output = Gaussian_smoother(data_multinomial_probit, parameters, exact_filter_sampler_output)

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




save(data_multinomial_probit, file    = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_data_multinomial_probit.Rdata"))
save(EP_filter_sampler_output, file    = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_EP_filter.Rdata"))
save(ADF_filter_sampler_output, file   = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_ADF_filter.Rdata"))
save(exact_filter_sampler_output, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_exact_filter.Rdata"))

EP_filter_sampler_output$log_likelihood
ADF_filter_sampler_output$log_likelihood
exact_filter_sampler_output$log_likelihood

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
    facet_col = "Sliding window"
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
