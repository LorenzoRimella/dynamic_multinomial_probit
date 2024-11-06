library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(gridExtra)

source("R_directory/model.R")
source("R_directory/multinomial_probit_sampling_filter.R")

# COMMENTS:
# - maybe consider the covariates of the expedia dataset directly
# - plot the quantiles
# - compute the predictive either using SUN conjugacy or Monte Carlo
# - closing condition: max |r_i - r_i^new| (or any matrix)

# Next questions:
# - prediction and log-likelihood in EP?


L = 2 # number of choices
Time= 50 # number of time steps

p = 2 # number of subject dependent covariate

N_max = 10 # maximum number of customers
lambda_M = 3
hyperparameters = list(lambda_M = lambda_M, N_max = N_max)

X_cov = matrix(runif(N_max*p, -0.5, 0.5), nrow = N_max, ncol = p)

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = 0.1*diag(beta_size)

G = 0.3*diag(beta_size)
G[2,1] = -0.3
G[1,2] = -0.3

W = diag(beta_size)

Sigma = matrix(0.1, ncol = L, nrow = L) + 0.3*diag(L)

parameters = list(a_0 = a_0, P_0 = P_0, G = G, W = W, Sigma = Sigma)

data_multinomial_probit = simulator(X_cov, parameters, hyperparameters, Time)

# Convert Beta matrix to long format
df_beta      = as.data.frame(data_multinomial_probit$Beta)
df_beta$Time = seq_len(Time+1)  # Add index column for x-axis

df_beta_long = df_beta %>%
  gather(Column, Value, -Time) %>%
  mutate(Type = "Beta")  # Add a type column to distinguish Beta

time_steps <- seq_along(data_multinomial_probit$Y)  # Time steps (1, 2, 3, ...)
data_list <- lapply(seq_along(data_multinomial_probit$Y), 
                    function(i) {data.frame(table(data_multinomial_probit$Y[[i]])) %>%
                        mutate(time = i)})

# Combine all time steps into one data frame
df_y <- do.call(rbind, data_list)
names(df_y) <- c("Category", "Count", "Time")  # Rename columns
df_y$Time <- as.factor(df_y$Time) 
df_y$Category = factor(df_y$Category, levels = c("0", "1", "2", "3", "4"))

# Plotting the distribution over time using stacked bars
plot_1 = ggplot(df_y, aes(x = Time, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  scale_x_discrete(breaks = c(1, seq(0, 100, by = 10)[2:11])) +
  labs(x = "Time", y = "Count", fill = "Category") +
  theme_minimal() 
# theme(
#   strip.placement = "outside",
#   strip.text.y.left = element_text(angle = 0),  # Rotate strip text
#   legend.position = "none"
# )

# Plot using ggplot2 with facet_wrap
plot_2 =ggplot(df_beta_long, aes(x = Time, y = Value, color = Type)) +
  geom_line(alpha = 0.8) +
  scale_color_manual(values = c("blue", "red")) +  # Assign colors to Beta and Y
  #facet_wrap(~ Type, scales = "free_y", strip.position = "left", ncol = 1) +
  labs(x = "Time", y = "") +
  scale_x_discrete(breaks = c(1, seq(0, 100, by = 10)[2:11])) +
  theme_minimal() #+
# theme(
#   strip.placement = "outside",
#   legend.position = "none",
#   axis.title.y = element_blank()  # Remove the default y-axis title
# )

grid.arrange(plot_1, plot_2, ncol = 1)

ADF_filter_sampler_output = ADF_sampling_filter(data_multinomial_probit, parameters)

source("R_directory/multinomial_probit_filter_exact.R")

filter_exact_output = exact_filter(data_multinomial_probit, parameters, 15)

source("R_directory/multinomial_probit_filter_EP.R")

EP_filter_output = EP_sampling_filter(data_multinomial_probit, parameters, EP_iterations =1)

color_list <- c("blue", "black", "red", "green", "brown")

# Create a data frame for plotting
data_list <- list()
for (j in 1:4) {
  i =1
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = ADF_filter_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = ADF_filter_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_filter_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_filter_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i],
    facet_row = j,
    facet_col = "ADF"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = EP_filter_output$a_history[, 1*(j-1) + i],
    y_upper = EP_filter_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(EP_filter_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = EP_filter_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(EP_filter_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+1],
    facet_row = j,
    facet_col = "EP"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = filter_exact_output$a_history[, 1*(j-1) + i],
    y_upper = filter_exact_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(filter_exact_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = filter_exact_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(filter_exact_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
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


# SMOOTHING

source("R_directory/multinomial_probit_gaussian_smoothing.R")

ADF_smoother_sampler_output = Gaussian_smoother(parameters, ADF_filter_sampler_output)
smoother_exact_output       = Gaussian_smoother(parameters, filter_exact_output)
EP_smoother_output          = Gaussian_smoother(parameters, EP_filter_output)

color_list <- c("blue", "black", "red", "green", "brown")

# Create a data frame for plotting
data_smooth_list <- list()
for (j in 1:4) {
  i =1
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y3 = ADF_filter_sampler_output$a_history[, 1*(j-1) + i],
    y2 = ADF_smoother_sampler_output$a_history[, 1*(j-1) + i],
    y_upper = ADF_smoother_sampler_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_smoother_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_smoother_sampler_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_smoother_sampler_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i],
    facet_row = j,
    facet_col = "ADF"
  )
  data_smooth_list[[length(data_smooth_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y3 = EP_filter_output$a_history[, 1*(j-1) + i],
    y2 = EP_smoother_output$a_history[, 1*(j-1) + i],
    y_upper = EP_smoother_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(EP_smoother_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = EP_smoother_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(EP_smoother_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+1],
    facet_row = j,
    facet_col = "EP"
  )
  data_smooth_list[[length(data_smooth_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y3 = filter_exact_output$a_history[, 1*(j-1) + i],
    y2 = smoother_exact_output$a_history[, 1*(j-1) + i],
    y_upper = smoother_exact_output$a_history[, 1*(j-1) + i] + 1.96 * sqrt(smoother_exact_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = smoother_exact_output$a_history[, 1*(j-1) + i] - 1.96 * sqrt(smoother_exact_output$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+1],
    facet_row = j,
    facet_col = "Exact"
  )
  data_smooth_list[[length(data_smooth_list) + 1]] <- df_filter
  
}

plot_smooth_data <- bind_rows(data_smooth_list)

# Create the ggplot
ggplot(plot_smooth_data, aes(x = x)) +
  geom_line(aes(y = y1, color = color), alpha = 0.5) +
  geom_line(aes(y = y3, color = color), linetype = "dashed", alpha = 1.0) +
  geom_line(aes(y = y2, color = color), alpha = 1.0) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = color), alpha = 0.5) +
  facet_grid(facet_row ~ facet_col) +
  scale_color_identity() +
  scale_fill_identity() +
  theme_minimal() +
  theme( panel.spacing = unit(1, "lines"))+
  labs(x = "Time", y = "Beta")    
