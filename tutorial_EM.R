library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(gridExtra)

source("R_directory/model.R")
source("R_directory/multinomial_probit_sampling_filter.R")
source("R_directory/multinomial_probit_filter_exact.R")
source("R_directory/multinomial_probit_filter_EP.R")
source("R_directory/multinomial_probit_gaussian_smoothing.R")
source("R_directory/multinomial_probit_EM.R")

L = 2 # number of choices
Time= 50 # number of time steps

p = 2 # number of subject dependent covariate

N_max = 10 # maximum number of customers
lambda_M = 5
hyperparameters = list(lambda_M = lambda_M, N_max = N_max)

X_cov = matrix(runif(N_max*p, -0.5, 0.5), nrow = N_max, ncol = p)

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = 0.1*diag(beta_size)

G = 0.3*diag(beta_size)
G[2,1] = -0.3
G[1,2] = -0.3

W = diag(beta_size)/2

Sigma = matrix(0.0, ncol = L, nrow = L) + 0.1*diag(L)

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

ADF_filter_output_exact = ADF_sampling_filter(data_multinomial_probit, parameters)

ADF_filter_output_exact$log_likelihood

initial_parameters = list(a_0 = a_0, P_0 = P_0, G = 0.8*diag(beta_size), W = 2*diag(beta_size), Sigma = Sigma)
# initial_parameters = list(a_0 = a_0, P_0 = P_0, G = 0.8*diag(beta_size), W = W, Sigma = Sigma)

EM_parameters = multinomial_probit_EM(data_multinomial_probit, ADF_sampling_filter, initial_parameters, 100)

ADF_filter_output_EM = ADF_sampling_filter(data_multinomial_probit, EM_parameters)

ADF_filter_output_EM$log_likelihood

color_list <- c("blue", "black", "red", "green", "brown")

# Create a data frame for plotting
data_list <- list()
for (j in 1:4) {
  i =1
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = ADF_filter_output_exact$a_history[, 1*(j-1) + i],
    y_upper = ADF_filter_output_exact$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_filter_output_exact$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_filter_output_exact$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_filter_output_exact$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i],
    facet_row = j,
    facet_col = "ADF exact"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = ADF_filter_output_EM$a_history[, 1*(j-1) + i],
    y_upper = ADF_filter_output_EM$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_filter_output_EM$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_filter_output_EM$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_filter_output_EM$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+1],
    facet_row = j,
    facet_col = "ADF EM"
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
