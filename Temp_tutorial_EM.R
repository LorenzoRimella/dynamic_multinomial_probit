rm(list=ls())
# Set the directory to the path of the current script
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(gridExtra)

# Directions
# Try to use EP to exclude the (i,t) observation this should be on the lines of using r_-i and Q_-i
# make the code usable with ADF and EP

# set.seed(5)

source("R_directory/dynamic_multinomial_probit/model.R")
source("R_directory/dynamic_multinomial_probit/multinomial_probit_sampling_filter.R")
source("R_directory/dynamic_multinomial_probit/multinomial_probit_filter_exact.R")
source("R_directory/dynamic_multinomial_probit/multinomial_probit_filter_EP.R")
source("R_directory/dynamic_multinomial_probit/multinomial_probit_gaussian_smoothing.R")
source("R_directory/dynamic_multinomial_probit/Temp_multinomial_probit_EM.R")

L = 2 # number of choices
Time= 20 # number of time steps

p = 2 # number of subject dependent covariate

N_max = 10 # maximum number of customers
lambda_M = 1

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
Sigma[cbind(1:(L - 1), 2:L)] = 0.1
Sigma = Sigma +t(Sigma)

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
df_y$Category = factor(df_y$Category)#, levels = c("0", "1", "2", "3", "4"))

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

plot_3 =ggplot(df_beta_long, aes(x = Time, y = Value, color = Column)) +
  geom_line(alpha = 0.8) +
  labs(x = "Time", y = "Beta") +
  theme(legend.position="none")

grid.arrange(plot_1, plot_3, ncol = 1)

ADF_filter_output_exact = ADF_sampling_filter(data_multinomial_probit, parameters)

ADF_filter_output_exact$log_likelihood

EP_sampling_filter_output_exact = EP_sampling_filter(data_multinomial_probit, parameters)

EP_sampling_filter_output_exact$log_likelihood

# EM algorithm
EM_iterations = 100 
learn_G = FALSE
learn_W = FALSE
learn_Sigma = TRUE

if((learn_G & learn_W) & !(learn_Sigma))
{
  initial_parameters = list(a_0 = a_0, P_0 = P_0, 
                            G = 0.1*diag(beta_size), 
                            W = 0.1*diag(dim(parameters$W)[1]), 
                            Sigma = Sigma) #diag(L))
}

if(learn_Sigma & !(learn_G) & !(learn_W))
{
  initial_parameters = list(a_0 = a_0, P_0 = P_0, 
                            G = G, 
                            W = W, 
                            Sigma = diag(L))
}

if(learn_Sigma & learn_G & learn_W)
{
  initial_parameters = list(a_0 = a_0, P_0 = P_0, 
                            G = 0.1*diag(beta_size), 
                            W = 0.1*diag(dim(parameters$W)[1]), 
                            Sigma = diag(L))
}

EM_output = multinomial_probit_EM(data_multinomial_probit, EP_sampling_filter, initial_parameters, 
                                  EM_iterations,
                                  learn_G, learn_W, learn_Sigma)

# Create a data frame for the line plot
df <- data.frame(
  iterations = 1:(EM_iterations+1),
  log_likelihood = EM_output$log_likelihood_history
)

# Create the plot
ggplot() +
  geom_line(data = df, aes(x = iterations, y = log_likelihood, color = "EM Log-Likelihood"), size = 1) + # Line plot
  geom_hline(aes(yintercept = EP_sampling_filter_output_exact$log_likelihood, color = "EP Log-Likelihood"), 
             linetype = "dashed", size = 1) + # Horizontal line
  labs(
    x = "Iterations",  # X-axis label
    y = "Log-Likelihood",  # Y-axis label
    color = "Legend"  # Legend title
  ) +
  scale_color_manual(
    values = c("EM Log-Likelihood" = "blue", "EP Log-Likelihood" = "red"), # Custom colors
    labels = c("EM Log-Likelihood" = "EM Log-Likelihood", "EP Log-Likelihood" = "EP Exact Output") # Custom labels
  ) +
  theme_minimal()  # Minimal theme

# Extract the matrix dimensions
n <- dim(EM_output$G_history)[2] # Number of rows/columns in G

# Reshape data for ggplot
df_list <- lapply(1:n, function(i) {
  lapply(1:n, function(j) {
    data.frame(
      iteration = 1:100,
      G_value = EM_output$G_history[1:100, i, j],
      i = i,
      j = j,
      true_value = parameters$G[i, j]
    )
  })
})
df <- do.call(rbind, do.call(rbind, df_list))

# Create the plot
ggplot(df, aes(x = iteration, y = G_value)) +
  geom_line(color = "blue", size = 1) +                           # Line plot for G history
  geom_hline(aes(yintercept = true_value),                       # Horizontal line for true value
             color = "red", linetype = "dashed", size = 1) +
  facet_grid(i ~ j, labeller = label_both) +                     # Subfigures for each (i, j)
  labs(
    x = "Iterations",                                            # X-axis label
    y = "G Value",                                               # Y-axis label
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8)                          # Adjust facet label size
  )


# Extract the matrix dimensions
n <- dim(EM_output$W_history)[2] # Number of rows/columns in W

# Reshape data for ggplot
df_list <- lapply(1:n, function(i) {
  lapply(1:n, function(j) {
    data.frame(
      iteration = 1:100,
      W_value = EM_output$W_history[1:100, i, j],
      i = i,
      j = j,
      true_value = parameters$W[i, j]
    )
  })
})

df <- do.call(rbind, do.call(rbind, df_list))

# Create the plot
ggplot(df, aes(x = iteration, y = W_value)) +
  geom_line(color = "blue", size = 1) +                           # Line plot for G history
  geom_hline(aes(yintercept = true_value),                       # Horizontal line for true value
             color = "red", linetype = "dashed", size = 1) +
  facet_grid(i ~ j, labeller = label_both) +                     # Subfigures for each (i, j)
  labs(
    x = "Iterations",                                            # X-axis label
    y = "W Value",                                               # Y-axis label
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8)                          # Adjust facet label size
  )



# Extract the matrix dimensions
n <- dim(EM_output$Sigma_history)[2] # Number of rows/columns in Sigma

# Reshape data for ggplot
df_list <- lapply(1:n, function(i) {
  lapply(1:n, function(j) {
    data.frame(
      iteration = 1:100,
      Sigma_value = EM_output$Sigma_history[1:100, i, j],
      i = i,
      j = j,
      true_value = parameters$Sigma[i, j]
    )
  })
})

df <- do.call(rbind, do.call(rbind, df_list))

# Create the plot
ggplot(df, aes(x = iteration, y = Sigma_value)) +
  geom_line(color = "blue", size = 1) +                           # Line plot for G history
  geom_hline(aes(yintercept = true_value),                       # Horizontal line for true value
             color = "red", linetype = "dashed", size = 1) +
  facet_grid(i ~ j, labeller = label_both) +                     # Subfigures for each (i, j)
  labs(
    x = "Iterations",                                            # X-axis label
    y = "Sigma Value",                                               # Y-axis label
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8)                          # Adjust facet label size
  )






ADF_filter_output_EM = ADF_sampling_filter(data_multinomial_probit, EM_parameters, 30)

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

