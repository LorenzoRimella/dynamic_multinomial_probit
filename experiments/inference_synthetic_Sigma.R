rm(list=ls())
# Set the directory to the path of the current script
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../scripts/utils.R")
source("../scripts/model.R")
source("../scripts/multinomial_probit_filter_EP.R")
source("../scripts/multinomial_probit_filter_SMC.R")
source("../scripts/multinomial_probit_filter_ADF.R")
source("../scripts/multinomial_probit_filter_exact.R")
source("../scripts/multinomial_probit_gaussian_smoothing.R")
source("../scripts/multinomial_probit_EM_Sigma.R")

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
Time= 5 # number of time steps

p_X = 2 # number of subject dependent covariates
p_Y = 1 # number of choice dependent covariates
p = p_X*L + p_Y

expected_m_t = 20

name_init="inference_synthetic_ADF_trace_p_"

output_name = sim_name(name_init, L, Time, p, expected_m_t)
name_simulation = output_name$name
directory_name  = paste("../data/", output_name$name_for_dir, sep ="")

# Create the directory if needed
if (!dir.exists(directory_name)) {
  dir.create(directory_name)
}

set.seed(345)

# GENERATE PARAMETERS

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = 0.5*diag(beta_size)

G = 0.1*diag(beta_size)

W = 0.05*diag(beta_size)

for (i in 1:100){
  P_0 = G%*%P_0%*%t(G) + W
}

Sigma = 0.2*diag(L)
Sigma[1,1] = 0.5
for(i in 2:L)
{
  sign = sample(c(-1, 1), 1)
  Sigma[i,i-1] = sign*0.1
  Sigma[i-1,i] = sign*0.1
}

target_trace = 5*L

invSigma = solve(Sigma)
invSigma = target_trace*invSigma/(sum(diag((invSigma))))
Sigma_update = solve(invSigma)

parameters = list(a_0 = a_0, P_0 = P_0, G = G, W = W, Sigma = Sigma_update)

# GENERATE DATA
# Trace=1

hyperparameters = list(lambda_M = expected_m_t, N_max = 4*expected_m_t, Y_max = 20, prediction_bool= FALSE)

X_cov = matrix(runif(hyperparameters$N_max*p_X, -0.5, 0.5), nrow = hyperparameters$N_max, ncol = p_X)
Y_cov = matrix(runif(hyperparameters$Y_max*p_Y, -0.5, 0.5), nrow = hyperparameters$Y_max, ncol = p_Y)

data_multinomial_probit = simulator(X_cov, Y_cov, parameters, hyperparameters, Time)

# Check mean variance of the fixed effect
Sigma_cumul = parameters$P_0
sum_var     = array(0, dim= c(L,L))
elem = 0 
for(t in 1:(Time))
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

# True likelihood

ADF_filter_output_truth   = ADF_sampling_filter(data_multinomial_probit, parameters)
ADF_smoother_output_truth = Gaussian_smoother(data_multinomial_probit, parameters, ADF_filter_output_truth)

# EM algorithm
EM_iterations = 1000
eps_sample_size = 10000

initial_parameters = list(a_0 = parameters$a_0, P_0 = parameters$P_0, 
                            G = parameters$G, 
                            W = parameters$W, 
                            Sigma = diag(L))

EM_output_ADF = multinomial_probit_EM_invSigma(data_multinomial_probit, 
                                               ADF_sampling_filter, "std",
                                               initial_parameters, target_trace,
                                               EM_iterations, eps_sample_size)

par(mfrow = c(3, 3), mar = c(2, 2, 2, 1), oma = c(0, 0, 2, 0))  # optional: adjust margins

# Plot each matrix element
for (i in 1:3) {
  for (j in 1:3) {
    
    balanced   = EM_output_ADF$Sigma_history[1:1000,i,j]
    
    balanced   = balanced
  
    # Define y-limits to fit both
    ylims <- range(balanced)
    
    # Plot
    plot(balanced, type = "l", col = "red", lwd = 2,
         ylim = ylims, main = paste0("(", i, ",", j, ")"),
         xlab = "", ylab = "")
  }
}

save(data_multinomial_probit, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_data_multinomial_probit.Rdata"))
save(ADF_filter_output_truth, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_ADF_truth.Rdata"))
save(EM_output_ADF,           file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_ADF_EM.Rdata"))

EM_parameters = list(a_0 = a_0, P_0 = P_0, 
                     G = G, 
                     W = W, 
                     Sigma = EM_output_ADF$Sigma_history[EM_iterations+1,,])

ADF_filter_output_EM   = ADF_sampling_filter(data_multinomial_probit, EM_parameters)
ADF_smoother_output_EM = Gaussian_smoother(data_multinomial_probit, EM_parameters, ADF_filter_output_EM)

df_EM_ADF = data.frame(type ="EM ADF", 
                   log_likelihood_distance    = ADF_filter_output_EM$log_likelihood-ADF_filter_output_truth$log_likelihood,
                   filter_distance_from_truth = mean_wasserstein_sampling(ADF_filter_output_EM$sampling_history,
                                                                          ADF_filter_output_truth$sampling_history),
                   smoothing_distance_from_truth = mean_wasserstein_gaussian(ADF_smoother_output_EM, 
                                                                             ADF_smoother_output_truth))

# Create a data frame for the line plot
df <- data.frame(
  iterations = 1:(EM_iterations+1),
  log_likelihood = EM_output_ADF$log_likelihood_history
)

# Create the plot
ggplot() +
  geom_line(data = df, aes(x = iterations, y = log_likelihood, color = "EM Log-Likelihood"), size = 1) + # Line plot
  geom_hline(aes(yintercept = ADF_filter_output_truth$log_likelihood, color = "ADF Log-Likelihood"), 
             linetype = "dashed", size = 1) + # Horizontal line
  labs(
    x = "Iterations",  # X-axis label
    y = "Log-Likelihood",  # Y-axis label
    color = "Legend"  # Legend title
  ) +
  scale_color_manual(
    values = c("EM Log-Likelihood" = "blue", "ADF Log-Likelihood" = "red"), # Custom colors
    labels = c("EM Log-Likelihood" = "EM Log-Likelihood", "ADF Log-Likelihood" = "EP Exact Output") # Custom labels
  ) +
  theme_minimal()  # Minimal theme


color_list <- c("blue", "black", "red", "green", "brown")

# PLOT FILTERING
# Create a data frame for plotting
data_list <- list()
for (j in 1:4) {
  i =1
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = ADF_filter_output_EM$a_history[, 1*(j-1) + i],
    y_upper = ADF_filter_output_EM$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_filter_output_EM$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_filter_output_EM$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_filter_output_EM$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i],
    facet_row = j,
    facet_col = "ADF"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = ADF_filter_output_truth$a_history[, 1*(j-1) + i],
    y_upper = ADF_filter_output_truth$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_filter_output_truth$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_filter_output_truth$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_filter_output_truth$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i+1],
    facet_row = j,
    facet_col = "True parameters"
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
    y2 = ADF_smoother_output_EM$a_history[, 1*(j-1) + i],
    y_upper = ADF_smoother_output_EM$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_smoother_output_EM$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_smoother_output_EM$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_smoother_output_EM$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    color = color_list[i],
    facet_row = j,
    facet_col = "EP"
  )
  data_list[[length(data_list) + 1]] <- df_filter
  
  df_filter <- data.frame(
    x = seq(0, Time, length.out = Time+1),
    y1 = data_multinomial_probit$Beta[, 1*(j-1) + i],
    y2 = ADF_smoother_output_truth$a_history[, 1*(j-1) + i],
    y_upper = ADF_smoother_output_truth$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_smoother_output_truth$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_smoother_output_truth$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_smoother_output_truth$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
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


# # Extract the matrix dimensions
# n <- 3 #dim(EM_output_ADF$G_history)[2] # Number of rows/columns in G
# 
# # Reshape data for ggplot
# df_list <- lapply(1:n, function(i) {
#   lapply(1:n, function(j) {
#     data.frame(
#       iteration = 1:(EM_iterations+1),
#       G_value = EM_output_ADF$G_history[1:(EM_iterations+1), i, j],
#       i = i,
#       j = j,
#       true_value = parameters$G[i, j]
#     )
#   })
# })
# df <- do.call(rbind, do.call(rbind, df_list))
# 
# # Create the plot
# ggplot(df, aes(x = iteration, y = G_value)) +
#   geom_line(color = "blue", size = 1) +                           # Line plot for G history
#   geom_hline(aes(yintercept = true_value),                       # Horizontal line for true value
#              color = "red", linetype = "dashed", size = 1) +
#   facet_grid(i ~ j, labeller = label_both) +                     # Subfigures for each (i, j)
#   labs(
#     x = "Iterations",                                            # X-axis label
#     y = "G Value",                                               # Y-axis label
#   ) +
#   theme_minimal() +
#   theme(
#     strip.text = element_text(size = 8)                          # Adjust facet label size
#   )
# 
# # Extract the matrix dimensions
# n <- 3 # dim(EM_output$Sigma_history)[2] # Number of rows/columns in Sigma
# 
# # Reshape data for ggplot
# df_list <- lapply(1:n, function(i) {
#   lapply(1:n, function(j) {
#     data.frame(
#       iteration = 1:(EM_iterations+1),
#       Sigma_value = EM_output_ADF$Sigma_history[1:(EM_iterations+1), i, j],
#       i = i,
#       j = j,
#       true_value = parameters$Sigma[i, j]
#     )
#   })
# })
# 
# df <- do.call(rbind, do.call(rbind, df_list))
# 
# # Create the plot
# ggplot(df, aes(x = iteration, y = Sigma_value)) +
#   geom_line(color = "blue", size = 1) +                           # Line plot for G history
#   geom_hline(aes(yintercept = true_value),                       # Horizontal line for true value
#              color = "red", linetype = "dashed", size = 1) +
#   facet_grid(i ~ j, labeller = label_both) +                     # Subfigures for each (i, j)
#   labs(
#     x = "Iterations",                                            # X-axis label
#     y = "Sigma Value",                                               # Y-axis label
#   ) +
#   theme_minimal() +
#   theme(
#     strip.text = element_text(size = 8)                          # Adjust facet label size
#   )
