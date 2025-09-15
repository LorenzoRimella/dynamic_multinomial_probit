rm(list=ls())
# Set the directory to the path of the current script
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# If we change L it does not run
# Check rango X matrix e se min autovalore di X^t X non è troppo piccolo
# Standardizzare le variabili

library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(gridExtra)

load("dataCarbo.RData")

# shorten X and Y (to be removed)
TimeEnd <- 10
X <- head(X,TimeEnd)
Y <- head(Y,TimeEnd)

L = 5

source("../scripts/model.R")
source("../scripts/utils.R")
source("../scripts/multinomial_probit_filter_ADF.R")
source("../scripts/multinomial_probit_filter_exact.R")
source("../scripts/multinomial_probit_filter_EP.R")
source("../scripts/multinomial_probit_gaussian_smoothing.R")
source("../expedia_data/Temp_expedia_EM.R")
source("../scripts/multinomial_probit_EM_Sigma.R")

Time= length(Y)

p = 1
L = 5

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = 3*diag(beta_size)

c = 0.9
G = c*diag(beta_size)
W = (1-c^2)*P_0
Sigma = diag(L)

data_multinomial_probit = list(X = X, Y = Y)

time_steps <- seq_along(data_multinomial_probit$Y)  # Time steps (1, 2, 3, ...)
data_list <- lapply(seq_along(data_multinomial_probit$Y), 
                    function(i) {data.frame(table(data_multinomial_probit$Y[[i]])) %>%
                        mutate(time = i)})

# Combine all time steps into one data frame
df_y <- do.call(rbind, data_list)
names(df_y) <- c("Category", "Count", "Time")  # Rename columns

df_y$Time <- as.factor(df_y$Time) 
df_y$Category <- factor(df_y$Category, levels = as.character(sort(as.numeric(levels(df_y$Category)))))

ggplot(df_y, aes(x = Time, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  scale_x_discrete(breaks = c(1, seq(0, Time, by = L))) +
  labs(x = "Time", y = "Count", fill = "Category") +
  theme_minimal() 

# Parameters
beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = 1/beta_size*diag(beta_size)

c = 0.8
G = c*P_0
W = (1-c^2)*P_0

load("~/R_directory/dynamic_multinomial_probit/dunnhumby_Carbo-Loading/CarboLoading_EM_output_T30_200.RData")

EM_iterations = 200

Sigma_estimated = EM_output$Sigma_history[EM_iterations,,]

plot(EM_output$Sigma_history[,1,1], type = "l")
for(j in 2:10)
{
  lines(EM_output$Sigma_history[,1,j])
}
for(i in 2:10)
{
  for(j in i:10)
  {
    lines(EM_output$Sigma_history[,i,j])
  }
}

plot(EM_output$log_likelihood_history, type = "l")

parameters = list(a_0 = a_0,
                  P_0 = P_0,
                  G   = G,
                  W   = W,
                  Sigma = Sigma_estimated)

# Filter
EP_max_iterations = 5
MC_samples = 10000
eps_sample_size = 10000
tol = 10e-3

EP_sampling_filter_to_use <- function(data, parameters)
{
  return(EP_sampling_filter(data, parameters, EP_max_iterations, MC_samples, tol))
}

EP_filter_output_EM = EP_sampling_filter_to_use(data_multinomial_probit, parameters)


EP_output = EP_filter_output_EM
data = data_multinomial_probit

predictive = EP_predictive(data, parameters, EP_output, MC_samples)


cross_entropy_1 = 0 
cross_entropy_2 = 0 
mean_prob     = rep(0, L+1)
hist_prob = as.vector(table(unlist(data_multinomial_probit$Y)+1)/length(unlist(data_multinomial_probit$Y)))

for(t in 1:100){
  
  y_t     = data_multinomial_probit$Y[[t]]+1
  indexes = 1:dim(ADF_filter_output_EM$prediction_prob_history[[t]])[2]
  
  prob_t = ADF_filter_output_EM$prediction_prob_history[[t]][cbind(y_t, indexes)]
  
  cross_entropy_1 = cross_entropy_1 +mean(-log(prob_t))
  cross_entropy_2 = cross_entropy_2 +mean(-log(hist_prob[y_t]))
  
  mean_prob = mean_prob + rowMeans(ADF_filter_output_EM$prediction_prob_history[[t]])
}

cross_entropy_1/100
cross_entropy_2/100
-log(1/11)

mean_prob/100



























save(ADF_filter_output_EM, file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_ADF_out.Rdata"))



plot(EM_output_EP$Sigma_history[,1,1], type = "l")
for(j in 2:10)
{
  lines(EM_output_EP$Sigma_history[,1,j])
}
for(i in 2:10)
{
  for(j in i:10)
  {
    lines(EM_output_EP$Sigma_history[,i,j])
  }
}

plot(EM_output_EP$log_likelihood_history[100:200])

save(EM_output_EP,           file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_EP_EM.Rdata"))

Trace_EM = sum(diag(EM_output_ADF$Sigma_history[EM_iterations+1,,]))
EM_parameters = list(a_0 = a_0, P_0 = P_0/Trace_EM, 
                     G = EM_output_ADF$G_history[EM_iterations+1,,], 
                     W = W/Trace_EM, 
                     Sigma = EM_output_ADF$Sigma_history[EM_iterations+1,,]/Trace_EM)

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
             linetype = "dashed", linewidth = 1) + # Horizontal line
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
