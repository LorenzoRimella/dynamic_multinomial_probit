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
source("../DMP_cineca/scripts/multinomial_probit_EM_Sigma.R")

library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(transport)
library(expm) # for sqrtm (matrix square root)
library(future.apply)

library(dplyr)
library(tidyr)

# LOAD DATA
data_X = read.csv("../expedia_data/data_X.csv")
#data_X$prop_location_score1 = as.double(scale(data_X$prop_location_score1))
#data_X$log_price_usd = as.double(scale(data_X$log_price_usd))
#data_X$prop_starrating = as.double(scale(data_X$prop_starrating))

data_Y = read.csv("../expedia_data/data_Y.csv")

covariates_names = c( "prop_brand_bool",   
                      "prop_location_score1",
                      "promotion_flag",
                      "log_price_usd",
                      "prop_starrating")

L = 10
t_0  = min(data_X$t)
Time = max(data_X$t)

end   = 50 #Time - t_0 + 1

X = list()
Y = list()

set.seed(5)

percent_search = 0.05

current_t = t_0
for(i in 1:end)
{
  current_X = subset(data_X,    t == current_t)
  current_Y = subset(data_Y, t == current_t)
  searches  = unique(current_X$srch_id)
  
  n_searches = floor(percent_search*length(searches))
  
  if(n_searches<2){n_searches = 2}
  
  searches = sample(searches, n_searches, replace = FALSE, prob = NULL)
  
  X_array = array(NA, dim = c(length(searches), length(covariates_names)+L, L))
  Y_array = array(NA, dim = c(length(searches)))
  index_s=1
  for(s in searches)
  {
    #print(s)
    current_X_search = subset(current_X, srch_id == s)
    current_X_search = current_X_search[order(current_X_search$position), ]
    current_X_search = current_X_search[, covariates_names]
    
    current_Y_search = subset(current_Y, srch_id == s)
    current_Y_search = current_Y_search$booking_pos_sum
    
    X_array[index_s,,] = rbind(t(as.matrix(current_X_search)), diag(rep(1, L)))
    Y_array[index_s]   = current_Y_search
    index_s = index_s +1
  }
  
  X[[i]] = X_array
  Y[[i]] = Y_array
  
  current_t = current_t+1
}

# NAME SIMULATION
Time= length(X)

p = dim(X[[1]])[2]-L # number of subject dependent covariate

name_init="inference_expedia_EP_bigP0bigSigma_p_"

expected_m_t = mean(sapply(seq(1,length(X)), function(i) return(dim(X[[i]])[1])))

output_name = sim_name(name_init, L, Time, p, expected_m_t)
name_simulation = output_name$name
directory_name  = paste("../data/", output_name$name_for_dir, sep ="")

# Create the directory if needed
if (!dir.exists(directory_name)) {
  dir.create(directory_name)
}

# CREATE AND PLOT DATA

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
  scale_x_discrete(breaks = c(1, seq(0, Time, by = 10))) +
  labs(x = "Time", y = "Count", fill = "Category") +
  theme_minimal() 

# Parameters
beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = diag(beta_size)

alpha = 0.5
G = alpha*diag(beta_size)

W = (1-alpha*alpha)*P_0

# the invariant distr is N(0,P_0)

Sigma = diag(L)

initial_parameters = list(a_0 = a_0,
                          P_0 = P_0,
                          G   = G,
                          W   = W,
                          Sigma = Sigma)

# Check mean variance of the fixed effect
Sigma_cumul = initial_parameters$P_0
sum_var     = array(0, dim= c(L,L))
elem = 0 
for(t in 1:10)
{
  Sigma_cumul = initial_parameters$G%*%Sigma_cumul%*%t(initial_parameters$G) + initial_parameters$W
  
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

# Filter
EP_max_iterations = 5
MC_samples = 10000
eps_sample_size = 10000
tol = 10e-3

target_trace = L

EP_sampling_filter_to_use <- function(data, parameters)
{
  return(EP_sampling_filter(data, parameters, EP_max_iterations, MC_samples, tol))
}

# EM algorithm
EM_iterations = 200

# EM_output_EP = multinomial_probit_EM_Sigma_G( data_multinomial_probit, 
#                                             EP_sampling_filter_to_use, type = "nonstd",
#                                             initial_parameters,
#                                             EM_iterations, eps_sample_size, 
#                                             log = c(TRUE, name_simulation))
# 

EM_output_EP = multinomial_probit_EM_Sigma( data_multinomial_probit,
                                            EP_sampling_filter_to_use, type = "nonstd",
                                            initial_parameters, target_trace,
                                            EM_iterations, eps_sample_size,
                                            log = c(TRUE, name_simulation))

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

save(EM_output_EP,           file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_EP_EM_LL_lower.Rdata"))

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
