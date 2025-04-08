source("R_directory/dynamic_multinomial_probit/scripts/model.R")
source("R_directory/dynamic_multinomial_probit/scripts/multinomial_probit_filter_EP.R")
source("R_directory/dynamic_multinomial_probit/scripts/multinomial_probit_filter_ADF.R")
source("R_directory/dynamic_multinomial_probit/scripts/multinomial_probit_filter_exact.R")
source("R_directory/dynamic_multinomial_probit/scripts/multinomial_probit_gaussian_smoothing.R")
source("R_directory/dynamic_multinomial_probit/scripts/multinomial_probit_EM.R")

library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

name_simulation = "inference_synthetic_ADF"
directory_name  = "R_directory/dynamic_multinomial_probit/data/inference_synthetic_ADF/"

# Create the directory if needed
if (!dir.exists(directory_name)) {
  dir.create(directory_name)
}

set.seed(123)

# GENERATE PARAMETERS

L = 5 # number of choices
Time= 20 # number of time steps

p = 10 # number of subject dependent covariate

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

hyperparameters = list(lambda_M = 10, N_max = 30, prediction_bool= FALSE)

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

# True likelihood

ADF_filter_output_exact = ADF_sampling_filter(data_multinomial_probit, parameters)

ident_parameters = list(a_0 = a_0, P_0 = 100*P_0, G = G, W = 100*W, Sigma = 100*Sigma)
ident_ADF_filter_output_exact = ADF_sampling_filter(data_multinomial_probit, ident_parameters)


ident_parameters = list(a_0 = a_0, P_0 = 100*P_0, G = G, W = 100*W, Sigma = 100*Sigma)
ident_ADF_filter_output_exact = ADF_sampling_filter(data_multinomial_probit, ident_parameters)

# EM algorithm
EM_iterations = 20
learn_G = TRUE
learn_W = FALSE
learn_Sigma = TRUE

if((learn_G & learn_W) & !(learn_Sigma))
{
  initial_parameters = list(a_0 = a_0, P_0 = P_0, 
                            G = diag(beta_size), 
                            W = diag(dim(parameters$W)[1]), 
                            Sigma = Sigma)
}

if(learn_Sigma & !(learn_G) & !(learn_W))
{
  initial_parameters = list(a_0 = a_0, P_0 = P_0, 
                            G = G, 
                            W = W, 
                            Sigma = diag(L))
}

if(learn_Sigma & learn_G & !(learn_W))
{
  initial_parameters = list(a_0 = a_0, P_0 = P_0, 
                            G = diag(beta_size), 
                            W = W, 
                            Sigma = diag(L))
}

if(learn_Sigma & learn_G & learn_W)
{
  initial_parameters = list(a_0 = a_0, P_0 = P_0, 
                            G = diag(beta_size), 
                            W = diag(dim(parameters$W)[1]), 
                            Sigma = diag(L))
}

EM_output_ADF = multinomial_probit_EM(data_multinomial_probit, ADF_sampling_filter, initial_parameters, 
                                  EM_iterations,
                                  learn_G, learn_W, learn_Sigma, type = "ADF")

save(data_multinomial_probit, file    = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_data_multinomial_probit.Rdata"))
save(EM_output_ADF, file    = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_EM_output.Rdata"))

ADF_filter_output_exact = ADF_sampling_filter(data_multinomial_probit, parameters)
ADF_smoother_output_exact = Gaussian_smoother(data_multinomial_probit, parameters, ADF_filter_output_exact)

EM_parameters = list(a_0 = a_0, P_0 = P_0, 
                     G = EM_output_ADF$G_history[EM_iterations+1,,], 
                     W = W, 
                     Sigma = EM_output_ADF$Sigma_history[EM_iterations+1,,])
                     
ADF_filter_output_EM = ADF_sampling_filter(data_multinomial_probit, EM_parameters)
ADF_smoother_output_EM = Gaussian_smoother(data_multinomial_probit, EM_parameters, ADF_filter_output_EM)

# Create a data frame for the line plot
df <- data.frame(
  iterations = 1:(EM_iterations+1),
  log_likelihood = EM_output_ADF$log_likelihood_history
)

# Create the plot
ggplot() +
  geom_line(data = df, aes(x = iterations, y = log_likelihood, color = "EM Log-Likelihood"), size = 1) + # Line plot
  geom_hline(aes(yintercept = ADF_filter_output_exact$log_likelihood, color = "ADF Log-Likelihood"), 
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

# Extract the matrix dimensions
n <- 3 #dim(EM_output_ADF$G_history)[2] # Number of rows/columns in G

# Reshape data for ggplot
df_list <- lapply(1:n, function(i) {
  lapply(1:n, function(j) {
    data.frame(
      iteration = 1:(EM_iterations+1),
      G_value = EM_output_ADF$G_history[1:(EM_iterations+1), i, j],
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
n <- 3 # dim(EM_output$Sigma_history)[2] # Number of rows/columns in Sigma

# Reshape data for ggplot
df_list <- lapply(1:n, function(i) {
  lapply(1:n, function(j) {
    data.frame(
      iteration = 1:(EM_iterations+1),
      Sigma_value = EM_output_ADF$Sigma_history[1:(EM_iterations+1), i, j],
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
    y2 = ADF_filter_output_exact$a_history[, 1*(j-1) + i],
    y_upper = ADF_filter_output_exact$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_filter_output_exact$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_filter_output_exact$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_filter_output_exact$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
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
for (j in 5:10) {
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
    y2 = ADF_smoother_output_exact$a_history[, 1*(j-1) + i],
    y_upper = ADF_smoother_output_exact$a_history[, 1*(j-1) + i] + 1.96 * sqrt(ADF_smoother_output_exact$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
    y_lower = ADF_smoother_output_exact$a_history[, 1*(j-1) + i] - 1.96 * sqrt(ADF_smoother_output_exact$P_history[, 1*(j-1) + i, 1*(j-1) + i]),
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
