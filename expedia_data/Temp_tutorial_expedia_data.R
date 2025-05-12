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

data_X = read.csv("data_X.csv")
data_Y = read.csv("data_Y.csv")

# covariates_names = c(  "prop_brand_bool",   "prop_location_score1",
#                        "promotion_flag", "log_price_usd",
#                        "prop_starrating_0", "prop_starrating_1", "prop_starrating_2",
#                        "prop_starrating_3", "prop_starrating_4", "prop_starrating_5")

covariates_names = c(  "prop_brand_bool",   "prop_location_score1",
                       "promotion_flag", "log_price_usd",
                       "prop_starrating")

L = 10

t_0  = min(data_X$t)
Time = max(data_X$t)

end   = 60 #Time - t_0 + 1

X = list()
Y = list()

set.seed(5)

data_X$log_price_usd= 0.5*scale(data_X$log_price_usd)
data_X$prop_starrating= 0.5*scale(data_X$prop_starrating)
data_X$prop_location_score1= 0.5*scale(data_X$prop_location_score1)

current_t = t_0
for(i in 1:end)
{
  #print(i)
  current_X = subset(data_X,    t == current_t)
  current_Y = subset(data_Y, t == current_t)
  searches  = unique(current_X$srch_id)
  
  n_searches = floor(0.05*length(searches))
  
  # if(n_searches<2){n_searches = 2}
  
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


set.seed(5)

source("../scripts/model.R")
source("../scripts/multinomial_probit_filter_ADF.R")
source("../scripts/multinomial_probit_filter_exact.R")
source("../scripts/multinomial_probit_filter_EP.R")
source("../scripts/multinomial_probit_gaussian_smoothing.R")
source("Temp_expedia_EM.R")
source("../scripts/multinomial_probit_EM.R")

Time= length(Y)

p = 5 #10 # number of subject dependent covariate

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = 0.1*diag(beta_size)

G = 0.9*diag(beta_size)
W = 0.05*diag(beta_size)
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

# EM algorithm
EM_iterations = 100 
learn_G = F
learn_W = F
learn_Sigma = TRUE

initial_parameters = list(a_0 = a_0,
                            P_0 = P_0,
                            G = G,
                            W = W,
                            Sigma = diag(L))

# EM_output = multinomial_probit_EM(data_multinomial_probit, EP_sampling_filter, initial_parameters, 
#                                   EM_iterations,
#                                   learn_G, learn_W, learn_Sigma)
# 
# 
# 
# initial_parameters = list(a_0 = matrix(rep(0, beta_size), ncol = 1),
#                           P_0 = diag(beta_size),
#                           G = EM_output$G_history[5,,],
#                           W = EM_output$W_history[5,,],
#                           Sigma = EM_output$Sigma_history[5,,])

# EM_output = multinomial_probit_EM(data_multinomial_probit, EP_sampling_filter, initial_parameters, 
#                                   EM_iterations,
#                                   learn_G, learn_W, learn_Sigma)

EM_output <- multinomial_probit_EM_constrained(data_multinomial_probit, EP_sampling_filter, initial_parameters,
                                               EM_iterations,
                                               learn_G, learn_W, learn_Sigma, type = "nonstd",
                                               log = c(FALSE, NA),
                                               target_trace = L)


# as.matrix(read.csv("R_directory/dynamic_multinomial_probit/expedia_data/G_history.csv"))









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

