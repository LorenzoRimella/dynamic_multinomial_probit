library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(gridExtra)

source("R_directory/model.R")
source("R_directory/multinomial_probit_sampling_filter.R")


L = 4 # number of choices
Time= 50 # number of time steps

p = 5 # number of subject dependent covariate

N_max = 100 # maximum number of customers
lambda_M = 40
hyperparameters = list(lambda_M = lambda_M, N_max = N_max)

X_cov = matrix(runif(N_max*p, -0.5, 0.5), nrow = N_max, ncol = p)

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = 0.1*diag(beta_size)

G = 0.5*diag(beta_size)
G[2,1] = -0.3
G[1,2] = -0.3
G[4,3] = 0.3
G[3,4] = 0.3

W = 0.1*diag(beta_size)

n_factors = 3

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

color_list <- c("blue", "black", "red", "green", "brown")

# Create a data frame for plotting
data_list <- list()
for (j in 1:3) {
  for (i in 1:3) {
    df_filter <- data.frame(
      x = seq(0, Time, length.out = Time+1),
      y1 = data_multinomial_probit$Beta[, 3*(j-1) + i],
      y2 = ADF_filter_sampler_output$a_history[, 3*(j-1) + i],
      y_upper = ADF_filter_sampler_output$a_history[, 3*(j-1) + i] + 1.96 * sqrt(ADF_filter_sampler_output$P_history[, 3*(j-1) + i, 3*(j-1) + i]),
      y_lower = ADF_filter_sampler_output$a_history[, 3*(j-1) + i] - 1.96 * sqrt(ADF_filter_sampler_output$P_history[, 3*(j-1) + i, 3*(j-1) + i]),
      color = color_list[i],
      facet_row = j,
      facet_col = i
    )
    data_list[[length(data_list) + 1]] <- df_filter
  }
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

ADF_prediction_prob = cbind(1-apply(ADF_filter_sampler_output$prediction_prob_history, 1, sum), 
                            ADF_filter_sampler_output$prediction_prob_history)
prediction_prob = cbind(1-apply(data_multinomial_probit$prediction_prob, 1, sum), 
                        data_multinomial_probit$prediction_prob)

prediction_prob_long <- melt(prediction_prob, varnames = c("Time", "Category"), value.name = "Predicted_Prob")
prediction_prob_long$Category = as.factor(prediction_prob_long$Category-1)
prediction_prob_long$Time <- as.numeric(prediction_prob_long$Time)

ADF_prediction_prob_long <- melt(ADF_prediction_prob, varnames = c("Time", "Category"), value.name = "ADF_Predicted_Prob")
ADF_prediction_prob_long$Category = as.factor(ADF_prediction_prob_long$Category-1)
ADF_prediction_prob_long$Time <- as.numeric(ADF_prediction_prob_long$Time)

# Merge the original dataframe df_y with prediction_prob_long by Time and Category
df_merged <- merge(prediction_prob_long, ADF_prediction_prob_long, by = c("Time", "Category"))

# Now, plot the Mean and Predicted_Prob together
ggplot(df_merged, aes(x = Time, group = Category, color = Category)) +
  geom_line(aes(y = Predicted_Prob), linetype = "solid") +  # Solid line for Mean
  geom_line(aes(y = ADF_Predicted_Prob), linetype = "dashed") +  # Dashed line for Predicted Probabilities
  geom_hline(yintercept = 0.2, linetype = "dotted", color = "black") +
  scale_x_discrete(breaks = c(1, seq(0, 100, by = 10)[2:11])) +
  facet_grid(Category ~ ., scales = "free_y") +  # Facet by Category (one plot per row)
  labs(x = "Time", y = "Expected scaled counts") +
  theme_minimal()

# L1 distance from random classifier
time_seq = seq(1, dim(prediction_prob)[1])
ADF_l1 = apply(abs(prediction_prob-ADF_prediction_prob), 1, sum)
print(mean(ADF_l1))
#0.5443504
random_l1 = apply(abs(prediction_prob-0.2), 1, sum)
print(mean(random_l1))
#0.6615634
ADF_l2 = sqrt(apply((prediction_prob-ADF_prediction_prob)^2, 1, sum))
print(mean(ADF_l2))
#0.298156
random_l2 = sqrt(apply((prediction_prob-0.2)^2, 1, sum))
print(mean(random_l2))
#0.3537693

df = data.frame(Time = time_seq, ADF_l1 = ADF_l1, random_l1 = random_l1, ADF_l2 = ADF_l2, random_l2 = random_l2)

df_long <- df %>%
  pivot_longer(cols = -Time, names_to = "Metric", values_to = "Value") %>%
  mutate(Group = ifelse(Metric %in% c("ADF_l1", "random_l1"), "ADF_l2", "random_l2"))

# Plot using ggplot with facet_wrap
ggplot(df_long, aes(x = Time, y = Value, color = Metric, group = Metric)) +
  geom_line() +
  facet_wrap(~ Group, ncol = 1, scales = "free_y") +  # Facet by Group, one plot below the other
  labs(x = "Time", y = "Value", color = "Lines") +
  theme_minimal()

