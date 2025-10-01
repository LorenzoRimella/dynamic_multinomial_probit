library(dplyr)
library(tidyr)

library(ggplot2)


set.seed(678)

# SET UP THE SIMULATION
L = 3 # number of choices
Time= 10 # number of time steps

p_X = 2 # number of subject dependent covariates
p_Y = 3 # number of choice dependent covariates
p = p_X*L + p_Y

expected_m_t = 15

# GENERATE PARAMETERS

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = diag(beta_size)

G = 0.5*diag(beta_size)
for(i in 2:beta_size)
{
  sign = sample(c(-1, 1), 1)
  G[i,i-1] = sign*0.25
  G[i-1,i] = sign*0.25
}

W = 0.05*diag(beta_size)

Sigma = 0.2*diag(L)
for(i in 2:L)
{
  sign = sample(c(-1, 1), 1)
  Sigma[i,i-1] = sign*0.1
  Sigma[i-1,i] = sign*0.1
}

invSigma = solve(Sigma)

target_trace = sum(diag(invSigma))

load("DMP_cineca/data/trace_constrained_repeat_EM_p_9_mt_15_L_3_Time_10/trace_constrained_repeat_EM_p_9_mt_15_L_3_Time_10_Sigma_repeated_sample.Rdata")
load("DMP_cineca/data/trace_constrained_repeat_EM_p_9_mt_15_L_3_Time_10/trace_constrained_repeat_EM_p_9_mt_15_L_3_Time_10_log_likelihood_true.Rdata")
load("DMP_cineca/data/trace_constrained_repeat_EM_p_9_mt_15_L_3_Time_10/trace_constrained_repeat_EM_p_9_mt_15_L_3_Time_10_log_likelihood_history.Rdata")

EM_iterations = 200
n_samples = 50

df <- as.data.frame(log_likelihood_history - t(t(log_likelihood_true))%*%array(1, dim = c(1, EM_iterations+1)))
df$id <- factor(1:n_samples)  # Add an ID for each row

# Optional: give time-based column names
colnames(df)[1:(EM_iterations+1)] <- paste0("t", 1:(EM_iterations+1))

df_long <- pivot_longer(df, 
                        cols = starts_with("t"), 
                        names_to = "time", 
                        values_to = "value")

# Convert 'time' to numeric (optional, for ordering)
df_long$time <- as.numeric(sub("t", "", df_long$time))

ggplot(df_long, aes(x = time, y = value, group = id, color = id)) +
  geom_line() +
  labs(x = "Time", y = "Value", title = "N Trajectories Over T Time Points") +
  theme_minimal()

matrix_normalization <- function(Sigma)
{
  trace_Sigma = sum(diag(Sigma))
  
  return(L*Sigma/trace_Sigma)
}

mean_Sigma = Reduce("+", lapply(seq(1, n_samples), function(i) Sigma_repeated_sample[i,,]))/n_samples
sd_Sigma   = sqrt(Reduce("+", lapply(seq(1, n_samples), function(i) (Sigma_repeated_sample[i,,] - mean_Sigma)^2))/n_samples)


parameters = list(a_0 = a_0, P_0 = P_0, G = G, W = W, Sigma = Sigma)

# Get row and column indices for each element
positions <- expand.grid(row = 1:L, col = 1:L)

# Flatten the matrices
vec1 <- as.vector(mean_Sigma)
vec2 <- as.vector(sd_Sigma)
vec3 <- as.vector(parameters$Sigma)

# Create the formatted table
result <- data.frame(
  element = paste(positions$row, positions$col, sep = ","),
  EM_mean_sd = paste(round(vec1,2), "+/-", round(vec2,2)),
  Sigma_value = round(vec3,2)
)

print(result)

# optional: filter only lower diagonal part of Sigma
indeces <- c(1,2,5,3,6,9)
vec1 <- vec1[indeces]
vec2 <- vec2[indeces]
vec3 <- vec3[indeces]

# plot estimated values and real values
lb <- vec1 - vec2
ub <- vec1 + vec2

df <- data.frame(
  lb = lb,
  ub = ub,
  mSigma = vec1,
  trueSigma = vec3,
  id = 1:length(lb)  # add an ID to distinguish each interval
)

# Create the labels as expressions like expression(Sigma[x[1]])
y_labels <- c(expression(Sigma[11]),
              expression(Sigma[21]),
              expression(Sigma[31]),
              expression(Sigma[12]),
              expression(Sigma[22]),
              expression(Sigma[32]),
              expression(Sigma[13]),
              expression(Sigma[23]),
              expression(Sigma[33]))

# optional
y_labels <- y_labels[indeces]

y_pos <- 1:length(vec1)
points_df <- data.frame(
  x = c(vec1, vec3),
  y = rep(y_pos, 2),
  group = factor(rep(c("Estimated mean", "True value"), each = length(y_pos)))
)

ggplot(df) +
  geom_segment(aes(x = lb, xend = ub, y = id, yend = id), size = 1.5, color = "gray") +
  geom_point(data = points_df, aes(x = x, y = y, color = group), size = 4) +
  
  # Legend & aesthetics
  scale_color_manual(values = c("Estimated mean" = "blue", "True value" = "red")) +
  theme_minimal() +
  labs(
    x = "Value",
    y = "",
    title = expression("True vs estimated values of " * Sigma),
    color = "Group"
  )+
  theme(
    plot.title = element_text(hjust = 0.5)
  )+
  scale_y_continuous(breaks = df$id, labels = y_labels)



# plot the two matrices
library(reshape2)

# For mat1
df1 <- melt(parameters$Sigma)
colnames(df1) <- c("row", "col", "value")
df1$matrix <- "True"

df2 <- melt(mean_Sigma)
colnames(df2) <- c("row", "col", "value")
df2$matrix <- "Estimated mean"

combined_df <- rbind(df1, df2)
combined_df$matrix <- factor(combined_df$matrix, levels = c("True", "Estimated mean"))

ggplot(combined_df, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_y_reverse() +
  facet_wrap(~ matrix) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = expression(Sigma*" true vs estimated heatmaps"),
    x = "Column",
    y = "Row",
    fill = "Value"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

