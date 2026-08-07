library(dplyr)
library(tidyr)

library(ggplot2)


set.seed(678)

# Simulation setting
L = 3
Time = 10

p_X = 2
p_Y = 3
p = p_X*L + p_Y

expected_m_t = 15

beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = diag(beta_size)

alpha = 0.5
G = alpha*diag(beta_size)

W = (1-alpha*alpha)*P_0

Sigma = matrix(c(1.6, -0.1, 0,
                 -0.1, 0.2, 0.1,
                 0, 0.1, 1.2),
               nrow = L, byrow = TRUE)

Sigma_time_mean = array(NA, c(3, 3, 3))
Sigma_time_sd = array(NA, c(3, 3, 3))
Sigma_time_quantile = array(NA, c(3, 2, 3, 3))

n_samples = 48

listing = c(10, 30, 60)

for(i in 1:3)
{
  Time = listing[i]
  string = paste("data/Sigma_recovery_p_9_mt_15_L_3_Time_", Time, sep = "")
  string = paste(string, "/Sigma_recovery_p_9_mt_15_L_3_Time_" ,sep = "")
  string = paste(string, Time ,sep = "")
  
  load(paste(string, "_Sigma_repeated_sample.Rdata",sep = ""))
  
  mean_Sigma = Reduce("+", lapply(seq(1, n_samples), function(i) Sigma_repeated_sample[i,,]))/n_samples
  sd_Sigma   = sqrt(Reduce("+", lapply(seq(1, n_samples), function(i) (Sigma_repeated_sample[i,,] - mean_Sigma)^2))/n_samples)
  
  Sigma_time_mean[i,,] =(mean_Sigma)
  Sigma_time_sd[i,,] = (sd_Sigma)
  
  Sigma_time_quantile[i,,,] = apply(Sigma_repeated_sample, c(2, 3), quantile, probs = c(0.025, 0.975))
  
  
}

positions <- expand.grid(row = 1:L, col = 1:L)

vec1_10 <- as.vector(Sigma_time_mean[1,,])
vec2_10 <- as.vector(Sigma_time_sd[1,,])

vec1_30 <- as.vector(Sigma_time_mean[2,,])
vec2_30 <- as.vector(Sigma_time_sd[2,,])

vec1_60 <- as.vector(Sigma_time_mean[3,,])
vec2_60 <- as.vector(Sigma_time_sd[3,,])

vec3 <- as.vector(Sigma)

result <- data.frame(
  element = paste(positions$row, positions$col, sep = ","),
  EM_10_mean_sd = paste(round(vec1_10,2), "+/-", round(vec2_10,2)),
  EM_30_mean_sd = paste(round(vec1_30,2), "+/-", round(vec2_30,2)),
  EM_60_mean_sd = paste(round(vec1_60,2), "+/-", round(vec2_60,2)),
  Sigma_value = round(vec3,2)
)


library(xtable)

print(result)
print(xtable(result), include.rownames = FALSE)

indeces <- c(1,2,5,3,6,9)

mean_10 = vec1_10[indeces]
mean_30 = vec1_30[indeces]
mean_60 = vec1_60[indeces]

vec1_10 = as.vector(Sigma_time_quantile[1,1,,])
vec2_10 = as.vector(Sigma_time_quantile[1,2,,])

vec1_30 = as.vector(Sigma_time_quantile[2,1,,])
vec2_30 = as.vector(Sigma_time_quantile[2,2,,])

vec1_60 = as.vector(Sigma_time_quantile[3,1,,])
vec2_60 = as.vector(Sigma_time_quantile[3,2,,])

vec1_10 <- vec1_10[indeces]
vec2_10 <- vec2_10[indeces]

vec1_30 <- vec1_30[indeces]
vec2_30 <- vec2_30[indeces]

vec1_60 <- vec1_60[indeces]
vec2_60 <- vec2_60[indeces]

vec3 <- vec3[indeces]

lb_10 <- vec1_10
ub_10 <- vec2_10

lb_30 <- vec1_30
ub_30 <- vec2_30

lb_60 <- vec1_60
ub_60 <- vec2_60

df <- data.frame(
  mean_10 = mean_10,
  lb_10 = lb_10,
  ub_10 = ub_10,
  mean_30 = mean_30,
  lb_30 = lb_30,
  ub_30 = ub_30,
  mean_60 = mean_60,
  lb_60 = lb_60,
  ub_60 = ub_60,
  trueSigma = vec3,
  id = seq_along(lb_10)
)

y_labels <- c(
  expression(Sigma[11]),
  expression(Sigma[21]),
  expression(Sigma[22]),
  expression(Sigma[31]),
  expression(Sigma[32]),
  expression(Sigma[33])
)

df_long <- df %>%
  pivot_longer(
    cols = matches("^(mean|lb|ub)_"),
    names_to = c(".value", "level"),
    names_pattern = "(mean|lb|ub)_(\\d+)"
  )

points_df <- df_long %>%
  select(id, level, trueSigma, mean) %>%
  pivot_longer(
    cols = c(trueSigma, mean),
    names_to = "group",
    values_to = "x"
  ) %>%
  mutate(
    y = id + as.numeric(level)/100 - 0.2
  )

ggem = ggplot(df_long) +
  geom_segment(
    aes(
      x = lb, xend = ub,
      y = id + as.numeric(level)/100 - 0.2,
      yend = id + as.numeric(level)/100 - 0.2,
      color = level
    ),
    linewidth = 1.7
  ) +
  
  geom_point(
    data = points_df,
    aes(x = x, y = y, color = group),
    size = 3.4
  ) +
  
  scale_color_manual(
    name = "",
    values = c(
      "10" = "gray80",
      "30" = "gray60",
      "60" = "gray40",
      "mean" = "blue",
      "trueSigma" = "red"
    ),
    breaks = c("10", "30", "60", "mean", "trueSigma"),
    labels = c("T = 10", "T = 30", "T = 60",
               "Estimated mean", "True value")
  ) +
  
  theme_minimal(base_size = 17) +
  labs(
    x = "Value",
    y = ""
  ) +
  theme(
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 15),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.35, color = "grey88"),
    plot.margin = margin(8, 10, 8, 8)
  ) +
  scale_y_continuous(breaks = df$id, labels = y_labels)


ggsave("figures/em_repeated.png", plot = ggem, width = 9, height = 5.4,
       dpi = 300, units = "in", bg = "white")
ggsave("figures/em_repeated.pdf", plot = ggem, width = 9, height = 5.4,
       units = "in", bg = "white")
