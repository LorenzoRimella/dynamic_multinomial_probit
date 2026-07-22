rm(list = ls())

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  getwd()
}
setwd(script_dir)

library(ggplot2)
library(patchwork)

dir.create("figures", showWarnings = FALSE)

theme_dmnp <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(5.5, 6, 5.5, 6)
    )
}

save_plot <- function(plot, stem, width, height) {
  ggsave(
    filename = file.path("figures", paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
  ggsave(
    filename = file.path("figures", paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

# -----------------------------------------------------------------------------
# Figure 1: EM history and estimated correlation matrix
# -----------------------------------------------------------------------------

result_file <- file.path(
  "data",
  "inference_expedia_longEM_alpha_0.99_EP_p_3_mt_14.64_L_10_Time_50",
  "inference_expedia_longEM_alpha_0.99_EP_p_3_mt_14.64_L_10_Time_50_EM.Rdata"
)
if (!file.exists(result_file)) {
  stop("Expedia EM output not found: ", result_file)
}

result_env <- new.env(parent = emptyenv())
load(result_file, envir = result_env)
em_output <- result_env$EM_output_EP

loglik <- as.numeric(em_output$log_likelihood_history)
iteration <- seq_along(loglik) - 1L
window <- 25L
moving_average <- vapply(seq_along(loglik), function(i) {
  mean(loglik[max(1L, i - window + 1L):i], na.rm = TRUE)
}, numeric(1))

loglik_df <- data.frame(
  iteration = iteration,
  loglik = loglik,
  moving_average = moving_average
)

p_loglik <- ggplot(loglik_df, aes(iteration, loglik)) +
  geom_line(color = "grey70", linewidth = 0.3) +
  geom_line(aes(y = moving_average), color = "#0072B2", linewidth = 0.75) +
  scale_x_continuous(breaks = seq(0, 500, by = 100), expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "EM iteration", y = "Log marginal likelihood") +
  theme_dmnp()

Sigma_final <- em_output$Sigma_history[dim(em_output$Sigma_history)[1], , ]
correlation_final <- cov2cor(Sigma_final)
correlation_df <- expand.grid(
  row = seq_len(nrow(correlation_final)),
  column = seq_len(ncol(correlation_final))
)
correlation_df$correlation <- as.vector(correlation_final)

p_correlation <- ggplot(
  correlation_df,
  aes(column, row, fill = correlation)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Correlation"
  ) +
  scale_x_continuous(breaks = 1:10, expand = c(0, 0)) +
  scale_y_reverse(breaks = 1:10, expand = c(0, 0)) +
  coord_fixed() +
  labs(x = "Position", y = "Position") +
  theme_dmnp() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.key.height = grid::unit(0.34, "in"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5)
  )

training_plot <- p_loglik + p_correlation +
  plot_layout(widths = c(1.55, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_plot(training_plot, "expedia_training_Gio", width = 7.2, height = 3.15)

# -----------------------------------------------------------------------------
# Figure 2: EP filtering coefficients over the full observation period
# -----------------------------------------------------------------------------

source(file.path("scripts", "utils.R"))
source(file.path("scripts", "multinomial_probit_filter_ADF.R"))
source(file.path("scripts", "multinomial_probit_filter_EP.R"))

data_X <- read.csv(file.path("data", "data_X10_rand.csv"))
data_Y <- read.csv(file.path("data", "data_Y10_rand.csv"))

covariate_names <- c(
  "prop_location_score1",
  "log_price_usd",
  "prop_starrating"
)
for (variable in covariate_names) {
  data_X[[variable]] <- as.double(scale(data_X[[variable]]))
}

L <- 10L
time_ids <- sort(unique(data_X$t))
if (!identical(time_ids, seq.int(min(time_ids), max(time_ids)))) {
  stop("The Expedia time index is not consecutive.")
}

X <- vector("list", length(time_ids))
Y <- vector("list", length(time_ids))

# This reproduces the ordering used by the original analysis script.
set.seed(678)
for (i in seq_along(time_ids)) {
  current_X <- data_X[data_X$t == time_ids[i], , drop = FALSE]
  current_Y <- data_Y[data_Y$t == time_ids[i], , drop = FALSE]
  searches <- unique(current_X$srch_id)
  searches <- sample(searches, length(searches), replace = FALSE)

  X_array <- array(NA_real_, dim = c(length(searches), length(covariate_names) + L, L))
  Y_array <- numeric(length(searches))

  for (j in seq_along(searches)) {
    search_X <- current_X[current_X$srch_id == searches[j], , drop = FALSE]
    search_X <- search_X[order(search_X$position), covariate_names, drop = FALSE]
    search_Y <- current_Y$booking_pos_sum[current_Y$srch_id == searches[j]]

    if (nrow(search_X) != L || length(search_Y) != 1L) {
      stop("Unexpected Expedia search structure for srch_id = ", searches[j])
    }

    X_array[j, , ] <- rbind(t(as.matrix(search_X)), diag(L))
    Y_array[j] <- search_Y
  }

  X[[i]] <- X_array
  Y[[i]] <- Y_array
}

beta_size <- length(covariate_names) + L
alpha <- 0.99
parameters <- list(
  a_0 = matrix(0, nrow = beta_size, ncol = 1),
  P_0 = diag(beta_size),
  G = alpha * diag(beta_size),
  W = (1 - alpha^2) * diag(beta_size),
  Sigma = Sigma_final
)

filter_output <- EP_sampling_filter(
  data = list(X = X, Y = Y),
  parameters = parameters,
  EP_max_iterations = 11,
  MC_samples = 10000,
  tol = 1e-3
)

# Drop the prior stored in row 1; the remaining rows correspond to Mar 1--Jun 30.
a_hist <- filter_output$a_history[-1, , drop = FALSE]
P_hist <- filter_output$P_history[-1, , , drop = FALSE]
dates <- as.Date("2013-03-01") + time_ids

coefficient_labels <- c(
  "Location score",
  "Log price",
  "Star rating",
  paste("Position", 1:10)
)

sd_hist <- t(vapply(
  seq_len(nrow(P_hist)),
  function(i) sqrt(pmax(diag(P_hist[i, , ]), 0)),
  numeric(beta_size)
))

coefficient_df <- data.frame(
  date = rep(dates, times = beta_size),
  coefficient = factor(
    rep(coefficient_labels, each = length(dates)),
    levels = coefficient_labels
  ),
  mean = as.vector(a_hist),
  sd = as.vector(sd_hist)
)
coefficient_df$lower <- coefficient_df$mean - 1.96 * coefficient_df$sd
coefficient_df$upper <- coefficient_df$mean + 1.96 * coefficient_df$sd

month_breaks <- seq(as.Date("2013-03-01"), as.Date("2013-06-01"), by = "month")
attribute_df <- coefficient_df[coefficient_df$coefficient %in% coefficient_labels[1:3], ]
position_df <- coefficient_df[coefficient_df$coefficient %in% coefficient_labels[4:13], ]

attribute_palette <- c(
  "Location score" = "#009E73",
  "Log price" = "#D55E00",
  "Star rating" = "#0072B2"
)

p_attributes <- ggplot(
  attribute_df,
  aes(date, mean, color = coefficient, fill = coefficient)
) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.16, color = NA) +
  geom_line(linewidth = 0.55) +
  facet_wrap(~coefficient, nrow = 1) +
  scale_color_manual(values = attribute_palette, guide = "none") +
  scale_fill_manual(values = attribute_palette, guide = "none") +
  scale_x_date(
    limits = range(dates),
    breaks = month_breaks,
    date_labels = "%b",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(x = NULL, y = expression(beta[t])) +
  theme_dmnp()

p_positions <- ggplot(position_df, aes(date, mean)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#8DA0CB", alpha = 0.25) +
  geom_line(color = "#2F4B7C", linewidth = 0.5) +
  facet_wrap(~coefficient, ncol = 5) +
  scale_x_date(
    limits = range(dates),
    breaks = month_breaks,
    date_labels = "%b",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(x = "Date", y = expression(beta[t])) +
  theme_dmnp(base_size = 8.5) +
  theme(axis.text.x = element_text(size = 7.5))

coefficient_plot <- p_attributes / p_positions +
  plot_layout(heights = c(1, 1.75)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_plot(coefficient_plot, "expedia_coefficients_Gio", width = 7.2, height = 5.15)
