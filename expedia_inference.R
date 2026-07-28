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

source(file.path("scripts", "utils.R"))
source(file.path("scripts", "multinomial_probit_filter_ADF.R"))
source(file.path("scripts", "multinomial_probit_filter_EP.R"))
source(file.path("scripts", "multinomial_probit_gaussian_smoothing.R"))
source(file.path("scripts", "multinomial_probit_EM_Sigma.R"))

dir.create("figures", showWarnings = FALSE)

# Common Expedia configuration. The first 50 days are used to estimate Sigma;
# filtering is subsequently performed over the complete observation period.
alpha <- 0.99
em_time <- 50L
em_iterations <- 500L
EP_max_iterations <- 11L
MC_samples <- 10000L
eps_sample_size <- 10000L
tol <- 1e-3
L <- 10L

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
# Prepare the Expedia data once for both EM estimation and full-period filtering
# -----------------------------------------------------------------------------

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

time_ids <- sort(unique(data_X$t))
if (!identical(time_ids, seq.int(min(time_ids), max(time_ids)))) {
  stop("The Expedia time index is not consecutive.")
}
if (length(time_ids) < em_time) {
  stop("Fewer than ", em_time, " Expedia days are available.")
}

X <- vector("list", length(time_ids))
Y <- vector("list", length(time_ids))

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

data_full <- list(X = X, Y = Y)
data_em <- list(
  X = X[seq_len(em_time)],
  Y = Y[seq_len(em_time)]
)

beta_size <- length(covariate_names) + L
initial_parameters <- list(
  a_0 = matrix(0, nrow = beta_size, ncol = 1),
  P_0 = diag(beta_size),
  G = alpha * diag(beta_size),
  W = (1 - alpha^2) * diag(beta_size),
  Sigma = diag(L)
)

# -----------------------------------------------------------------------------
# Figure 1: EM history and estimated covariance matrix
# -----------------------------------------------------------------------------

expected_m_t <- mean(lengths(data_em$Y))
name_init <- paste0("inference_expedia_longEM_alpha_", alpha, "_EP_p_")
output_name <- sim_name(
  name_init = name_init,
  L = L,
  Time = em_time,
  p = length(covariate_names),
  expected_m_t = expected_m_t
)
name_simulation <- output_name$name
directory_name <- file.path("data", sub("/$", "", output_name$name_for_dir))
dir.create(directory_name, recursive = TRUE, showWarnings = FALSE)
result_file <- file.path(directory_name, paste0(name_simulation, "_EM.Rdata"))

if (file.exists(result_file)) {
  result_env <- new.env(parent = emptyenv())
  load(result_file, envir = result_env)
  if (!exists("EM_output_EP", envir = result_env, inherits = FALSE)) {
    stop("The saved Expedia file does not contain EM_output_EP: ", result_file)
  }
  EM_output_EP <- result_env$EM_output_EP
  if (length(EM_output_EP$log_likelihood_history) != em_iterations + 1L) {
    stop(
      "The saved Expedia EM run has ",
      length(EM_output_EP$log_likelihood_history) - 1L,
      " iterations instead of ", em_iterations,
      ". Move or remove it before rerunning: ", result_file
    )
  }
  message("Using the completed alpha = ", alpha, " EM run: ", result_file)
} else {
  EP_sampling_filter_to_use <- function(data, parameters) {
    EP_sampling_filter(
      data = data,
      parameters = parameters,
      EP_max_iterations = EP_max_iterations,
      MC_samples = MC_samples,
      tol = tol
    )
  }

  set.seed(678)
  EM_output_EP <- multinomial_probit_EM_Sigma(
    data = data_em,
    filter = EP_sampling_filter_to_use,
    type = "nonstd",
    initial_parameters = initial_parameters,
    target_trace = L,
    EM_iterations = em_iterations,
    eps_sample_size = eps_sample_size,
    log = list(TRUE, c(paste0(directory_name, .Platform$file.sep), name_simulation))
  )
  save(EM_output_EP, file = result_file)
}

em_output <- EM_output_EP

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
  scale_x_continuous(
    breaks = seq(0, em_iterations, by = 100),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(x = "EM iteration", y = "Log marginal likelihood") +
  theme_dmnp()

Sigma_final <- em_output$Sigma_history[dim(em_output$Sigma_history)[1], , ]
Sigma_df <- expand.grid(
  row = seq_len(nrow(Sigma_final)),
  column = seq_len(ncol(Sigma_final))
)
Sigma_df$value <- as.vector(Sigma_final)
Sigma_limit <- max(abs(Sigma_final))

p_sigma <- ggplot(
  Sigma_df,
  aes(column, row, fill = value)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-Sigma_limit, Sigma_limit),
    name = "Value"
  ) +
  scale_x_continuous(breaks = 1:10, expand = c(0, 0)) +
  scale_y_reverse(breaks = 1:10, expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title = expression(hat(Sigma)),
    x = "Position",
    y = "Position"
  ) +
  theme_dmnp(base_size = 8.5) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.key.height = grid::unit(0.3, "in"),
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7)
)

training_plot <- p_loglik + p_sigma +
  plot_layout(widths = c(1.55, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 10))

save_plot(
  training_plot,
  "expedia_training",
  width = 7.2,
  height = 3.15
)

if (identical(Sys.getenv("EXPEDIA_TRAINING_ONLY"), "1")) {
  quit(save = "no", status = 0)
}

# -----------------------------------------------------------------------------
# Figure 2: EP filtering coefficients over the full observation period
# -----------------------------------------------------------------------------

parameters <- list(
  a_0 = initial_parameters$a_0,
  P_0 = initial_parameters$P_0,
  G = initial_parameters$G,
  W = initial_parameters$W,
  Sigma = Sigma_final
)

set.seed(679)
filter_output <- EP_sampling_filter(
  data = data_full,
  parameters = parameters,
  EP_max_iterations = EP_max_iterations,
  MC_samples = MC_samples,
  tol = tol
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

save_plot(
  coefficient_plot,
  "expedia_coefficients",
  width = 7.2,
  height = 5.15
)
