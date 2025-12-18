library(dplyr)
library(tidyr)

library(xtable)


load("DMP_cineca/data/parallel_filt_smoot_lik_p_9_mt_10_L_3_Time_10/parallel_filt_smoot_lik_p_9_mt_10_L_3_Time_10_data.Rdata")

# Compute mean and standard deviation for each type
summary_table <- data_frame_output %>%
  group_by(type) %>%
  summarise(across(
    where(is.numeric),
    list(mean = ~mean(.), sd = ~sd(.)),
    .names = "{.col}_{.fn}"
  )) %>%
  mutate(across(
    ends_with("_mean"),
    ~paste0(round(.x, 3), " ± ", round(get(sub("_mean$", "_sd", cur_column())), 3)),
    .names = "{.col}_formatted"
  )) %>%
  dplyr::select(type, ends_with("_formatted")) %>%
  rename_with(~sub("_mean_formatted$", "", .x))

# Print table
print(summary_table, row.names = FALSE)
xtable(t(summary_table))

load("DMP_cineca/data/full_parallel_filt_smoot_lik_p_9_mt_10_L_3_Time_50/full_parallel_filt_smoot_lik_p_9_mt_10_L_3_Time_50_data.Rdata")

# Compute mean and standard deviation for each type
summary_table <- data_frame_output %>%
  group_by(type) %>%
  summarise(across(
    where(is.numeric),
    list(mean = ~mean(.), sd = ~sd(.)),
    .names = "{.col}_{.fn}"
  )) %>%
  mutate(across(
    ends_with("_mean"),
    ~paste0(round(.x, 3), " ± ", round(get(sub("_mean$", "_sd", cur_column())), 3)),
    .names = "{.col}_formatted"
  )) %>%
  dplyr::select(type, ends_with("_formatted")) %>%
  rename_with(~sub("_mean_formatted$", "", .x))

# Print table
print(summary_table, row.names = FALSE)
xtable(t(summary_table))

load("DMP_cineca/data/full_parallel_filt_smoot_lik_p_9_mt_50_L_3_Time_50/full_parallel_filt_smoot_lik_p_9_mt_50_L_3_Time_50_data.Rdata")

# Compute mean and standard deviation for each type
summary_table <- data_frame_output %>%
  group_by(type) %>%
  summarise(across(
    where(is.numeric),
    list(mean = ~mean(.), sd = ~sd(.)),
    .names = "{.col}_{.fn}"
  )) %>%
  mutate(across(
    ends_with("_mean"),
    ~paste0(round(.x, 3), " ± ", round(get(sub("_mean$", "_sd", cur_column())), 3)),
    .names = "{.col}_formatted"
  )) %>%
  dplyr::select(type, ends_with("_formatted")) %>%
  rename_with(~sub("_mean_formatted$", "", .x))

# Print table
print(summary_table, row.names = FALSE)
xtable(t(summary_table))