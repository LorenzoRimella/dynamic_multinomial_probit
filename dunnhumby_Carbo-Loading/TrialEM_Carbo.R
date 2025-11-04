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

# a_0 = matrix(rep(0, beta_size), ncol = 1)
# P_0 = 3*diag(beta_size)
# 
# G = 0.9*diag(beta_size)
# W = 0.05*diag(beta_size)

a_0 = matrix(rep(0, beta_size), ncol = 1)

P_0 =  diag(beta_size) # 1/(beta_size)*diag(beta_size) # 3*diag(beta_size) # 1/(beta_size)*diag(beta_size)

c = 0.8

G = c*diag(dim(P_0)[1])

W = (1-c^2)*P_0 # 0.01


Sigma = matrix(0.5,L,L)
diag(Sigma) = 1

# Sigma = diag(L)

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
                          Sigma = Sigma)

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

set.seed(5)

output_name = sim_name(paste("CarboLoading_c_",c,"P_0_",round(P_0[1,1],2), sep = ""), L, TimeEnd, beta_size, 100)
name_simulation = output_name$name
directory_name  = paste("../data/", output_name$name_for_dir, sep ="")

# Create the directory if needed
if (!dir.exists(directory_name)) {
  dir.create(directory_name)
}

cat(P_0[1,1])


# EM_output <- multinomial_probit_EM_invSigma(data = data_multinomial_probit,
#                                             filter = EP_sampling_filter,
#                                             type = "nonstd",
#                                             initial_parameters = initial_parameters,
#                                             target_trace = L,
#                                             EM_iterations = EM_iterations,
#                                             eps_sample_size = 1e4,
#                                             log = list(TRUE, c(directory_name, name_simulation)))


EM_output <- multinomial_probit_EM_Sigma(data = data_multinomial_probit,
                                            filter = EP_sampling_filter,
                                            type = "nonstd",
                                            initial_parameters = initial_parameters,
                                            target_trace = L,
                                            EM_iterations = EM_iterations,
                                            eps_sample_size = 1e4,
                                            log = list(TRUE, c(directory_name, name_simulation)))

# save(EM_output,P_0,c, file = "CarboLoading_EM_output_T10_500_P0_05_c_09_W_005.RData")