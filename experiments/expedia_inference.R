rm(list=ls())
# Set the directory to the path of the current script

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path)) {
  setwd(dirname(script_path))
  message("Working directory set to: ", getwd())
} else {
  message("Script path not detected. Working directory not changed.")
}

source("../scripts/utils.R")
source("../scripts/model.R")
source("../scripts/multinomial_probit_filter_EP.R")
source("../scripts/multinomial_probit_filter_SMC.R")
source("../scripts/multinomial_probit_filter_ADF.R")
source("../scripts/multinomial_probit_filter_exact.R")
source("../scripts/multinomial_probit_gaussian_smoothing.R")
source("../scripts/multinomial_probit_EM_Sigma.R")

library(tidyr)
library(dplyr)
library(transport)
library(expm) # for sqrtm (matrix square root)

library(future)
library(future.apply)

library(dplyr)
library(tidyr)

# LOAD DATA
data_X = read.csv("../data/data_X10_rand.csv")
data_X$prop_location_score1 = as.double(scale(data_X$prop_location_score1))
data_X$log_price_usd = as.double(scale(data_X$log_price_usd))
data_X$prop_starrating = as.double(scale(data_X$prop_starrating))

data_Y = read.csv("../data/data_Y10_rand.csv")

covariates_names = c( #"prop_brand_bool",   
                      "prop_location_score1",
                      #"promotion_flag",
                      "log_price_usd",
                      "prop_starrating")

L = 10
t_0  = min(data_X$t)
Time = max(data_X$t)

end   = 50 #Time - t_0 + 1

X = list()
Y = list()

set.seed(678)

percent_search = 1

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

alpha_num = 0.99

p = dim(X[[1]])[2]-L # number of subject dependent covariate

name_init = paste("inference_expedia_alpha_", alpha_num, sep="")
name_init = paste(name_init, "_EP_p_", sep="")

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

# Parameters
beta_size = p + L

a_0 = matrix(rep(0, beta_size), ncol = 1)
P_0 = diag(beta_size)

G = alpha_num*diag(beta_size)

W = (1-alpha_num*alpha_num)*P_0

# the invariant distr is N(0,P_0)

Sigma = diag(L)

initial_parameters = list(a_0 = a_0,
                          P_0 = P_0,
                          G   = G,
                          W   = W,
                          Sigma = Sigma)

# Filter
EP_max_iterations = 11
MC_samples = 10000
eps_sample_size = 10000
tol = 1e-3

target_trace = L

EP_sampling_filter_to_use <- function(data, parameters)
{
  return(EP_sampling_filter(data, parameters, EP_max_iterations, MC_samples, tol))
}

# EM algorithm
EM_iterations = 100

EM_output_EP = multinomial_probit_EM_Sigma( data_multinomial_probit,
                                            EP_sampling_filter_to_use, type = "nonstd",
                                            initial_parameters, target_trace,
                                            EM_iterations, eps_sample_size,
                                            log = c(TRUE, list(c(directory_name, name_simulation))))

save(EM_output_EP,           file = paste(paste(directory_name, sep = "", name_simulation), sep = "", "_EM.Rdata"))
