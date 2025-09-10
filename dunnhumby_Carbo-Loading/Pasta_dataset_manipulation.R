rm(list = ls())
# load the library
library(dplyr)
library(stringr)


# Load a CSV file
data_products <- read.csv("dunnhumby - Carbo-Loading CSV/dh_product_lookup.csv")
data_transactions <- read.csv("dunnhumby - Carbo-Loading CSV/dh_transactions.csv")

# View the first few rows
head(data_products)
head(data_transactions)

table(data_products$commodity)

pasta_subset <- data_products[data_products$commodity=="pasta",]

# Select only the columns you want from df2 before joining
pasta_subset <- pasta_subset %>% select(upc, brand, product_size)

# Convert product size to ounces
convert_to_oz <- function(x) {
  # standardize
  x <- str_trim(toupper(x))
  
  if (str_detect(x, "LB")) {
    # extract number before "LB"
    lbs <- as.numeric(str_extract(x, "\\d+\\.?\\d*"))
    return(lbs * 16)
  } else if (str_detect(x, "Z")) {
    # extract number before "OZ"
    oz <- as.numeric(str_extract(x, "\\d+\\.?\\d*"))
    return(oz)
  } else if (str_detect(x, "OUNCE")) {
    oz <- as.numeric(str_extract(x, "\\d+\\.?\\d*"))
    return(oz)
  } else {
    return(NA) # not a valid weight
  }
}

pasta_subset$product_size <- sapply(pasta_subset$product_size, convert_to_oz)

pasta_subset <- na.omit(pasta_subset)

# Keep only transactions related to pasta purchases
# Filter rows where product_id is in the vector
filtered_transactions <- data_transactions %>%
  filter(upc %in% pasta_subset$upc)

# Exclude coupons
filtered_transactions <- filtered_transactions %>%
  filter(coupon == 0)

# Exclude negative prices
filtered_transactions <- filtered_transactions %>%
  filter(dollar_sales > 0)

# # Take subset of transactions with only one
# filtered_transactions <- filtered_transactions %>%
#   filter(units == 1)

# Join datasets
merged_data <- filtered_transactions %>%
  left_join(pasta_subset, by = "upc")

# Exclude private labels
merged_data <- merged_data %>%
  filter(!(brand %in% c("Private Label", "Private Label Premium", "Private Label Value")))

# take only second year
merged_data <- merged_data %>%
  filter(day > 364)


print(merged_data)

sort(table(merged_data$brand), decreasing = T)

L <- 5 # numer of top brands to consider
sum(sort(table(merged_data$brand), decreasing = T)[1:(L+1)])

# cumulative sums
cumsum(sort(table(merged_data$brand), decreasing = T)/sum(table(merged_data$brand)))

# # obtain total sales by brand
# # Group by brand and sum sales
# brand_sales <- merged_data %>%
#   group_by(brand) %>%
#   summarise(total_sales = sum(dollar_sales), .groups = "drop")
# 
# print(brand_sales, n=60)
# 
# sorted_brand_sales <- brand_sales %>% 
#   arrange(desc(total_sales))
# 
# print(sorted_brand_sales, n=60)

# take the names of the L categories



# # add baseline category (not used here, where we use Barilla as reference class)
# categories <- sorted_brand_sales$brand[1:L]
# 
# merged_data_L <- merged_data
# merged_data_L$brand[!(merged_data_L$brand %in% categories)] <- "L0"

# filter only observations for the 6 top-selling brands (covering 92% of sales)
categories <- names(sort(table(merged_data$brand), decreasing = T)[1:(L+1)])
merged_data_L <- merged_data[merged_data$brand %in% categories,]


# compute price per ounce
merged_data_L$log_price_per_ounce <- log(merged_data_L$dollar_sales/(merged_data_L$units*merged_data_L$product_size))


options(max.print = 100)
A <- table(merged_data_L$store,merged_data_L$brand)

# Get row names where all entries > 0
valid_rows <- rownames(A)[apply(A, 1, function(x) all(x > 0))]

# merged_data_L <- merged_data_L %>%
#   filter(store %in% valid_rows)
# 
# dim(merged_data_L)
# 
# table(merged_data_L$day)



storesToKeep <- names(sort(rowSums(A[valid_rows,]), decreasing = T))[1:10] # keep the 10 most-selling stores

merged_data_L <- merged_data_L %>%
  filter(store %in% storesToKeep)

merged_data_L$t <- merged_data_L$day
table(merged_data_L$t)


##########################
# CREATE DATASETS
##########################
merged_data_L$transID <- 1:nrow(merged_data_L)
merged_data_L$t       <- merged_data_L$t - min(merged_data_L$t) + 1
merged_data_L$brand   <- as.integer(factor(merged_data_L$brand, levels = categories)) - 1


####### TO DO: IMPROVE THIS PART
merged_data_L$t[merged_data_L$t>160] <- merged_data_L$t[merged_data_L$t>160] - 1 # 161 was missing
##################

# create tidy dataset storing average log-prices per ounce of each (day, store, brand) tuple

# 1) daily averages where data exist
daily_avg <- merged_data_L %>%
  group_by(store, brand, t) %>%
  summarise(avg_log_oz = mean(log_price_per_ounce), .groups = "drop")

# 2) build the full grid of (store, brand, t) you want to cover
all_t      <-  1:length(unique(merged_data_L$t))
all_pairs  <- merged_data_L %>% distinct(store, brand)
panel <- tidyr::expand_grid(t = all_t, all_pairs)

# 3) join and nearest-neighbor fill within each (store, brand) ## TO DOUBLE CHECK
nearest_fill <- function(t, y) {
  # y has NAs at missing days
  if (all(is.na(y))) return(rep(NA_real_, length(y)))   # no data for this store-brand
  o <- order(t); t <- t[o]; y <- y[o]                   # ensure sorted by t
  known <- which(!is.na(y))
  kt <- t[known]; ky <- y[known]
  # for each t, pick value from closest known day (tie -> earlier day)
  picked <- vapply(t, function(tt) ky[ which.min(abs(kt - tt)) ], numeric(1))
  # keep original where present
  out <- ifelse(is.na(y), picked, y)
  out[order(o)]  # restore original order (not strictly needed here)
}

avg_log_prices_t_store_brand <- panel %>%
  left_join(daily_avg, by = c("t","store","brand")) %>%
  group_by(store, brand) %>%
  mutate(avg_log_oz_filled = nearest_fill(t, avg_log_oz)) %>%
  ungroup()

sum(is.na(avg_log_prices_t_store_brand$avg_log_oz_filled))


Tmin  = min(merged_data_L$t)
Tmax  = max(merged_data_L$t)

data_Y         <- merged_data_L[,c("transID","t","brand")]
data_X         <- merged_data_L[,c("transID","t","log_price_per_ounce","brand","store")]

covariates_names <- c("log_price_per_ounce")

X = list()
Y = list()

for(current_t in Tmin:Tmax)
{
  current_X <- subset(data_X, t == current_t)
  current_Y <- subset(data_Y, t == current_t)
  
  n_t <- nrow(current_X)
  
  # create design matrix for time t
  X_array <- array(NA, dim = c(n_t, length(covariates_names)+L, L))

  for (i in 1:n_t) {
    choice_i                     <- current_Y$brand[i]
    store_i                      <- current_X$store[i]
    log_price_per_ounce_choice_i <- current_X$log_price_per_ounce[i]
    
    # create design matrix for individual i of time n_t
    X_i <- double(L) # this is ok only if we have only one covariate (choice specific)
    
    if(choice_i == 0){
      log_price_per_ounce_baseline <- current_X$log_price_per_ounce[i]
    } else{
      log_price_per_ounce_baseline <- avg_log_prices_t_store_brand[avg_log_prices_t_store_brand$t == current_t &
                                                                     avg_log_prices_t_store_brand$store == store_i &
                                                                     avg_log_prices_t_store_brand$brand == 0, ]$avg_log_oz_filled
    }
    
    for (choice in 1:L) {
      if(choice == choice_i){
        log_price_per_ounce_i <- log_price_per_ounce_choice_i
      } else{
        log_price_per_ounce_i <- avg_log_prices_t_store_brand[avg_log_prices_t_store_brand$t == current_t &
                                                                avg_log_prices_t_store_brand$store == store_i &
                                                                avg_log_prices_t_store_brand$brand == choice_i, ]$avg_log_oz_filled
      }
      
      X_i[choice] <- log_price_per_ounce_i - log_price_per_ounce_baseline
    }
    
    X_array[i,,] = rbind(X_i,diag(L))
    rbind(X_i,diag(L))
    
  }
  
  X[[current_t]] = X_array
  Y[[current_t]] = current_Y$brand
  
}


save(X,Y,file = "dataCarbo.RData")
