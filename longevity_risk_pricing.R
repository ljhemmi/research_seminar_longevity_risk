#####################################################################
# Title: 
# Research Seminar - Population Ageing
#
# Authors: 
# - Yannick Cadonau
# - Claudio Hess
# - Lucas Jan Hemmi
#
# Description:
# this script models the pricing of longevity risk for a 
# hypothetical insurance
#
# Notation:
# variables with M (F) correspond to Male (Female) population data        
#       
#####################################################################

# 0. Setup ----

## 0.0 Libraries ----

library(svd)
library(patchwork)
library(tidyverse)
library(readxl)
library(writexl)
library(parallel)
library(zoo)

# clear environment
rm(list=ls())

getwd()


## 0.1 General Functions ----

# function to convert mx (central death rate) to qx (probability of death)
# according to two different approaches
mx_to_qx <- function(m_x, approximation = TRUE) {
  if (approximation) {
    # Approximation formula
    q_x <- m_x / (1 + 0.5 * m_x)
  } else {
    # Exact formula
    q_x <- 1 - exp(-m_x)
  }
  return(q_x)
}

# function to replace mx (central death rates) of value = zero in the mortality data
# this ensures that the natural logarithm can be calculated
# mode = "average": takes the average of the left and right adjacent cell of m_x which are 0
#     this creates smoothing over the year (horizontal axis)
# mode = "add": add a small constant (defined by "small_number") to the 0 to avoid
#     taking the logarithm of 0 but 0+small_number instead.
replace_zeros <- function(df, mode = "average", small_number = 1e-10) {
  df <- as.data.frame(df)  # Ensure that it's a data frame
  
  replace_zero_with_average <- function(row) {
    n <- length(row)
    for (i in 2:(n - 1)) {
      if (row[i] == 0) {
        left <- row[i - 1]
        right <- row[i + 1]
        if (left != 0 | right != 0) {
          row[i] <- mean(c(left, right), na.rm = TRUE)
        }
      }
    }
    # Handling edge cases for the first and last elements
    if (n > 1 && row[1] == 0 && row[2] != 0) {
      row[1] <- row[2]
    }
    if (n > 1 && row[n] == 0 && row[n - 1] != 0) {
      row[n] <- row[n - 1]
    }
    return(row)
  }
  
  if (mode == "average") {
    df <- t(apply(df, 1, replace_zero_with_average))  # Apply the function across rows
    df <- as.data.frame(df)  # Convert back to data frame
  } else if (mode == "add") {
    df[df == 0] <- small_number
  } else {
    stop("Invalid mode. Use 'average' or 'add'.")
  }
  
  return(df)
}


# function that generates mx-tables based on k_t, a_x, & b_x
generate_mx_table <- function(k_t, a_x, b_x) {
  # Ensure k_t, a_x, and b_x are vectors of correct lengths
  if(length(a_x) != length(b_x)) {
    stop("a_x and b_x must have the same length.")
  }
  
  # Initialize an empty matrix to store the mx values
  mx_table <- matrix(nrow = length(a_x), ncol = length(k_t))
  
  # Calculate the mx values for each age group and each k_t
  for (i in 1:length(a_x)) {
    for (j in 1:length(k_t)) {
      mx_table[i, j] <- exp(a_x[i] + b_x[i] * k_t[j])
    }
  }
  
  # Assign row names and column names for clarity
  rownames(mx_table) <- paste("AgeGroup", 1:length(a_x))
  colnames(mx_table) <- paste("Year", seq_along(k_t))
  
  return(mx_table)
}


## 0.2 Simulation Functions ----

# Function to automatically derive a Lee-Carter Model given a mx Table
mx_to_LC_model <- function(mx_input){
  # Preparing Matrix and fitting a_x
  # get the ln(m_x) matrix
  log_mx <- log(mx_input[,-1])
  
  # get the mean log-rate for each age-group (over all t) (avg(ln(m_x)))
  mean_log_mx <- rowMeans(log_mx) 
  
  # subtract this average from each age-group for each year
  for (i in 1:nrow(log_mx)){
    log_mx[i,] <- log_mx[i,] - mean_log_mx[i]
  }
  
  # final estimates of a_x (the mean over time of the log mortality rate for the age groups)
  a_x <- mean_log_mx
  
  # 2.1 fitting Lee-Carter Model (first step estimation)
  # Singular Value Decomposition to fit b_x & k_t
  d <- svd(as.matrix(log_mx), 
           nv=nrow(log_mx), 
           nu=ncol(log_mx))
  # due to the sign ambiguity we can rescale b_x and k_t (the -1 cancel eachother in the full formula)
  b_x <- d$u[,1]*(-1)
  
  # k_t calculation based on Girosi & King (2007, p. 4)
  k_t <- as.vector(b_x %*% as.matrix(log_mx))
  
  # consolidate the outputs of the LC model
  LC_output <- list(a_x = a_x, b_x = b_x, k_t = k_t, age_categories = mx_input$Age, years = as.numeric(colnames(mx_input[,-1])))
  return(LC_output)
}


# function that performs the monte carlo simulations of the k_t process
# returns as output a nested list
LC_to_simulation_optimized <- function(input_LC, # takes the output of mx_to_LC_model as input
                                       end_year_T = 2050, # up until which year to simulate
                                       n_simulations = 10, # number of monte carlo iterations
                                       drift_type = "full", # if set to "full" all years are considered for drift (e.g. "trend" calculation), alternative "window" 
                                       drift_window = 20, # if drift_type = "window", define how many years to look back for non-full period drift
                                       sigma_type = "full", # # if set to "full" all years are considered for sigma (e.g. "volatility" calculation), alternative "window" 
                                       sigma_window = 20, # if drift_type = "window", define how many years to look back for non-full period sigma
                                       alpha_VaR_ES = 0.01) { # define the alpha percentile for the VaR / ES calculation
  t_0 <- min(input_LC$years) # define the start year t0
  t_n <- max(input_LC$years) # define the last available historic year tn
  years_range <- t_0:end_year_T # define the total year range
  n_years <- length(years_range) # define the total number of years (historic + simulated)
  
  # drift calculation
  if (drift_type == "full") {
    drift <- (tail(input_LC$k_t, 1) - head(input_LC$k_t, 1)) / (length(input_LC$k_t) - 1)
  } else if (drift_type == "window") {
    drift <- (tail(input_LC$k_t, 1) - input_LC$k_t[length(input_LC$k_t) - drift_window + 1]) / (drift_window - 1)
  } else {
    stop("Invalid drift_type: Use 'full' or 'window'.")
  }
  
  
  # Sigma calculation
  if (sigma_type == "full") {
    sigma <- sum((diff(input_LC$k_t) - drift)^2) / (length(input_LC$k_t) - 1)
  } else if (sigma_type == "window") {
    if (length(input_LC$k_t) >= sigma_window) {
      k_t_window <- tail(input_LC$k_t, sigma_window)
      sigma <- sum((diff(k_t_window) - drift)^2) / (sigma_window - 1)
    } else {
      sigma <- sum((diff(input_LC$k_t) - drift)^2) / (length(input_LC$k_t) - 1)
      warning("Not enough data for the specified sigma window. Falling back to 'full' calculation.")
    }
  } else {
    stop("Invalid sigma_type: Use 'full' or 'window'.")
  }
  
  # preparing for parallel processing to speed up the function
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, varlist = c("input_LC", "end_year_T", "n_simulations", "drift_type", "drift_window", "sigma_type","sigma_window","drift", "sigma", "t_0", "t_n", "n_years"), envir = environment())
  
  results <- parLapply(cl, 1:n_simulations, function(sim) {
    set.seed(sim)
    sim_k_t <- numeric(n_years)
    sim_k_t[1:length(input_LC$k_t)] <- input_LC$k_t
    
    for (i in (length(input_LC$k_t) + 1):n_years) {
      sim_k_t[i] <- sim_k_t[i - 1] + drift + rnorm(1, mean = 0, sd = sqrt(sigma))
    }
    
    simulated_mx <- matrix(exp(input_LC$a_x + outer(input_LC$b_x, sim_k_t[(t_n - t_0 + 2):(end_year_T - t_0 + 1)])), nrow = length(input_LC$age_categories), ncol = (end_year_T - t_n))
    colnames(simulated_mx) <- (t_n + 1):end_year_T
    rownames(simulated_mx) <- input_LC$age_categories
    
    return(list(kt = sim_k_t, mx = simulated_mx))
  })
  
  stopCluster(cl)
  
  # Prepare plot_data matrix
  sim_kt_matrix <- matrix(nrow = n_years, ncol = n_simulations)
  for (i in 1:n_simulations) {
    sim_kt_matrix[, i] <- results[[i]]$kt
  }
  
  # Include observed k_t values for the period they cover
  observed_kt_extended <- rep(NA, n_years)
  observed_kt_extended[1:length(input_LC$k_t)] <- input_LC$k_t
  
  plot_data <- cbind(Observed = observed_kt_extended, sim_kt_matrix)
  rownames(plot_data) <- years_range
  colnames(plot_data)[2:ncol(plot_data)] <- c(1:n_simulations)
  
  stats <- list(
    p05 = apply(plot_data[, -1], 1, function(x) quantile(x, probs = 0.05, na.rm = TRUE)),
    p25 = apply(plot_data[, -1], 1, function(x) quantile(x, probs = 0.25, na.rm = TRUE)),
    mean = rowMeans(plot_data[, -1], na.rm = TRUE),
    median = apply(plot_data[, -1], 1, function(x) median(x, na.rm = TRUE)),
    p75 = apply(plot_data[, -1], 1, function(x) quantile(x, probs = 0.75, na.rm = TRUE)),
    p95 = apply(plot_data[, -1], 1, function(x) quantile(x, probs = 0.95, na.rm = TRUE)),
    VaR = apply(plot_data[, -1], 1, function(x) quantile(x, probs = alpha_VaR_ES)),
    ES = apply(plot_data[, -1], 1, function(x) quantile(x, probs = (alpha_VaR_ES/2))),
    min = apply(plot_data[, -1], 1, function(x) min(x))
  )
  
  return(list(simulations = results, statistics = stats, plot_data = plot_data))
}

# combined function of modeling & pricing calculations
pricing_calculation <- function(t_0_year = 2023, # starting year of pricing modeling
                                maturity = 10, # maturity of the S-forward contract
                                l_0 = 10, # initial cohort size
                                x_0 = 65, # initial age of cohort
                                mx_gender = mxM, # selected mx Table (either "mxM" or "mxF")
                                CoC = 0.06, # cost of capital requirements of SSt
                                rf_input = SST_input_data, # Risk free rate input
                                input_n_simulations = 100, # number of Monte carlo iterations
                                input_drift_type = "window", # window approach to drift calculation (alternative is "full")
                                input_drift_window = 20, # window size if input_drift_type = "window"
                                input_sigma_type = "window", # window approacht to variance calculation (alternative is "full")
                                input_sigma_window = 20, # window size if input_sigma_type = "window"
                                input_alpha_VaR_ES = 0.01){ # quantile to get VaR (and half of this to get ES)
  
  # Maturity date t = T
  T_maturity = t_0_year + maturity 
  
  # selects the historic risk-free rates available at the start of t_0_year
  # e.g. if t_0_year = 2019, we take the historic rates from December 2018
  rf_data <- rf_input %>%
    mutate(year = year(Date), month = month(Date)) %>% 
    filter (year == (t_0_year -1) & month == 12) %>% 
    select(-c("year","month")) %>% 
    # Use gather to convert to long format
    gather(key = "maturity", value = "rf", -Date) %>%
    # Convert maturity time to numeric
    mutate(maturity = as.numeric(maturity))
  
  # get the interpolated yield curve to calculate forwards at each year
  interpolated_data <- rf_data %>%
    complete(maturity = full_seq(maturity, period = 1)) %>%
    mutate(rf = if_else(is.na(rf), spline(rf_data$maturity, rf_data$rf, xout = maturity)$y, rf)) %>%
    arrange(maturity) %>%
    mutate(forward_rate = case_when(
      maturity == 1 ~ rf, # For the first row, forward rate is the same as the spot rate
      TRUE ~ ((1 + rf)^(maturity) / lag((1 + rf))^(maturity-1)) - 1 # get the forward rates
    )) %>%
    mutate(d_i_iplus1 = 1 / (1 + forward_rate)) # Calculate the discount factor
  
  
  # fit LC based on input mx table
  lc_output_mx <- mx_to_LC_model(mx_gender)
  
  # get simulated k_t based on LC fit
  sim_output_mx <- LC_to_simulation_optimized(input_LC = lc_output_mx, 
                                              end_year_T = T_maturity, 
                                              n_simulations = input_n_simulations,
                                              drift_type = input_drift_type,
                                              drift_window = input_drift_window,
                                              sigma_type = input_sigma_type,
                                              sigma_window = input_sigma_window,
                                              alpha_VaR_ES = input_alpha_VaR_ES)
  
  # median/mean k_t timeseries for the full forecast horizon (including historic) used as Best Estimate
  k_t_BE <- sim_output_mx$statistics$median
  names(k_t_BE) <- NULL
  
  mx_BE <- generate_mx_table(k_t = k_t_BE, a_x = lc_output_mx$a_x, b_x = lc_output_mx$b_x)
  colnames(mx_BE) <- c(min(lc_output_mx$years):T_maturity) # rename the columns to years
  rownames(mx_BE) <- lc_output_mx$age_categories
  
  # turn the mx table into a px table
  px_BE <- apply(mx_BE, c(1, 2), function(x) exp(-x))
  
  # VaR k_t timeseries for the full forecast horizon (including historic) used as worst case for insurer
  k_t_VaR <- sim_output_mx$statistics$VaR
  names(k_t_VaR) <- NULL
  
  mx_VaR <- generate_mx_table(k_t = k_t_VaR, a_x = lc_output_mx$a_x, b_x = lc_output_mx$b_x)
  colnames(mx_VaR) <- c(min(lc_output_mx$years):T_maturity) # rename the columns to years
  rownames(mx_VaR) <- lc_output_mx$age_categories
  
  # turn the mx table into a px table
  px_VaR <- apply(mx_VaR, c(1, 2), function(x) exp(-x))
  
  
  # ES k_t timeseries for the full forecast horizon (including historic) used as worst case for insurer
  k_t_ES <- sim_output_mx$statistics$ES
  names(k_t_ES) <- NULL
  
  mx_ES <- generate_mx_table(k_t = k_t_ES, a_x = lc_output_mx$a_x, b_x = lc_output_mx$b_x)
  colnames(mx_ES) <- c(min(lc_output_mx$years):T_maturity) # rename the columns to years
  rownames(mx_ES) <- lc_output_mx$age_categories
  
  # turn the mx table into a px table
  px_ES <- apply(mx_ES, c(1, 2), function(x) exp(-x))
  
  
  # get the relevant probabilities given the maturity of the endowment and the cohort start age
  gt_px_BE <- diag(px_BE[(x_0+1):(x_0+maturity+1),as.character(t_0_year:T_maturity)])
  gt_px_VaR <- diag(px_VaR[(x_0+1):(x_0+maturity+1),as.character(t_0_year:T_maturity)])
  gt_px_ES <- diag(px_ES[(x_0+1):(x_0+maturity+1),as.character(t_0_year:T_maturity)])
  
  
  # create helper tibble for calculation
  helper_px_BE <- tibble(t = 0:maturity, year = t_0_year:T_maturity, px = gt_px_BE)
  helper_px_VaR <- tibble(t = 0:maturity, year = t_0_year:T_maturity, px = gt_px_VaR)
  helper_px_ES <- tibble(t = 0:maturity, year = t_0_year:T_maturity, px = gt_px_ES)
  
  
  ### 4.2.3 Calculation using the premium and spread formula ----
  
  # get the Tp^x,0 = BE of probability in t = 0 that individual is still alive at t = T
  
  # BE_T_px_0 <- prod(helper_px_BE$px)
  BE_T_px_0 <- prod(helper_px_BE[helper_px_BE$t %in% c(0:(maturity-1)),]$px)
  
  # get the sum of the probabilities
  sum_px <- 0
  for (i in 0:(maturity-1)){
    # BE_i_px_0 <- prod(helper_px_BE[helper_px_BE$t %in% 0:i,]$px)
    # ES_Tminus1_px_i <- prod(helper_px_ES[helper_px_ES$t %in% i:maturity,]$px) 
    # discount_i_iplus1 <- interpolated_data[interpolated_data$maturity == (i+1),]$d_i_iplus1
    
    BE_i_px_0 <- prod(helper_px_BE[helper_px_BE$t %in% c(0:(i-1)),]$px)
    ES_Tminus1_px_i <- prod(helper_px_ES[helper_px_ES$t %in% c(i:(maturity-1)),]$px) 
    discount_i_iplus1 <- interpolated_data[interpolated_data$maturity == (i+1),]$d_i_iplus1
    
    sum_px = sum_px + BE_i_px_0*ES_Tminus1_px_i*discount_i_iplus1
  }
  
  # get the sum of the discount factors
  sum_df <- 0
  for (i in 0:(maturity-1)){
    sum_df = sum_df + interpolated_data[interpolated_data$maturity == (i+1),]$d_i_iplus1
  }
  
  # get the final maximum price (pi) and maximum spread
  pi_max <- CoC*((sum_px/BE_T_px_0)-sum_df)
  spread_max <- log((CoC*sum_px+BE_T_px_0*(1-CoC*sum_df))/(BE_T_px_0))*(1/maturity)
  
  pricing_output <- list(pi = pi_max, spread = spread_max)
  return(pricing_output)
}


# 1. Loading Data ----

# this first data processing step is crucial to ensure that the starting mx
# table is correct for further processing through the manual LC process or the functions
# year_range needs to be in line with the forecasting start (e.g. we want 
# to look at historical data from 1935 until 2018 and start forecasting in 2019)

# selecting date-range (depends on data availability)
year_range <- seq.default(from = 1935,to = 2018, by = 1)

## 1.1 Data from Human Mortality Database ----

# select country (countries need to be stored in 00_Data folder structure)
selected_country <- "Switzerland"

# load complete HMD data (1x1 life tables)
# for male population
HMD_male <- read.table(paste0("./00_Data/00_Mortality Data/",selected_country,"/00_HMD/mltper_1x1.txt"), skip = 2, header=TRUE)
HMD_male <- HMD_male %>% 
  mutate(Age = as.numeric(str_extract(HMD_male$Age, "[0-9]+"))) %>% 
  filter(Age != 110) # since mx of 110 is always 1

# for female population
HMD_female <- read.table(paste0("./00_Data/00_Mortality Data/",selected_country,"/00_HMD/fltper_1x1.txt"), skip = 2, header=TRUE)
HMD_female <- HMD_female %>% 
  mutate(Age = as.numeric(str_extract(HMD_female$Age, "[0-9]+"))) %>% 
  filter(Age != 110) # since mx of 110 is always 1


# select mx (central death rates) for males and females
mxM <- HMD_male %>%
  select(c("Age","Year","mx")) %>% 
  pivot_wider(names_from = "Year", values_from = "mx") %>% 
  select(c("Age",as.character(year_range)))

mxF <- HMD_female %>%
  select(c("Age","Year","mx")) %>% 
  pivot_wider(names_from = "Year", values_from = "mx") %>% 
  select(c("Age",as.character(year_range)))


# slight data cleaning since some mx are equal to 0
which(colSums(mxM[,-1] == 0) != 0)
which(colSums(mxF[,-1] == 0) != 0)

# ensure that mx is not exactly equal to zero since this would lead to
# -infinity after taking the log 

# apply the replace_zeros() function to the dataframe 
# (when mode is "add" you can define the small_number added to the 0s)
mxM[,-1] <- replace_zeros(mxM[,-1], mode = "average", small_number = 1e-10)
mxF[,-1] <- replace_zeros(mxF[,-1], mode = "average", small_number = 1e-10)


# note: In sections 2 & 3 of the code the steps of Lee Carter Fitting and Monte-
#       Carlo simulation of k_t are performed in manual steps instead of using the
#       the functions defined in section 0.1 & 0.2. This is done to create 
#       visualizations that are used in the research paper. However, to recreate
#       the table in Chapter 4 of the research paper, one could also skip the code
#       of sections 2 & 3 and only run from section 4.



# 2. Lee-Carter (Manual Calculation for plots) ----

## 2.1 fitting Lee-Carter Model (first step estimation) ----

# the calculations here are performed outside of the main function mc_to_LC_model()
# to ensure that we can create some nice data visualizations

### 2.1.1 Preparing Matrix and fitting a_x ----

# get the ln(m_x) matrix
log_mxM <- log(mxM[,-1])
log_mxF <- log(mxF[,-1])

# get the mean log-rate for each age-group (over all t) (avg(ln(m_x)))
mean_log_mxM <- rowMeans(log_mxM) 
mean_log_mxF <- rowMeans(log_mxF) 

# subtract this average from each age-group for each year
for (i in 1:nrow(log_mxM)){
  log_mxM[i,] <- log_mxM[i,] - mean_log_mxM[i]
}

for (i in 1:nrow(log_mxF)){
  log_mxF[i,] <- log_mxF[i,] - mean_log_mxF[i]
}

# final estimates of a_x (the mean over time of the log mortality rate for the age groups)
a_x_M <- mean_log_mxM
a_x_F <- mean_log_mxF

# plotting a_x
a_x_plot_data <- tibble(Age = c(0:(nrow(log_mxF)-1)), Male_a_x = a_x_M, Female_a_x = a_x_F)

ggplot(aes(x = Age), data = a_x_plot_data) +
  geom_line(aes(y = Male_a_x, color = "Male"))+
  geom_line(aes(y = Female_a_x, color = "Female"))+
  theme_minimal()+
  theme(
    panel.background = element_blank(), # Make panel background transparent
    panel.grid.major.x = element_blank(), # Remove major x grid lines
    panel.grid.minor.x = element_blank(), # Remove minor x grid lines
    panel.grid.major.y = element_line(color = "grey"), # Customize or keep major y grid lines
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.8) # Add black box around the plot
  )+
  labs(x = "age-category", y = "ax", color = "Gender")+
  scale_color_manual(values = c("Male" = "#003049", "Female" = "#d62828"))


### 2.1.2 Singular Value Decomposition to fit b_x & k_t ----

# for Male Population
d_M <- svd(as.matrix(log_mxM), 
         nv=nrow(log_mxM), 
         nu=ncol(log_mxM))
# due to the sign ambiguity we can rescale b_x and k_t (the -1 cancel eachother)
b_x_M <- d_M$u[,1]*(-1)

# k_t calculation based on Girosi & King (2007, p. 4)
k_t_M <- as.vector(b_x_M %*% as.matrix(log_mxM))


# for Female Population
d_F <- svd(as.matrix(log_mxF), 
           nv=nrow(log_mxF), 
           nu=ncol(log_mxF))
# due to the sign ambiguity we can rescale b_x and k_t (the -1 cancel eachother)
b_x_F <- d_F$u[,1]*(-1)


# k_t calculation based on Girosi & King (2007, p. 4)
k_t_F <- as.vector(b_x_F %*% as.matrix(log_mxF))


# check if the constraint holds
# unity constraint of b_x^2 (needs to sum to 1)
sum(b_x_M^2)
sum(b_x_F^2)

# zero constraint of k_t (needs to sum to 0)
sum(k_t_M)
sum(k_t_F)


# graphical comparison
# b_x
b_x_plot_data <- tibble(Age = c(0:(nrow(log_mxF)-1)), Male_b_x = b_x_M, Female_b_x = b_x_F)

ggplot(aes(x = Age), data = b_x_plot_data) +
  geom_line(aes(y = Male_b_x, color = "Male"))+
  geom_line(aes(y = Female_b_x, color = "Female"))+
  theme_minimal()+
  theme(
    panel.background = element_blank(), # Make panel background transparent
    panel.grid.major.x = element_blank(), # Remove major x grid lines
    panel.grid.minor.x = element_blank(), # Remove minor x grid lines
    panel.grid.major.y = element_line(color = "grey"), # Customize or keep major y grid lines
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.8) # Add black box around the plot
  )+
  labs(x = "age-category", y = "bx", color = "Gender")+
  scale_color_manual(values = c("Male" = "#003049", "Female" = "#d62828"))

# k_t
k_t_plot_data <- tibble(Year = year_range, Male_k_t = k_t_M, Female_k_t = k_t_F)

ggplot(aes(x = Year), data = k_t_plot_data) +
  geom_line(aes(y = Male_k_t, color = "Male"))+
  geom_line(aes(y = Female_k_t, color = "Female"))+
  geom_vline(xintercept = 1995, linetype = "dashed", color = "#003049")+
  geom_vline(xintercept = 1955, linetype = "dashed", color = "#d62828")+
  labs(x = "year", y = "kt", color = "Gender")+
  theme(
    panel.background = element_blank(), # Make panel background transparent
    panel.grid.major.x = element_blank(), # Remove major x grid lines
    panel.grid.minor.x = element_blank(), # Remove minor x grid lines
    panel.grid.major.y = element_line(color = "grey"), # Customize or keep major y grid lines
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.8) # Add black box around the plot
  )+
  labs(x = "year", y = "kt", color = "Gender")+
  scale_color_manual(values = c("Male" = "#003049", "Female" = "#d62828"))


# add sigma plots
# Males
var_plot_data <- data.frame(year = year_range, k_t = k_t_M)
var_plot_data$diff_k_t <- c(NA, diff(var_plot_data$k_t))  # NA for the first difference

# Calculate rolling variance of first differences
window_sizes <- c(5, 10, 20) # Define the window sizes
for (w in window_sizes) {
  roll_var_name <- paste0('rolling_', w, 'y_var')
  var_plot_data[[roll_var_name]] <- rollapply(var_plot_data$diff_k_t, width = w, FUN = sd, fill = NA, align = 'right', na.rm = TRUE)
}

# Transform the data to long format for plotting
df_long <- var_plot_data %>%
  gather(key = "type", value = "var", starts_with("roll"), na.rm = TRUE)

# Create the plot for the rolling variance of the first differences
varplot_male <-ggplot(df_long, aes(x = year, y = var, group = type, color = type)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Year", y = "Variance",title = "Male") +
  scale_color_manual(values = c("rolling_5y_var" = "#00798c", "rolling_10y_var" = "#d1495b", "rolling_20y_var" = "#edae49")) +
  theme(legend.title = element_blank())+
  theme(
    panel.background = element_blank(), # Make panel background transparent
    panel.grid.major.x = element_blank(), # Remove major x grid lines
    panel.grid.minor.x = element_blank(), # Remove minor x grid lines
    panel.grid.major.y = element_line(color = "grey"), # Customize or keep major y grid lines
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.8) # Add black box around the plot
  )+
  scale_y_continuous(limits = c(0.1, 1))


# Females
var_plot_data <- data.frame(year = year_range, k_t = k_t_F)
var_plot_data$diff_k_t <- c(NA, diff(var_plot_data$k_t))  # NA for the first difference

# Calculate rolling variance of first differences
window_sizes <- c(5, 10, 20) # Define the window sizes
for (w in window_sizes) {
  roll_var_name <- paste0('rolling_', w, 'y_var')
  var_plot_data[[roll_var_name]] <- rollapply(var_plot_data$diff_k_t, width = w, FUN = sd, fill = NA, align = 'right', na.rm = TRUE)
}

# Transform the data to long format for plotting
df_long <- var_plot_data %>%
  gather(key = "type", value = "var", starts_with("roll"), na.rm = TRUE)

# Create the plot for the rolling variance of the first differences
varplot_female <- ggplot(df_long, aes(x = year, y = var, group = type, color = type)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Year", y = "", title = "Female") +
  scale_color_manual(values = c("rolling_5y_var" = "#00798c", "rolling_10y_var" = "#d1495b", "rolling_20y_var" = "#edae49")) +
  theme(legend.title = element_blank())+
  theme(
    panel.background = element_blank(), # Make panel background transparent
    panel.grid.major.x = element_blank(), # Remove major x grid lines
    panel.grid.minor.x = element_blank(), # Remove minor x grid lines
    panel.grid.major.y = element_line(color = "grey"), # Customize or keep major y grid lines
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.8) # Add black box around the plot
  )+
  scale_y_continuous(limits = c(0.1, 1))


(varplot_male + varplot_female) + 
  plot_layout(guides = 'collect') & # Collect and share the legend
  theme(legend.position = "bottom")



# 3. Forecasting using the fitted Lee-Carter model ----

## 3.1 Parameter Estimation ----

# general simulation assumptions

n_simulations = 100 # number of iterations of the monte carlo simulation

t_0 = head(year_range,1) # first observation
t_n = tail(year_range,1) # last date for which we have info (e.g. end of year_range)
end_year_T = 2060 # projection limit


### 3.1.1 drift term ----
# calculated as average drift = (k_T - k_1)/(T-1) (Girosi & King, 2007, p. 5)
drift_M <- (k_t_M[length(k_t_M)] - k_t_M[1]) / (length(k_t_M) - 1)
drift_F <- (k_t_F[length(k_t_F)] - k_t_F[1]) / (length(k_t_F) - 1)

# define the starting point for a shorter range
start_short_drift_M <- 1995
start_short_drift_F <- 1955

drift_M_short_sample <- (k_t_M[length(k_t_M)] - k_t_M[start_short_drift_M-t_0+1]) / (length(k_t_M) - (start_short_drift_M-t_0) - 1)
drift_F_short_sample <- (k_t_F[length(k_t_F)] - k_t_F[start_short_drift_F-t_0+1]) / (length(k_t_F) - (start_short_drift_F-t_0) - 1)


### 3.1.2 volatility parameter ----
# according to Girosi & King (2007, p. 5)
# 1/(T-1)*Sum over t=1 until T-1 ((k_t+1 - k_t - theta^hat)^2), where theta is the drift
sigma_M <- sum((diff(tail(k_t_M,t_n-start_short_drift_M)) - drift_M_short_sample)^2) / (t_n-start_short_drift_M)
sigma_F <- sum((diff(tail(k_t_F,t_n-start_short_drift_F)) - drift_F_short_sample)^2) / (t_n-start_short_drift_F)


### 3.1.3 Check Correlation between Male and Female mortality rates
cor(k_t_M,k_t_F)


### 3.1.4 Random Walk Forecast ----

#### 3.1.4.1 Setup of Simulation tables ----

# simulation tables
full_sim_k_t_M <- tibble(year= seq(t_0,end_year_T,1))
full_sim_k_t_F <- tibble(year= seq(t_0,end_year_T,1))

# pre-allocate a simulation table where the columns represent the iterations of 
# monte-carlo simulations (first rows represent the existing k_t data)
k_t_M_proj <- c(k_t_M, numeric(end_year_T - t_n))

for (i in 1:(n_simulations)){
  full_sim_k_t_M[,i+1] <- k_t_M_proj
  colnames(full_sim_k_t_M)[i+1] <- i
}


k_t_F_proj <- c(k_t_F, numeric(end_year_T - t_n))

for (i in 1:(n_simulations)){
  full_sim_k_t_F[,i+1] <- k_t_F_proj
  colnames(full_sim_k_t_F)[i+1] <- i
}


#### 3.1.4.2 Simulation of k_t ----

# ensure that the correct drift is chosen
for (j in 2:(n_simulations+1)){ # for the number of simulations we want to perform
  for (i in (t_n-t_0+2):(end_year_T-t_0+1)){ # for all rows where we want to simulate
    full_sim_k_t_M[i,j] <- full_sim_k_t_M[i-1,j] + drift_M_short_sample + rnorm(1,mean=0, sd=sqrt(sigma_M))
  }
}


for (j in 2:(n_simulations+1)){ # for the number of simulations we want to perform
  for (i in (t_n-t_0+2):(end_year_T-t_0+1)){ # for all rows where we want to simulate
    full_sim_k_t_F[i,j] <- full_sim_k_t_F[i-1,j] + drift_F_short_sample + rnorm(1,mean=0, sd=sqrt(sigma_F))
  }
}


#### 3.1.4.3 creation of stochastic mortality tables ----

# recreate stochastic death Tables (GT / PT)
# for each simulation we will have a table

all_simulated_mx_M <- list() # pre-allocate list to store the stochastic GTs
all_simulated_mx_F <- list() # pre-allocate list to store the stochastic GTs


all_simulated_px_M <- list() # pre-allocate list to store the stochastic GTs
all_simulated_px_F <- list() # pre-allocate list to store the stochastic GTs



# for males
for (n in 1:n_simulations){ # for all simulations
  
  # for each simulation select the generated k_t vector
  simulated_k_t_M <- full_sim_k_t_M[which(full_sim_k_t_M$year == (t_n+1)):which(full_sim_k_t_M$year == end_year_T),c(1,n+1)]
  
  # pre-allocate a new mxM matrix that will be filled in for each age-group and each year
  simulated_mxM <- tibble(Age = mxM$Age)
  
  
  for (t in simulated_k_t_M$year){ # for each year in a given stochastic GT
    
    temp_sim_k_t <- simulated_k_t_M[which(simulated_k_t_M$year == t),2][[1]]
    
    temp_mx_year_vector <- exp(a_x_M + b_x_M*temp_sim_k_t)
    
    simulated_mxM <- cbind(simulated_mxM, temp_mx_year_vector)
    
    names(simulated_mxM)[ncol(simulated_mxM)] <- t
    
  }
  
  
  # conversion from mxM to pxm (qx = 1-exp(-mx) & px = 1-qx --> px = exp(-mx))
  simulated_pxM <- simulated_mxM
  simulated_pxM[,-1] <- apply(simulated_pxM[,-1],c(1,2), function(x) exp(-x))
  
  all_simulated_mx_M[[n]] <- simulated_mxM # store the simulated GT to the table
  all_simulated_px_M[[n]] <- simulated_pxM # store the simulated GT to the table
  
}

# for females
for (n in 1:n_simulations){ # for all simulations
  
  # for each simulation select the generated k_t vector
  simulated_k_t_F <- full_sim_k_t_F[which(full_sim_k_t_F$year == (t_n+1)):which(full_sim_k_t_F$year == end_year_T),c(1,n+1)]
  
  # pre-allocate a new mxF matrix that will be filled in for each age-group and each year
  simulated_mxF <- tibble(Age = mxF$Age)
  
  
  for (t in simulated_k_t_F$year){ # for each year in a given stochastic GT
    
    temp_sim_k_t <- simulated_k_t_F[which(simulated_k_t_F$year == t),2][[1]]
    
    temp_mx_year_vector <- exp(a_x_F + b_x_F*temp_sim_k_t)
    
    simulated_mxF <- cbind(simulated_mxF, temp_mx_year_vector)
    
    names(simulated_mxF)[ncol(simulated_mxF)] <- t
    
  }
  
  
  # conversion from mxF to pxm (qx = 1-exp(-mx) & px = 1-qx --> px = exp(-mx))
  simulated_pxF <- simulated_mxF
  simulated_pxF[,-1] <- apply(simulated_pxF[,-1],c(1,2), function(x) exp(-x))
  
  all_simulated_mx_F[[n]] <- simulated_mxF # store the simulated GT to the table
  all_simulated_px_F[[n]] <- simulated_pxF # store the simulated GT to the table
  
}

#### 3.1.4.3 visualizing the simulated k_t paths ----

# for males
# some summary statistics
# mean
full_sim_k_t_M$mean <- apply(full_sim_k_t_M[,2:ncol(full_sim_k_t_M)], 1, mean)
# 5% percentile
full_sim_k_t_M$perc_5<- apply(full_sim_k_t_M[,2:ncol(full_sim_k_t_M)], 1, function(x) quantile(x, probs = 0.05))
# 50% percentile (median)
full_sim_k_t_M$perc_50 <- apply(full_sim_k_t_M[,2:ncol(full_sim_k_t_M)], 1, function(x) quantile(x, probs = 0.50))
# 95% percentile
full_sim_k_t_M$perc_95 <- apply(full_sim_k_t_M[,2:ncol(full_sim_k_t_M)], 1, function(x) quantile(x, probs = 0.95))

# change to a long format for ggplot
k_t_M_simulation_data_long <- pivot_longer(full_sim_k_t_M, cols = -year, names_to = "observation", values_to = "value")

# Plot the data using ggplot2
highlighted_obs <- c("mean", "perc_5", "perc_50", "perc_95") # define 

ggplot(data = k_t_M_simulation_data_long, aes(x = year, y = value, color = observation)) +
  geom_line(data = k_t_M_simulation_data_long %>% filter(year <= t_n), aes(color = observation), alpha = 0.2, linewidth = 0.5) +  # Draw lighter lines before t_n
  geom_line(data = k_t_M_simulation_data_long %>% filter(year >= t_n, !observation %in% highlighted_obs), aes(color = observation), alpha = 0.2, linewidth = 0.5) +  # Draw lighter lines for the rest after t_n
  geom_line(data = k_t_M_simulation_data_long %>% filter(year >= t_n, observation %in% highlighted_obs), aes(color = observation), alpha = 1, linewidth = 1.2) +  # Draw bolder lines for "mean" and "perc_x" after t_n
  geom_vline(xintercept = t_n, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(x = "Year",
       y = "kt")+
  theme(
    legend.position = "none",
    panel.background = element_blank(), # Make panel background transparent
    panel.grid.major.x = element_blank(), # Remove major x grid lines
    panel.grid.minor.x = element_blank(), # Remove minor x grid lines
    panel.grid.major.y = element_line(color = "grey"), # Customize or keep major y grid lines
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.8) # Add black box around the plot
  )


# 4. Simulation of a hypothetical insurer's longevity risks ----

# background: 
# the hypothetical insurer is exposed to longevity risk through a portfolio of
# pure endowments that pay policyholders a lump sum (fixed amount of 1 monetary unit)
# on the condition that at time t = T_final the policy holder is still alive.

# The insurer is located in Switzerland and subject to Switzerland's 
# Swiss Solvency Test (SST) regulation and thus must have a Market Value Margin (MVM) on 
# its Solvency Capital Requirements (SCR) that arise from the risks on its balance sheet.

# main assumptions (based on Levantesi & Menzietti (2017)):
# no counter-party risk introduced by the S-forwards
# no basis risk since we assume same reference population for insurer and for the forward


## 4.1 model set-up ----

# start and end of the endowment contract
t_0_year = 2019 # ensure that this matches with the last available historical year in the Lee-Carter calibration
maturity = 5 # years until maturity
T_maturity = t_0_year + maturity # Maturity date t = T

# information on cohort
l_0 = 10 # initial cohort size
x_0 = 65 # initial age of the members of the cohort
mx_gender = mxM # chose "mxM" or "mxF" to look at male or female cohort

# financial data
CoC = 0.06 # Cost of Capital on SCR based on SST assumption


### 4.1.1 Discount rates (based on FINMA SST Input Data): ----
# source: https://www.finma.ch/de/~/media/finma/dokumente/dokumentencenter/myfinma/2ueberwachung/sst/sst-inputdaten.xlsx?sc_lang=de&hash=01307EAB5A49D579BA21929885209FDD
SST_input_data <- read_excel("./00_Data/01_SST/SST Inputdaten.xlsx",
                             sheet = "History_SNB",
                             range = "B4:O257")
colnames(SST_input_data)[1] <- "Date"
SST_input_data[,-1] <- SST_input_data[,-1]/100 # turn the rates into decimal values

# selects the historic risk-free rates available at the start of t_0_year
# e.g. if t_0_year = 2019, we take the historic rates from December 2018
rf_data <- SST_input_data %>%
  mutate(year = year(Date), month = month(Date)) %>% 
  filter (year == (t_0_year -1) & month == 12) %>% 
  select(-c("year","month")) %>% 
  # Use gather to convert to long format
  gather(key = "maturity", value = "rf", -Date) %>%
  # Convert maturity time to numeric
  mutate(maturity = as.numeric(maturity))

# get the interpolated yield curve to calculate forwards at each year
interpolated_data <- rf_data %>%
  complete(maturity = full_seq(maturity, period = 1)) %>%
  mutate(rf = if_else(is.na(rf), spline(rf_data$maturity, rf_data$rf, xout = maturity)$y, rf)) %>%
  arrange(maturity) %>%
  mutate(forward_rate = case_when(
    maturity == 1 ~ rf, # For the first row, forward rate is the same as the spot rate
    TRUE ~ ((1 + rf)^(maturity) / lag((1 + rf))^(maturity-1)) - 1 # get the forward rates
  )) %>%
  mutate(d_i_iplus1 = 1 / (1 + forward_rate)) # Calculate the discount factor


# plot the interpolated spot and forward rates curve
ggplot(interpolated_data, aes(x = maturity)) +
  geom_point(aes(y = rf, color = 'Spot Rate')) +
  geom_point(aes(y = forward_rate, color = 'Forward Rate')) +
  geom_line(aes(x = maturity, y = rf, color = 'Spot Rate')) +
  geom_line(aes(x = maturity, y = forward_rate, color = 'Forward Rate')) +
  theme_minimal() +
  theme(
    panel.background = element_blank(), # Make panel background transparent
    panel.grid.major.x = element_blank(), # Remove major x grid lines
    panel.grid.minor.x = element_blank(), # Remove minor x grid lines
    panel.grid.major.y = element_line(color = "grey"), # Customize or keep major y grid lines
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.8) # Add black box around the plot
  )+
  labs(title = t_0_year, x = "Maturity", y = "Rate") +
  scale_color_manual(values = c('Spot Rate' = '#345995', 'Forward Rate' = '#eac435')) +
  guides(color = guide_legend(title = ""))


## 4.2 Pricing  ----

### 4.2.1 Forecasting k_t using the functions ----

# fit LC based on input mx table
lc_output_mx <- mx_to_LC_model(mx_gender)

# get simulated k_t based on LC fit
sim_output_mx <- LC_to_simulation_optimized(input_LC = lc_output_mx, 
                                            end_year_T = T_maturity, 
                                            n_simulations = 2000,
                                            drift_type = "window",
                                            drift_window = 20,
                                            sigma_type = "full",
                                            sigma_window = 20,
                                            alpha_VaR_ES = 0.01)


# median/mean k_t timeseries for the full forecast horizon (including historic) used as Best Estimate
k_t_BE <- sim_output_mx$statistics$median
names(k_t_BE) <- NULL

mx_BE <- generate_mx_table(k_t = k_t_BE, a_x = lc_output_mx$a_x, b_x = lc_output_mx$b_x)
colnames(mx_BE) <- c(min(lc_output_mx$years):T_maturity) # rename the columns to years
rownames(mx_BE) <- lc_output_mx$age_categories

# turn the mx table into a px table
px_BE <- apply(mx_BE, c(1, 2), function(x) exp(-x))


# VaR k_t timeseries for the full forecast horizon (including historic) used as worst case for insurer
k_t_VaR <- sim_output_mx$statistics$VaR
names(k_t_VaR) <- NULL

mx_VaR <- generate_mx_table(k_t = k_t_VaR, a_x = lc_output_mx$a_x, b_x = lc_output_mx$b_x)
colnames(mx_VaR) <- c(min(lc_output_mx$years):T_maturity) # rename the columns to years
rownames(mx_VaR) <- lc_output_mx$age_categories

# turn the mx table into a px table
px_VaR <- apply(mx_VaR, c(1, 2), function(x) exp(-x))


# ES k_t timeseries for the full forecast horizon (including historic) used as worst case for insurer
k_t_ES <- sim_output_mx$statistics$ES
names(k_t_ES) <- NULL

mx_ES <- generate_mx_table(k_t = k_t_ES, a_x = lc_output_mx$a_x, b_x = lc_output_mx$b_x)
colnames(mx_ES) <- c(min(lc_output_mx$years):T_maturity) # rename the columns to years
rownames(mx_ES) <- lc_output_mx$age_categories

# turn the mx table into a px table
px_ES <- apply(mx_ES, c(1, 2), function(x) exp(-x))


### 4.2.2 Get the survival probabilities ----

# get the relevant probabilities given the maturity of the endowment and the cohort start age
gt_px_BE <- diag(px_BE[(x_0+1):(x_0+maturity+1),as.character(t_0_year:T_maturity)])
gt_px_VaR <- diag(px_VaR[(x_0+1):(x_0+maturity+1),as.character(t_0_year:T_maturity)])
gt_px_ES <- diag(px_ES[(x_0+1):(x_0+maturity+1),as.character(t_0_year:T_maturity)])


# create helper tibble for calculation
helper_px_BE <- tibble(t = 0:maturity, year = t_0_year:T_maturity, px = gt_px_BE)
helper_px_VaR <- tibble(t = 0:maturity, year = t_0_year:T_maturity, px = gt_px_VaR)
helper_px_ES <- tibble(t = 0:maturity, year = t_0_year:T_maturity, px = gt_px_ES)


### 4.2.3 Calculation using the premium and spread formula ----

# get the Tp^x,0 = BE of probability in t = 0 that individual is still alive at t = T
BE_T_px_0 <- prod(helper_px_BE[helper_px_BE$t %in% c(0:(maturity-1)),]$px)

# get the sum of the probabilities
sum_px <- 0
for (i in 0:(maturity-1)){
  BE_i_px_0 <- prod(helper_px_BE[helper_px_BE$t %in% c(0:(i-1)),]$px)
  ES_Tminus1_px_i <- prod(helper_px_ES[helper_px_ES$t %in% c(i:(maturity-1)),]$px) 
  discount_i_iplus1 <- interpolated_data[interpolated_data$maturity == (i+1),]$d_i_iplus1
  sum_px = sum_px + BE_i_px_0*ES_Tminus1_px_i*discount_i_iplus1
}

# get the sum of the discount factors
sum_df <- 0
for (i in 0:(maturity-1)){
  sum_df = sum_df + interpolated_data[interpolated_data$maturity == (i+1),]$d_i_iplus1
}


# get the final maximum price (pi) and maximum spread
pi_max <- CoC*((sum_px/BE_T_px_0)-sum_df)
spread_max <- log((CoC*sum_px+BE_T_px_0*(1-CoC*sum_df))/(BE_T_px_0))*(1/maturity)



### 4.2.4 Pricing for different maturities and starting ages ----
# define the dimensions of the sensitivity pricing output table
maturity_input = seq(from = 5, to = 25, by = 5)
age_input = seq(from = 65, to = 80, by = 5)

pricing_sensitivity <- matrix(nrow = length(maturity_input), 
                              ncol = length(age_input))

colnames(pricing_sensitivity) <- age_input
rownames(pricing_sensitivity) <- maturity_input

pricing_sensitivity_spread <- pricing_sensitivity
pricing_sensitivity_premium <- pricing_sensitivity

# loop through all defined maturities of the S-forward and starting ages for 
# a given start year
start_years_sensitivity_table = c(2008,2019,2023)
collected_sensitivity_table_male <- list()
collected_sensitivity_table_female <- list()

# for males
# for each start year
for (t_0 in start_years_sensitivity_table){
  # for each maturity
  for (i in 1:length(maturity_input)){
    maturity_input_i <- maturity_input[i]
    # for each starting age
    for(j in 1:length(age_input)){
      age_input_i = age_input[j]
      
      pricing_output <- pricing_calculation(t_0_year = t_0, 
                                            maturity = maturity_input_i,
                                            l_0 = 10,
                                            x_0 = age_input_i,
                                            mx_gender = if (t_0 >= 2019) {
                                              mxM
                                            } else {
                                              mxM[c("Age",as.character(1935:(t_0-1)))]
                                            },
                                            CoC = 0.06,
                                            rf_input = SST_input_data,
                                            input_n_simulations = 2000,
                                            input_drift_type = "window",
                                            input_drift_window = t_0-1995,
                                            input_sigma_type = "window",
                                            input_sigma_window = t_0-1995,
                                            input_alpha_VaR_ES = 0.01)
      
      pricing_sensitivity_spread[i,j] <- pricing_output$spread
      pricing_sensitivity_premium[i,j] <- pricing_output$pi
    }
  }
  collected_sensitivity_table_male[[as.character(t_0)]] <- pricing_sensitivity_spread
}


# for females
# for each start year
for (t_0 in start_years_sensitivity_table){
  # for each maturity
  for (i in 1:length(maturity_input)){
    maturity_input_i <- maturity_input[i]
    # for each starting age
    for(j in 1:length(age_input)){
      age_input_i = age_input[j]
      
      pricing_output <- pricing_calculation(t_0_year = t_0, 
                                            maturity = maturity_input_i,
                                            l_0 = 10,
                                            x_0 = age_input_i,
                                            mx_gender = if (t_0 >= 2019) {
                                              mxF
                                            } else {
                                              mxF[c("Age",as.character(1935:(t_0-1)))]
                                            },
                                            CoC = 0.06,
                                            rf_input = SST_input_data,
                                            input_n_simulations = 2000,
                                            input_drift_type = "window",
                                            input_drift_window = t_0-1955,
                                            input_sigma_type = "window",
                                            input_sigma_window = t_0-1955,
                                            input_alpha_VaR_ES = 0.01)
      
      pricing_sensitivity_spread[i,j] <- pricing_output$spread
      pricing_sensitivity_premium[i,j] <- pricing_output$pi
    }
  }
  collected_sensitivity_table_female[[as.character(t_0)]] <- pricing_sensitivity_spread
}


collected_sensitivity_table_male
collected_sensitivity_table_female

# uncomment this step if you want to extract the sensitivity tables to excel
# Convert each matrix in the list to a data frame

# collected_sensitivity_table_male_df <- lapply(collected_sensitivity_table_male, as.data.frame)
# collected_sensitivity_table_female_df <- lapply(collected_sensitivity_table_female, as.data.frame)

# Write the list of data frames to an Excel file
# write_xlsx(collected_sensitivity_table_male_df, "male.xlsx")
# write_xlsx(collected_sensitivity_table_female_df, "female.xlsx")
