### This script file will allow the user to run the SAM model to 
### evaluate controls on T (or WUE?) for sites in the NMEG. 

### Install and load packages and functions ###

# Run this if necessary packages are not installed
if(!"rjags" %in% installed.packages()) {
  install.packages("rjags")
}
if(!"mcmcplots" %in% installed.packages()) {
  install.packages("mcmcplots")
}
if(!"tidyverse" %in% installed.packages()) {
  install.packages("tidyverse")
}
if(!"postjags" %in% installed.packages()) {
  devtools::install_github("fellmk/PostJAGS/postjags")
}
if(!"coda4dummies" %in% installed.packages()) {
  devtools::install_github("egreich/coda4dummies")
}

# Load packages
library(rjags)
load.module('dic')
library(mcmcplots)
library(tidyverse)
library(postjags)
library(coda4dummies)


# Load self-made functions
source("./scripts/SAM_initialize_function.R")
source("./scripts/SAM_function.R")

### Read in data ###
# Read in the covariate file of interest. This should be a file that has all 
# covariates for all time periods to ensure that 6-months of precipitation is 
# supplied to the model - not just growing season precipitation
dataIN_T = read.csv("./clean_data/d_B_mpj.csv") %>%
  filter(year %in% c(2017,2018,2019,2020)) # shorten df for testing
dataIN_WUE = read.csv("./clean_data/d_B_wue_mpj.csv") %>%
  filter(year %in% c(2017,2018,2019,2020))

dataIN <- left_join(dataIN_T, dataIN_WUE)

### Get initials ###
# Run the below function to get initial estimates
# Note: If you've already run the model, you do not need to do this
# Note: Do not run this function using an HPC, 
#unless you comment out the mcmcplot function in the function code
get_SAM_inits(dataIN, "mpj")

### Run the model ###
output <- SAM_WUE(dataIN, "mpj")



### Look at stuff  ###
load("./output_coda/zc_mpj.RData")
mcmcplot(zc1)
sum_tab <- coda.fast(chains=3, burn.in=0, thin=1, coda=zc1)

# Function to extract posterior means, and 2.5 and 97.5 CI quantiles
# takes a list of variable names and the coda summary
# all variables in the list MUST have the same length posterior outputs
# Similar version called coda_to_rows available in the coda4dummies package via devtools::install_github("egreich/coda4dummies")
get_coda_rows_to_cols <- function(var_list, coda_sum){
  
  sum_tb <- coda_sum[["statistics"]]
  quan_tb <- coda_sum[["quantiles"]]
  
  voi_list <- list()
  column_names <- list()
  j = 1
  
  for(i in c(1:length(var_list))){
    
    searchterm <- paste("^", var_list[i], "\\[", sep = "")
    
    # Check if there is more than instance of the variable or not
    if(length(grep(searchterm, row.names(sum_tb))) == 0){ # if we find nothing
      searchterm <- paste("^", var_list[i], sep = "") # check if there is only one instance and correct the search term
      if(length(grep(searchterm, row.names(sum_tb))) == 0){ # if we still find nothing
        print(paste("Warning: ", var_list[i], " not found in coda summary output", sep = ""))
        next
      }
    }
    
    voi_list[[j]] <- sum_tb[grep(searchterm, row.names(sum_tb)),1]
    voi_list[[j+1]] <- quan_tb[grep(searchterm, row.names(quan_tb)),1]
    voi_list[[j+2]] <- quan_tb[grep(searchterm, row.names(quan_tb)),5]
    
    column_names[[j]] <- paste("", var_list[i], sep = "")
    column_names[[j+1]] <- paste("cred2.5_", var_list[i], sep = "")
    column_names[[j+2]] <- paste("cred97.5_", var_list[i], sep = "")
    
    j = j + 3
    
  }
  suppressMessages(df <- dplyr::bind_cols(voi_list))
  colnames(df) <- column_names
  return(df)
  
}

varlist <- c("wT","wV","wSs")
weights <- get_coda_rows_to_cols(varlist, summary(zc1))
weights <- weights %>%
  rowid_to_column("lag")
varlist <- c("wP")
weightsP <- get_coda_rows_to_cols(varlist, summary(zc1))
weightsP <- weightsP %>%
  rowid_to_column("lag")
varlist <- c("wP.weekly")
weightsPweekly <- get_coda_rows_to_cols(varlist, summary(zc1))
weightsPweekly <- weightsPweekly  %>%
  rowid_to_column("lag")
varlist <- c("wP.monthly")
weightsPmonthly <- get_coda_rows_to_cols(varlist, summary(zc1))
weightsPmonthly <- weightsPmonthly %>%
  rowid_to_column("lag")

weights <- weights %>%
  pivot_longer(cols = c("wT",
                        "wV",
                        "wSs"), names_to = "var")

weightsP <- weightsP %>%
  pivot_longer(cols = c("wP"), names_to = "var")

weightsPweekly <- weightsPweekly %>%
  pivot_longer(cols = c("wP.weekly"), names_to = "var")

weightsPmonthly <- weightsPmonthly  %>%
  pivot_longer(cols = c("wP.monthly"), names_to = "var")

weights <- full_join(weights, weightsP, by = c("lag","var","value"))
weights <- full_join(weights, weightsPweekly, by = c("lag","var","value"))
weights <- full_join(weights, weightsPmonthly, by = c("lag","var","value"))

# pivot_longer ci2.5
temp <- c()
var_list <- c("wT","wV","wSs", "wP$", "wP.weekly", "wP.monthly")
col_list <-c("cred2.5_wT","cred2.5_wV","cred2.5_wSs", "cred2.5_wP$", "cred2.5_wP.weekly", "cred2.5_wP.monthly")
for( i in c(1:length(var_list))){
  for( j in c(1:nrow(weights))){
    
    # if the var name in the df row matches i
    # then assign the correct sd value to temp
    if(grepl(var_list[i],weights$var[j])){
      temp2 <- select(weights,contains(col_list[i]))[j,]
      temp[[j]] <- as.numeric(temp2)
    }
    
  }
}
weights_longer <- weights %>% 
  mutate(ci2.5 = temp) %>%
  select(-contains(col_list))


# pivot_longer ci97.5
temp <- c()
var_list <- c("wT","wV","wSs", "wP$", "wP.weekly", "wP.monthly")
col_list <-c("cred97.5_wT","cred97.5_wV","cred97.5_wSs", "cred97.5_wP$", "cred97.5_wP.weekly", "cred97.5_wP.monthly")
for( i in c(1:length(var_list))){
  for( j in c(1:nrow(weights_longer))){
    
    # if the var name in the df row matches i
    # then assign the correct sd value to temp
    if(grepl(var_list[i],weights_longer$var[j])){
      temp2 <- select(weights_longer,contains(col_list[i]))[j,]
      temp[[j]] <- as.numeric(temp2)
    }
    
  }
}
weights_longer <- weights_longer %>% 
  mutate(ci97.5 = temp) %>%
  select(-contains(col_list))


# fix P

weights_longer$ci2.5[22:30] <- weights_longer$cred2.5_wP[22:30]
weights_longer$ci97.5[22:30] <-weights_longer$cred97.5_wP[22:30]

weights_longer$ci2.5<-as.numeric(weights_longer$ci2.5)
weights_longer$ci97.5 <-as.numeric(weights_longer$ci97.5)


# Graph
library(ggforce)
weights_longer %>%
ggplot(aes(x=lag, y=value)) + 
  geom_errorbar(aes(ymax = ci97.5, ymin = ci2.5), position = "dodge") +
  geom_point(pch=21, size = 3) +
  facet_col("var") +
  labs(title = NULL, y = "weight", x = NULL)+
  theme_bw() +
  theme(legend.position = "none",
        legend.text=element_text(size=14),
        text = element_text(size=14),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))


#####################################################################
#####################################################################
#####################################################################
#####################################################################


################################################### Notes:
#adapted from heili lowman on kelp project 2022
# t_max <- function(pointID, start_date, end_date){
#   
#   dat <- monthly_t %>%
#     filter(MeasurementID == pointID) %>%
#     filter(ClimateYear == as.numeric(str_sub(MeasurementID, -4, length(MeasurementID)))) %>%
#     filter(Month >= start_date & Month <= end_date)
#   
#   return(max(dat$Tmax, na.rm = TRUE))
#   
# }
# 
# # output vector y which will be filled in by the for loop
# prod1$Tmax <- rep(NA, length(prod1$Nest_ID)) # create an empty vector of NAs
# 
# for(i in 1:length(prod1$Nest_ID)) {
#   prod1$Tmax[i] <- t_max(pointID = prod1$Nest_ID[i], 
#                          start_date = prod1$start_date[i],
#                          end_date = prod1$end_date[i]) 
#}