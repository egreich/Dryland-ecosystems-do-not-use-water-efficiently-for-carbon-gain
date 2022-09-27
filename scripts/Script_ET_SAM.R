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

# Load packages
library(rjags) # install rjags if necessary
library(mcmcplots)
library(tidyverse)
library(postjags)
load.module('dic')

# Load self-made functions
source("./scripts/SAM_initialize_function.R")
source("./scripts/SAM_function.R")

### Read in data ###
# Read in the covariate file of interest. This should be a file that has all 
# covariates for all time periods to ensure that 6-months of precipitation is 
# supplied to the model - not just growing season precipitation
dataIN = read.csv("./output_dfs/d_B_mpj.csv") %>%
  filter(water_year %in% c(2017,2018,2019,2020)) # shorten df for testing
#dataIN = read.csv("./output_dfs/d_B_wue_mpj.csv")

### Get initials ###
# Run the below function to get initial estimates
# Note: If you've already run the model, you do not need to do this
# Note: Do not run this function using an HPC, 
#unless you comment out the mcmcplot function in the function code
get_SAM_inits(dataIN, "mpj")

### Run the model ###
output <- SAM_WUE(dataIN, "mpj")

load("./output_coda/zc_mpj.RData")
mcmcplot(zc1)
sum_tab <- coda.fast(chains=3, burn.in=0, thin=1, coda=zc1)

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