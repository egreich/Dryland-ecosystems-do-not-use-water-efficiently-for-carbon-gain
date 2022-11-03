#!/usr/bin/env Rscript

### This script file will allow the user to run the SAM model to 
### evaluate controls on T (or WUE?) for sites in the NMEG. 

# Set run params
args<-commandArgs(TRUE)
print(args)
print("chain:")
(chain <- as.numeric(args[1]))
print("site:")
(site <- as.numeric(args[2]))
print("seed:")
(SEED <- as.numeric(args[3]))

# Set defined R seed
set.seed(SEED, kind = NULL, normal.kind = NULL)
# Generate "random" seed for jags
JAGS.seed<-ceiling(runif(1,1,10000000))

# key for which site corresponds to which index
if(site == 1){
  key = "seg"
} else if(site == 2){
  key = "ses"
} else if(site == 3){
  key = "wjs"
} else if(site == 4){
  key = "mpj"
} else if(site == 5){
  key = "vcp"
} else if(site == 6){
  key = "vcm1"
} else if(site == 7){
  key = "vcm2"
} else if(site == 8){
  key = "vcs"
}

# Load packages
library(rjags)
load.module('dic')
library(mcmcplots)
library(tidyverse)


# Load self-made functions
source("./scripts/SAM_initialize_function.R")
source("./scripts/SAM_function.R")
source("./scripts/coda_functions.R")

### Read in data ###
# Read in the covariate file of interest. This should be a file that has all 
# covariates for all time periods to ensure that 6-months of precipitation is 
# supplied to the model - not just growing season precipitation
#load(paste("./clean_data/dataIN_",key,".RData",sep=""))
#load(paste("./clean_data/dataIN_wue_",key,".RData",sep=""))
#load(paste("./clean_data/dataIN_gpp_",key,".RData",sep=""))
if(key == "vcm1"){
  
  dataIN_T = read.csv(paste("./clean_data/d_B_","vcm",".csv",sep=""))
  dataIN_WUE = read.csv(paste("./clean_data/d_B_wue_","vcm",".csv",sep=""))
  
  dataIN <- left_join(dataIN_T, dataIN_WUE)
  
  dataIN <- dataIN %>%
    filter(year < 2013)
  
}

if(key == "vcm2"){
  
  dataIN_T = read.csv(paste("./clean_data/d_B_","vcm",".csv",sep=""))
  dataIN_WUE = read.csv(paste("./clean_data/d_B_wue_","vcm",".csv",sep=""))
  
  dataIN <- left_join(dataIN_T, dataIN_WUE)
  
  dataIN <- dataIN %>%
    filter(year > 2013)
  
}

if(key %in% c("seg", "ses", "wjs", "mpj", "vcp", "vcs")){
  
  dataIN_T = read.csv(paste("./clean_data/d_B_",key,".csv",sep=""))
  dataIN_WUE = read.csv(paste("./clean_data/d_B_wue_",key,".csv",sep=""))
  
  dataIN <- left_join(dataIN_T, dataIN_WUE)
  
}

### Get initials ###
# Run the below function to get initial estimates
# Note: If you've already run the model, you do not need to do this
# Note: Do not run this function using an HPC, 
#unless you comment out the mcmcplot function in the function code
get_SAM_inits(dataIN, key, chain)

### Run the model ###
SAM_WUE(dataIN, key, chain)

