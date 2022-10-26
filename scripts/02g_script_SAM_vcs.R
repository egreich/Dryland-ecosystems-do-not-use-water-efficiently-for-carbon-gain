#!/usr/bin/env Rscript

### This script file will allow the user to run the SAM model to 
### evaluate controls on T (or WUE?) for sites in the NMEG. 


# Load packages
library(rjags)
load.module('dic')
library(mcmcplots)
library(tidyverse)
#library(postjags)
#library(coda4dummies)


# Load self-made functions
source("./scripts/SAM_initialize_function.R")
source("./scripts/SAM_function.R")
source("./scripts/coda_functions.R")

### Read in data ###
# Read in the covariate file of interest. This should be a file that has all 
# covariates for all time periods to ensure that 6-months of precipitation is 
# supplied to the model - not just growing season precipitation
dataIN_T = read.csv("./clean_data/d_B_vcs.csv")
dataIN_WUE = read.csv("./clean_data/d_B_wue_vcs.csv")

dataIN <- left_join(dataIN_T, dataIN_WUE)

### Get initials ###
# Run the below function to get initial estimates
# Note: If you've already run the model, you do not need to do this
# Note: Do not run this function using an HPC, 
#unless you comment out the mcmcplot function in the function code
get_SAM_inits(dataIN, "vcs")

### Run the model ###
SAM_WUE(dataIN, "vcs")

### diagnostics

#raftery.diag(zc1)
#gelman.diag(zc1)

