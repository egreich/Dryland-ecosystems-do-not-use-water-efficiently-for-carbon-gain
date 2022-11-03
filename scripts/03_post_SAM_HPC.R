#!/usr/bin/env Rscript

### This script file will re-combine 3 chains for all sites
### after running the SAM model on an HPC


# Load packages
library(rjags)
load.module('dic')
library(tidyverse)


# Load self-made functions
source("./scripts/post_SAM_HPC_function.R")
source("./scripts/coda_functions.R")

### recombine chains, summarize as a dataframe ###
post_SAM(1)
post_SAM(2)
post_SAM(3)
post_SAM(4)
post_SAM(5)
post_SAM(6)
post_SAM(7)
post_SAM(8)



