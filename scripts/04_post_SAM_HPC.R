#!/usr/bin/env Rscript

### This script file will re-combine 3 chains for all sites
### after running the SAM model on an HPC

# Set run params
args<-commandArgs(TRUE)
print(args)
print("site:")
(site <- as.numeric(args[1]))
print("model version:")
(modelv <- as.numeric(args[2]))


# Load packages
library(rjags)
load.module('dic')
library(tidyverse)
library(mcmcplots)


# Load self-made functions
source("./scripts/post_SAM_HPC_function.R")
source("./scripts/check_convergence_function.R")
source("./scripts/coda_functions.R")

### recombine chains, summarize as a dataframe ###
#post_SAM(site,modelv)

############################  start temp
s=site
# Create necessary folders if they do not already exist
if(!file.exists("output_dfs")){dir.create("output_dfs")}

# Get the chain and site info for each slurm job:
slurm <- read.csv("Slurm_jobs_3chains.csv")
chains = slurm$chainENDs
sites = slurm$siteEND
modelvs = slurm$modelvEND

# key for which site corresponds to which site ID
if(s == 1){
  key = "seg"
} else if(s == 2){
  key = "ses"
} else if(s == 3){
  key = "wjs"
} else if(s == 4){
  key = "mpj"
} else if(s == 5){
  key = "vcp"
} else if(s == 6){
  key = "vcm1"
} else if(s == 7){
  key = "vcm2"
} else if(s == 8){
  key = "vcs"
}
print("about to define filenames")
# Define filenames
zcfilename <- paste("./output_coda/coda_all_", key, "_v", modelv, ".RData", sep = "")
zcxfilename <- paste("./output_coda/coda_all_x_", key, "_v", modelv, ".RData", sep = "")
table.dicfilename <- paste("./output_dfs/table.dic_", key,"_v", modelv, ".csv", sep = "")
dffilename <- paste("./output_dfs/df_sum_", key, "_v", modelv, ".csv", sep = "")
print("filenames defined")
# Combine chains into one mcmc list (coda object) for each site:


# make an empty mcmc list
coda_all <- mcmc.list()
coda_all_x <- mcmc.list()
# fill list with single-chain coda:
r = c()
x = list()

# Pick off chain indices for each site
r = which(sites == s & modelvs == modelv)

for(c in 1:length(r)){
  # Load single-chain coda objects
  load(paste("./output_coda/zc_", c,"_", key, "_v", modelv, ".RData", sep = ""))
  load(paste("./output_coda/zc_x_", c,"_", key, "_v", modelv, ".RData", sep = ""))
  # Create coda objects with all chains:
  coda_all[[c]] <- as.mcmc(zc1[[1]])
  coda_all_x[[c]] <- as.mcmc(zc1x[[1]])
}
print("about to save coda_all")
# Save coda values:
save(coda_all, file = zcfilename)
save(coda_all_x, file = zcxfilename)
print("coda_all saved")

### Create summary data frames that are compatible with tidyverse
df_sum <- coda.fast(chains=3, burn.in=0, thin=1, coda=coda_all) # Summarizing chains via Mike Fell's code
df_sum_x <- coda.fast(chains=3, burn.in=0, thin=1, coda=coda_all_x)
df_sum <- rbind(df_sum,df_sum_x)
df_sum <- rownames_to_column(df_sum, "var")
df_sum <- df_sum %>% # make index column
  mutate(ID = sub('.*\\[(.*)\\]', '\\1', df_sum$var))
df_sum <- df_sum %>% # separate index column into 1st and 2nd dimension
  mutate(ID1 = sub('(.*)\\,.*', '\\1', df_sum$ID),
         ID2 = sub('.*\\,(.*)', '\\1', df_sum$ID))
df_sum$ID2 <- ifelse(!grepl(',', df_sum$ID), NA, df_sum$ID2) # get rid of ID2 if there's no 2nd dimension
df_sum$ID1 <- ifelse(!grepl('[^[:alpha:]]', df_sum$ID), 1, df_sum$ID1) # make ID1=1 if there is only 1 instance
df_sum <- df_sum %>% 
  mutate(var = sub('(.*)\\[.*', '\\1', df_sum$var)) # get rid of numbers in var col
df_sum <- df_sum %>% 
  select("var","ID1","ID2","mean","median","sd","pc2.5","pc97.5") %>% #reorder columns, drop ID
  mutate(site = key)

write.csv(df_sum, dffilename)
print("df_sum saved")

### DIC

# Load and name jags.model objects based on chains
chain = 1 # for chain 1
jagsfilename <- paste("./output_coda/jagsmodel_", chain,"_", key,"_v", modelv, ".csv", sep = "")
load(jagsfilename)
jm1 <- jm1.b # define jags.model names based on chain number
chain = 2 # for chain 2
jagsfilename <- paste("./output_coda/jagsmodel_", chain,"_", key,"_v", modelv, ".csv", sep = "")
load(jagsfilename)
jm2 <- jm1.b # define jags.model names based on chain number
chain = 3 # for chain 3
jagsfilename <- paste("./output_coda/jagsmodel_", chain,"_", key,"_v", modelv, ".csv", sep = "")
load(jagsfilename)
jm3 <- jm1.b # define jags.model names based on chain number

# make a new jm object to hold all chains
new_jm <- jm1
# combining initials
environment(new_jm[["ptr"]])[["init.values"]][[2]] <- environment(jm2[["ptr"]])[["init.values"]][[1]]
environment(new_jm[["ptr"]])[["inits"]][[2]] <- environment(jm2[["ptr"]])[["inits"]][[1]]
environment(new_jm[["ptr"]])[["init.values"]][[3]] <- environment(jm3[["ptr"]])[["init.values"]][[1]]
environment(new_jm[["ptr"]])[["inits"]][[3]] <- environment(jm3[["ptr"]])[["inits"]][[1]]
environment(new_jm[["ptr"]])[["n.chains"]] <- 3
# adding data
environment(new_jm[["ptr"]])[["model.state"]][[2]] <- environment(jm2[["ptr"]])[["model.state"]][[1]]
environment(new_jm[["ptr"]])[["model.state"]][[3]] <- environment(jm3[["ptr"]])[["model.state"]][[1]]

# recompile
new_jm$recompile()

# set DIC monitor
dic1<- dic.samples(new_jm, n.iter=3000)

# DIC = the model fit + the penalty (effect of number of parameters)
DIC1 = sum(dic1$deviance)+sum(dic1$penalty)
pD1 = sum(dic1$penalty)

# Set up a table with comparison statistics
table.dic <- data.frame(site = key, modelv = modelv, DIC = DIC1, pD = pD1)

write.csv(table.dic, file = table.dicfilename)  # save for model comparisons

##################### end temp


### create mcmc plots, save raftery tables ###
check_convergence(site,modelv)



