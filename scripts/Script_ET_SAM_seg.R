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
dataIN_T = read.csv("./clean_data/d_B_seg.csv") %>%
  filter(year %in% c(2017,2018,2019,2020)) # shorten df for testing
dataIN_WUE = read.csv("./clean_data/d_B_wue_seg.csv") %>%
  filter(year %in% c(2017,2018,2019,2020))

dataIN <- left_join(dataIN_T, dataIN_WUE)

### Get initials ###
# Run the below function to get initial estimates
# Note: If you've already run the model, you do not need to do this
# Note: Do not run this function using an HPC, 
#unless you comment out the mcmcplot function in the function code
get_SAM_inits(dataIN, "seg")

### Run the model ###
output <- SAM_WUE(dataIN, "seg")



### Look at stuff  ###
load("./output_coda/zc_seg.RData")
mcmcplot(zc1)
#sum_tab <- coda.fast(chains=3, burn.in=0, thin=1, coda=zc1)


# Function to extract posterior means, and 2.5 and 97.5 CI quantiles
# takes a list of variable names and the coda summary
# all variables in the list MUST have the same length posterior outputs
coda_pivot_longer <- function(var_list, coda_sum, colnams = NULL){
  
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
        searchterm <- paste(utils::glob2rx(var_list[i]), sep = "") # check if user is using *
        if(length(grep(searchterm, row.names(sum_tb))) == 0){ # if we still find nothing
          print(paste("Warning: ", var_list[i], " not found in coda summary output", sep = ""))
          next
        }
      }
    }
    
    voi_list[[j]] <- sum_tb[grep(searchterm, row.names(sum_tb)),1]
    voi_list[[j+1]] <- quan_tb[grep(searchterm, row.names(quan_tb)),1]
    voi_list[[j+2]] <- quan_tb[grep(searchterm, row.names(quan_tb)),5]
    
    if(is.null(colnams)){
      
      column_names[[j]] <- paste("B_", var_list[i], sep = "")
      column_names[[j+1]] <- paste("ci2.5_", var_list[i], sep = "")
      column_names[[j+2]] <- paste("ci97.5_", var_list[i], sep = "")
      
    }
    
    if(!is.null(colnams)){
      
      column_names[[j]] <- paste(colnams[i], sep = "")
      column_names[[j+1]] <- paste("ci2.5_", colnams[i], sep = "")
      column_names[[j+2]] <- paste("ci97.5_", colnams[i], sep = "")
      
    }
    
    j = j + 3
    
  }
  suppressMessages(df <- dplyr::bind_cols(voi_list))
  colnames(df) <- column_names
  
  # pivot longer posteriors
  if(!is.null(colnams)){
    var_list <- colnams
  }
  df <- df %>%
    pivot_longer(cols = c(var_list), names_to = "var")
  
  # pivot_longer ci2.5
  temp <- c()
  
  col_list <- unlist(column_names[grep("2.5",column_names)])
  for( i in c(1:length(var_list))){
    for( j in c(1:nrow(df))){
      
      # if the var name in the df row matches i
      # then assign the correct ci value to temp
      if(grepl(var_list[i],df$var[j])){
        temp2 <- select(df, ends_with(col_list[i]))[j,]
        temp[[j]] <- as.numeric(temp2)
      }
      
    }
  }
  df_longer <- df %>%
    mutate(ci2.5 = temp) %>%
    select(-contains(col_list))
  
  # pivot_longer ci97.5
  temp <- c()
  if(!is.null(colnams)){
    var_list <- colnams
  }
  col_list <- unlist(column_names[grep("97.5",column_names)])
  for( i in c(1:length(var_list))){
    for( j in c(1:nrow(df_longer))){
      
      # if the var name in the df row matches i
      # then assign the correct ci value to temp
      if(grepl(var_list[i],df_longer$var[j])){
        temp2 <- select(df_longer, ends_with(col_list[i]))[j,]
        temp[[j]] <- as.numeric(temp2)
      }
      
    }
  }
  df_longer <- df_longer %>%
    mutate(ci97.5 = temp) %>%
    select(-contains(col_list))
  
  
  return(df_longer)
  
}

sum <- summary(zc1)
varlist <- c("wT","wV","wSs")
weights <- coda_pivot_longer(varlist, sum)
weights <- weights %>%
  rowid_to_column("lag")
varlist <- c("wP")
weightsP <- coda_pivot_longer(varlist, sum)
weightsP <- weightsP %>%
  rowid_to_column("lag")
varlist <- c("wP.weekly")
weightsPweekly <- coda_pivot_longer(varlist, sum)
weightsPweekly <- weightsPweekly  %>%
  rowid_to_column("lag")
varlist <- c("wP.monthly")
weightsPmonthly <- coda_pivot_longer(varlist, sum)
weightsPmonthly <- weightsPmonthly %>%
  rowid_to_column("lag")


weights_longer <- rbind(weights, weightsP, weightsPweekly, weightsPmonthly)


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
