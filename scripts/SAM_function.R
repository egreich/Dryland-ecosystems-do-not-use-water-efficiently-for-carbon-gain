#!/usr/bin/env Rscript

SAM_WUE <- function(dataIN, key, chain=NULL){
  
  # Create necessary folders if they do not already exist
  if(!file.exists("output_coda")) { dir.create("output_coda")}
  if(!file.exists("output_dfs")) { dir.create("output_dfs")}
  if(!file.exists("models/inits")) { dir.create("models/inits")}
  
  # Define filenames
  # If not running on an HPC
  if(is.null(chain)){
  initfilename <- paste("./models/inits/inits_", key, ".RData", sep = "")
  zcfilename <- paste("./output_coda/coda_all_", key, ".RData", sep = "")
  dffilename <- paste("./output_dfs/df_sum_", key, ".csv", sep = "")
  }
  # If running on an HPC, we will run a postscript later to combine chains
  if(!is.null(chain)){
    initfilename <- paste("./models/inits/inits_", chain,"_", key, ".RData", sep = "")
    zcfilename <- paste("./output_coda/zc_", chain,"_", key, ".RData", sep = "")
  }
  
  # Make the response variable file of interest. This file includes 
  # growing-season response (Y) variables of interest along with indices 
  # that link the Y variables back to covariates of interest.
  # NOTE: Kym did this to only test growing season response variables (4/1 to 10/31 for each year)
  YIN = dataIN %>%
    select(c("date","year","month","day","B_WUE.pred")) %>% # select response variable of interest
    #mutate(WUE_test = GPP/B_T) %>% # test WUE
    #mutate(WUE_test =ifelse(is.na(WUE_test), 0, WUE_test)) %>%
    rowid_to_column("dayind") # dayind: index to link the growing season Y variables back to the appropriate 
  # row in the covariate data set
  
  # Growing season starts and stops for each year
  Gstart = YIN %>%
    filter(month == 4 & day == 1) %>%
    select(dayind)
  
  Gstop = YIN %>%
    filter(month == 10 & day == 31) %>%
    select(dayind)
  
  # Filter YIN to just be growing season (April-Oct)
  YIN = YIN %>%
    filter(month %in% c(4,5,6,7,8,9,10))
  
  # Choose column in YIN that provides indices linking response variables with covariates
  Yday = YIN$dayind
  
  # Change Y to the column for the response variable of interest (T or WUE) 
  Y  = YIN$B_WUE.pred
  
  # jIND file provides indices to calculate interactions between covariates
  # Basically a matrix version of X2, defined in SAM_initialize_function
  # X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,3], X1[,1]*X1[,5], X1[,1]*X1[,6], X1[,2]*X1[,3], X1[,2]*X1[,5], X1[,2]*X1[,6], 
  # X1[,3]*X1[,5], X1[,3]*X1[,6], X1[,5]*X1[,6])
  jIND <- data.frame(j = c(1:10),
                     ID1 = c(1,1,1,1,2,2,2,3,3,5),
                     ID2 = c(2,3,5,6,3,5,6,5,6,6))
  
  # Choose the starting index. This is an index for a row in the Y data file. 
  # The value in the indexed row should be greater than 1 to accommodate 
  # calculation of antecedent values.
  Nstart = Gstart$dayind[1]
  Nstart2 = Gstart$dayind[3] # starting at index 3 allows us to test years into the past
  # Choose the ending index. 
  Nend   = nrow(YIN)  
  
  # Prepare data for JAGS -- covariates are scaled
  data = list(Nstart = Nstart,
              Nend = Nend,
              Nlag = 11,
              NlagP = 9, 
              Nparms = 6, # Nparms is the number of driving variables included to calculate main effects
              Yday = YIN$dayind, # Choose column in YIN that provides indices linking response variables with covariates
              ID1 = jIND[,2], 
              ID2 = jIND[,3],
              jlength = nrow(jIND),
              Y = Y,
              VPD = as.vector(scale(dataIN$VPD,center=TRUE,scale=TRUE)), # scale function takes vector of values, centers and scales by SD
              Tair = as.vector(scale(dataIN$Tair,center=TRUE,scale=TRUE)),
              P = as.vector(scale(dataIN$P,center=TRUE,scale=TRUE)),
              PAR = as.vector(scale(dataIN$PPFD_IN,center=TRUE,scale=TRUE)),
              Sshall = as.vector(scale(dataIN$S,center=TRUE,scale=TRUE)),
              Sdeep = as.vector(scale(dataIN$Sdeep,center=TRUE,scale=TRUE)),
              C1 = c(0, 1, 2, 4, 6, 13, 20, 27, 55, 83, 111), #stop times for covariates
              C2 = c(0, 1, 2, 3, 5, 7, 14, 21, 28, 56, 84), #start times for covariates
              P1 = c(6, 13, 20, 27, 55, 83, 111, 139, 167), #stop times for precip
              P2 = c(0, 7, 14, 21, 28, 56, 84, 112, 140)) #start times for precip
  
  # If running the vcp site, run the split model to account for data gaps
  if(key == "vcp"){
    Nsplitstart = YIN %>%
      filter(year == 2012 & month == 10 & day == 31) %>%
      select(dayind)
    Nsplitend = YIN %>%
      filter(year == 2014 & month == 4 & day == 1) %>%
      select(dayind)
    Nstart2split = YIN %>%
      filter(year == 2016 & month == 4 & day == 1) %>% # give two years of room
      select(dayind)
    # Prepare data for JAGS -- covariates are scaled
    data = list(Nstart = Nstart,
                Nsplitstart = Nsplitstart,
                Nstart2split = Nstart2split,
                Nsplitend = Nsplitend,
                Nend = Nend,
                Nlag = 11,
                NlagP = 9, 
                Nparms = 6, # Nparms is the number of driving variables included to calculate main effects
                Yday = YIN$dayind, # Choose column in YIN that provides indices linking response variables with covariates
                ID1 = jIND[,2], 
                ID2 = jIND[,3],
                jlength = nrow(jIND),
                Y = Y,
                VPD = as.vector(scale(dataIN$VPD,center=TRUE,scale=TRUE)), # scale function takes vector of values, centers and scales by SD
                Tair = as.vector(scale(dataIN$Tair,center=TRUE,scale=TRUE)),
                P = as.vector(scale(dataIN$P,center=TRUE,scale=TRUE)),
                PAR = as.vector(scale(dataIN$PPFD_IN,center=TRUE,scale=TRUE)),
                Sshall = as.vector(scale(dataIN$S,center=TRUE,scale=TRUE)),
                Sdeep = as.vector(scale(dataIN$Sdeep,center=TRUE,scale=TRUE)),
                # Set stop and start indices for time into the past
                # The amount of time should accumulate into greater blocks as we move further into the past
                C1 = c(0, 1, 2, 4, 6, 13, 20, 27, 55, 83, 111), #stop times for covariates
                C2 = c(0, 1, 2, 3, 5, 7, 14, 21, 28, 56, 84), #start times for covariates
                P1 = c(6, 13, 20, 27, 55, 83, 111, 139, 167), #stop times for precip
                P2 = c(0, 7, 14, 21, 28, 56, 84, 112, 140)) #start times for precip
  }
  

  # Load initial values from previous run
  load(initfilename)
 
   #####################################################################
  # Part 2: Initialize JAGS Model
  n.adapt = 500 # adjust this number (and n.iter) as appropriate 
  n.chains = 3
  # If running on an HPC, make n.chains=1
  if(!is.null(chain)){
    n.chains = 1
  }
  model.name <- ifelse(key != "vcp", "./models/Model_SAM_ETpart.R", "./models/Model_SAM_ETpart_split.R")
  
  start<-proc.time() # set start time
  jm1.b=jags.model(model.name,
                   data=data,
                   n.chains=n.chains,
                   n.adapt=n.adapt,
                   inits = saved.state[[2]])
  end<-proc.time()
  elapsed<- (end-start)/60
  print("jags.model done running; minutes to completion:")
  print(elapsed[3])
  
  #####################################################################
  # Part 3: Run coda.samples with JAGS model
  
  # Choose the parameters to monitor. For this analysis, we allowed the model to 
  # converge while monitoring variables included in "zc1"
  # below before monitoring dYdX (zc1dYdX), Xant (zc1X), or Y (zc1Y)
  
  
  n.iter = 25000
  thin = 5
  
  # parameters to track
  params = c("Y", "Y.rep","beta0","beta1","beta1a",
             "beta2","deviance","dYdX","wP","wSd","wSs","wT","wV",
             "sig")
  
  start<-proc.time()
  zc1 = coda.samples(jm1.b,variable.names=params,
                     n.iter=n.iter,thin = thin)
  end<-proc.time()
  elapsed<- (end-start)/(60*60)
  print("coda.samples done running; hours to completion:")
  print(elapsed[3])
  
  save(zc1, file = zcfilename)  # save the model output for graphs
  
  #####################################################################
  # Part 4: Save coda summary
  
  # Only run this part (making a df) if not running on an HPC
  if(is.null(chain)){
    
  # Summarizing chains via Mike Fell's code
  df_sum <- coda.fast(chains=3, burn.in=0, thin=1, coda=zc1)
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
  }
  
  # Save output
  # sumzc <- summary(zc1)
  # sumstats <- sumzc[["statistics"]] # save coda summary in table form
  # write.csv(sumstats, file = sumfilename)
  # quanstats <- sumzc[["quantiles"]] # save coda quantiles in table form
  # write.csv(quanstats, file = quanfilename)
  
  #####################################################################
  # Part 5: Save inits for future runs
  
  # inits to save
  init_names = c("beta0","beta1","beta1a","beta2", "sig")

  # find which variables in the coda object to remove
  remove_vars = get_remove_index(init_names, params)
  
  #extract final iteration to reinitialize model if needed
  newinits<-initfind(zc1, OpenBUGS = F)
  #remove non-root node variables
  saved.state <- removevars(initsin = newinits, variables=remove_vars) # remove non-variable nodes
  #check both items in list
  save(saved.state, file=initfilename)
  
  
  #####################################################################

}