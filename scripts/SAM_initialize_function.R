#!/usr/bin/env Rscript

get_SAM_inits <- function(dataIN, key, chain=NULL){
  
  # Create necessary folders if they do not already exist
  if(!file.exists("models/inits")) { dir.create("models/inits")}
  
  # Define filenames
  initfilename <- paste("./models/inits/inits_", key, ".RData", sep = "")
  # If running on an HPC, save the initials by chain number
  if(!is.null(chain)){
    initfilename <- paste("./models/inits/inits_", chain,"_",key, ".RData", sep = "")
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
  # Basically a matrix version of X2, defined later
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
              Nstart2 = Nstart2,
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
                Nstart2 = Nstart2,
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
  
  # Initial values are estimated using a linear model. As in the data list (above), 
  # covariates are centered and standardized. Replace name of covariate in quotes
  # with the appropriate column number in the dataIN file.
  X1  = cbind(data$VPD[Yday[Nstart:Nend]], #1
              data$Tair[Yday[Nstart:Nend]], #2
              data$P[Yday[Nstart:Nend]], #3
              data$PAR[Yday[Nstart:Nend]], #4
              data$Sshall[Yday[Nstart:Nend]], #5
              data$Sdeep[Yday[Nstart:Nend]]) #6
  
  # Notes: The code below is indexed numerically, which you will have to pay attention to as you change covariates of interest
  # Squared terms calculated for VPD and Tair
  X1a = cbind(X1[,1]^2, X1[,2]^2) 
  # Put all covariates together;
  # Interactions incorporated into linear model used to estimate initial values
  X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,3], X1[,1]*X1[,5], X1[,1]*X1[,6], X1[,2]*X1[,3], X1[,2]*X1[,5], X1[,2]*X1[,6], 
              X1[,3]*X1[,5], X1[,3]*X1[,6], X1[,5]*X1[,6])
  # Fit simple linear model
  fit <- lm(Y[Nstart:Nend] ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + X1[,6] + # main effects
              X1a[,1] + X1a[,2] + # squared
              X2[,1] + X2[,2] + X2[,3] + X2[,4] + X2[,5] + X2[,6] + X2[,7] + X2[,8] + X2[,9] + X2[,10]) # interactions
  # Extract coefficient estimates:
  beta0  = fit$coefficients[1] # the intercept
  beta1  = fit$coefficients[2:7] # main effects
  beta1a = fit$coefficients[8:9] # squared effects
  
  beta2 <- matrix(data = 0, nrow = 10, ncol = 1)
  beta2[1:4,1]   = as.numeric(fit$coefficients[10:13]) # X1[,1] interactions (VPD with Tair, P, and soil moistures)
  beta2[5:7,1]   = as.numeric(fit$coefficients[14:16]) # X1[,2] interactions (Tair with P and soil moistures)
  beta2[8:9,1]   = as.numeric(fit$coefficients[17:18]) # X1[,3] interactions (P and soil moistures)
  beta2[10,1]   = as.numeric(fit$coefficients[19]) # X1[,5]*X1[,6] interaction (soil moistures)
  
  # Create initials based on the above estimates:
  inits = list(list(beta0 = beta0, beta1 = beta1, beta1a = beta1a, beta2 = beta2, sig = 1), #pink
               list(beta0 = beta0/10, beta1 = beta1/10, beta1a = beta1a/10, beta2 = beta2/10, sig = 1/2), #green
               list(beta0 = beta0*10, beta1 = beta1*10, beta1a = beta1a*10, beta2 = beta2*10, sig = 2)) #blue
  
  # If running on an HPC, make inits 1 chain corresponding to the chain number
  if(!is.null(chain)){
    inits = inits[[chain]]
  }
  
  #####################################################################
  # Part 2: Initialize JAGS Model
  n.adapt = 500
  n.chains = 3
  # If running on an HPC, make n.chains=1
  if(!is.null(chain)){
    n.chains = 1
  }
  # If running the vcp site, run the split model to account for data gaps
  model.name <- ifelse(key != "vcp", "./models/Model_SAM_ETpart.R", "./models/Model_SAM_ETpart_split.R")
  
  start<-proc.time()
  
  jm1.b=jags.model(model.name,
                   data=data,
                   n.chains=n.chains,
                   n.adapt=n.adapt,
                   inits = inits)
  end<-proc.time()
  elapsed<- (end-start)/60
  print("jags.model done running; minutes to completion:")
  print(elapsed[3])
  
  #####################################################################
  #Part 3: Run coda.samples with JAGS model
  
  # Choose the parameters to monitor. For this analysis, we allowed the model to 
  # converge while monitoring variables included in "zc1"
  # below before monitoring dYdX (zc1dYdX), Xant (zc1X), or Y (zc1Y)
  
  n.iter = 1000
  thin = 1
  
  # parameters to track
  params = c("deviance","beta0","beta1","beta1a",
             "beta2", "sig")

  start<-proc.time()
  zc1 = coda.samples(jm1.b,variable.names=params,
                     n.iter=n.iter,thin = thin)
  end<-proc.time()
  elapsed<- (end-start)/(60*60)
  print("coda.samples done running; hours to completion:")
  print(elapsed[3])
  
  #####################################################################
  # Part 4: Check diagnostics
  
  #plotting to visualize chains, diagnose convergence issues, etc
  #mcmcplot(zc1)
  
  #check convergence via gelman diagnostics, psrf should be <1.2
  #gel<-gelman.diag(zc1, multivariate = F)
  #print(gel)
  
  #check how much to run
  #raft<-raftery.diag(zc1)
  #raft<-maxraft(chains=3,coda=zc1) #find the min number of iterations needed per chain
  #print(raft)
  
  #####################################################################
  # Part 5: Save inits for future runs
  
  # inits to save
  #init_names = c("beta0","beta1","beta1a","beta2", "sig")
  
  # find which variables in the coda object to remove
  #remove_vars = get_remove_index(init_names, params)
  
  #extract final iteration to reinitialize model
  newinits<-initfind(zc1, OpenBUGS = F)
  #newinits[[1]]
  
  
  #saved.state <- removevars(initsin = newinits, variables = remove_vars) # remove non-root nodes
  saved.state <- newinits # if nothing needs to be removed
  
  #check items in list
  #saved.state[[1]]

  # Save initials
  save(saved.state, file=initfilename) 
  
}