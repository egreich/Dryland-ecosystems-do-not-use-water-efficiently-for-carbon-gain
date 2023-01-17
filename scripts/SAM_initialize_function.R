#!/usr/bin/env Rscript

get_SAM_inits <- function(dataIN, dataIN_wue, key, modelv){
  
  # Create necessary folders if they do not already exist
  if(!file.exists("models/inits")) { dir.create("models/inits")}
  
  # Define filenames for each chain
  initfilename1 <- paste("./models/inits/inits_1_", key, ".RData", sep = "")
  initfilename2 <- paste("./models/inits/inits_2_", key, ".RData", sep = "")
  initfilename3 <- paste("./models/inits/inits_3_", key, ".RData", sep = "")
  
  # End for ETpart model (all days)
  N = nrow(dataIN)
  
  # Make the response variable file of interest. This file includes 
  # growing-season response (Y) variables of interest along with indices 
  # that link the Y variables back to covariates of interest.
  # NOTE: Kym did this to only test growing season response variables (4/1 to 10/31 for each year)
  YIN = dataIN %>%
    mutate(Tperiod = ifelse(Season = "Spring", 1, NA)) %>%
    mutate(Tperiod = ifelse(Season = "Summer", 2, Tperiod)) %>%
    mutate(Tperiod = ifelse(Season = "Fall", 3, Tperiod)) %>%
    mutate(Tperiod = ifelse(Season = "Winter", 4, Tperiod)) %>%
    #select(c("date","year","month","day","B_WUE.pred")) %>% # select response variable of interest
    rowid_to_column("dayind") # dayind: index to link the growing season Y variables back to the appropriate 
  # row in the covariate data set
  
  # Growing season starts and stops for each year
  Gstart = YIN %>%
    filter(month == 4 & day == 1) %>%
    select(dayind)
  
  
  # If running the burned site, start later to accommodate January data start
  if(key == "vcm2"){
    Gstart = YIN %>%
    filter(month == 7 & day == 1) %>%
    select(dayind)}
  
  Gstop = YIN %>%
    filter(month == 10 & day == 31) %>%
    select(dayind)
  
  # Filter YIN to just be growing season (April-Oct)
  #YIN = YIN %>%
  #  filter(month %in% c(4,5,6,7,8,9,10))
  
  # Change Y to the column for the response variable of interest (T or WUE) 
  #Y  = YIN$WUE
  Yday = YIN$dayind
  
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
  # Choose the ending index. 
  Nend   = nrow(YIN)
  
  # Prepare data for JAGS -- covariates are scaled
  data = list(N = N, # Number of rows
              Nperiods = 4, # Number of seasons
              Tperiod = YIN$Tperiod,
              Esnow = dataIN$Esnow,
              Tsoil = dataIN$Tsoil,
              S = dataIN$S,
              S.min = min(dataIN$S),
              ET = dataIN$ET,
              GPP = dataIN$GPP,
              ws = dataIN$ws,
              conv.fact = 0.0864 * 0.408,
              rho = dataIN$pair,
              Ri = dataIN$Ri,
              rah_unstable = dataIN$rah, # in the data, rah is just rah_unstable when appropriate. We are letting rah_stable vary (stochastic) so that's why we read this in.
              Cp = 1000,
              pi = 3.14159265359,
              Rwv = 461.52, # specific gas constant for water vapor in J/(kg*K)
              sres = 0.15*dataIN$fclay[1], # residual soil moisture
              ssat = 0.489 - (0.126*dataIN$fsand[1]),
              psisat = -10*exp(1.88-1.31*dataIN$fsand[1]), # parameterized air entry pressure, in mm of water
              rssmin = 50, # minimum soil resistance in s/m
              g  = 9.80665,     # acceleration due to gravity in m/s^2
              gamma = dataIN$gamma,
              #rah = dataIN$rah,
              e.sat = dataIN$es,
              e.a = dataIN$ea,
              Z = dataIN$Z[1], # reference height (m)
              #fclay = dataIN$fclay[1],
              fc = dataIN$fc[1],
              bch = 2.91 + 15.9*dataIN$fclay[1], # The Clapp and Hornberger parameter est. as in Cosby et al. [1984]
              LAI = dataIN$LAI_mod,
              P.unscaled = dataIN$P,
              
              # data for SAM model
              Nstart = Nstart,
              Nend = Nend,
              Nlag = 11,
              NlagP = 9, 
              Nparms = 6, # Nparms is the number of driving variables included to calculate main effects
              Yday = YIN$dayind, # Choose column in YIN that provides indices linking response variables with covariates
              ID1 = jIND[,2], 
              ID2 = jIND[,3],
              jlength = nrow(jIND),
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
    # Prepare data for JAGS -- covariates are scaled
    data = list(N = N, # Number of rows
                Nperiods = 4, # Number of seasons
                Tperiod = YIN$Tperiod,
                Esnow = dataIN$Esnow,
                Tsoil = dataIN$Tsoil,
                S = dataIN$S,
                S.min = min(dataIN$S),
                ET = dataIN$ET,
                GPP = dataIN$GPP,
                ws = dataIN$ws,
                conv.fact = 0.0864 * 0.408,
                rho = dataIN$pair,
                Ri = dataIN$Ri,
                rah_unstable = dataIN$rah, # in the data, rah is just rah_unstable when appropriate. We are letting rah_stable vary (stochastic) so that's why we read this in.
                Cp = 1000,
                pi = 3.14159265359,
                Rwv = 461.52, # specific gas constant for water vapor in J/(kg*K)
                sres = 0.15*dataIN$fclay[1], # residual soil moisture
                ssat = 0.489 - (0.126*dataIN$fsand[1]),
                psisat = -10*exp(1.88-1.31*dataIN$fsand[1]), # parameterized air entry pressure, in mm of water
                rssmin = 50, # minimum soil resistance in s/m
                g  = 9.80665,     # acceleration due to gravity in m/s^2
                gamma = dataIN$gamma,
                #rah = dataIN$rah,
                e.sat = dataIN$es,
                e.a = dataIN$ea,
                Z = dataIN$Z[1], # reference height (m)
                #fclay = dataIN$fclay[1],
                fc = dataIN$fc[1],
                bch = 2.91 + 15.9*dataIN$fclay[1], # The Clapp and Hornberger parameter est. as in Cosby et al. [1984]
                LAI = dataIN$LAI_mod,
                P.unscaled = dataIN$P,
                
                # data for SAM model
                Nstart = Nstart,
                Nsplitstart = Nsplitstart,
                #Nstart2split = Nstart2split,
                Nsplitend = Nsplitend,
                Nend = Nend,
                Nlag = 11,
                NlagP = 9, 
                Nparms = 6, # Nparms is the number of driving variables included to calculate main effects
                Yday = YIN$dayind, # Choose column in YIN that provides indices linking response variables with covariates
                ID1 = jIND[,2], 
                ID2 = jIND[,3],
                jlength = nrow(jIND),
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
  
  
  # Get initials for WUE precision based on the log of previous model versions
  tau.WUE = 1/(sd(log(dataIN_wue$B_WUE.pred))**2)
  # predicted WUE based on previous model version, to use as a response variable in simple regressions
  Y = log(dataIN_wue$B_WUE.pred)
  
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
  inits = list(list(beta0 = beta0, beta1 = beta1, beta1a = beta1a, beta2 = beta2, tau.ET = 1/(sd(dataIN$ET)**2), tau.log.WUE = tau.WUE), #pink
               list(beta0 = beta0/10, beta1 = beta1/10, beta1a = beta1a/10, beta2 = beta2/10, tau.ET =  (1/(sd(dataIN$ET)**2))/10, tau.log.WUE =  tau.WUE/3), #green
               list(beta0 = beta0*10, beta1 = beta1*10, beta1a = beta1a*10, beta2 = beta2*10, tau.ET =  (1/(sd(dataIN$ET)**2))*10, tau.log.WUE =  tau.WUE*3)) #blue
  
  #####################################################################
  # Part 2: Initialize JAGS Model
  n.adapt = 500
  n.chains = 3
  
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
  
  n.iter = 500
  thin = 1
  
  # parameters to track
  params = c("deviance","beta0","beta1","beta1a",
             "beta2", "tau.ET", "tau.log.WUE")

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
  #init_names = c("beta0","beta1","beta1a","beta2", "tau.ET", "tau.log.WUE")
  
  #extract final iteration to reinitialize model
  newinits<-initfind(zc1, OpenBUGS = F)
  #newinits[[1]] # tracked parameters, will not include deviance
  
  # find which variables in the coda object to remove
  #remove_vars = get_remove_index(init_names, newinits[[1]])
  

    #saved.state <- removevars(initsin = newinits, variables = remove_vars) # remove non-root nodes
    saved.state <- newinits # if nothing needs to be removed
  
    # Save inits based on chains with lowest deviance
    dev_col <- which(colnames(zc1[[1]]) == "deviance")
    dev1<- mean(zc1[[1]][,dev_col])
    dev2<- mean(zc1[[2]][,dev_col])
    dev3<- mean(zc1[[3]][,dev_col])
    dev_min <- min(dev1, dev2, dev3)
    if(dev1 == dev_min){
      devin = 1
    } else if(dev2 == dev_min){
      devin = 2
    } else if(dev3 == dev_min){
      devin = 3
    }
    
    chain1 = saved.state[[2]][[devin]] # Best (low dev) initials for chain 1
    chain2 = lapply(saved.state[[2]][[devin]],"*",2)
    chain3 = lapply(saved.state[[2]][[devin]],"/",2)
    
    saved.state <- chain1
    save(saved.state, file=initfilename1)  # Save new initials
    saved.state <- chain2
    save(saved.state, file=initfilename2)  # Save new initials
    saved.state <- chain3
    save(saved.state, file=initfilename3)  # Save new initials
  
}


