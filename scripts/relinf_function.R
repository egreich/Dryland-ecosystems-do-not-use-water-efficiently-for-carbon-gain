#!/usr/bin/env Rscript

relinf_WUE <- function(dataIN, key, modelv, voi, newinits, lowdev=F, post_only=F){
  
  # Create necessary folders if they do not already exist
  if(!file.exists("output_relinf")) { dir.create("output_relinf")}
  if(!file.exists("output_relinf/output_coda")) { dir.create("output_relinf/output_coda")}
  if(!file.exists("output_relinf/output_dfs")) { dir.create("output_relinf/output_dfs")}
  if(!file.exists("output_relinf/inits")) { dir.create("output_relinf/inits")}
  if(!file.exists("output_relinf/convergence")) { dir.create("output_relinf/convergence")}
  if(!file.exists(paste("output_relinf/convergence/", key, "_v", modelv, "_voi", voi, sep = ""))) { dir.create(paste("output_relinf/convergence/", key, "_v", modelv, "_voi", voi, sep = ""))}
  
  # Define filenames
  initfilename <- paste("./output_relinf/inits/inits_", key,"_v", modelv, "_voi", voi, ".RData", sep = "")
  zcfilename <- paste("./output_relinf/output_coda/coda_all_", key,"_v", modelv, "_voi", voi, ".RData", sep = "")
  dffilename <- paste("./output_relinf/output_dfs/df_sum_", key,"_v", modelv, "_voi", voi, ".csv", sep = "")
  
  
  # End for ETpart model (all days)
  N = nrow(dataIN)
  
  # Make the response variable file of interest. This file includes 
  # growing-season response (Y) variables of interest along with indices 
  # that link the Y variables back to covariates of interest.
  # NOTE: Kym did this to only test growing season response variables (4/1 to 10/31 for each year)
  YIN = dataIN %>%
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
  
  # Filter YIN to just be growing season (April-Oct)
  #YIN = YIN %>%
   # filter(month %in% c(4,5,6,7,8,9,10))
  
  # Get rid of pre-April starts (Oct and Nov)
  YIN = YIN[Gstart$dayind[1]:nrow(YIN),]
  # Filter out winter
  YIN <- YIN %>% filter(Season != "Winter")
  
  # jIND file provides indices to calculate interactions between covariates
  # Basically a matrix version of X2, defined in SAM_initialize_function and below
  jIND <- data.frame(j = c(1:15),
                     ID1 = c(1,1,1,1,2,2,2,3,3,5,4,4,4,4,4),
                     ID2 = c(2,3,5,6,3,5,6,5,6,6,1,2,3,5,6))
  
  # Choose the starting index. This is an index for a row in the Y data file. 
  # The value in the indexed row should be greater than 1 to accommodate 
  # calculation of antecedent values.
  Nstart = Gstart$dayind[1]
  # Choose the ending index. 
  Nend   = nrow(YIN)
  
  
  C1 = c(0, 1, 2, 4, 6, 13, 20) #stop times for covariates
  C2 = c(0, 1, 2, 3, 5, 7, 14) #start times for covariates
  P1 = c(6, 13, 20, 27, 55, 83, 111, 139, 167) #stop times for precip
  P2 = c(0, 7, 14, 21, 28, 56, 84, 112, 140) #start times for precip
  
  if(modelv==7){
    C1 = c(0, 1, 2, 3, 4, 5, 6) #stop times for covariates
    C2 = c(0, 1, 2, 3, 4, 5, 6) #start times for covariates
  }

  Nlag = length(C1)
  NlagP = length(P1)
  
  
  # Prepare data for JAGS -- covariates are scaled
  data = list(N = Nend, # Number of rows, same as Nend, I just changed it here so I don't have to change the model code again (!)
              Nperiods = 3, # Number of seasons
              Esnow = YIN$Esnow,
              Tsoil = YIN$Tsoil,
              S = YIN$S,
              S.min = min(YIN$S),
              ET = YIN$ET,
              GPP = YIN$GPP,
              ws = YIN$ws,
              conv.fact = (60*60*24)/((2.501 - 0.00237*YIN$Tair)*10^6), # latent heat of vaporization
              rho = YIN$pair,
              Ri = YIN$Ri,
              rah_unstable = YIN$rah, # in the data, rah is just rah_unstable when appropriate. We are letting rah_stable vary (stochastic) so that's why we read this in.
              Cp = 1000,
              pi = 3.14159265359,
              Rwv = 461.52, # specific gas constant for water vapor in J/(kg*K)
              sres = 0.15*YIN$fclay[1], # residual soil moisture
              ssat = 0.489 - (0.126*YIN$fsand[1]),
              psisat = -10*exp(1.88-1.31*YIN$fsand[1]), # parameterized air entry pressure, in mm of water
              rssmin = 50, # minimum soil resistance in s/m
              g  = 9.80665,     # acceleration due to gravity in m/s^2
              gamma = YIN$gamma,
              e.sat = YIN$es,
              e.a = YIN$ea,
              Z = YIN$Z[1], # reference height (m)
              fc = YIN$fc[1],
              bch = 2.91 + 15.9*YIN$fclay[1], # The Clapp and Hornberger parameter est. as in Cosby et al. [1984]
              LAI = YIN$LAI_mod,
              P.unscaled = YIN$P,
    
              # data for SAM model
              Nstart = Nstart,
              Nend = Nend,
              Nlag = Nlag,
              NlagP = NlagP, 
              Nparms = 6, # Nparms is the number of driving variables included to calculate main effects
              Yday = YIN$dayind, # Choose column in YIN that provides indices linking response variables with covariates
              ID1 = jIND[,2], 
              ID2 = jIND[,3],
              jlength = nrow(jIND),
              #Y = Y,
              Astar = 1, # Astar is the standard deviation parameter
              VPD = as.vector(scale(dataIN$VPD,center=TRUE,scale=TRUE)), # scale function takes vector of values, centers and scales by SD
              Tair = as.vector(scale(dataIN$Tair,center=TRUE,scale=TRUE)),
              P = as.vector(scale(dataIN$P,center=TRUE,scale=TRUE)),
              PAR = as.vector(scale(dataIN$PPFD_IN,center=TRUE,scale=TRUE)),
              Sshall = as.vector(scale(dataIN$S,center=TRUE,scale=TRUE)),
              Sdeep = as.vector(scale(dataIN$Sdeep,center=TRUE,scale=TRUE)),
              C1 = C1,
              C2 = C2,
              P1 = P1,
              P2 = P2)

  
  # If running the vcp site, run the split model to account for data gaps
  if(key == "vcp"){
    Nsplitstart = YIN %>%
      filter(year == 2012 & month == 10 & day == 31) %>%
      select(dayind)
    Nsplitend = YIN %>%
      filter(year == 2014 & month == 4 & day == 1) %>%
      select(dayind)
    
    data$Nsplitstart <- Nsplitstart
    data$Nsplitend <- Nsplitend
  }
  
  
  if(modelv == 9){
    # Change Y to the column for the response variable of interest (GPP or ET) 
    Y  = YIN$GPP
    
    # First 30 in list are for ET partitioning, remove them
    data[2:29] <- NULL
    data$Y <- Y
  }
  


  ################# initials code
  # Get initials for WUE precision based on the log of GPP/ET
  WUE.unlogged <- YIN$GPP/YIN$ET
  Y = log(WUE.unlogged)
  if(modelv==9){
    Y = YIN$GPP
  }
  if(modelv==8){
    Y = WUE.unlogged
  }
  for(i in 1:length(Y)){
    if(is.infinite(Y[i])){
      Y[i] <- 0
    }
  }
  
  sig.Y = sd(Y, na.rm = T)**2
  tau.Y = 1/(sd(Y, na.rm = T)**2)
  
  # Initial values are estimated using a linear model. As in the data list (above), 
  # covariates are centered and standardized. Replace name of covariate in quotes
  # with the appropriate column number in the dataIN file.
  X1  = cbind(data$VPD[YIN$dayind[Nstart:Nend]], #1
              data$Tair[YIN$dayind[Nstart:Nend]], #2
              data$P[YIN$dayind[Nstart:Nend]], #3
              data$PAR[YIN$dayind[Nstart:Nend]], #4
              data$Sshall[YIN$dayind[Nstart:Nend]], #5
              data$Sdeep[YIN$dayind[Nstart:Nend]]) #6
  
  # Notes: The code below is indexed numerically, which you will have to pay attention to as you change covariates of interest
  # Squared terms calculated for VPD and Tair
  X1a = cbind(X1[,1]^2, X1[,2]^2, X1[,3]^2, X1[,4]^2, X1[,5]^2, X1[,6]^2) 
  # Put all covariates together;
  # Interactions incorporated into linear model used to estimate initial values
  X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,3], X1[,1]*X1[,5], X1[,1]*X1[,6], X1[,2]*X1[,3], X1[,2]*X1[,5], X1[,2]*X1[,6], 
              X1[,3]*X1[,5], X1[,3]*X1[,6], X1[,5]*X1[,6], 
              X1[,4]*X1[,1], X1[,4]*X1[,2], X1[,4]*X1[,3], X1[,4]*X1[,5], X1[,4]*X1[,6])
  # Fit simple linear model
  fit <- lm(Y[Nstart:Nend] ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + X1[,6] + # main effects
              X1a[,1] + X1a[,2] + X1a[,3] + X1a[,4] + X1a[,5] + X1a[,6] + # squared
              X2[,1] + X2[,2] + X2[,3] + X2[,4] + X2[,5] + X2[,6] + X2[,7] + X2[,8] + X2[,9] + X2[,10] + X2[,11] + X2[,12] + X2[,13] + X2[,14] + X2[,15]) # interactions
  # Extract coefficient estimates:
  beta0  = fit$coefficients[1] # the intercept
  beta1  = fit$coefficients[2:7] # main effects
  beta1a = fit$coefficients[8:13] # squared effects
  # interactions, only written like this because we added on interactions later, and this helps us keep track of indexing
  beta2 <- matrix(data = 0, nrow = 15, ncol = 1)
  beta2[1:4,1]   = as.numeric(fit$coefficients[14:17]) # X1[,1] interactions (VPD with Tair, P, and soil moistures)
  beta2[5:7,1]   = as.numeric(fit$coefficients[18:20]) # X1[,2] interactions (Tair with P and soil moistures)
  beta2[8:9,1]   = as.numeric(fit$coefficients[21:22]) # X1[,3] interactions (P and soil moistures)
  beta2[10,1]   = as.numeric(fit$coefficients[23]) # X1[,5]*X1[,6] interaction (soil moistures)
  beta2[11:15,1]   = as.numeric(fit$coefficients[24:28]) # X1[,4] interactions (PAR) with all variables
  
  # Create initials based on the above estimates:
  # Note: These are only used when we don't already have initials saved from a previous run
  # For our models that have main effects vary by season, beta1 will have to be adjusted to be a matrix (an initial for Spring, Summer, Fall)
  # Because we ran the models that vary by season later in the modeling process, we adjusted the beta1 from previous model runs, and did not account for season here
  # init_names = c("beta0","beta1","beta1a","beta2", "tau.ET", "tau.log.WUE")
  # if(modelv==8){ # theoretical, simple autoregressive model
  #   init_names = c("beta0","beta1","beta1a","beta2", "tau.ET", "sig.WUE")
  # }
  # if(modelv==9){ # simple GPP-only model
  #   init_names = c("beta0","beta1","beta1a","beta2", "sig")
  # }
  inits = list(list(beta0 = beta0, beta1 = beta1, beta1a = beta1a, beta2 = beta2, tau.ET = 1/(sd(dataIN$ET)**2), tau.log.WUE = tau.Y), #pink
               list(beta0 = beta0/10, beta1 = beta1/10, beta1a = beta1a/10, beta2 = beta2/10, tau.ET =  (1/(sd(dataIN$ET)**2))/10, tau.log.WUE =  tau.Y/3), #green
               list(beta0 = beta0*10, beta1 = beta1*10, beta1a = beta1a*10, beta2 = beta2*10, tau.ET =  (1/(sd(dataIN$ET)**2))*10, tau.log.WUE =  tau.Y*3)) #blue
  if(modelv==8){
    inits = list(list(beta0 = beta0, beta1 = beta1, beta1a = beta1a, beta2 = beta2, tau.ET = 1/(sd(dataIN$ET)**2), sig.WUE = sig.Y), #pink
                 list(beta0 = beta0/10, beta1 = beta1/10, beta1a = beta1a/10, beta2 = beta2/10, tau.ET =  (1/(sd(dataIN$ET)**2))/10, sig.WUE = sig.Y/10), #green
                 list(beta0 = beta0*10, beta1 = beta1*10, beta1a = beta1a*10, beta2 = beta2*10, tau.ET =  (1/(sd(dataIN$ET)**2))*10, sig.WUE = sig.Y*10)) #blue
  }
  if(modelv==9){
    inits = list(list(beta0 = beta0, beta1 = beta1, beta1a = beta1a, beta2 = beta2, sig = sig.Y), #pink
                 list(beta0 = beta0/10, beta1 = beta1/10, beta1a = beta1a/10, beta2 = beta2/10, sig =  sig.Y/10), #green
                 list(beta0 = beta0*10, beta1 = beta1*10, beta1a = beta1a*10, beta2 = beta2*10, sig =  sig.Y*10)) #blue
  }
  
  # Initial values: from saved state (if we already ran the model, we want to start where we left off)
  if(newinits==F){
    if(file.exists(initfilename)){
      load(initfilename)
      
      if(lowdev == T){
        if(length(saved.state)==3){ # temp
          saved.state <- lowdevrestart(saved.state, vary_by = 2)
        }else if(file.exists(initlowfilename)){ # temp
          load(initlowfilename) # initlow object is just the lowest dev chain number, 1,2, or 3
          saved.state[[3]] <- initlow
          saved.state <- lowdevrestart(saved.state, vary_by = 2)
        }
      }
      
      initslist <- saved.state[[2]]
      
    }else if(!file.exists(initfilename)){
      initslist <- inits
    }
  }else if(newinits==T){
    initslist <- inits
  }
  
  
  
  # Remove voi from data to negate in model
  if(voi==1){ # VPD #1
    data$VPD[1:length(data$VPD)] <- 0
  }else if(voi==2){ # Tair #2
    data$Tair[1:length(data$Tair)] <- 0
  }else if(voi==3){ # P #3
    data$P[1:length(data$P)] <- 0
  }else if(voi==4){ # PAR #4
    data$PAR[1:length(data$PAR)] <- 0
  }else if(voi==5){ # Sshall #5
    data$Sshall[1:length(data$Sshall)] <- 0
  }else if(voi==6){ # Sdeep #6
    data$Sdeep[1:length(data$Sdeep)] <- 0
  }


  
  #####################################################################
  # Part 2: Run JAGS Model (jagsui initializes and runs in one step)
  
  if(post_only==F){
    
  model.name <- ifelse(key != "vcp", paste("./models/Model_v", modelv, ".R", sep=""), paste("./models/Model_split_v", modelv, ".R", sep=""))
  if(modelv==7){ # if running model 7, use the code for 3
    model.name <- ifelse(key != "vcp", paste("./models/Model_v", 3, ".R", sep=""), paste("./models/Model_split_v", 3, ".R", sep=""))
  }
  
  # Choose the parameters to monitor.
  
  # parameters to track
  params = c("deviance", # for 3 or 7
             "beta0","beta1","beta1a","beta2",
             "tau.ET", "tau.log.WUE",
             "ldx","dx", "R2", "Dsum", "beta0_p_temp", "beta1_p_temp", "beta1a_p_temp", "beta2_p_temp")
  
  if(modelv %in% c(1)){
    params = c("deviance",
               "beta0","beta1","beta1a","beta2",
               "tau.ET", "tau.log.WUE",
               "ldx","dx", "R2", "Dsum", "beta0_p_temp", "beta1_p_temp", "beta1a_p_temp", "beta2_p_temp")
  }
  
  if(modelv==8){
    params = c("deviance",
               "beta0","beta1","beta1a","beta2",
               "tau.ET", "sig.WUE",
               "ldx","dx", "R2", "Dsum", "beta0_p_temp", "beta1_p_temp", "beta1a_p_temp", "beta2_p_temp")
  }
  if(modelv==9){
    params = c("deviance",
               "beta0","beta1","beta1a","beta2",
               "tau",
               "ldx","dx", "R2", "Dsum", "beta0_p_temp", "beta1_p_temp", "beta1a_p_temp", "beta2_p_temp")
  }

  start<-proc.time()
  jagsui <- jags(data = data,
                 inits = initslist,
                 model.file = model.name,
                 parameters.to.save = params,
                 n.chains = 3,
                 n.adapt = 500,
                 n.thin = 10,
                 n.iter = 50000,
                 parallel = TRUE)
  end<-proc.time()
  elapsed<- (end-start)/(60*60)
  print("jags done running; hours to completion:")
  print(elapsed[3])
  
  #####################################################################
  # Part 3: Save coda summary

  save(jagsui, file = zcfilename)  # save the model output for graphs
  } # end post_only==F
  
  if(post_only==T){
    load(zcfilename) # named jagsui
  }
  
  jm_coda <- jagsui$samples # convert to coda form to work with postjags functions
  
  df_sum <- dumsum(jagsui, type = "jagsUI") # organize the coda object as a dataframe
  write.csv(df_sum, file = dffilename)
  
  #####################################################################
  # Part 4: Save inits for future runs
  
  # inits to save
  init_names = names(initslist[[1]])

  # create a saved.state object with initials for next run
  # saved.state[[3]] is the chain number with lowest deviance
  saved.state <- keepvars(codaobj = jm_coda, to_keep = init_names, paramlist = params, type="jagsUI")
  
  # save
  save(saved.state, file=initfilename)
  
  
  #####################################################################
  
  # check diagnostics
  mcmcplot(jm_coda,
           random = 15,
           dir = paste("./output_relinf/convergence/", key, "_v", modelv, "_voi", voi, sep = ""),
           filename = paste("MCMC_", key, "_v", modelv, "_voi", voi, sep = ""))

}