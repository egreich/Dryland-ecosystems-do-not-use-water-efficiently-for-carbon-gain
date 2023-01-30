#!/usr/bin/env Rscript

SAM_WUE <- function(dataIN, key, modelv, chain=NULL){
  
  # Create necessary folders if they do not already exist
  if(!file.exists("output_coda")) { dir.create("output_coda")}
  if(!file.exists("output_dfs")) { dir.create("output_dfs")}
  if(!file.exists("models/inits")) { dir.create("models/inits")}
  
  # Define filenames
  # If not running on an HPC
  if(is.null(chain)){
  initfilename <- paste("./models/inits/inits_", key,"_v", modelv, ".RData", sep = "")
  zcfilename <- paste("./output_coda/coda_all_", key,"_v", modelv, ".RData", sep = "")
  zcxfilename <- paste("./output_coda/coda_all_x_", key,"_v", modelv, ".RData", sep = "")
  dffilename <- paste("./output_dfs/df_sum_", key,"_v", modelv, ".csv", sep = "")
  }
  # If running on an HPC, we will run a postscript later to combine chains
  if(!is.null(chain)){
    initfilename <- paste("./models/inits/inits_", chain,"_", key,"_v", modelv, ".RData", sep = "")
    zcfilename <- paste("./output_coda/zc_", chain,"_", key,"_v", modelv, ".RData", sep = "")
    zcxfilename <- paste("./output_coda/zc_x_", chain,"_", key,"_v", modelv, ".RData", sep = "")
    jagsfilename <- paste("./output_coda/jagsmodel_", chain,"_", key,"_v", modelv, ".csv", sep = "")
  }
  
  
  # End for ETpart model (all days)
  N = nrow(dataIN)
  
  # Make the response variable file of interest. This file includes 
  # growing-season response (Y) variables of interest along with indices 
  # that link the Y variables back to covariates of interest.
  # NOTE: Kym did this to only test growing season response variables (4/1 to 10/31 for each year)
  YIN = dataIN %>%
    mutate(Tperiod = ifelse(Season == "Spring", 1, NA)) %>%
    mutate(Tperiod = ifelse(Season == "Summer", 2, Tperiod)) %>%
    mutate(Tperiod = ifelse(Season == "Fall", 3, Tperiod)) %>%
    mutate(Tperiod = ifelse(Season == "Winter", 4, Tperiod)) %>%
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
  
  # Filter YIN to just be growing season (April-Oct)
  #YIN = YIN %>%
   # filter(month %in% c(4,5,6,7,8,9,10))
  
  # Change Y to the column for the response variable of interest (T or WUE) 
  Y  = YIN$WUE
  
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
  
  
  C1 = c(0, 1, 2, 4, 6, 13, 20) #stop times for covariates
  C2 = c(0, 1, 2, 3, 5, 7, 14) #start times for covariates
  P1 = c(6, 13, 20, 27, 55, 83, 111, 139, 167) #stop times for precip
  P2 = c(0, 7, 14, 21, 28, 56, 84, 112, 140) #start times for precip

  Nlag = length(C1)
  NlagP = length(P1)
  
  
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
              Nlag = Nlag,
              NlagP = NlagP, 
              Nparms = 6, # Nparms is the number of driving variables included to calculate main effects
              Yday = YIN$dayind, # Choose column in YIN that provides indices linking response variables with covariates
              ID1 = jIND[,2], 
              ID2 = jIND[,3],
              jlength = nrow(jIND),
              #Y = Y,
              #Yoff = Yoff,
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
    #Nstart2split = YIN %>%
     # filter(year == 2016 & month == 4 & day == 1) %>% # give two years of room
      #select(dayind)
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
                Nlag = Nlag,
                NlagP = NlagP,
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
                C1 = C1,
                C2 = C2,
                P1 = P1,
                P2 = P2)
  }
  

  # Load initial values from previous run
  initfilenametemp <- paste("./models/inits/inits_", chain,"_", key, ".RData", sep = "") # temp, only for first run
  load(initfilenametemp) # temp, change to initfilename after first run
  if(file.exists(initfilename)){
    load(initfilename)
  }
 
   #####################################################################
  # Part 2: Initialize JAGS Model
  n.adapt = 500 # adjust this number (and n.iter) as appropriate
  
  # Define the correct dimension of the saved state
  if(length(saved.state)!=2){
  inits = saved.state #[[2]]
  # if(modelv==7){
  #   saved.state[["sig.WUE"]] <- runif(1, 5, 15) # temp for model 6
  #   inits = saved.state
  #}
  }
  # If we ran the model already, and saved the initials with their names
  if(length(saved.state)==2){
    inits = saved.state[[2]]
    # if(modelv==7){
    #   saved.state[["initials"]][[1]][["sig.WUE"]] <- runif(1, 5, 15) # temp for model 6
    #   inits = saved.state[[2]]
    # }
    # if(modelv %in% c(2,5)){ # temp for model 2,5 - to have main effects vary by season
    #   beta1inits <- saved.state[["initials"]][[1]][["beta1"]]
    #   saved.state[["initials"]][[1]][["beta1"]] <- matrix(data = 0, nrow = 6, ncol = 4)
    #   saved.state[["initials"]][[1]][["beta1"]][,1] <- beta1inits
    #   saved.state[["initials"]][[1]][["beta1"]][,2] <- beta1inits
    #   saved.state[["initials"]][[1]][["beta1"]][,3] <- beta1inits
    #   saved.state[["initials"]][[1]][["beta1"]][,4] <- beta1inits
    #   inits = saved.state[[2]]
    # }
  }
  # If running on an HPC, make n.chains=1
  if(!is.null(chain)){
    n.chains = 1
  } else if(is.null(chain)){
    n.chains = 3
  }
  model.name <- ifelse(key != "vcp", paste("./models/Model_v", modelv, ".R", sep=""), paste("./models/Model_split_v", modelv, ".R", sep=""))
  
  start<-proc.time() # set start time
  jm1.b=jags.model(model.name,
                   data=data,
                   n.chains=n.chains,
                   n.adapt=n.adapt,
                   inits = inits)
  end<-proc.time()
  elapsed<- (end-start)/60
  print("jags.model done running; minutes to completion:")
  print(elapsed[3])
  
  # Save jags.model object to calculate DIC and pD later
  save(jm1.b, file=jagsfilename)
  
  #####################################################################
  # Part 3: Run coda.samples with JAGS model
  
  # Choose the parameters to monitor. For this analysis, we allowed the model to 
  # converge while monitoring variables included in "zc1"
  # below before monitoring dYdX (zc1dYdX), Xant (zc1X), or Y (zc1Y)
  
  
  n.iter = 30000
  thin = 6
  
  # parameters to track
  params = c("deviance",
             "beta0","beta1","beta1a","beta2",#"dYdX",
             "wP","wSd","wSs","wT","wV",
             "tau.ET", "tau.log.WUE",
             "ET", "E.model", "ET.pred", "T.pred", "T.ratio",
             "WUE.pred")
  
  if(modelv %in% c(1,2)){
    params = c("deviance",
               "beta0","beta1","beta1a","beta2",#"dYdX",
               "wP",
               "tau.ET", "tau.log.WUE",
               "ET", "E.model", "ET.pred", "T.pred", "T.ratio",
               "WUE.pred")
  }
  
  if(modelv==7){
    params = c("deviance",
               "beta0","beta1","beta1a","beta2",
               "tau.ET", "sig.WUE",
               "ET", "E.model", "ET.pred", "T.pred", "T.ratio",
               "WUE.pred")
  }
  
  start<-proc.time()
  zc1 = coda.samples(jm1.b,variable.names=params,
                     n.iter=n.iter,thin = thin)
  end<-proc.time()
  elapsed<- (end-start)/(60*60)
  print("coda.samples done running; hours to completion:")
  print(elapsed[3])
  
  # New code object to monitor quantities for Dinf, WAIC, and R2
  zc1x = coda.samples(jm1.b, variable.names = c("Dsum","ldx","dx","R2","ET.rep"), n.iter = 3000)
  
  save(zc1, file = zcfilename)  # save the model output for graphs
  save(zc1x, file = zcxfilename)  # save the model output for graphs
  
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
  init_names = c("beta0","beta1","beta1a","beta2", "tau.ET", "tau.log.WUE")
  if(modelv==7){
    init_names = c("beta0","beta1","beta1a","beta2", "tau.ET", "sig.WUE")
  }

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