get_SAM_inits <- function(dataIN, key){
  
  # Create necessary folders if they do not already exist
  if(!file.exists("models/inits")) { dir.create("models/inits")}
  
  # Define filenames
  initfilename <- paste("./models/inits/inits_", key, ".RData", sep = "")
  
  # Make the response variable file of interest. This file includes 
  # growing-season response (Y) variables of interest along with indices 
  # that link the Y variables back to covariates of interest.
  # NOTE: Kym did this to only test growing season response variables (4/1 to 10/31 for each year)
  YIN = dataIN %>%
    select(c("date","year","month","day","B_T", "GPP")) %>% # select response variable of interest
    mutate(WUE_test = GPP/B_T) %>% # test WUE
    mutate(WUE_test =ifelse(is.na(WUE_test), 0, WUE_test)) %>%
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
  Y  = YIN$WUE_test
  
  # jIND file provides indices to calculate interactions between covariates
  # Basically a matrix version of X2, defined later
  #X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,3], X1[,1]*X1[,5], X1[,2]*X1[,3], X1[,2]*X1[,5], X1[,3]*X1[,5]) 
  jIND <- data.frame(j = c(1:6),
                     ID1 = c(1,1,1,2,2,3),
                     ID2 = c(2,3,5,3,5,5))
  
  # Choose the starting index. This is an index for a row in the Y data file. 
  # The value in the indexed row should be greater than 1 to accommodate 
  # calculation of antecedent values.
  Nstart = Gstart$dayind[1] 
  # Choose the ending index. 
  Nend   = nrow(YIN)  
  
  # Prepare data for JAGS -- covariates are scaled
  data = list(Nstart = Nstart, 
              Nend = Nend, 
              Nlag = 7, 
              NlagP = 9, 
              Nparms = 5, # Nparms is the number of driving variables included to calculate main effects
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
              P1 = c(6, 13, 20, 27, 55, 83, 111, 139, 167), #stop times
              P2 = c(0, 7, 14, 21, 28, 56, 84, 112, 140)) #start times
  
  # Initial values are estimated using a linear model. As in the data list (above), 
  # covariates are centered and standardized. Replace name of covariate in quotes
  # with the appropriate column number in the dataIN file.
  X1  = cbind(data$VPD[Yday[Nstart:Nend]], #1
              data$Tair[Yday[Nstart:Nend]], #2
              data$P[Yday[Nstart:Nend]], #3
              data$PAR[Yday[Nstart:Nend]], #4
              data$Sshall[Yday[Nstart:Nend]]) #5
  
  # Notes: The code below is indexed numerically, which you will have to pay attention to as you change covariates of interest
  # Squared terms calculated for VPD and Tair
  X1a = cbind(X1[,1]^2, X1[,2]^2) 
  # Put all covariates together;
  # Interactions incorporated into linear model used to estimate initial values
  X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,3], X1[,1]*X1[,5], X1[,2]*X1[,3], X1[,2]*X1[,5], X1[,3]*X1[,5]) 
  # Fit simple linear model
  fit <- lm(Y[Nstart:Nend] ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + X1a[,1] + 
              X1a[,2] + X2[,1] + X2[,2]  + X2[,3]  + X2[,4]  + X2[,5]  + X2[,6])
  # Extract coefficient estimates:
  beta0  = fit$coefficients[1] # the intercept
  beta1  = fit$coefficients[2:6] # main effects
  beta1a = fit$coefficients[7:8] # squared effects
  
  beta2 <- matrix(data = 0, nrow = 6, ncol = 1)
  beta2[1:3,1]   = as.numeric(fit$coefficients[9:11]) # X1[,1] interactions
  beta2[4:5,1]   = as.numeric(fit$coefficients[12:13]) # X1[,2] interactions
  beta2[6,1]   = as.numeric(fit$coefficients[14]) # X1[,3] interactions
  
  # Create initials based on the above estimates:
  inits = list(list(beta0 = beta0, beta1 = beta1, beta1a = beta1a, beta2 = beta2, sig = 1), #pink
               list(beta0 = beta0/10, beta1 = beta1/10, beta1a = beta1a/10, beta2 = beta2/10, sig = 1/2), #green
               list(beta0 = beta0*10, beta1 = beta1*10, beta1a = beta1a*10, beta2 = beta2*10, sig = 2)) #blue
  
  #####################################################################
  # Part 2: Initialize JAGS Model
  n.adapt = 500 # adjust this number (and n.iter) as appropriate 
  n.iter = 1000
  n.chains = 3
  
  start<-proc.time()
  jm1.b=jags.model("./models/Model_SAM_ETpart.R",
                   data=data,
                   n.chains=n.chains,
                   n.adapt=n.adapt,
                   inits = inits)
  end<-proc.time()
  elapsed<- end-start
  print(elapsed[3])
  
  #####################################################################
  #Part 3: Run coda.samples with JAGS model
  
  # Choose the parameters to monitor. For this analysis, we allowed the model to 
  # converge while monitoring variables included in "zc1"
  # below before monitoring dYdX (zc1dYdX), Xant (zc1X), or Y (zc1Y)
  
  #n.iter = 40000
  n.iter = 1000
  #thin = 40
  thin = 1
  
  # parameters to track
  params = c("deviance","beta0","beta1","beta1a",
             "beta2", "wT","wV","wP","wSs",
             "wP.weekly","wP.monthly", "sig")

  zc1 = coda.samples(jm1.b,variable.names=params,
                     n.iter=n.iter,thin = thin)
  
  #####################################################################
  # Part 4: Check convergence
  
  #plotting to visualize chains, diagnose convergence issues, etc
  mcmcplot(zc1)
  
  #check convergence
  gel<-gelman.diag(zc1, multivariate = F)
  print(gel)
  
  #####################################################################
  # Part 5: Save inits for future runs
  
  # inits to save
  init_names = c("beta0","beta1","beta1a","beta2", "sig")
  
  # variables to remove
  get_remove_index <- function(to_keep, list){
    out_list <- c()
    for(j in c(1:length(list))){
      if(list[j] %in% to_keep){
        out_list[j] = NA
      } else{
        out_list[j] = j
      }
    }
    out_list <- out_list[!is.na(out_list)]
    out_list
  }
  
  remove_vars = get_remove_index(init_names, params)
  
  #extract final iteration to reinitialize model if needed
  newinits<-initfind(zc1, OpenBUGS = F)
  #newinits[[1]]
  #remove non-root node variables
  saved.state <- removevars(initsin = newinits, variables=remove_vars) # remove non-variable nodes
  #check both items in list
  #saved.state[[1]]
  save(saved.state, file=initfilename) 
  
}