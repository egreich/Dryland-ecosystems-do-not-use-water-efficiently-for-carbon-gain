SAM_WUE <- function(dataIN, key){
  
  # Create necessary folders if they do not already exist
  if(!file.exists("output_coda")) { dir.create("output_coda")}
  if(!file.exists("output_dfs")) { dir.create("output_dfs")}
  if(!file.exists("models/inits")) { dir.create("models/inits")}
  
  # Define filenames
  initfilename <- paste("./models/inits/inits_", key, ".RData", sep = "")
  zcfilename <- paste("./output_coda/zc_", key, ".RData", sep = "")
  sumfilename <- paste("./output_coda/sum_", key, ".csv", sep = "")
  quanfilename <- paste("./output_coda/quan_", key, ".csv", sep = "")
  
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
              P1 = c(6, 13, 20, 27, 55, 83, 111, 139, 167), #stop times for precip
              P2 = c(0, 7, 14, 21, 28, 56, 84, 112, 140)) #start times for precip
  

  # Load initial values from previous run
  load(initfilename)
 
   #####################################################################
  # Part 2: Initialize JAGS Model
  n.adapt = 500 # adjust this number (and n.iter) as appropriate 
  n.iter = 10000
  n.chains = 3
  
  start<-proc.time() # set start time
  jm1.b=jags.model("./models/Model_SAM_ETpart.R",
                   data=data,
                   n.chains=n.chains,
                   n.adapt=n.adapt,
                   inits = saved.state[[2]])
  end<-proc.time()
  elapsed<- end-start
  print(elapsed[3])
  
  #####################################################################
  # Part 3: Run coda.samples with JAGS model
  
  # Choose the parameters to monitor. For this analysis, we allowed the model to 
  # converge while monitoring variables included in "zc1"
  # below before monitoring dYdX (zc1dYdX), Xant (zc1X), or Y (zc1Y)
  
  #n.iter = 40000
  n.iter = 10000
  #thin = 40
  thin = 1
  
  # parameters to track
  params = c("deviance","beta0","beta1","beta1a",
             "beta2", "wT","wV","wP","wSs",
             "wP.weekly","wP.monthly", "sig")
  
  zc1 = coda.samples(jm1.b,variable.names=params,
                     n.iter=n.iter,thin = thin)
  
  save(zc1, file = zcfilename)  # save the model output for graphs
  
  #####################################################################
  # Part 4: Save coda summary
  
  # Summarizing chains via Mike Fell's code
  sum_tab <- coda.fast(chains=3, burn.in=0, thin=1, coda=zc1)
  
  # Save output
  sumzc <- summary(zc1)
  sumstats <- sumzc[["statistics"]] # save coda summary in table form
  write.csv(sumstats, file = sumfilename)
  quanstats <- sumzc[["quantiles"]] # save coda quantiles in table form
  write.csv(quanstats, file = quanfilename)
  
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
  
  
  
  sum_tab # This will return the summary as the output in row form
}