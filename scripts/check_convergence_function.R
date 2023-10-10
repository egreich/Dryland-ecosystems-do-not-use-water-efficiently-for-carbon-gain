check_convergence<- function(s,modelv){

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
  
  # Define filenames
  foldername <- paste("models/convergence/", key, "_v", modelv, sep = "")
  dirname <- paste("./models/convergence/", key, "_v", modelv, sep = "")
  filename <- paste("MCMC_", key, "_v", modelv, sep = "")
  zcfilename <- paste("./output_coda/coda_all_", key, "_v", modelv, ".RData", sep = "")
  
  # Create necessary folders if they do not already exist
  if(!file.exists("models/convergence")) { dir.create("models/convergence")}
  if(!file.exists(foldername)) { dir.create(foldername)}

  load(zcfilename) # called jagsui
  
  jm_coda <- jagsui$samples 
  
  mcmcplot(jm_coda,
           random = 15,
           dir = dirname,
           filename = filename)
  

}