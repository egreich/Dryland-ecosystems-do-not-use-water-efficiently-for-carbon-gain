#!/usr/bin/env Rscript

### Helper functions
# Many of these functions are on github and are edited here in minor ways here
# Functions by Michael Fell (PostJAGS package): https://github.com/fellmk/PostJAGS
# Functions by Emma Reich and Biplabendu (Billu) Das (coda4dummies package): https://github.com/egreich/coda4dummies

library(purrr)
# Create a "not in" function using negate from the purrr package
`%nin%` <- negate(`%in%`)


dumsum <- function(jagsobj, type, col.names = NULL){
  
  # Create a "not in" function using negate from the purrr package
  `%nin%` <- purrr::negate(`%in%`)
  
  if(type %nin% c("rjags", "jagsUI")){
    paste("Please indicate whether this is a rjags or jagsUI samples object")
  }
  
  
  if(type == "jagsUI"){
    
    if(is.null(jagsobj$samples)){
      print("This is not a full jagsUI object. If this is just the 'samples' from jagsUI, set 'type' to rjags.")
    }
    
    jagsui <- jagsobj
    jm_coda <- jagsui$samples # convert to coda form to work with postjags functions
    
    # Organize the coda object as a dataframe
    df_sum <- coda.fast(jm_coda)
    df_sum <- tibble::rownames_to_column(df_sum, "var")
    df_sum <- df_sum %>% # make index column
      mutate(ID = sub('.*\\[(.*)\\]', '\\1', df_sum$var))
    df_sum$ID <- ifelse(grepl('[[:alpha:]]', df_sum$ID), 1, df_sum$ID) # make ID=1 if there is only 1 instance
    
    # make a lists of list of indices
    IDlist <- strsplit(df_sum$ID, ",") #temp
    
    # get number of dims in ID
    counter <- 1
    for(i in 1:length(IDlist)){
      counter <- ifelse(length(IDlist[[i]])>counter, length(IDlist[[i]]), counter)
    }
    
    # create a character vector of column names based on max dim
    new_columns <- list()
    for(i in 1:counter){
      new_columns[[i]] <- paste("ID",i, sep="")
    }
    new_columns <- as.character(new_columns)
    
    # for each dimension, create a new column with the correct ID
    df_sum <- df_sum %>%
      tidyr::separate(ID,new_columns,sep=",")
    
    df_sum <- df_sum %>%
      mutate(var = sub('(.*)\\[.*', '\\1', df_sum$var)) # get rid of numbers in var col
    
    if(exists("c")){ # remove c if overwritten
      rm(c)
    }
    
    # Fix length error with dYdX
      # extract overlap0
    overlap0 = do.call(c, jagsui$overlap0)
    overlap0 = lapply(list, function(x) overlap0[!is.na(overlap0)]) # remove NAs
    overlap0 = as.data.frame(overlap0)
    colnames(overlap0)[1] <- "overlap0"
    overlap0 <- tibble::rownames_to_column(overlap0, "varnum")
      # extract gel, but don't remove NAs
    gel = do.call(c, jagsui$Rhat) 
    gel = as.data.frame(gel)
    colnames(gel)[1] <- "gel"
    gel <- tibble::rownames_to_column(gel, "varnum")
    
    # join overlap0 and gel, using overlap0 as a guide for what parameters will match the coda object
    uiinfo <- left_join(overlap0, gel, by = "varnum")
    
    df_mod <- df_sum %>%
      select("var",starts_with("ID"),"mean","median","pc2.5","pc97.5") %>% #reorder columns, drop ID
      mutate(overlap0 = uiinfo$overlap0, gel = uiinfo$gel)
    
    df_mod[,2:(counter+1)] <- lapply(2:(counter+1), function(x) as.numeric(df_mod[[x]])) # make appropraite columns numeric
    
    if(!is.null(col.names)){
      colnames(df_mod)[2:(counter+1)] <- col.names
    }
    
    return(df_mod)
  }
  
  if(type == "rjags"){
    jm_coda <- jagsobj
    
    # Organize the coda object as a dataframe
    df_sum <- coda.fast(jm_coda)
    df_sum <- tibble::rownames_to_column(df_sum, "var")
    df_sum <- df_sum %>% # make index column
      mutate(ID = sub('.*\\[(.*)\\]', '\\1', df_sum$var))
    df_sum$ID <- ifelse(grepl('[[:alpha:]]', df_sum$ID), 1, df_sum$ID) # make ID=1 if there is only 1 instance
    
    # make a lists of list of indices
    IDlist <- strsplit(df_sum$ID, ",") #temp
    
    # get number of dims in ID
    counter <- 1
    for(i in 1:length(IDlist)){
      counter <- ifelse(length(IDlist[[i]])>counter, length(IDlist[[i]]), counter)
    }
    
    # create a character vector of column names based on max dim
    new_columns <- list()
    for(i in 1:counter){
      new_columns[[i]] <- paste("ID",i, sep="")
    }
    new_columns <- as.character(new_columns)
    
    # for each dimension, create a new column with the correct ID
    suppressWarnings({
      df_sum <- df_sum %>%
        tidyr::separate(ID,new_columns,sep=",")
    })
    
    df_sum <- df_sum %>%
      mutate(var = sub('(.*)\\[.*', '\\1', df_sum$var)) # get rid of numbers in var col
    
    df_mod <- df_sum %>%
      select("var",starts_with("ID"),"mean","median","pc2.5","pc97.5") #%>% #reorder columns, drop ID
    #mutate(overlap0 = do.call(c, jagsui$overlap0), gel = do.call(c, jagsui$Rhat))
    
    df_mod[,2:(counter+1)] <- lapply(2:(counter+1), function(x) as.numeric(df_mod[[x]])) # make appropraite columns numeric
    
    if(!is.null(col.names)){
      colnames(df_mod)[2:(counter+1)] <- col.names
    }
    
    return(df_mod)
  }
  
}

lowdevrestart <- function(saved_state, vary_by = 10){
  
  initlow <- saved_state[[3]] # initlow is just the lowest dev chain number
  
  # take chain with lowest deviance, and make remaining chains vary around it
  saved_state[[2]][[1]] = saved_state[[2]][[initlow]] # Best (low dev) initials for chain 1
  saved_state[[2]][[2]] = lapply(saved_state[[2]][[initlow]],"*",vary_by)
  saved_state[[2]][[3]] = lapply(saved_state[[2]][[initlow]],"/",vary_by)
  
  return(saved_state)
  
}

findlowdev <- function(codaobj){
  
  if(length(codaobj)>3){
    codaobj <- codaobj$samples
  }
  
  jm_coda <- codaobj
  # Save inits based on chains with lowest deviance
  dev_col <- which(colnames(jm_coda[[1]]) == "deviance")
  dev1<- mean(jm_coda[[1]][,dev_col])
  dev2<- mean(jm_coda[[2]][,dev_col])
  dev3<- mean(jm_coda[[3]][,dev_col])
  dev_min <- min(dev1, dev2, dev3)
  if(dev1 == dev_min){
    devin = 1
  } else if(dev2 == dev_min){
    devin = 2
  } else if(dev3 == dev_min){
    devin = 3
  }
  
  initlow <- devin
  print(paste("chain with lowest deviance: ", initlow, sep=""))
  
  
  return(initlow) #returns the number of the lowest deviance chain
  
}

keepvars <- function(codaobj, to_keep, paramlist, type){
  
  if(!is.null(codaobj$samples)){
    codaobj <- codaobj$samples
  }
  
  # Create a "not in" function using negate from the purrr package
  `%nin%` <- purrr::negate(`%in%`)
  
  remove_vars <- get_remove_index(to_keep, paramlist, type)
  
  newinits <- initfind(codaobj, OpenBUGS = FALSE)
  saved_state <- removevars(initsin = newinits,
                            variables = remove_vars)
  
  initlow <- findlowdev(codaobj) # find the lowest deviance chain
  saved_state[[3]] <- initlow
  names(saved_state[[3]]) <- "lowdevchain"
  
  return(saved_state)
  
}

dateconnect <- function(dfobj, datevect, datename = "date", identifier = NULL, varlist){
  
  dflistj <- list() # create an empty list for output
  for(j in 1:length(varlist)){
    dfobj2 <- dfobj %>%
      filter(var == varlist[j])
    
    if(!is.null(identifier)){ # if the variables we want to connect to dates have the same names but are indexed
      
      IDlist <- unique(as.vector(dfobj2[,identifier]))
      endloop <- length(IDlist)
      
      dflisti <- list() # create an empty list
      for( i in c(1:endloop)){
        
        dfobj3 <- dfobj2 %>%
          filter(ID2==as.numeric(IDlist[i]))
        
        #dfobj3 <- dfobj3[-c(1,2),] # temp for testing because I changed the processing btwn runs to include PAR
        
        dfobj3[,datename] <- datevect
        dflisti[[i]] <- dfobj3
      }
      dflistj[[j]] <- bind_rows(dflisti)
      
    }else{ # if the variables we want to connect to dates are just named differently
      dfobj3 <- dfobj2
      dfobj3[,datename] <- datevect
      dflistj[[j]] <- dfobj3
    }
  }
  
  dfout <- bind_rows(dflistj)
  
  if(datename %in% colnames(dfobj)){ # if the column name already exists
    
    dfout[,"tempcol"] <- dfout[,datename] # define a temp column name
    dfout <- dfout[,!names(dfout)==datename] # remove the datename column, or else it messes with joining
    
    dfoutout <- left_join(dfobj, dfout)
    
    # when the datename column has an NA, fill it with the tempcol
    for(k in 1:nrow(dfoutout)){
      if(is.na(dfoutout[k,datename])){
        dfoutout[k,datename] <- dfoutout$tempcol[k]
      }
    }
    
    dfoutout <- dfoutout %>% # delete the temp column
      select(-tempcol)
    
  }else{
    dfoutout <- left_join(dfobj, dfout)
  }
  
  return(dfoutout)
  
}



# function that finds the index number for variables from rjags or jagsUI output
# intended to work with Michael Fell's "removevars" function
# Inputs: 
# to_keep: string of variable names
# list: list of all parameters tracked
# type: rjags or jagsUI
# Output: list of index values for all variables NOT included in the input
get_remove_index <- function(to_keep, list, type){
  
  if(type %nin% c("rjags", "jagUI")){
    paste("Please indicate whether this is a rjags or jagsUI samples object")
  }
  
  if(type == "rjags"){
    list <- list[list != "deviance"] # remove deviance
    list <- sort(list, method = "radix")
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
  
  if(type == "jagsUI"){
    list <- list[list != "deviance"] # remove deviance
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
  
}

###############################################################################
#
# A function that does pairwise comparisons among defined outputs from JAGS
#
###############################################################################

# A function that does pairwise comparisons among defined outputs from JAGS
# requires multcompView, pacman
library(multcompView)
library(pacman)
bd.check_similarity <- function(data, params, which.params, m, l, u, anchor=NULL, ids) {
  pacman::p_load(tidyverse)
  foo <- data;
  if(!is.null(anchor)) {
    all.letters <- list()
    for(i in 1:length(which.params)) {
      bar1 <-
        foo %>% 
        select({{params}},{{m}},{{l}},{{u}},{{anchor}},{{ids}})
      bar1 <- bar1[which(bar1[1]==which.params[i]),]
      colnames(bar1) <- c("p","m","l","u","anchor","id")
      bar2 <-
        bar1 %>% 
        inner_join(bar1, by=c("anchor")) %>% 
        mutate(c = ifelse(m.x==m.y, "ns",
                          ifelse(m.y>=l.x & m.y<=u.x | m.x>=l.y & m.x<=u.y, 
                                 "ns","sig"))) %>% 
        group_by(anchor, id.x, id.y) %>% 
        mutate(ids = paste(sort(c(id.x,id.y)),collapse = "-")) %>% 
        ungroup() %>% 
        select(anchor, ids , c) %>% 
        distinct() %>% 
        mutate(id2=ids) %>% 
        separate(id2, c("id1","id2"), "-")
      
      which.anchor <- unique(bar1$anchor)
      which.ids <- unique(bar1$id)
      l.letters <- list()
      for (j in 1:length(which.anchor)) {
        dif3 <- bar2 %>% filter(anchor==which.anchor[j]) %>% mutate(c=ifelse(c=="ns",F,T)) %>% pull(c)
        names(dif3) <- bar2 %>% filter(anchor==which.anchor[j]) %>% pull(ids)
        dif3L <- multcompView::multcompLetters(dif3)
        l.letters[[j]] <-
          data.frame(id = which.ids,
                     anchor = which.anchor[j],
                     params=which.params[i],
                     letters = dif3L$Letters)
        colnames(l.letters[[j]]) <- c(ids, anchor, params, "letters")
      }
      all.letters[[i]] <- do.call(rbind,l.letters)
    }
    return(do.call(rbind,all.letters))
    
  } else {
    writeLines("Please specify an anchor column.")
  }  
}




# maxraft <- function(chains, burn.in=0, coda)
# Finds the maximum number of MCMC iterations needed 
# across all chains and monitored quantities in coda, and 
# divides by the number of chains.
maxraft <- function(chains, burn.in=0, skip.rows=0,coda){
  
  #codaraft <- window(coda,start=burn.in)
  codaraft = coda
  
  raft<-raftery.diag(coda)
  test<-c()
  for(i in 1:chains){
    # For each chain, grab the maximum number of iterations needed (N) a
    # across all quantities monitored:
    ifelse(length(skip.rows)>1,test[i]<-max(raft[[i]][[2]][-skip.rows,2], na.rm = TRUE),
           test[i]<-max(raft[[i]][[2]][,2], na.rm = TRUE))
  }
  # find max number of iterations required among the chains, 
  # then divided by (# chains) to get iterations per chain:
  max(test)/chains
}

###############################################################################
# Written for work on the Multicomp project
# Updated by Michael Fell 9/10/2018
#   Added an option for an arbitrary function
#   Added more informative error messages
###############################################################################
# A function to summarize output from a JAGS or OpenBUGS model.
coda.fast <- function(coda=NULL, thin=1, FUN=NULL, colname = "optfun", ...){
  
  if(is.null(coda)){
    message("No coda object provided. Summarizing nothing is too philosophical")
    message("a task for this function.")
    stop()
  }
  
  # Get the number of chains
  chains <- length(coda)
  
  codal <- length(coda[[1]][,1])
  
  # Combine chains
  Ftab <- numeric()
  for(i in 1:chains){
    Ftab <- rbind(Ftab, coda[[i]][(0:(codal-1))%%thin==0,])
  }
  
  # mean, sd, 95% CrI table
  pred <- matrix(nrow=dim(Ftab)[2], ncol=5)
  colnames(pred)<-c("mean", "median", "sd","pc2.5","pc97.5")
  
  # Fill table with stats
  pred[,1] <- colMeans(Ftab) #calculate predicted mean RW's
  pred[,2] <- apply(X=Ftab, MARGIN=2, FUN=median, na.rm=TRUE)
  pred[,3] <- apply(X=Ftab, MARGIN=2, FUN=sd,na.rm=TRUE) #calculate sd, apply column wise
  pred[,4] <- apply(X=Ftab, MARGIN=2, FUN=quantile,probs=c(0.025),na.rm=TRUE) 
  pred[,5] <- apply(X=Ftab, MARGIN=2, FUN=quantile,probs=c(0.975),na.rm=TRUE)
  
  pred <- data.frame(pred)
  if(length(rownames(pred)) == length(colnames(coda[[1]]))){
    rownames(pred) <- colnames(coda[[1]])
  }else{
    message("Could not determine row (variable) names from coda.")
  }
  
  # Optional Function
  if(!is.null(FUN))
  {
    placeholder <- tryCatch(
      {
        out <- apply(X=Ftab, MARGIN=2, FUN=FUN, na.rm=TRUE, ...)
        out <- as.matrix(out)
        if(ncol(out) == nrow(pred)){
          out <- t(out)
        }
        
        pred <- cbind(pred, out)
        colnames(pred) <- c("mean", "median", "sd","pc2.5","pc97.5", colname)
      },
      error=function(cond){
        message(paste0("A problem led to an error executing the optional function."))
        message("The result without the added function will be returned.")
        message("Here is the original error:")
        message(cond)
      },
      warning=function(cond){
        message("A warning occurred executing the optional function.")
        message("The result without the added function will be returned.")
        message("Here is the original warning:")
        message(cond)
      },
      finally={
        return(pred)
      }
    )
  }
  
  # Return the summary values
  return(pred)
}

# A function to find initial values for a JAGS or OpenBUGS model.
# Output:
# The output from this function is a list containing two elements. The first
# contains the names of the variables and their indicies. These are useful 
# when using removevars to remove variables that don't need initial values
# in JAGS. The second element contains a list of initial values (this is a 
# list of lists).


initfind <- function(coda, iteration=0, OpenBUGS=FALSE){
  mcmcin <- coda # TODO change all mcmcin to coda in the future MKF 11/27/18
  # If mcmc.list convert to mcmc
  if(is.mcmc.list(mcmcin)==TRUE){
    mcmcin <- mcmc(data=mcmcin, thin=1)
  }
  
  # Get the number of chains
  n.chains <- length(mcmcin)
  
  # get variable names from a list
  var.names <- colnames(mcmcin[[1]])
  var.dims <- dim(mcmcin[[1]])
  if(iteration==0){
    iteration <- var.dims[1]
  }
  
  if(sum(var.names=="deviance")>0){
    var.names <- var.names[-which(var.names=="deviance")]
    var.dims[2] <- var.dims[2]-1 # correct for removing deviance
  }
  
  # Get names and dimension of each variable since the output is a table
  var.names2 <- apply(X=as.matrix(var.names), MARGIN=c(1), FUN=strsplit, split="\\x5B", perl=TRUE)
  var.names2 <- lapply(X=var.names2, FUN=unlist)
  var.names2 <- lapply(X=var.names2, FUN=gsub, pattern="\\x5D", replacement="", perl=TRUE)
  
  # Create a table of names and dimensions
  # Column 1 is the variable me column 2 has the dimensions
  var.info <- matrix(nrow=length(var.names), ncol=3)
  for(i in 1:length(var.names2)){
    if(length(var.names2[[i]]) > 1){
      var.info[i,] <- c(var.names2[[i]], var.names[i])
    }else if(length(var.names2[[i]]) == 1){
      var.info[i,] <- c(var.names2[[i]], 1, var.names[i])
      #print(i)
      #print(var.names2[[i]])
    }else{
      stop("A variable name has incorrect dimensions for parsing.") 
    }
  }
  
  # Get variable names
  unique.names <- unique(var.info[,1])
  initsoutall <- list()
  
  
  for(k in 1:n.chains){
    initsout <- list()
    for(i in 1:length(unique.names)){
      sel <- which(var.info[,1]==unique.names[i])
      #sel2 <- grep(pattern=paste0("^",unique.names[i],"\\x5B"), x=var.names)
      
      # Make sure the above selections worked
      #if(length(sel) != length(sel2)){
      #  stop("Error matching unique variable names with MCMC output")  
      #}
      name.sel <- var.info[sel,3]
      index <- apply(X=as.matrix(var.info[sel,2]), MARGIN=1, FUN=strsplit, split=",", perl=TRUE)
      index <- lapply(X=index, FUN=unlist)
      index <- matrix(data=as.numeric(unlist(index)), nrow=length(index), ncol=length(index[[1]]), byrow=TRUE)
      
      # There are possibly easier ways to do this but they make more assumptions
      dims <- as.numeric(apply(X=index, MARGIN=2, FUN=max))
      variable <- array(data=NA, dim=dims)
      
      # Fill the new variable with the correct values
      for(j in 1:dim(index)[1]){
        # The output into mcmc objects lists names in the order R stacks them
        # in arrays so the single index for the variable references the 
        # correct array location.
        variable[j] <- mcmcin[[k]][iteration, which(colnames(mcmcin[[k]])==name.sel[j])]
      }
      
      # Use dims to produce a new array to store the data
      initsout[[i]] <- variable
    } # End of variable loop
    names(initsout) <- unique.names
    initsoutall[[k]] <- initsout
  } # End of chain loop
  
  listout <- list(unique.names, initsoutall)
  names(listout) <- c("variables", "initials")
  
  # Account for OpenBUGS by outputing 1 dimensional arrays as vectors.
  if(OpenBUGS==TRUE){
    for(i in 1:n.chains){
      for(j in 1:length(listout[[2]][[i]])){
        if(length(dim(listout[[2]][[i]][[j]]))==1){
          listout[[2]][[i]][[j]] <- as.vector(listout[[2]][[i]][[j]])
        }
      }
    }
  }
  
  return(listout)
} # End of function


###############################################################################
#
# Removes specific variables from the initial values
#
###############################################################################

# A function to remove variables that don't need initial values in JAGS.

removevars <- function(initsin, variables){
  n.chains <- length(initsin[[2]])
  n.vars <- 1:length(initsin[[1]])
  n.vars <- n.vars[-variables]
  
  var.names <- initsin[[1]][n.vars]
  
  new.inits <- list()
  for(i in 1:n.chains){
    chain.inits <- list()
    for(k in 1:length(n.vars)){
      chain.inits[[k]] <- initsin[[2]][[i]][[n.vars[k]]] 
    } # End of k loop
    names(chain.inits) <- var.names
    new.inits[[i]] <- chain.inits
  } # End of i loop
  
  output <- list(var.names, new.inits)
  names(output) <- c("variables", "initials")
  
  return(output)
  
} # End of function


