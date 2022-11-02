model{
  
  ################################ split 1 ################################
  for(i in Nstart:Nsplitstart){
    # Likelihood of observed Y values (defined in script file)    
    Y[i] ~ dnorm(mu[i],tau)  
    # resid[i] <- Y[i] - mu[i] # explore residuals to see if distribution is appropriate
    # Replicated data for evaluating model fit
    Y.rep[i] ~ dnorm(mu[i],tau) 
    # Regression (mean) model
    mu[i] <- beta0 + main.effects[i] + squared.terms[i] + interactions[i] 
    
    # Define parts involving main effects, quadratic effects, and 2-way interactions
    main.effects[i]  <- sum(X.effect[,i])
    interactions[i]  <- sum(sum.XX.int[,i])
    squared.terms[i] <- sum(X2.effect[,i])
    
    # Define components of main.effects, squared.terms, and interactions:
    # Main effect parts:
    for(j in 1:Nparms){
      X.effect[j,i]<- beta1[j]*X[j,i]
    }
    
    # Squared terms:
    for(j in 1:2){
      X2.effect[j,i] <- beta1a[j]*pow(X[j,i],2)
    }
    
    # Two-way interaction terms:
    for(j in 1:jlength){
      XX.int[j,i] <- beta2[j,1]*X[ID1[j],i]*X[ID2[j],i] # check if beta2 needs the ,1 (not needed)
      sum.XX.int[j,i] <- sum(XX.int[j,i])
    }
    
    # Creating antecedent covariates; 
    # This matrix of values will end up being used to calculate 
    # the parts involving main effects, interactions, and squared terms
    # in the regression model:
    X[1,i] <- VPDant[i]       ## Also included as squared term
    X[2,i] <- TAant[i]        ## Also included as squared term
    X[3,i] <- PPTant[i]
    X[4,i] <- PAR[Yday[i]]    ## not included in interactions
    X[5,i] <- Sshall_ant[i]
    X[6,i] <- Sdeep_ant[i]
    
    # Computed antecedent values. 
    # PAR is assumed to instantaneously affect ecosystem fluxes, 
    # so no antecedent term calculated
    VPDant[i]     <- sum(VPDtemp[i,]) # summing over all lagged j's
    TAant[i]      <- sum(Tairtemp[i,])
    PPTant[i]     <- sum(PPTtemp[i,])
    Sshall_ant[i] <- sum(Sshalltemp[i,])
    Sdeep_ant[i] <- sum(Sdeeptemp[i,])
    
    # Intermediate weighted values of covariates with influence over 
    # flux over the past few days (or months for ppt)
    for(j in 1:Nlag){
      VPDtemp[i,j]    <- wV[j]*VPD[Yday[i]-j+1] # linking two datasets together
      Tairtemp[i,j]     <- wT[j]*Tair[Yday[i]-j+1]
      Sshalltemp[i,j] <- wSs[j]*Sshall[Yday[i]-j+1]
      Sdeeptemp[i,j] <- wSd[j]*Sdeep[Yday[i]-j+1]
    }
    # For precip (ppt):
    for(j in 1:NlagP){
      PPTtemp[i,j] <- wP[j]*P_temp[i,j] # P_temp is total precip during that block period
      P_temp[i,j] <- sum(P[(Yday[i]-P1[j]):(Yday[i]-P2[j])])
    }
    
    # Calculate net sensitivities (derivative) -- derived quantities
    
    dYdVPD[i] <- beta1[1] + 2*beta1a[1]*VPDant[i] + beta2[1,1]*TAant[i] + beta2[2,1]*PPTant[i] + beta2[3,1]*Sshall_ant[i] + beta2[4,1]*Sdeep_ant[i]
    dYdT[i]   <- beta1[2] + 2*beta1a[2]*TAant[i] + beta2[1,1]*VPDant[i] + beta2[5,1]*PPTant[i] + beta2[6,1]*Sshall_ant[i] + beta2[7,1]*Sdeep_ant[i]
    dYdP[i]   <- beta1[3] + beta2[2,1]*VPDant[i] + beta2[5,1]*TAant[i] + beta2[8,1]*Sshall_ant[i] + beta2[9,1]*Sdeep_ant[i]
    dYdSs[i]  <- beta1[4] + beta2[3,1]*VPDant[i] + beta2[6,1]*TAant[i] + beta2[8,1]*PPTant[i] + beta2[10,1]*Sdeep_ant[i]
    dYdSd[i]  <- beta1[5] + beta2[4,1]*VPDant[i] + beta2[7,1]*TAant[i] + beta2[9,1]*PPTant[i] + beta2[10,1]*Sshall_ant[i]
    
    # Put all net sensitivities into one array, for easy monitoring
    dYdX[i,1] <- dYdVPD[i]
    dYdX[i,2] <- dYdT[i]
    dYdX[i,3] <- dYdP[i]
    dYdX[i,4] <- dYdSs[i]
    dYdX[i,5] <- dYdSd[i]
  }  

  ################################ split 2 ################################
  for(i in Nsplitend:Nend){
    # Likelihood of observed Y values (defined in script file)    
    Y[i] ~ dnorm(mu[i],tau)  
    # resid[i] <- Y[i] - mu[i] # explore residuals to see if distribution is appropriate
    # Replicated data for evaluating model fit
    Y.rep[i] ~ dnorm(mu[i],tau) 
    # Regression (mean) model
    mu[i] <- beta0 + main.effects[i] + squared.terms[i] + interactions[i] 
    
    # Define parts involving main effects, quadratic effects, and 2-way interactions
    main.effects[i]  <- sum(X.effect[,i])
    interactions[i]  <- sum(sum.XX.int[,i])
    squared.terms[i] <- sum(X2.effect[,i])
    
    # Define components of main.effects, squared.terms, and interactions:
    # Main effect parts:
    for(j in 1:Nparms){
      X.effect[j,i]<- beta1[j]*X[j,i]
    }
    
    # Squared terms:
    for(j in 1:2){
      X2.effect[j,i] <- beta1a[j]*pow(X[j,i],2)
    }
    
    # Two-way interaction terms:
    for(j in 1:jlength){
      XX.int[j,i] <- beta2[j,1]*X[ID1[j],i]*X[ID2[j],i] # check if beta2 needs the ,1 (not needed)
      sum.XX.int[j,i] <- sum(XX.int[j,i])
    }
    
    # Creating antecedent covariates; 
    # This matrix of values will end up being used to calculate 
    # the parts involving main effects, interactions, and squared terms
    # in the regression model:
    X[1,i] <- VPDant[i]       ## Also included as squared term
    X[2,i] <- TAant[i]        ## Also included as squared term
    X[3,i] <- PPTant[i]
    X[4,i] <- PAR[Yday[i]]    ## not included in interactions
    X[5,i] <- Sshall_ant[i]
    X[6,i] <- Sdeep_ant[i]
    
    # Computed antecedent values. 
    # PAR is assumed to instantaneously affect ecosystem fluxes, 
    # so no antecedent term calculated
    VPDant[i]     <- sum(VPDtemp[i,]) # summing over all lagged j's
    TAant[i]      <- sum(Tairtemp[i,])
    PPTant[i]     <- sum(PPTtemp[i,])
    Sshall_ant[i] <- sum(Sshalltemp[i,])
    Sdeep_ant[i] <- sum(Sdeeptemp[i,])
    
    # Intermediate weighted values of covariates with influence over 
    # flux over the past few days (or months for ppt)
    for(j in 1:Nlag){
      VPDtemp[i,j]    <- wV[j]*VPD[Yday[i]-j+1] # linking two datasets together
      Tairtemp[i,j]     <- wT[j]*Tair[Yday[i]-j+1]
      Sshalltemp[i,j] <- wSs[j]*Sshall[Yday[i]-j+1]
      Sdeeptemp[i,j] <- wSd[j]*Sdeep[Yday[i]-j+1]
    }
    # For precip (ppt):
    for(j in 1:NlagP){
      PPTtemp[i,j] <- wP[j]*P_temp[i,j] # P_temp is total precip during that block period
      P_temp[i,j] <- sum(P[(Yday[i]-P1[j]):(Yday[i]-P2[j])])
    }
    
    # Calculate net sensitivities (derivative) -- derived quantities
    
    dYdVPD[i] <- beta1[1] + 2*beta1a[1]*VPDant[i] + beta2[1,1]*TAant[i] + beta2[2,1]*PPTant[i] + beta2[3,1]*Sshall_ant[i] + beta2[4,1]*Sdeep_ant[i]
    dYdT[i]   <- beta1[2] + 2*beta1a[2]*TAant[i] + beta2[1,1]*VPDant[i] + beta2[5,1]*PPTant[i] + beta2[6,1]*Sshall_ant[i] + beta2[7,1]*Sdeep_ant[i]
    dYdP[i]   <- beta1[3] + beta2[2,1]*VPDant[i] + beta2[5,1]*TAant[i] + beta2[8,1]*Sshall_ant[i] + beta2[9,1]*Sdeep_ant[i]
    dYdSs[i]  <- beta1[4] + beta2[3,1]*VPDant[i] + beta2[6,1]*TAant[i] + beta2[8,1]*PPTant[i] + beta2[10,1]*Sdeep_ant[i]
    dYdSd[i]  <- beta1[5] + beta2[4,1]*VPDant[i] + beta2[7,1]*TAant[i] + beta2[9,1]*PPTant[i] + beta2[10,1]*Sshall_ant[i]
    
    # Put all net sensitivities into one array, for easy monitoring
    dYdX[i,1] <- dYdVPD[i]
    dYdX[i,2] <- dYdT[i]
    dYdX[i,3] <- dYdP[i]
    dYdX[i,4] <- dYdSs[i]
    dYdX[i,5] <- dYdSd[i]
  }  
  
  
  # Relatively non-informative priors for regression parameters:
  
  # Overall intercept:
  beta0 ~ dnorm(0,0.00001)
  
  # Main effects:
  for(j in 1:Nparms){
    beta1[j] ~ dnorm(0,0.00001)
  }
  
  # Quadratic effects
  for(j in 1:2){
    beta1a[j] ~ dnorm(0,0.00001)
  }
  
  # Two-way interaction effects:
  for(j in 1:jlength){
    beta2[j,1] ~ dnorm(0,0.00001)
  }
  
  # Priors for importance weights for each covariate, "delta" or gamma "trick"
  # for imposing Dirichlet(1) priors for the weigths:  
  for(j in 1:Nlag){
    # Priors for unnormalized weights
    dV[j]    ~ dgamma(1,1)
    dT[j]    ~ dgamma(1,1)
    dSs[j]   ~ dgamma(1,1)
    dSd[j]   ~ dgamma(1,1)
    
    # Compute normalized weights:
    wV[j]    <- dV[j]/sum(dV[])
    wT[j]    <- dT[j]/sum(dT[])
    wSs[j]   <- dSs[j]/sum(dSs[])
    wSd[j]   <- dSd[j]/sum(dSd[])
  }
  
  # Priors for importance weights for precipitation:
  for(j in 1:(NlagP-1)){
    dP[j] ~ dgamma(1,1)
    wP[j+1] <- dP[j]/sum(dP[])  
  }
  # Precipitation in first time step i assumed to have no effect - 
  # moisture effects incorporated via soil water of current week
  wP[1] <- 0    
  
  # Rearrange precip weights into weights at the weekly and monthly scales.
  for(j in 1:4){
    wP.weekly[j] <- wP[j]
  }
  
  for(j in 1:6){
    wP.monthly[j] <- equals(j,1)*sum(wP[1:4]) + (1-equals(j,1))*wP[j+3] # a weight that sums over the first 4 weeks, then gets the monthly weights from before that 
  }
  
  #Prior for standard deviation in data likelihood 
  tau <- pow(sig,-2)
  sig ~ dunif(0,1000)
}