
model{
  for(i in Nstart:Nend){
    # Likelihood for ET data
    ET[i] ~ dnorm(ET.pred[i], tau.ET) # the 'true' ET should be somewhere around the predicted ET, with some precision
    ET.rep[i] ~ dnorm(ET.pred[i], tau.ET)
    ET.pred[i] <- E.model[i] + T.pred[i] # for evaluating model fit, not connected to ET data
    
    # Compute log predictive density for WAIC
    ldx[i] <- logdensity.norm(ET[i], ET.rep[i], tau.ET) # for WAIC
    # Compute predictive density for WAIC:
    dx[i] <- exp(ldx[i]) # for WAIC
    
    # Compute squared difference for calculating posterior predictive loss
    sqdiff[i] <- pow(ET.rep[i]-ET[i],2)

    # Soil evaporation process-based model

    # Calculate aerodynamic resistance to heat transfer (rah)- allow k to vary btwn 0.35-0.42
    # 0.001 is Z0m (m), the momentum soil roughness, [Yang et al., 2008; Stefan et al., 2015]
    # Notes: look for redundant code, and simplify to make it faster
    rah0[i] <- (1/((vk.pred**2)*ws[i]))*((log(Z/0.001))**2) # s/m
    #rah_unstable[i] <- rah0[i]/((1 + Ri[i])**0.75) # aerodynamic resistance to heat transfer s/m, estimated as in Choudhury et al. [1986]
    rah_stable[i] <- ((1 - 0.75*Ri[i])*((log(Z/0.001))**2))/((vk.pred**2)*ws[i])
    rah[i] <- ifelse( Ri[i] > 0, rah_unstable[i], rah_stable[i])

    # calculate parameters for alternative alpha and Bowen
    psi[i] <- psisat.pred * (S[i]/ssat.pred)**(-bch.pred)

    # calculate soil resistance (rss)- allow Arss and Brss to vary
    rss1[i] <- ((fc.pred - sres.pred)/(S[i] - sres.pred)) * rssmin # soil resistance for S > sres
    rss2[i] <- exp(8.206 - 4.255*(S[i]/fc.pred)) # another soil resistance derivation
    rss[i] <- ifelse(S[i] > sres.pred, rss1[i], rss2[i]) # the resistance to the diffusion of vapor in large soil pores

    # Calculate alpha and Bowen ratios
    alpha[i] <- exp((psi[i]*g)/(10000 * Rwv * Tsoil[i])) # Kelvin equation
    bowen[i] <- ifelse(S[i] > fc.pred, 1, (0.5 - 0.5*cos((S[i] * pi)/fc.pred))**2) # wetness function

    # Soil evaporation from E1 (CLM 4.5)
    LE4.5[i] <- ifelse(Tsoil[i] >= 0, bowen[i]*((rho[i]*Cp)/gamma[i])*((alpha[i]*(e.sat[i] - e.a[i]))/rah[i]), 0)
    Esoil4.5[i] <- conv.fact[i]*LE4.5[i]
    # Soil evaporation from E2 (CLM 3.5)
    LE3.5[i] <- ifelse(Tsoil[i] >= 0, ((rho[i]*Cp)/gamma[i])*((alpha[i]*(e.sat[i] - e.a[i]))/(rah[i] + rss[i])), 0)
    Esoil3.5[i] <- conv.fact[i]*LE3.5[i]

    # e.scalar ~ dunif(0,1) - maybe try this if the fit still isn't that great
    E.model[i] <- p*Esoil4.5[i] + (1-p)*Esoil3.5[i] + Eint[i] + Esnow[i]
    #E.model[i] <- Esoil3.5[i] + Eint[i] + Esnow[i]

    # Predicted transpiration, based on assumption that transpiration is proportional
    # to GPP (following Scott and Biederman); slope = 1/WUE,
    # and WUE is informed by Ecostress WUE precision
    # Note the proportionality term (slope) varies temporally, according to defined
    # time "blocks".
    T.pred[i] <- (1/WUE.pred[i])*GPP[i]
    # T/ET ratio
    ET.int[i] <- ifelse(ET.pred[i]==0, 0.0000000000000001, ET.pred[i]) # intermediate calculated to ensure the denominator is not 0
    T.ratio[i] <- ifelse(ET.pred[i]==0, 0, T.pred[i]/ET.int[i])
    
    # Likelihood (? not really) of predicted WUE values
    # WUE.log when combining with ETpart model, transform WUE to log scale
    WUE.pred[i] <- exp(log.WUE[i])
    #log.WUE[i] ~ dnorm(mu.log.WUE[i], tau.log.WUE)
    # Regression (mean) model
    log.WUE[i] <- beta0 + main.effects[i] + squared.terms[i] + interactions[i] 
    
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
    for(j in 1:Nparms){
      X2.effect[j,i] <- beta1a[j]*pow(X[j,i],2)
    }
    
    # Two-way interaction terms:
    for(j in 1:jlength){
      XX.int[j,i] <- beta2[j,1]*X[ID1[j],i]*X[ID2[j],i] # check if beta2 needs the ,1 (not needed)
      sum.XX.int[j,i] <- sum(XX.int[j,i])
    }
    
    # This matrix of values will end up being used to calculate 
    # the parts involving main effects, interactions, and squared terms
    # in the regression model:
    X[1,i] <- VPD[Yday[i]]       ## Also included as squared term
    X[2,i] <- Tair[Yday[i]]        ## Also included as squared term
    X[3,i] <- PPTant[i] # P is the only antecedent term
    X[4,i] <- PAR[Yday[i]]
    X[5,i] <- Sshall[Yday[i]]
    X[6,i] <- Sdeep[Yday[i]]
    
    # Computed antecedent values. 
    PPTant[i] <- sum(PPTtemp[i,])
    
    # For precip (ppt):
    for(j in 1:NlagP){
      PPTtemp[i,j] <- wP[j]*P_temp[i,j] # P_temp is total precip during that block period
      P_temp[i,j] <- sum(P[(Yday[i]-P1[j]):(Yday[i]-P2[j])])
    }
    
    # Calculate net sensitivities (derivative) -- derived quantities
    
    dYdVPD[i] <- beta1[1] + 2*beta1a[1]*VPD[Yday[i]] + beta2[1,1]*Tair[Yday[i]] + beta2[2,1]*PPTant[i] + beta2[3,1]*Sshall[Yday[i]] + beta2[4,1]*Sdeep[Yday[i]]+ beta2[11,1]*PAR[Yday[i]]
    dYdT[i]   <- beta1[2] + 2*beta1a[2]*Tair[Yday[i]] + beta2[1,1]*VPD[Yday[i]] + beta2[5,1]*PPTant[i] + beta2[6,1]*Sshall[Yday[i]] + beta2[7,1]*Sdeep[Yday[i]]+ beta2[12,1]*PAR[Yday[i]]
    dYdP[i]   <- beta1[3] + 2*beta1a[2]*PPTant[i] + beta2[2,1]*VPD[Yday[i]] + beta2[5,1]*Tair[Yday[i]] + beta2[8,1]*Sshall[Yday[i]] + beta2[9,1]*Sdeep[Yday[i]]+ beta2[13,1]*PAR[Yday[i]]
    dYdPAR[i]   <- beta1[4] + 2*beta1a[2]*PAR[Yday[i]] + beta2[11,1]*VPD[Yday[i]] + beta2[12,1]*Tair[Yday[i]] + beta2[13,1]*PPTant[i] + beta2[14,1]*Sshall[Yday[i]] + beta2[15,1]*Sdeep[Yday[i]]
    dYdSs[i]  <- beta1[5] + 2*beta1a[2]*Sshall[Yday[i]] + beta2[3,1]*VPD[Yday[i]] + beta2[6,1]*Tair[Yday[i]] + beta2[8,1]*PPTant[i] + beta2[10,1]*Sdeep[Yday[i]]+ beta2[14,1]*PAR[Yday[i]]
    dYdSd[i]  <- beta1[6] + 2*beta1a[2]*Sdeep[Yday[i]] + beta2[4,1]*VPD[Yday[i]] + beta2[7,1]*Tair[Yday[i]] + beta2[9,1]*PPTant[i] + beta2[10,1]*Sshall[Yday[i]] + beta2[15,1]*PAR[Yday[i]]
    
    # Put all net sensitivities into one array, for easy monitoring
    dYdX[i,1] <- dYdVPD[i]
    dYdX[i,2] <- dYdT[i]
    dYdX[i,3] <- dYdP[i]
    dYdX[i,4] <- dYdPAR[i]
    dYdX[i,5] <- dYdSs[i]
    dYdX[i,6] <- dYdSd[i]
  }

  # Intercepted E
  Eint[1] <- (P[1])*(1 - exp(-k.pred*(LAI[1])))
  for(i in 2:Nend){
    Eint[i] <- (P.unscaled[i] + P.unscaled[i-1])*(1 - exp(-k.pred*(LAI[i])))
  }

  # given p (proportion "contribution" pf Esoil4.5 to estimated Esoil) a prior, or set = 0.5 in data list
  p ~ dunif(0,1)


  # for(i in 1:(Nstart-1)){
  #   # This assumes that the precision for the "predicted" or "true" WUE of the site varies
  #   # around the WUE derived from ECOSTRESS for that area, with some uncertainty
  #   WUE.pred[i] ~ dnorm(WUE.pred[i], tau.WUE)T(0,) # WUE.pred[i-1], tau.ecostress T(0,)
  # }
  
  # Relatively non-informative priors for regression parameters:
  
  # Overall intercept:
  beta0 ~ dnorm(0,0.00001)
  
  # Main effects:
  for(j in 1:Nparms){
    beta1[j] ~ dnorm(0,0.00001)
  }
  
  # Quadratic effects
  for(j in 1:Nparms){
    beta1a[j] ~ dnorm(0,0.00001)
  }
  
  # Two-way interaction effects:
  for(j in 1:jlength){
    beta2[j,1] ~ dnorm(0,0.00001)
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
  #tau <- pow(sig,-2)
  #sig ~ dunif(0,1000)
  
  # sum across days
  Dsum <- sum(sqdiff[Nstart:Nend])
  
  # Compute quantities for calculating Bayesian R2
  var.pred <- pow(sd(ET.pred[Nstart:Nend]),2)
  var.resid <- 1/tau.ET
  R2 <- var.pred/(var.pred + var.resid)
  
  # Priors for ET and WUE:
  tau.ET ~ dgamma(0.1,0.1) # since this is associated with the data model for ET.
  sig.ET <- 1/sqrt(tau.ET)
  tau.log.WUE ~ dgamma(0.01,0.01)
  #sig.WUE ~ dunif(0,20)
  #tau.WUE <- pow(sig.WUE,-2)
  #sig.ecostress ~ dunif(0,10)
  #tau.ecostress <- pow(sig.ecostress,-2)
  
  # Priors for stochastic parameters for E equations
  vk.pred ~ dunif(0.35, 0.42) # the von Karman constant is usually 0.40
  bch.pred ~ dnorm(bch, 0.40)T(0,) # Clapp and Hornberger parameter
  k.pred ~ dnorm(0.5, 10)T(0,) # k = 0.5 # decay function k for intercepted E
  fc.pred ~ dnorm(fc, 200)T(S.min,1) # upper limit 1
  # pick a precision that's less precise
  # Look at clay across all sites, maybe base precision off that, make more flexible
  # decrease precision, increase variance, but check with sites
  sres.pred ~ dnorm(sres, 80000)T(0,S.min)# residual soil moisture
  # Plot ssat to see what the range is
  ssat.pred ~ dnorm(ssat, 50)T(0,) # soil moisture at saturation
  #check this equation- it gives negative numbers when fsand is above .38, which is weird (should it be positive???)
  psisat.pred ~ dnorm(psisat, 0.015)T(-1000,0)  # parameterized air entry pressure, in mm of water, check equation and limits
}