
model{
  for(i in Nstart:N){
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
    Esoil4.5[i] <- conv.fact*LE4.5[i]
    # Soil evaporation from E2 (CLM 3.5)
    LE3.5[i] <- ifelse(Tsoil[i] >= 0, ((rho[i]*Cp)/gamma[i])*((alpha[i]*(e.sat[i] - e.a[i]))/(rah[i] + rss[i])), 0)
    Esoil3.5[i] <- conv.fact*LE3.5[i]

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
    
  }
    
  # Intercepted E
  Eint[1] <- (P[1])*(1 - exp(-k.pred*(LAI[1])))
  for(i in 2:N){
    Eint[i] <- (P.unscaled[i] + P.unscaled[i-1])*(1 - exp(-k.pred*(LAI[i])))
    
    # This assumes that the precision for the "predicted" or "true" WUE of the site varies
    # around the WUE from the days around it, with some uncertainty
    WUE.pred[i] ~ dnorm(WUE.pred[i-1], tau.WUE)T(0,)
  }

  # given p (proportion "contribution" pf Esoil4.5 to estimated Esoil) a prior, or set = 0.5 in data list
  p ~ dunif(0,1)


  # Priors for "initial conditions" for WUE.
  WUE.pred[1] ~ dunif(0,30)
  
  #Prior for standard deviation in data likelihood 
  #tau <- pow(sig,-2)
  #sig ~ dunif(0,1000)
  
  # sum across days
  Dsum <- sum(sqdiff[Nstart:N])
  
  # Compute quantities for calculating Bayesian R2
  var.pred <- pow(sd(ET.pred[Nstart:N]),2)
  var.resid <- 1/tau.ET
  R2 <- var.pred/(var.pred + var.resid)
  
  # Priors for ET and WUE:
  tau.ET ~ dgamma(0.1,0.1) # since this is associated with the data model for ET.
  sig.ET <- 1/sqrt(tau.ET)
  #tau.log.WUE ~ dgamma(0.01,0.01)
  sig.WUE ~ dunif(0,20)
  tau.WUE <- pow(sig.WUE,-2)
  
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