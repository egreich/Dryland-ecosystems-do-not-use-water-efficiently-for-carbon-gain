# Calculate ant variables from coda samples, for each day in each site's timeseries
for(d in 1:Ndays){
  for(i in 1:Niters){
    # loop through coda iterations (use 3000-4000 samples, total), grab sampled weigths (w's), 
    # use observed daily covariate data:
    Tant[d,i] <- sum(w1[i]*Temp[d] + w2[i]*Temp[d-1]+w3[i]*Temp[d-2]...)
    # Do for other ant variables...
  }
  # Get posterior means for each day:
  Tant.daily[d] = mean(Tant[d,])
}
# Compute overall mean for each ant variable; describes the "average" conditions
# for each variable.
Tant.mean = mean(Tant.daily[])

# Compute lower and upper 2.5th percentiles
Tant.lower = quantile(Tant.daily,probs=c(0.025))
Tant.upper = quantile(Tant.daily,probs=c(0.975))



## Calculate net sensitivites as functions of covariates.
## Simulate values for the ant covariate that span range realisitic values
Tant.vals = seq(from=Tant.lower,to=Tant.upper,length = 10)  # Do for other variables
Sshall.ant.vals = seq(from=Sshall.ant.lower,to=Sshall.ant.upper,length = 10)  # Do for other variables

for(t in 1:10){
  for(s in 1:10){
    for(i in 1:Niter){
    dYdVPD[t,s,i] <- beta1[1,i] + 2*beta1a[1,i]*VPDant.mean + beta2[1,1,i]*TAant.vals[t] + beta2[2,1,i]*PPTant.mean + 
      beta2[3,1,i]*Sshall.ant.vals[s] + beta2[4,1,i]*Sdeep_ant.mean + beta2[11,1,i]*PAR_ant.mean
    }
  }
    
}
dYdVPD[i] <- beta1[1] + 2*beta1a[1]*VPDant[i] + beta2[1,1]*TAant[i] + beta2[2,1]*PPTant[i] + beta2[3,1]*Sshall_ant[i] + beta2[4,1]*Sdeep_ant[i] + beta2[11,1]*PAR_ant[i]