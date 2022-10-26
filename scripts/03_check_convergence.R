#!/usr/bin/env Rscript

library(mcmcplots)


params = c("beta0","beta1","beta1a",
           "beta2","deviance","dYdX","wP","wSd","wSs","wT","wV",
           "wP.monthly","wP.weekly",
           "sig")

load("output_coda/zc_seg.RData")
zc_seg <- zc1

mcmcplot(zc_seg,
         random = 15,
         dir = "./models/convergence",
         filename = "MCMC_seg") 


load("output_coda/zc_ses.RData")
zc_ses <- zc1

mcmcplot(zc_ses,
         random = 15,
         dir = "./models/convergence",
         filename = "MCMC_ses") 


load("output_coda/zc_wjs.RData")
zc_wjs <- zc1

mcmcplot(zc_wjs,
         random = 15,
         dir = "./models/convergence",
         filename = "MCMC_wjs") 


load("output_coda/zc_mpj.RData")
zc_mpj <- zc1

mcmcplot(zc_mpj,
         random = 15,
         dir = "./models/convergence",
         filename = "MCMC_mpj") 


load("output_coda/zc_vcp1.RData")
zc_vcp1 <- zc1

mcmcplot(zc_vcp1,
         random = 15,
         dir = "./models/convergence",
         filename = "MCMC_vcp1") 

load("output_coda/zc_vcp2.RData")
zc_vcp2 <- zc1

mcmcplot(zc_vcp2,
         random = 15,
         dir = "./models/convergence",
         filename = "MCMC_vcp2") 


load("output_coda/zc_vcm1.RData")
zc_vcm1 <- zc1

mcmcplot(zc_vcm1,
         random = 15,
         dir = "./models/convergence",
         filename = "MCMC_vcm1") 


load("output_coda/zc_vcm2.RData")
zc_vcm2 <- zc1

mcmcplot(zc_vcm2,
         random = 15,
         dir = "./models/convergence",
         filename = "MCMC_vcm2") 


load("output_coda/zc_vcs.RData")
zc_vcs <- zc1

mcmcplot(zc_vcs,
         random = 15,
         dir = "./models/convergence",
         filename = "MCMC_vcs")  



