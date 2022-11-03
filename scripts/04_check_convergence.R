#!/usr/bin/env Rscript

library(mcmcplots)

# Create necessary folders if they do not already exist
if(!file.exists("models/convergence")) { dir.create("models/convergence")}
if(!file.exists("models/convergence/seg")) { dir.create("models/convergence/seg")}
if(!file.exists("models/convergence/ses")) { dir.create("models/convergence/ses")}
if(!file.exists("models/convergence/wjs")) { dir.create("models/convergence/wjs")}
if(!file.exists("models/convergence/mpj")) { dir.create("models/convergence/mpj")}
if(!file.exists("models/convergence/vcp")) { dir.create("models/convergence/vcp")}
if(!file.exists("models/convergence/vcm1")) { dir.create("models/convergence/vcm1")}
if(!file.exists("models/convergence/vcm2")) { dir.create("models/convergence/vcm2")}
if(!file.exists("models/convergence/vcs")) { dir.create("models/convergence/vcs")}


# MCMC plots

params = c("beta0","beta1","beta1a",
           "beta2","deviance","dYdX","wP","wSd","wSs","wT","wV",
           "sig")

load("output_coda/coda_all_seg.RData")
zc_seg <- coda_all

mcmcplot(zc_seg,
         random = 15,
         dir = "./models/convergence/seg",
         filename = "MCMC_seg") 


load("output_coda/coda_all_ses.RData")
zc_ses <- coda_all

mcmcplot(zc_ses,
         random = 15,
         dir = "./models/convergence/ses",
         filename = "MCMC_ses") 


load("output_coda/coda_all_wjs.RData")
zc_wjs <- coda_all

mcmcplot(zc_wjs,
         random = 15,
         dir = "./models/convergence/wjs",
         filename = "MCMC_wjs") 


load("output_coda/coda_all_mpj.RData")
zc_mpj <- coda_all

mcmcplot(zc_mpj,
         random = 15,
         dir = "./models/convergence/mpj",
         filename = "MCMC_mpj") 


load("output_coda/coda_all_vcp.RData")
zc_vcp <- coda_all

mcmcplot(zc_vcp,
         random = 15,
         dir = "./models/convergence/vcp",
         filename = "MCMC_vcp1") 

load("output_coda/coda_all_vcm1.RData")
zc_vcm1 <- coda_all

mcmcplot(zc_vcm1,
         random = 15,
         dir = "./models/convergence/vcm1",
         filename = "MCMC_vcm1") 


load("output_coda/coda_all_vcm2.RData")
zc_vcm2 <- coda_all

mcmcplot(zc_vcm2,
         random = 15,
         dir = "./models/convergence/vcm2",
         filename = "MCMC_vcm2") 


load("output_coda/coda_all_vcs.RData")
zc_vcs <- coda_all

mcmcplot(zc_vcs,
         random = 15,
         dir = "./models/convergence/vcs",
         filename = "MCMC_vcs") 

# Raftery diagnostic to determine number of samples needed:

source("./scripts/coda_functions.R")

maxraft1 <- maxraft(chains=3,coda=zc_seg)
maxraft2 <- maxraft(chains=3,coda=zc_ses)
maxraft3 <- maxraft(chains=3,coda=zc_wjs)
maxraft4 <- maxraft(chains=3,coda=zc_mpj)
maxraft5 <- maxraft(chains=3,coda=zc_vcp)
maxraft6 <- maxraft(chains=3,coda=zc_vcm1)
maxraft7 <- maxraft(chains=3,coda=zc_vcm2)
maxraft8 <- maxraft(chains=3,coda=zc_vcs)
raft.out = matrix(data = NA, nrow = 8, ncol = 2)
row.names(raft.out) = c("Seg", "Ses", "Wjs", "Mpj", "Vcp", "Vcm1", "Vcm2", "Vcs")
colnames(raft.out) = c("Total its", "Intercept its")
for(m in 1:8){
  raft.out[m,1] <- eval(parse(text = paste0("maxraft",m)))
}

save(raft.out, file="models/convergence/raft_outstats.csv")

# Gel




