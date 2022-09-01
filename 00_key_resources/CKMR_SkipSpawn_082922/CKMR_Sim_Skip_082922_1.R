library(ggplot2)
library(doParallel)
library(gridExtra)
theme_set(theme_bw())
#source("D:\\Google Drive\\Research\\R Projects\\CKMR Simulation\\CKMR_Func_Sim_Skip_082922_1.R")

out = NULL
nsim = 20
nthreads = 5 #Note that the loop will actually require nthreads * nchains threads (e.g., 10 * 3 = 30 threads total)

Ntarget = 4000 #target overall N
Prop_f = 0.5 #proportion of pop that's female
phi = 0.8 #survival rate
mat = 3 #age of maturity
RunYears = 75 #number of years to run it for
StartYear = 1950 #to transform into "real" years for the convenience of my brain
a = 2 #periodicity of skipped spawning (e.g., age=2 is every other year)
psi = 0.8

Nsample = 2000
SampYear = 50
nYears = 1

ModelLocation = "D:\\Google Drive\\Research\\Postdoc\\Meetings\\CKMRgrp_082422\\Models\\"

cl = makeCluster(nthreads)
registerDoParallel(cl)    
  
out = foreach(i = seq(nsim), .combine = rbind, .packages = c("tidyverse", "rjags", "jagsUI", "runjags", "postpack")) %dopar%
#for(i in seq(nsim))
{
  #cat(paste(i, "\r"))
  ind.mat = CKMR.sim.pop.skipspawn(Ntarget = Ntarget, #target overall N
                                   Prop_f = Prop_f, #proportion of pop that's female
                                   phi = phi, #survival rate
                                   mat = mat, #age of maturity
                                   RunYears = RunYears, #number of years to run it for
                                   StartYear = StartYear, #to transform into "real" years for the convenience of my brain
                                   a = a, #periodicity of skipped spawning (e.g., age=2 is every other year)
                                   psi = psi) #proportion of population that engages in skipped spawning
  
  samp.pairwise = CKMR.sim.samp.skipspawn(ind.mat = ind.mat, Nsample = Nsample, SampYear = SampYear, nYears = nYears)
  
  est.out = CKMR.est.skipspawn(samp.pairwise, psi.fix = psi, a.fix = a,   
                               jags.model_location = ModelLocation)
  
  Nf.truth = sum(sapply(strsplit(ind.mat[SampYear,], "_"), `[`, 3) == "F" & 
                   ((SampYear + StartYear + 0) - as.numeric(sapply(strsplit(ind.mat[SampYear,], "_"), `[`, 2)) >= (mat)), na.rm = T)

  Nm.truth = sum(sapply(strsplit(ind.mat[SampYear,], "_"), `[`, 3) == "M" & 
                    ((SampYear + StartYear + 0) - as.numeric(sapply(strsplit(ind.mat[SampYear,], "_"), `[`, 2)) >= (mat)), na.rm = T)

  temp = data.frame(N = i, 
             Nfmean = est.out[1,"mean"], Nf50 = est.out[1,"Q50"],Nf2.5 = est.out[1,"Q2.5"],Nf97.5 = est.out[1,"Q97.5"],
             NfTruth = Nf.truth,
             Nmmean = est.out[2,"mean"], Nm50 = est.out[2,"Q50"],Nm2.5 = est.out[2,"Q2.5"],Nm97.5 = est.out[2,"Q97.5"],
             NmTruth = Nm.truth,
             psimean = est.out[4,"mean"], psi50 = est.out[4,"Q50"],psi2.5 = est.out[4,"Q2.5"],psi97.5 = est.out[4,"Q97.5"],
             phimean = est.out[3,"mean"], phi50 = est.out[3,"Q50"],phi2.5 = est.out[3,"Q2.5"],phi97.5 = est.out[3,"Q97.5"])

  return(temp)
}

stopCluster(cl)

out$NfBias = (out$Nf50 - out$NfTruth) / out$NfTruth
out$NmBias = (out$Nm50 - out$NmTruth) / out$NmTruth
out$phiBias = (out$phi50 - phi) / phi
out$psiBias = (out$psi50 - psi) / psi

plot.bias = ggplot(out)+
  geom_violin(aes(x = "Nf", y = NfBias), fill = 'red', alpha = 0.5, trim = F)+
  geom_violin(aes(x = "Nm", y = NmBias), fill = 'blue', alpha = 0.5, trim = F)+
  geom_violin(aes(x = "phi", y = phiBias), fill = 'darkgreen', alpha = 0.5, trim = F)+
  geom_violin(aes(x = "psi", y = psiBias), fill = 'purple', alpha = 0.5, trim = F)+
  geom_hline(aes(yintercept = 0), lty = 2)+
  scale_x_discrete("")+
  scale_y_continuous("Bias")+
  coord_cartesian(ylim = c(-2,2))

plot.est.N = ggplot(out)+
  geom_pointrange(aes(x = "Nf", y = Nf50, ymin = Nf2.5, ymax = Nf97.5), 
                  position = position_jitter(width = 0.2, height = 0), alpha = 0.3)+
  geom_point(aes(x = "Nf", y = NfTruth), 
                  position = position_jitter(width = 0.2, height = 0), col = 'red')+
  geom_pointrange(aes(x = "Nm", y = Nm50, ymin = Nm2.5, ymax = Nm97.5), 
                  position = position_jitter(width = 0.2, height = 0), alpha = 0.3)+
  geom_point(aes(x = "Nm", y = NmTruth), 
             position = position_jitter(width = 0.2, height = 0), col = 'blue')+
  scale_x_discrete("")+
  scale_y_continuous("Value")+
  coord_cartesian(ylim = c(0, Ntarget))

plot.est.par = ggplot(out)+
  geom_pointrange(aes(x = "phi", y = phi50, ymin = phi2.5, ymax = phi97.5), 
                  position = position_jitter(width = 0.2, height = 0), alpha = 0.3)+
  geom_point(aes(x = "phi", y = phi), 
             position = position_jitter(width = 0.2, height = 0), col = 'darkgreen')+
  geom_pointrange(aes(x = "psi", y = psi50, ymin = psi2.5, ymax = psi97.5), 
                  position = position_jitter(width = 0.2, height = 0), alpha = 0.3)+
  geom_point(aes(x = "psi", y = psi), 
             position = position_jitter(width = 0.2, height = 0), col = 'purple')+
  scale_x_discrete("")+
  scale_y_continuous("Value")+
  coord_cartesian(ylim = c(0,1))
  
grid.arrange(plot.bias, plot.est.N, plot.est.par,
             layout_matrix = matrix(c(1,1,2,3), nrow = 2, byrow = T))


 # png(paste(ModelLocation,"BenModel.png"), height = 5, width = 8, units = 'in', res = 600)
 # grid.arrange(plot.bias, plot.est.N, plot.est.par,
 #              layout_matrix = matrix(c(1,1,2,3), nrow = 2, byrow = T))
 # dev.off()

