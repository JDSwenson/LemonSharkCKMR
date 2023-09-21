## Diagnostics
library(coda)
library(tidyverse)
library(R2jags)
library(postpack)

rm(list=ls())

#############################################################################
###############Load posterior samples - Objective 2##########################
#############################################################################

MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.output/"







################################################################################
############################## Read in files ###################################
################################################################################

#------------------------Population growth----------------------
#Slight population decline, narrow lambda prior: scenario_2.2.1
#All juvs
slight.pop.decline_all.juvs <- readRDS(paste0(MCMC_location, "CKMR_modelout_06Aug2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_2.2.1_annual.model_sample.all.juvenile.ages"))

#Just YOY
slight.pop.decline_target.YOY <- readRDS(paste0(MCMC_location, "CKMR_modelout_06Aug2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_2.2.1_annual.model_target.YOY"))

#all ages
slight.pop.decline_all.ages <- readRDS(paste0(MCMC_location, "CKMR_modelout_06Aug2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_2.2.1_annual.model_sample.ALL.ages"))


#Slight population increase, narrow lambda prior: scenario_2.2.2
#All juvs
slight.pop.increase_all.juvs <- readRDS(paste0(MCMC_location, "CKMR_modelout_06Aug2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_2.2.2_annual.model_sample.all.juvenile.ages"))

#Just YOY
slight.pop.increase_target.YOY <- readRDS(paste0(MCMC_location, "CKMR_modelout_06Aug2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_2.2.2_annual.model_target.YOY"))

#All ages
slight.pop.increase_all.ages <- readRDS(paste0(MCMC_location, "CKMR_modelout_06Aug2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_2.2.2_annual.model_sample.ALL.ages"))


#Severe population decrease - wide lambda prior: scenario_2.3.3
#All juvs
severe.pop.decrease_all.juvs <- readRDS(paste0(MCMC_location, "CKMR_modelout_06Aug2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_2.3.3_annual.model_sample.all.juvenile.ages"))

#Just YOY
severe.pop.decrease_target.YOY <- readRDS(paste0(MCMC_location, "CKMR_modelout_06Aug2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_2.3.3_annual.model_target.YOY"))

#All ages
severe.pop.decrease_all.ages <- readRDS(paste0(MCMC_location, "CKMR_modelout_06Aug2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_2.3.3_annual.model_sample.ALL.ages"))


#Stable population
#All juvs
stable.pop_all.juvs <- readRDS(paste0(MCMC_location, "CKMR_modelout_08Sep2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_2.2.4_annual.model_sample.all.juvenile.ages"))

#Just YOY
stable.pop_target.YOY <- readRDS(paste0(MCMC_location, "CKMR_modelout_08Sep2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_2.2.4_annual.model_target.YOY"))

#All ages
stable.pop_all.ages <- readRDS(paste0(MCMC_location, "CKMR_modelout_08Sep2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_2.2.4_annual.model_sample.ALL.ages"))


#------------------------Multiennial breeding----------------------
### 100% biennial breeders: scenario_3.1.2
#All juvs
biennial_all.juvs <- readRDS(paste0(MCMC_location, "CKMR_modelout_06Aug2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_3.1.2_multiennial.model_sample.all.juvenile.ages"))

#Just YOY
biennial_target.YOY <- readRDS(paste0(MCMC_location, "CKMR_modelout_06Aug2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_3.1.2_multiennial.model_target.YOY"))

#All ages
biennial_all.ages <- readRDS(paste0(MCMC_location, "CKMR_modelout_06Aug2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_3.1.2_multiennial.model_sample.ALL.ages"))

### 100% annual breeders: scenario_3.7.2
#All juvs
annual_all.juvs <- readRDS(paste0(MCMC_location, "CKMR_modelout_08Sep2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_3.7.2_multiennial.model_sample.all.juvenile.ages"))

#Just YOY
annual_target.YOY <- readRDS(paste0(MCMC_location, "CKMR_modelout_08Sep2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_3.7.2_multiennial.model_target.YOY"))

#All ages
annual_all.ages <- readRDS(paste0(MCMC_location, "CKMR_modelout_08Sep2023_Seeds2022.04.15_1.5_prop.sampled_thin20_draw40000_burn50000_scenario_3.7.2_multiennial.model_sample.ALL.ages"))


################################################################################
################## Cross correlation among parameters ##########################
################################################################################

#If running a single iteration stored in sims.list
crosscorr(sims.list.1[[1]]$samples) #T0
crosscorr(sims.list.1[[2]]$samples) #T0-10
crosscorr(sims.list.1[[3]]$samples) #present-5
crosscorr(sims.list.1[[4]]$samples) #present

#If running just a few iterations
for(cr in 1:length(slight.pop.decline_t0)){
  s2.2.1_cross.list.J[[cr]] <- crosscorr(slight.pop.decline_t0[[cr]]$samples)
  
  s2.2.2_cross.list.J[[cr]] <- crosscorr(slight.pop.decline_t0_10[[cr]]$samples)
  
  s2.3.3_cross.list.J[[cr]] <- crosscorr(slight.pop.decline_present_5[[cr]]$samples)
  
    s2.2.4_cross.list.J[[cr]] <- crosscorr(slight.pop.decline_present[[cr]]$samples)

}

(t0_cross.mean.J <- as.data.frame(Reduce("+", s2.2.1_cross.list.J) / length(s2.2.1_cross.list.J)))
  
(t0_10_cross.mean.J <- as.data.frame(Reduce("+", s2.2.2_cross.list.J) / length(s2.2.2_cross.list.J)))

(present_5_cross.mean.J <- as.data.frame(Reduce("+", s2.3.3_cross.list.J) / length(s2.3.3_cross.list.J)))

(present_cross.mean.J <- as.data.frame(Reduce("+", s2.2.4_cross.list.J) / length(s2.2.4_cross.list.J)))

t0_cross.mean.J #T0
t0_10_cross.mean.J #T0-10
present_5_cross.mean.J #present-5
present_cross.mean.J #present

#-----------------Population growth-------------------
cross.temp <- NULL
s2.2.1_cross.list.A = s2.2.1_cross.list.J = s2.2.1_cross.list.Y <- list()
s2.2.2_cross.list.A = s2.2.2_cross.list.J = s2.2.2_cross.list.Y <- list()
s2.3.3_cross.list.A = s2.3.3_cross.list.J = s2.3.3_cross.list.Y <- list()
s2.2.4_cross.list.A = s2.2.4_cross.list.J = s2.2.4_cross.list.Y <- list()


for(cr in 1:length(slight.pop.decline_all.juvs)){
  s2.2.1_cross.list.A[[cr]] <- crosscorr(slight.pop.decline_all.ages[[cr]]$samples)
  s2.2.1_cross.list.J[[cr]] <- crosscorr(slight.pop.decline_all.juvs[[cr]]$samples)
  s2.2.1_cross.list.Y[[cr]] <- crosscorr(slight.pop.decline_target.YOY[[cr]]$samples)
  
  s2.2.2_cross.list.A[[cr]] <- crosscorr(slight.pop.increase_all.ages[[cr]]$samples)
  s2.2.2_cross.list.J[[cr]] <- crosscorr(slight.pop.increase_all.juvs[[cr]]$samples)
  s2.2.2_cross.list.Y[[cr]] <- crosscorr(slight.pop.increase_target.YOY[[cr]]$samples)
  
  s2.3.3_cross.list.A[[cr]] <- crosscorr(severe.pop.decrease_all.ages[[cr]]$samples)
  s2.3.3_cross.list.J[[cr]] <- crosscorr(severe.pop.decrease_all.juvs[[cr]]$samples)
  s2.3.3_cross.list.Y[[cr]] <- crosscorr(severe.pop.decrease_target.YOY[[cr]]$samples)
  
  s2.2.4_cross.list.A[[cr]] <- crosscorr(stable.pop_all.ages[[cr]]$samples)
  s2.2.4_cross.list.J[[cr]] <- crosscorr(stable.pop_all.juvs[[cr]]$samples)
  s2.2.4_cross.list.Y[[cr]] <- crosscorr(stable.pop_target.YOY[[cr]]$samples)
  
  
}

##Mean cross-correlation values
(s2.2.1_cross.mean.A <- as.data.frame(Reduce("+", s2.2.1_cross.list.A) / length(s2.2.1_cross.list.A)) %>% 
  mutate(sampling.scheme = "sample all ages", scenario = "scenario_2.2.1") %>% 
  rownames_to_column(var = "parameter"))

(s2.2.1_cross.mean.J <- as.data.frame(Reduce("+", s2.2.1_cross.list.J) / length(s2.2.1_cross.list.J)) %>% 
    mutate(sampling.scheme = "sample all juveniles", scenario = "scenario_2.2.1") %>% 
  rownames_to_column(var = "parameter"))

(s2.2.1_cross.mean.Y <- as.data.frame(Reduce("+", s2.2.1_cross.list.Y) / length(s2.2.1_cross.list.Y)) %>% 
    mutate(sampling.scheme = "target YOY", scenario = "scenario_2.2.1") %>% 
  rownames_to_column(var = "parameter"))


(s2.2.2_cross.mean.A <- as.data.frame(Reduce("+", s2.2.2_cross.list.A) / length(s2.2.2_cross.list.A)) %>% 
    mutate(sampling.scheme = "sample all ages", scenario = "scenario_2.2.2") %>% 
  rownames_to_column(var = "parameter"))

(s2.2.2_cross.mean.J <- as.data.frame(Reduce("+", s2.2.2_cross.list.J) / length(s2.2.2_cross.list.J)) %>% 
    mutate(sampling.scheme = "sample all juveniles", scenario = "scenario_2.2.2") %>% 
    rownames_to_column(var = "parameter"))

(s2.2.2_cross.mean.Y <- as.data.frame(Reduce("+", s2.2.2_cross.list.Y) / length(s2.2.2_cross.list.Y)) %>% 
    mutate(sampling.scheme = "target YOY", scenario = "scenario_2.2.2") %>% 
    rownames_to_column(var = "parameter"))


(s2.3.3_cross.mean.A <- as.data.frame(Reduce("+", s2.3.3_cross.list.A) / length(s2.3.3_cross.list.A)) %>% 
    mutate(sampling.scheme = "sample all ages", scenario = "scenario_2.3.3") %>% 
    rownames_to_column(var = "parameter"))

(s2.3.3_cross.mean.J <- as.data.frame(Reduce("+", s2.3.3_cross.list.J) / length(s2.3.3_cross.list.J)) %>% 
    mutate(sampling.scheme = "sample all juveniles", scenario = "scenario_2.3.3") %>% 
    rownames_to_column(var = "parameter"))

(s2.3.3_cross.mean.Y <- as.data.frame(Reduce("+", s2.3.3_cross.list.Y) / length(s2.3.3_cross.list.Y)) %>% 
    mutate(sampling.scheme = "target YOY", scenario = "scenario_2.3.3") %>% 
    rownames_to_column(var = "parameter"))


(s2.2.4_cross.mean.A <- as.data.frame(Reduce("+", s2.2.4_cross.list.A) / length(s2.2.4_cross.list.A)) %>% 
    mutate(sampling.scheme = "sample all ages", scenario = "scenario_2.2.4") %>% 
    rownames_to_column(var = "parameter"))

(s2.2.4_cross.mean.J <- as.data.frame(Reduce("+", s2.2.4_cross.list.J) / length(s2.2.4_cross.list.J)) %>% 
    mutate(sampling.scheme = "sample all juveniles", scenario = "scenario_2.2.4") %>% 
    rownames_to_column(var = "parameter"))

(s2.2.4_cross.mean.Y <- as.data.frame(Reduce("+", s2.2.4_cross.list.Y) / length(s2.2.4_cross.list.Y)) %>% 
    mutate(sampling.scheme = "target YOY", scenario = "scenario_2.2.4") %>% 
    rownames_to_column(var = "parameter"))

pop.growth_all <- bind_rows(s2.2.1_cross.mean.A, s2.2.1_cross.mean.J, s2.2.1_cross.mean.Y, s2.2.2_cross.mean.A, s2.2.2_cross.mean.J, s2.2.2_cross.mean.Y, s2.3.3_cross.mean.A, s2.3.3_cross.mean.J, s2.3.3_cross.mean.Y, s2.2.4_cross.mean.A, s2.2.4_cross.mean.J, s2.2.4_cross.mean.Y)


#saveRDS(pop.growth_all, file = paste0(MCMC_location, "cross.corr_pop.growth_all_20230919"))


#-----------------Multiennial breeding-------------------
cross.temp <- NULL
s3.1.2_cross.list.A = s3.1.2_cross.list.J = s3.1.2_cross.list.Y <- list()
s3.7.2_cross.list.A = s3.7.2_cross.list.J = s3.7.2_cross.list.Y <- list()


for(cr in 1:length(biennial_all.juvs)){
  s3.1.2_cross.list.A[[cr]] <- crosscorr(biennial_all.ages[[cr]]$samples)
  s3.1.2_cross.list.J[[cr]] <- crosscorr(biennial_all.juvs[[cr]]$samples)
  s3.1.2_cross.list.Y[[cr]] <- crosscorr(biennial_target.YOY[[cr]]$samples)
  
  s3.7.2_cross.list.A[[cr]] <- crosscorr(annual_all.ages[[cr]]$samples)
  s3.7.2_cross.list.J[[cr]] <- crosscorr(annual_all.juvs[[cr]]$samples)
  s3.7.2_cross.list.Y[[cr]] <- crosscorr(annual_target.YOY[[cr]]$samples)
  
  
}

##Mean cross-correlation values
(s3.1.2_cross.mean.A <- as.data.frame(Reduce("+", s3.1.2_cross.list.A) / length(s3.1.2_cross.list.A)) %>% 
    mutate(sampling.scheme = "sample all ages", scenario = "scenario_3.1.2") %>% 
    rownames_to_column(var = "parameter"))

(s3.1.2_cross.mean.J <- as.data.frame(Reduce("+", s3.1.2_cross.list.J) / length(s3.1.2_cross.list.J)) %>% 
    mutate(sampling.scheme = "sample all ages", scenario = "scenario_3.1.2") %>% 
    rownames_to_column(var = "parameter"))

(s3.1.2_cross.mean.Y <- as.data.frame(Reduce("+", s3.1.2_cross.list.Y) / length(s3.1.2_cross.list.Y)) %>% 
    mutate(sampling.scheme = "sample all ages", scenario = "scenario_3.1.2") %>% 
    rownames_to_column(var = "parameter"))


(s3.7.2_cross.mean.A <- as.data.frame(Reduce("+", s3.7.2_cross.list.A) / length(s3.7.2_cross.list.A)) %>% 
    mutate(sampling.scheme = "sample all ages", scenario = "scenario_3.7.2") %>% 
    rownames_to_column(var = "parameter"))

(s3.7.2_cross.mean.J <- as.data.frame(Reduce("+", s3.7.2_cross.list.J) / length(s3.7.2_cross.list.J)) %>% 
    mutate(sampling.scheme = "sample all ages", scenario = "scenario_3.7.2") %>% 
    rownames_to_column(var = "parameter"))

(s3.7.2_cross.mean.Y <- as.data.frame(Reduce("+", s3.7.2_cross.list.Y) / length(s3.7.2_cross.list.Y)) %>% 
    mutate(sampling.scheme = "sample all ages", scenario = "scenario_3.7.2") %>% 
    rownames_to_column(var = "parameter"))


multiennial.model_all <- bind_rows(s3.1.2_cross.mean.A, s3.1.2_cross.mean.J, s3.1.2_cross.mean.Y, s3.7.2_cross.mean.A, s3.7.2_cross.mean.J, s3.7.2_cross.mean.Y)


#saveRDS(multiennial.model_all, file = paste0(MCMC_location, "cross.corr_multiennial.model_all_20230919"))


s3.3.2Y[[1]]
summary(s3.3.2Y[[1]])
plot(s3.3.2Y[[1]])
lattice::densityplot(s3.3.2Y[[1]])
s3.3.2Y[[1]]$DIC




pl=1 #Just so I don't get knocked all the way down each time I run code above this















#-----------------------------Other code that may be useful----------------------------#
jags_params <- c("Nfb", "psi", "Nm", "survival", "lambda") #Specify parameters


head(results)

#-----------------------------Trace plots/convergence---------------------------------------------------
#Specify parameters to plot
jags_params_4plot <- c("Nfb", "psi", "Nm", "surv", "lam") #Specify parameters

#Specify save location for pdf of plots
tracePlot.file <- paste0(mcmc_plots_location, "TracePlots_", date.of.simulation, "_", seeds, "_", purpose, "_",  "_allSampleSizes.pdf")
pdf(file = tracePlot.file)

#Loop through each list element (i.e. each iteration aka mcmc object) and save to pdf.
#In the pdf, the page will correspond to the iteration
for(i in 1:length(s.all)){
  diag_plots(s.all[[i]], jags_params_4plot, layout = "4x1")
}

dev.off()

#Check convergence - should be 0 (1.01 is strict; 1.05 would work as a more liberal threshold)
results %>% summarize(no.convergence = sum(Rhat >= 1.01))
results %>% dplyr::filter(Rhat >= 1.01) %>% View()

#-----------------------------Gelman & Rubin---------------------------------------------------
# The Gelman diagnostic calculates the potential scale reduction factor (PSRF) for each variable. The PSRF estimates a factor by which the scale of the distribution might be reduced if the simulations were run for an infinite number of iterations. As the number of iterations approaches infinity, the PSRF should decline to 1.  ... it's kind of like an ANOVA, where it compares the within-chain and between-chain variance. This is essentially the same as the rhat metric.

#Calculate gelman diagnostic for each iteration
gelman.df = gelman.temp <- NULL

for(g in 1:length(s.all)){
  
  gelman.temp <- data.frame(t(gelman.diag(s.all[[g]])[[1]])) %>% 
    rownames_to_column(var = "type") %>% 
    mutate(iteration = g)
  
  gelman.df <- rbind(gelman.df, gelman.temp)
}

#Save plots of gelman diagnostic
#Specify save location for pdf of plots
gelman.file <- paste0(mcmc_plots_location, "Gelman_", date.of.simulation, "_", seeds, "_", purpose, "_", "allSampleSizes.pdf")
pdf(file = gelman.file) #Open pdf file for plotting

#Loop through each list element (i.e. each iteration aka mcmc object) and save to pdf.
#In the pdf, the page will correspond to the iteration
for(i in 1:length(s.all)){
  gelman.plot(s.all[[i]])
}

dev.off() #Close pdf file


#-----------------------------Autocorrelation---------------------------------------------------
#Specify save location for pdf of plots
autocorr.file <- paste0(mcmc_plots_location, "Autocorrelation.plots_", date.of.simulation, "_", seeds, "_", purpose, "_", "allSampleSizes.pdf")
pdf(file = autocorr.file)

#Loop through each list element (i.e. each iteration aka mcmc object) and save to pdf.
#In the pdf, the page will correspond to the iteration
for(i in 1:length(s.all)){
  autocorr.plot(s.all[[i]])
}

dev.off()

#Visual check for one iteration
autocorr.diag(s3[[5]], lags = c(0, 1, 5, 10, 15, 20))

#-----------------------------Cross-correlation---------------------------------------------------
#Specify save location for pdf of plots
crosscorr.file <- paste0(mcmc_plots_location, "Cross_correlation.plots_", date.of.simulation, "_", seeds, "_", purpose, "_", "allSampleSizes.pdf")
pdf(file = crosscorr.file)

#Loop through each list element (i.e. each iteration aka mcmc object) and save to pdf.
#In the pdf, the page will correspond to the iteration
for(c in 1:length(s.all)){
  crosscorr.plot(s.all[[c]])
}

dev.off()

#Visual check for one iteration
crosscorr(s4[[10]])

#-----------------------------Effective chain length---------------------------------------------------
#Initialize dataframes
effectiveSize.df = EF.temp <- NULL
for(j in 1:length(s.all)){
  EF.temp <- data.frame(t(effectiveSize(s.all[[j]]))) %>% 
    mutate(iteration = j)
  effectiveSize.df <- rbind(effectiveSize.df, EF.temp)
}

#This will be double bc I'm running two chains
effectiveSize.df

#-----------------------------Geweke---------------------------------------------------
# Tests the means of the first and last parts of the chain and compares their z scores. We expect that the means will be the same; otherwise, our model may not be converging until a later sample is drawn.
# the test statistic is the standard z score for the equality of two means - value greater than 1.64 (p<0.1) or 1.96 (p<0.05) two-tailed
# We do NOT want to reject the null here, as that means the beginning and end of the chains are very different

geweke.df = geweke.temp.c1 = geweke.temp.c2 <- NULL
g.thresh <- 1.65 # Set threshold for geweke failure

# Loop through all chains and iterations and store geweke diagnostic results
for(w in 1:length(s3)){
  geweke.temp.c1 <- data.frame(t(geweke.diag(s3[[w]])[[1]][[1]])) %>% 
    mutate(iteration = w, chain = 1)
  geweke.temp.c2 <- data.frame(t(geweke.diag(s3[[w]])[[2]][[1]])) %>% 
    mutate(iteration = w, chain = 2)
  geweke.df <- rbind(geweke.df, geweke.temp.c1, geweke.temp.c2)
}

# Identify index for each instance where the geweke diagnostic failed
Nf.Gew.vec <- which(abs(geweke.df$Nf) > g.thresh)
Nm.Gew.vec <- which(abs(geweke.df$Nm) > g.thresh)
surv.Gew.vec <- which(abs(geweke.df$surv) > g.thresh)
total.Gew.vec <- unique(c(Nf.Gew.vec, Nm.Gew.vec, surv.Gew.vec))

#Save iteration number of each failure so can analyze results
geweke.failed = G.temp <- NULL
for(f in 1:length(total.Gew.vec)){
  g.index <- total.Gew.vec[f]
  G.temp <- geweke.df$iteration[g.index]
  geweke.failed <- c(geweke.failed, G.temp)
}
geweke.failed

#Summary
#N Females
geweke.df %>% summarize(`0.10` = sum(abs(Nf) >1.64))
geweke.df %>% summarize(`0.05` = sum(abs(Nf) >1.96))

#N Males
geweke.df %>% summarize(`0.10` = sum(abs(Nm) >1.64))
geweke.df %>% summarize(`0.05` = sum(abs(Nm) > 1.96))

#Survival
geweke.df %>% summarize(`0.10` = sum(abs(surv) > 1.64))
geweke.df %>% summarize(`0.05` = sum(abs(surv) > 1.96))

#Save plots of geweke diagnostic
#Specify save location for pdf of plots
geweke.file <- paste0(mcmc_plots_location, "Geweke_", date.of.simulation, "_", sim.samples.3, "_", seeds, ".pdf")
pdf(file = geweke.file) #Open pdf file for plotting

#Loop through each list element (i.e. each iteration aka mcmc object) and save to pdf.
#In the pdf, the page will correspond to the iteration
for(i in 1:length(s3)){
  geweke.plot(s3[[i]])
}

dev.off() #Close pdf file

#-----------------------------Heidelberger & Welch---------------------------------------------------
#Heidelberger and Welch’s convergence diagnostic
# The HW diagnostic examines whether the draws from the posterior come from a stationary distribution. 
# It iteratively removes proportions of the samples, and reports the iteration at which we should start the chain i.e. increase the burn-in period by the iteration reported here.
# We do NOT want to reject the null.

heidel.diag(s3[[1]])


#-----------------------------Raftery-Lewis---------------------------------------------------
# Raftery-Lewis test
# Will determine the appropriate additional burn-in and thinning rate
# M = number of iterations that should be discarded
# N = total number of iterations that should be run for each variable
# Nmin = the minimum number of iterations that should be run for each variable
# I = the increase in number of iterations needed to reach convergence.
raftery.diag(s3[[1]])


superdiag(s3[[1]])
##===============================================================================


#----------------Subsetting a long chain to dial in MCMC parameters--------------------
head(s1)
length(s1)

mcmc.end <- nrow(s1[[1]][[1]])
burn_in <- 20000
thin <- 15

#Subset for proposed burn in and thinning rate
#Sim samples 1
s1.subset <- NULL
for(l in 1:length(s1)){
    s1.subset[[l]] <- window(s1[[l]], 
                                             start = burn_in +1,
                                             end = mcmc.end, 
                                             thin = thin)
  }


#Sim samples 2
s2.subset <- NULL
for(l in 1:length(s2)){
  s2.subset[[l]] <- window(s2[[l]], 
                                           start = burn_in +1,
                                           end = mcmc.end, 
                                           thin = thin)
}

#Sim samples 3
s3.subset <- NULL
for(l in 1:length(s3)){
  s3.subset[[l]] <- window(s3[[l]], 
                                           start = burn_in +1,
                                           end = mcmc.end, 
                                           thin = thin)
}

s1 <- s1.subset
s2 <- s2.subset
s3 <- s3.subset




################################################################################
#################################### DIC #######################################
################################################################################
#Doesn't seem like DIC can tell the models apart, at least not when estimating for yr 90. Can't really tell for any of the sampling scenarios. 
#If we want to try to do some model comparison, then I may need to re-run everything (bc I only saved the models for yr 90) and potentially try WAIC. Need to be more organized about this next time.
#If wanting to reevaluate all of this to draw conclusions, skip ahead and focus on/start with the extreme scenario. That's the better example and the code is more streamlined.

#Scenario slight population change
(scenario_2.1_post.all <- list.files(path = objective_2_MCMC_location,
                                     recursive = FALSE,
                                     pattern = "*scenario_2.1_*",
                                     full.names = TRUE))

(scenario_2.2.1_post.all <- list.files(path = objective_2_MCMC_location,
                                       recursive = FALSE,
                                       pattern = "*scenario_2.2.1_*",
                                       full.names = TRUE))

(scenario_2.3.1_post.all <- list.files(path = objective_2_MCMC_location,
                                       recursive = FALSE,
                                       pattern = "*scenario_2.3.1_*",
                                       full.names = TRUE))

#Slight population change
s2.1.1Y <- readRDS(scenario_2.1_post.all[5])
s2.1.1J_1 <- readRDS(scenario_2.1_post.all[6])
s2.1.1J_2 <- readRDS(scenario_2.1_post.all[7])
s2.1.1J_2.subset <- s2.1.1J_2[501:1000] #First 500 list elements are blank
s2.1.1J <- append(s2.1.1J_1, s2.1.1J_2.subset)
s2.1.1A <- readRDS(scenario_2.1_post.all[4])

s2.2.1Y <- readRDS(scenario_2.2.1_post.all[3])
s2.2.1J_1 <- readRDS(scenario_2.2.1_post.all[5])
s2.2.1J_2 <- readRDS(scenario_2.2.1_post.all[6])
s2.2.1J_2.subset <- s2.2.1J_2[501:1000] #First 500 list elements are blank
s2.2.1J <- append(s2.2.1J_1, s2.2.1J_2.subset)

s2.3.1Y <- readRDS(scenario_2.3.1_post.all[1])
s2.3.1J_1 <- readRDS(scenario_2.3.1_post.all[9])
s2.3.1J_2 <- readRDS(scenario_2.3.1_post.all[10])
s2.3.1J_2.subset <- s2.3.1J_2[501:1000] #First 500 list elements are blank
s2.3.1J <- append(s2.3.1J_1, s2.3.1J_2.subset)
s2.3.1A <- readRDS(scenario_2.3.1_post.all[4])

#Specify file to extract DIC from
file <- s2.3.1Y

DIC.temp <- NULL
DIC.df <- NULL

for(i in 1:length(file)){
  DIC.temp <- tibble(DIC_2.3.1 = file[[i]]$DIC, iteration = i, 
                     sampling.scheme = "target YOY")
  
  DIC.df <- bind_rows(DIC.df, DIC.temp)
}

#DIC.YOY_2.1.1 <- DIC.df
#DIC.YOY_2.2.1 <- DIC.df
#DIC.YOY_2.3.1 <- DIC.df

#DIC.juvs_2.1.1 <- DIC.df
#DIC.juvs_2.2.1 <- DIC.df
#DIC.juvs_2.3.1 <- DIC.df

#DIC.adults_2.1.1 <- DIC.df
#DIC.adults_2.3.1 <- DIC.df

DIC.YOY_slight <- DIC.YOY_2.1.1 %>% left_join(DIC.YOY_2.2.1, by = c("iteration", "sampling.scheme")) %>% 
  left_join(DIC.YOY_2.3.1, by = c("iteration", "sampling.scheme"))

DIC.juvs_slight <- DIC.juvs_2.1.1 %>% left_join(DIC.juvs_2.2.1, by = c("iteration", "sampling.scheme")) %>% 
  left_join(DIC.juvs_2.3.1, by = c("iteration", "sampling.scheme"))

DIC.adults_slight <- DIC.adults_2.1.1 %>% 
  left_join(DIC.adults_2.3.1, by = c("iteration", "sampling.scheme"))



#Extreme population change
(scenario_2.2.2_post.all <- list.files(path = objective_2_MCMC_location,
                                       recursive = FALSE,
                                       pattern = "*scenario_2.2.2_*",
                                       full.names = TRUE))

(scenario_2.3.2_post.all <- list.files(path = objective_2_MCMC_location,
                                       recursive = FALSE,
                                       pattern = "*scenario_2.3.2_*",
                                       full.names = TRUE))

#Extreme population change
s2.1.2Y <- readRDS(scenario_2.1_post.all[3])
s2.1.2J <- readRDS(scenario_2.1_post.all[2])
s2.1.2A <- readRDS(scenario_2.1_post.all[1])

s2.2.2Y <- readRDS(scenario_2.2.2_post.all[3])
s2.2.2J <- readRDS(scenario_2.2.2_post.all[2])
s2.2.2A <- readRDS(scenario_2.2.2_post.all[1])

s2.3.2Y <- readRDS(scenario_2.3.2_post.all[3])
s2.3.2J <- readRDS(scenario_2.3.2_post.all[2])
s2.3.2A <- readRDS(scenario_2.3.2_post.all[1])

DIC.df.temp = DIC.temp.1 = DIC.temp.2 = DIC.temp.3 <- NULL
DIC.YOY_extreme = DIC.juvs_extreme = DIC.adults_extreme <- NULL

for(i in 1:500){
  #YOY
  DIC.temp.1 <- tibble(DIC_2.1.2 = s2.1.2Y[[i]]$DIC, iteration = i, 
                       sampling.scheme = "target YOY")
  
  DIC.temp.2 <- tibble(DIC_2.2.2 = s2.2.2Y[[i]]$DIC, iteration = i, 
                       sampling.scheme = "target YOY")
  
  DIC.temp.3 <- tibble(DIC_2.3.2 = s2.3.2Y[[i]]$DIC, iteration = i, 
                       sampling.scheme = "target YOY")
  
  DIC.df.temp <- DIC.temp.1 %>% left_join(DIC.temp.2, by = c("iteration", "sampling.scheme")) %>% 
    left_join(DIC.temp.3, by = c("iteration", "sampling.scheme"))
  
  DIC.YOY_extreme <- bind_rows(DIC.YOY_extreme, DIC.df.temp)
  
  #juvs
  DIC.temp.1 <- tibble(DIC_2.1.2 = s2.1.2J[[i]]$DIC, iteration = i, 
                       sampling.scheme = "sample all juvenile ages")
  
  DIC.temp.2 <- tibble(DIC_2.2.2 = s2.2.2J[[i]]$DIC, iteration = i, 
                       sampling.scheme = "sample all juvenile ages")
  
  DIC.temp.3 <- tibble(DIC_2.3.2 = s2.3.2J[[i]]$DIC, iteration = i, 
                       sampling.scheme = "sample all juvenile ages")
  
  DIC.df.temp <- DIC.temp.1 %>% left_join(DIC.temp.2, by = c("iteration", "sampling.scheme")) %>% 
    left_join(DIC.temp.3, by = c("iteration", "sampling.scheme"))
  
  DIC.juvs_extreme <- bind_rows(DIC.juvs_extreme, DIC.df.temp)
  
  
  #adults
  #juvs
  DIC.temp.1 <- tibble(DIC_2.1.2 = s2.1.2A[[i]]$DIC, iteration = i, 
                       sampling.scheme = "sample all juvenile ages")
  
  DIC.temp.2 <- tibble(DIC_2.2.2 = s2.2.2A[[i]]$DIC, iteration = i, 
                       sampling.scheme = "sample all juvenile ages")
  
  DIC.temp.3 <- tibble(DIC_2.3.2 = s2.3.2A[[i]]$DIC, iteration = i, 
                       sampling.scheme = "sample all juvenile ages")
  
  DIC.df.temp <- DIC.temp.1 %>% left_join(DIC.temp.2, by = c("iteration", "sampling.scheme")) %>% 
    left_join(DIC.temp.3, by = c("iteration", "sampling.scheme"))
  
  DIC.adults_extreme <- bind_rows(DIC.adults_extreme, DIC.df.temp)
  
}
