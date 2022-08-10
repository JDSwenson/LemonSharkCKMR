## Diagnostics
library(coda)
library(tidyverse)
library(R2jags)
library(postpack)

rm(list=ls())

#############################################################################
###############Load posterior samples - Objective 3##########################
#############################################################################
#Objective 3.1
(scenario_3.1.1_post.all <- list.files(path = objective_3_MCMC_location,
                                       recursive = FALSE,
                                       pattern = "*1.5prop.*scenario_3.1.1_*",
                                       full.names = TRUE))

s3.1.1A <- readRDS(scenario_3.1.1_post.all[1])
s3.1.1J <- readRDS(scenario_3.1.1_post.all[2])
s3.1.1Y <- readRDS(scenario_3.1.1_post.all[3])

#-----------------Cross correlation among parameters-------------------#
cross.temp <- NULL
s3.1.1_cross.list.A = s3.1.1_cross.list.J = s3.1.1_cross.list.Y <- list()

for(cr in 1:length(s3.1.1A)){
  s3.1.1_cross.list.A[[cr]] <- crosscorr(s3.1.1A[[cr]])
  
  s3.1.1_cross.list.J[[cr]] <- crosscorr(s3.1.1J[[cr]])
  
  s3.1.1_cross.list.Y[[cr]] <- crosscorr(s3.1.1Y[[cr]])
}

##Mean cross-correlation values
(s3.1.1_cross.mean.A <- Reduce("+", s3.1.1_cross.list.A) / length(s3.1.1_cross.list.A))

(s3.1.1_cross.mean.J <- Reduce("+", s3.1.1_cross.list.J) / length(s3.1.1_cross.list.J))

(s3.1.1_cross.mean.Y <- Reduce("+", s3.1.1_cross.list.Y) / length(s3.1.1_cross.list.Y))


#Mostly interested in comparing scenarios 3.2.2 and 3.3.2 (both have lambda; the first has 5% non-conformists and the second has 0 non-conformists)
#----------------------------Objective 3.2.2------------------------------#
(scenario_3.2.2_post.all <- list.files(path = objective_3_MCMC_location,
                                       recursive = FALSE,
                                       pattern = "*1.5prop.*scenario_3.2.2_*",
                                       full.names = TRUE))

s3.2.2A <- readRDS(scenario_3.2.2_post.all[1])
s3.2.2J <- readRDS(scenario_3.2.2_post.all[2])
s3.2.2Y <- readRDS(scenario_3.2.2_post.all[3])

#-----------------Cross correlation among parameters-------------------#
cross.temp <- NULL
s3.2.2_cross.list.A = s3.2.2_cross.list.J = s3.2.2_cross.list.Y <- list()

for(cr in 1:length(s3.2.2A)){
  s3.2.2_cross.list.A[[cr]] <- crosscorr(s3.2.2A[[cr]])
  
  s3.2.2_cross.list.J[[cr]] <- crosscorr(s3.2.2J[[cr]])
  
  s3.2.2_cross.list.Y[[cr]] <- crosscorr(s3.2.2Y[[cr]])
}

##Mean cross-correlation values
(s3.2.2_cross.mean.A <- Reduce("+", s3.2.2_cross.list.A) / length(s3.2.2_cross.list.A))

(s3.2.2_cross.mean.J <- Reduce("+", s3.2.2_cross.list.J) / length(s3.2.2_cross.list.J))

(s3.2.2_cross.mean.Y <- Reduce("+", s3.2.2_cross.list.Y) / length(s3.2.2_cross.list.Y))




#-----------------------Objective 3.3.2-----------------------------#
(scenario_3.3.2_post.all <- list.files(path = objective_3_MCMC_location,
                                       recursive = FALSE,
                                       pattern = "*1.5prop.*scenario_3.3.2_*",
                                       full.names = TRUE))

s3.3.2A <- readRDS(scenario_3.3.2_post.all[1])
s3.3.2J <- readRDS(scenario_3.3.2_post.all[2])
s3.3.2Y <- readRDS(scenario_3.3.2_post.all[3])


#-----------------Cross correlation among parameters-------------------#
cross.temp <- NULL
s3.3.2_cross.list.A = s3.3.2_cross.list.J = s3.3.2_cross.list.Y <- list()

for(cr in 1:length(s3.3.2A)){
  s3.3.2_cross.list.A[[cr]] <- crosscorr(s3.3.2A[[cr]])
  
  s3.3.2_cross.list.J[[cr]] <- crosscorr(s3.3.2J[[cr]])
  
  s3.3.2_cross.list.Y[[cr]] <- crosscorr(s3.3.2Y[[cr]])
}

##Mean cross-correlation values
(s3.3.2_cross.mean.A <- Reduce("+", s3.3.2_cross.list.A) / length(s3.3.2_cross.list.A))

(s3.3.2_cross.mean.J <- Reduce("+", s3.3.2_cross.list.J) / length(s3.3.2_cross.list.J))

(s3.3.2_cross.mean.Y <- Reduce("+", s3.3.2_cross.list.Y) / length(s3.3.2_cross.list.Y))



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

