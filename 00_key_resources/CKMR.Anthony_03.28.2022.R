####-------------ANTHONY SEVEQUE, MISSISSIPPI STATE UNIVERSITY, 2022--------------####
####-------------FORWARD IN TIME POPULATION SIMULATION AND KINSHIP ANALYSIS FOR CKMR USE IN LARGE TERRESTRIAL MAMMALS--------------####

####---------- Load packages and source functions ----------####

library(dplyr)
library(tidyr)
library(nimble)
library(MCMCvis)

source("population_simulation.R")
source("pairwise_comparison.R")

####---------- A1 Promiscuity / Random sampling ----------####

###----- Input simulation parameters -----###

# Initial population
init.prop.female <- 0.5 # Proportion of the initial population size that is female
repro.age <- 5 # Set age of reproductive maturity

# Note: I could start my populations at stable age structure, but elected not to do so to let each reproductive
# strategy reach its equilibrium "naturally"... this might be dumb and I might have to change that later, dunnow, we'll see.
# (because maybe) pop won't be in equilibrium with same survivals in different harvest strategies, for instance.

# Year 0 breeding
birth.sex.ratio <- c(0.5,0.5) # The probability that each baby is F:M

# All breeding years
num.years <- 50 # The number of years to run in the simulation
yoy.survival <- 0.91 # Young of year survival / Based on black bears den checks in Michigan's Upper Peninsula
juvenile.survival <- 0.705 # Juvenile survival - between 1 years old and repro.age / Based on harvest mortality in Michigan's Upper Peninsula (mean 0.82 survival) and decreased until stable pop
adult.survival <- 0.705 # Adult survival - after age of repro.age/ Based on harvest mortality in Michigan's Upper Peninsula (mean 0.82 survival) and decreased until stable pop
max.age <- 25 # Set the maximum age allowed in the simulation

# Sampling variables
harvest.size = 10 #Harvest n% of the previous year's population
sample_first_y = 45 # First year from simulation where we start genetic sampling
sample_subset <- seq(from=sample_first_y, to=num.years, by=5) # Years of simulation where we will keep samples for CKMR. "by=n" is the sampling frequency.

# Do science good
iterations <- 4 # Replication for each scenario
rseeds <- sample(1:100000, iterations) # One different seed for each iteration



###----- Run population simulation -----###

a1_ckmr_results <- data.frame()

for(iter in 1:iterations) {
  set.seed(rseeds[iter])
  
  init.pop.size <- rnorm(1, mean = 2000, sd = 500)
  
# Run the population simulation
  a1 <- a1.promiscuity.random (init.pop.size = init.pop.size,
                               init.prop.female = init.prop.female,
                               repro.age = repro.age, 
                               birth.sex.ratio = birth.sex.ratio, 
                               num.years = num.years,
                               yoy.survival = yoy.survival,
                               juvenile.survival = juvenile.survival,
                               adult.survival = adult.survival,
                               max.age = max.age, 
                               harvest.size = harvest.size)

#Outputs a list, where the first element is a list of various population parameters for each year, and the second is the population size for each year
a1.loopy.list <- a1 [[1]] #List of dataframes for each year of simulation
a1.sample.list <- a1 [[2]] #List of sampled invidivuals for each year of simulation
a1.pop.size <- a1 [[3]] #population parameters for each year of simulation


# Verify that your population is stable (including age distribution)  
par(mfrow=c(2,2))
plot(a1.pop.size$population_size, pch=19, ylim = c(0.9*min(a1.pop.size$population_size), 1.1*max(a1.pop.size$population_size))) #Plot population size through time
plot(a1.pop.size$Nm_preharvest, pch=19, ylim = c(0.9*min(a1.pop.size$Nm_preharvest), 1.1*max(a1.pop.size$Nm_preharvest))) #Plot Nm through time
plot(a1.pop.size$Nf_preharvest, pch=19, ylim = c(0.9*min(a1.pop.size$Nf_preharvest), 1.1*max(a1.pop.size$Nf_preharvest))) #Plot Nm through time

###----- Create pairwise matrices -----###

# Select your samples that will be used to run the CKMR.
# You can play with this to test different sampling frequency.

samples <- data.frame()

for(s in sample_subset){
  
    samples_year <- a1.sample.list[[s]]
    samples_year$sampling_year <- rep(s)
    
    samples <- rbind(samples, samples_year)
  }

# Sort dataframe by birth year, so that individual 1 is always older than individual 2 (because it comes first in the df)
samples <- samples[order(samples$birth.year),]

# Run build.pairwise function, that will create HSP/POP pairwise matrices.
pairwise <- build.pairwise(samples = samples)

# Separate list elements into different dataframes
pairwise.df_all <- pairwise[[1]]
pop_positives <- pairwise[[2]]
pop_mom_comps <- pairwise[[4]]
pop_dad_comps <- pairwise[[5]]
hsp_positives <- pairwise[[6]]
hsp_mom_comps <- pairwise[[8]]
hsp_dad_comps <- pairwise[[9]]

# For curiosity's sake: check how many kin pairs
MPOP <- sum(pop_mom_comps$yes)
FPOP <- sum(pop_dad_comps$yes)
MHSP <- sum(hsp_mom_comps$yes)
FHSP <- sum(hsp_dad_comps$yes)

###----- Build CKMR model -----###

#-Simple because adult fecundity is constant with sex and age-#
#-Population abundance is also considered relatively stable, so no need to include lambda-#

###---We write the model---###

POP_model <- nimbleCode({
  
  # priors
  Nf ~ dunif(0, 1e5) # Uninformative prior for female abundance
  Nm ~ dunif(0, 1e5) # Uninformative prior for male abundance
  Phi ~ dunif(0, 1) # Uninformative prior for adult survival
  

  # likelihoods
  for (i in 1:pop_mom_length){
  MPOP[i] ~ dbinom((1/Nf), pop_mom_all_comps[i]) # Sex specific POP CKMR model equation.
  }
  
  for (j in 1:pop_dad_length){
  FPOP[j] ~ dbinom((1/Nm), pop_dad_all_comps[j]) # Sex specific POP CKMR model equation
  }
  
  for (k in 1:hsp_mom_length){
  MHSP[k] ~ dbinom((Phi^(hsp_mom_ys_birth[k] - hsp_mom_os_birth[k]))/(Nf), hsp_mom_all_comps[k]) # Sex specific HSP CKMR model equation
  }
   
  for (l in 1:hsp_dad_length){
  FHSP[l] ~ dbinom((Phi^(hsp_dad_ys_birth[l] - hsp_dad_os_birth[l]))/(Nm), hsp_dad_all_comps[l]) # Sex specific HSP CKMR model equation
  }
   
  # How does Nimble "understand" that Nf/Nm for HSP is Ind 2 birth, and not Ind 1?
  
  })

###---We put the constants in a list---###

my.constants <- list(pop_mom_length = nrow(pop_mom_comps), # Number of cohort comparisons to loop over
                     pop_dad_length = nrow(pop_dad_comps),
                     hsp_mom_length = nrow(hsp_mom_comps),
                     hsp_dad_length = nrow(hsp_dad_comps)
                     )
                     
my.constants

###---We put the data in a list---###

my.data <- list(# Moms POP
                MPOP = pop_mom_comps$yes, # Positive maternal parent offspring
                pop_mom_all_comps = pop_mom_comps$all, # Total maternal pop comparisons
                
                # Dads POP
                FPOP = pop_dad_comps$yes,
                pop_dad_all_comps = pop_dad_comps$all,
                
                # Moms HSP
                MHSP = hsp_mom_comps$yes,
                hsp_mom_all_comps = hsp_mom_comps$all,
                hsp_mom_ys_birth = hsp_mom_comps$Ind_2_birth, # Younger sibling year of birth
                hsp_mom_os_birth = hsp_mom_comps$Ind_1_birth, # Older sibling year of birth
                
                # Dads HSP
                FHSP = hsp_dad_comps$yes,
                hsp_dad_all_comps = hsp_dad_comps$all,
                hsp_dad_ys_birth = hsp_dad_comps$Ind_2_birth,
                hsp_dad_os_birth = hsp_dad_comps$Ind_1_birth
                ) 
my.data

###---We write a function that generates initial values---###

#The idea here is to provide different initial values for each chain.

initial.values <- function() list(Nf = rnorm(1, mean = 2000, sd = 100),
                                  Nm = rnorm(1, mean = 2000, sd = 100), 
                                  Phi = runif(1,0,1)
                                  )
initial.values()

###---We specify the variables we would like to monitor---###

parameters.to.save <- c("Nf", "Nm")
parameters.to.save

###---We specify a few details about the chains---###

n.iter <- 50000
n.burnin <- 10000
n.chains <- 2
n.thin <- 10

###---And now, we run nimble---###

mcmc.output <- nimbleMCMC(code = POP_model,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          thin = n.thin,
                          nburnin = n.burnin,
                          nchains = n.chains)

###---Let's look at the results---###

MCMCsummary <- MCMCsummary(mcmc.output, round = 0)

MCMCsummary_Nf <- MCMCsummary[1,]
MCMCsummary_Nm <- MCMCsummary[2,]

a1.pop.size_filt <- a1.pop.size [(sample_first_y-(max.age - repro.age)):num.years,] #Retain values of interest (based on first sampling year & number of years between repro.age and max.age)

iter_results <- cbind.data.frame(rseeds[iter], MPOP, FPOP, MHSP, FHSP,
                                 
                                 Nf_real = mean(a1.pop.size_filt$Nf_preharvest),
                                 Nf_CKMR_mean = MCMCsummary_Nf$mean, Nf_sd = MCMCsummary_Nf$sd,
                                 `Nf_2.5%` = MCMCsummary_Nf$`2.5%`,`Nf_50%` = MCMCsummary_Nf$`50%`, `Nf_97.5%` = MCMCsummary_Nf$`97.5%`,
                                 Nf_Rhat = MCMCsummary_Nf$Rhat, Nf_n.eff = MCMCsummary_Nf$n.eff,
                                 
                                 Nm_real = mean(a1.pop.size_filt$Nm_preharvest),
                                 Nm_CKMR_mean = MCMCsummary_Nm$mean, Nm_sd = MCMCsummary_Nm$sd,
                                 `Nm_2.5%` = MCMCsummary_Nm$`2.5%`,`Nm_50%` = MCMCsummary_Nm$`50%`, `Nm_97.5%` = MCMCsummary_Nm$`97.5%`,
                                 Nm_Rhat = MCMCsummary_Nm$Rhat, Nm_n.eff = MCMCsummary_Nm$n.eff) 

a1_ckmr_results <- rbind (a1_ckmr_results, iter_results)

gc() # Releases memory from unused objects between loops, or R will come to an error

# Anthony: gc() may still not be enough for 100 iterations - memory keeps stacking up.
# Error in substitute(x) : invalid environment specified
# Random errors after 7 iterations again, seriously
# Error in startsWith(tmpname, "I(") : object 'prefix' not found
# After 5 iterations


MCMCtrace(mcmc.output, pdf = FALSE)

}

write.csv(a1_ckmr_results, "a1_ckmr_results_hsponly.csv", row.names = TRUE)

# Calculate coefficient of variation

# CV = (sd / mean) * 100

# But do you calculate it for each iteration and then average CV over a scenario, or do you calculate one CV for the entire scenario? (meaning you'll lose the info on sd for each model)
# Also, is it better to use sd or 95% CI? How do you average a distribution?

# And/Or calculate bias from "real" value of Nf/Nm. 
# Need to find a way to combine precision and accuracy based on results itself and on bias 

####---------- A2 Promiscuity / Sex bias ----------####

# Problem: only harvesting dead individuals means that a sex-biased harvest does not affect adult survival and sex-ratios. 
# It will affect the sampling, but not the population per se. So not really realistic... 

# Either: separate natural survival and harvest
# Or: different survival (natural and harvest combined) for males and females

####---------- B1 Serial monogamy / Random sampling ----------####

####---------- B2 Serial monogamy / Sex bias ----------####

####---------- C1 Polygyny / Random sampling ----------####

####---------- C2 Polygyny / Sex bias ----------####

####---------- Michigan black bear ----------####
