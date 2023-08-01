library(tidyverse) # safe to ignore conflicts with filter() and lag()
library(MASS)
library(popbio)
library(mpmtools)
library(ggpubr)
library(rjags)
library(R2jags)
library(jagsUI)
library(Rlab)
library(runjags)
library(postpack)
library(coda)
library(FSA)

######################### Specify common output prefixes #########################
results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
jags.model.prefix <- "CKMR.JAGS_"
mom.comps.prefix <- "comparisons/mom.comps"
dad.comps.prefix <- "comparisons/dad.comps"
samples.prefix <- "samples/samples.missassigned"

######################### Specify other common settings #########################
PopSim.location <- "G://My Drive/Personal_Drive/R/CKMR/Population.simulations/Nov2022/" #Location of population simulation output files
temp_location <- "~/R/working_directory/temp_results/" #Location to save temporary files in case run gets killed
jags.model_location <- "G://My Drive/Personal_Drive/R/CKMR/JAGS_models/" #Location of JAGS models
today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
date.of.simulation <- today
rseeds <- readRDS("rseeds_2022.04.15.rda")
outSeeds <- "Seeds2022.04.15"


adult.survival <- 0.825 # CHANGED FROM 0.825; Adult survival
repro.age <- 12 # set age of reproductive maturity
max.age <- 50 #set the maximum age allowed in the simulation
n_yrs <- 90 #Number of years the simulation was run for
inSeeds <- "Seeds2022.04.15" #Seeds used for population simulation
sample.props <- "all" #Either label this with the percent we want to target if just one (e.g., 1.5)) or if wanting to run over all sample proportions, set as "all"

#------------------------- Set input file locations -------------------------#
PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding

#------------------------- Set output file locations -------------------------# 
 MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.2_population.growth/Nov2022/Model.output/"
 results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.2_population.growth/Nov2022/Model.results/"
 outSeeds <- "Seeds2022.04.15"

#------------------------- Objective 2 common settings -------------------------#
 estimation.years <- c(n_yrs - 10, n_yrs - 5, n_yrs)
 fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
 model <- "annual.model" #For naming output files
 source("./02_Estimation.model/functions/Obj2.functions.R") #Changed name of script that includes pairwise comparison and other functions


#========================= Scenario 2.1 =========================
#jags_file = paste0(jags.model_location, "HS.PO_noLambda_annual_model.txt") #Annual JAGS model w/o lambda

#------------------------- Scenario 2.1.1: Small population decline; no lambda in model
# PopSim.lambda <- "lambda.slight.decrease"
# scenario<- "scenario_2.1.1" #For naming output files
# jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 2.1.2: Small population growth; no lambda in model
# PopSim.lambda <- "lambda.slight.increase"
# scenario<- "scenario_2.1.2" #For naming output files
# jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 2.1.3: Substantial population decline; no lambda in model
# PopSim.lambda <- "lambda.extreme" 
# scenario<- "scenario_2.1.3" #For naming output files
# jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")


#========================= Scenario 2.2 =========================
jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")

#------------------------- Scenario 2.2.1: Small population decline; lambda in model w/ tight prior
 PopSim.lambda <- "lambda.slight.decrease"
 scenario<- "scenario_2.2.1" #For naming output files
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 2.2.2: Small population growth; lambda in model w/ tight prior
# PopSim.lambda <- "lambda.slight.increase"
# scenario<- "scenario_2.2.2" #For naming output files
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 2.2.3: Substantial population decline; lambda in model w/ tight prior
# PopSim.lambda <- "lambda.extreme" 
# scenario<- "scenario_2.2.3" #For naming output files
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")


#========================= Scenario 2.3 =========================
#jags_file = paste0(jags.model_location, "HS.PO_wideLambda_annual_model.txt")

#------------------------- Scenario 2.3.1: Small population decline; lambda in model w/ wide prior
# PopSim.lambda <- "lambda.slight.decrease"
# scenario<- "scenario_2.3.1" #For naming output files
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 2.3.2: Small population growth; lambda in model w/ wide prior
# PopSim.lambda <- "lambda.slight.increase"
# scenario<- "scenario_2.3.2" #For naming output files
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 2.3.3: Substantial population decline; lambda in model w/ wide prior
# PopSim.lambda <- "lambda.extreme" 
# scenario<- "scenario_2.3.3" #For naming output files
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")

#========================== Sampling
#------------------------- Target YOY -------------------------#
# sampling.scheme <- "target.YOY" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
# date.of.PopSim <- "17Nov2022"
# HS.only <- "yes" #Do we only want to filter HS relationships?
# PO.only <- "no" #Do we only want to filter PO relationships?

#------------------------- Sample all juveniles -------------------------#
 sampling.scheme <- "sample.all.juvenile.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
 date.of.PopSim <- "17Nov2022"
 HS.only <- "yes" #Do we only want to filter HS relationships?
 PO.only <- "no" #Do we only want to filter PO relationships?

#------------------------- Sample all ages -------------------------#
# sampling.scheme <- "sample.ALL.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
# date.of.PopSim <- "17Nov2022"
# HS.only <- "no" #Do we only want to filter HS relationships?
# PO.only <- "no" #Do we only want to filter PO relationships?

########################## Read in sampling and other dataframes #########################
samples.df <- readRDS(file = paste0(PopSim.location, "sample.info_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))

samples.df %>% group_by(age.x) %>% summarize(n())

pop_size.df <- readRDS(file = paste0(PopSim.location, "pop.size_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))

truth.df <- readRDS(file = paste0(PopSim.location, "truth_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))

n_yrs <- max(pop_size.df$year)

#Run the below twice: once after loading files from population increase, and once after loading files from population decrease

######################### Lambda comparison #################################
#--------------------------Compare over 90 years-------------------------
yr0 <- 1
yrs <- 90

#Initialize df for saving values
N90_estimates <- NULL

for(i in 1:max(pop_size.df$iteration)){
  #Subset for iteration of interest
  iter <- pop_size.df %>% dplyr::filter(iteration == i)

    #Change calculation of lambda
  mean.lam_yr1 <- (iter$Total.adult.pop[yrs]/iter$Total.adult.pop[1])^(1/(yrs - yr0))
  
  indv.lambdas <- iter$adult.lambda
  
  #Initialize vectors
  indv.method_yr1 = avg.method_yr1 <- NULL
  
  for(j in 2:length(indv.lambdas)){
    
    #Store N0 and save as first element of each vector
    N0 <- iter$Total.adult.pop[1]
    indv.method_yr1[1] = avg.method_yr1[1] <- N0
    
    #Save vector of individual lambda values
    indv.lambdas_yr1 <- iter$adult.lambda
    
    #Calculate mean lambda -- removed after Liz suggested a different approach
#    mean.lam_yr1 <- mean(iter$adult.lambda, na.rm = TRUE)
    
    #Calculate abundance
    indv.method_yr1[j] <- round(indv.method_yr1[j-1] * indv.lambdas_yr1[j], 0)
    avg.method_yr1[j] <- round(avg.method_yr1[j-1] * mean.lam_yr1, 0)
    
  }
  
  combined.df <- tibble(iter$Total.adult.pop[90], indv.method_yr1[90], avg.method_yr1[90]) %>% mutate(iteration = iter$iteration[j],
                                                                                                          yr0 = yr0,
                                                                                                          lambda.yrs = yrs-yr0)
  
  N90_estimates <- bind_rows(N90_estimates, combined.df)

  }
  
colnames(N90_estimates)[1:3] <- c("true.adult.pop", "indv.method.pop", "avg.method.pop")

N90_estimates



#-----------------------Compare over 40 years-------------------------
yr0 <- 50
yrs <- 90

#Initialize df for saving values
N50_estimates <- NULL

for(i in 1:max(pop_size.df$iteration)){
  #Subset for iteration of interest
  iter <- pop_size.df %>% dplyr::filter(iteration == i)
  
  #Change calculation of lambda
  mean.lam_yr50 <- (iter$Total.adult.pop[yrs]/iter$Total.adult.pop[yr0])^(1/(yrs - yr0))
  
  indv.lambdas <- iter$adult.lambda
  
  #Initialize vectors
  indv.method_yr50.vec = avg.method_yr50.vec <- NULL
  
  #Store N0 and save as first element of each vector
  N0 <- iter$Total.adult.pop[yr0]
  indv.method_yr50.vec[1] = avg.method_yr50.vec[1] <- N0
  
  #Save vector of individual lambda values
  indv.lambdas_yr50 <- iter$adult.lambda

  indx <- 1
    
  for(j in yr0:yrs){
    #Calculate mean lambda
    #mean.lam_yr50 <- mean(iter$adult.lambda, na.rm = TRUE)
    
    indx <- indx + 1
    
    #Calculate abundance
    indv.method_yr50.vec[indx] <- round(indv.method_yr50.vec[indx-1] * indv.lambdas_yr50[j+1], 0)

    avg.method_yr50.vec[indx] <- round(avg.method_yr50.vec[indx-1] * mean.lam_yr50, 0)

  }
  
  combined.df <- tibble(iter$Total.adult.pop[yrs], indv.method_yr50.vec[yrs-yr0+1], avg.method_yr50.vec[yrs-yr0+1]) %>% mutate(iteration = iter$iteration[j],
                                                                                                          yr0 = yr0,
                                                                                                          lambda.yrs = yrs-yr0)
  
  N50_estimates <- bind_rows(N50_estimates, combined.df)
  
}

colnames(N50_estimates)[1:3] <- c("true.adult.pop", "indv.method.pop", "avg.method.pop")

N50_estimates

#-----------------------Compare over 20 years-------------------------
yr0 <- 70
yrs <- 90

#Initialize df for saving values
N70_estimates <- NULL

for(i in 1:max(pop_size.df$iteration)){
  #Subset for iteration of interest
  iter <- pop_size.df %>% dplyr::filter(iteration == i)
  
  #Change calculation of lambda
  mean.lam_yr70 <- (iter$Total.adult.pop[yrs]/iter$Total.adult.pop[yr0])^(1/(yrs - yr0))
  
  indv.lambdas <- iter$adult.lambda
  
  #Initialize vectors
  indv.method_yr70.vec = avg.method_yr70.vec <- NULL
  
  #Store N0 and save as first element of each vector
  N0 <- iter$Total.adult.pop[yr0]
  indv.method_yr70.vec[1] = avg.method_yr70.vec[1] <- N0
  
  #Save vector of individual lambda values
  indv.lambdas_yr70 <- iter$adult.lambda
  
  indx <- 1
  
  for(j in yr0:yrs){
    #Calculate mean lambda
    #mean.lam_yr70 <- mean(iter$adult.lambda, na.rm = TRUE)
    
    indx <- indx + 1
    
    #Calculate abundance
    indv.method_yr70.vec[indx] <- round(indv.method_yr70.vec[indx-1] * indv.lambdas_yr70[j+1], 0)
    
    avg.method_yr70.vec[indx] <- round(avg.method_yr70.vec[indx-1] * mean.lam_yr70, 0)
    
  }
  
  combined.df <- tibble(iter$Total.adult.pop[yrs], indv.method_yr70.vec[yrs-yr0+1], avg.method_yr70.vec[yrs-yr0+1]) %>% mutate(iteration = iter$iteration[j],
                                                                                                                               yr0 = yr0,
                                                                                                                               lambda.yrs = yrs-yr0)
  
  N70_estimates <- bind_rows(N70_estimates, combined.df)
  
}

colnames(N70_estimates)[1:3] <- c("true.adult.pop", "indv.method.pop", "avg.method.pop")

N70_estimates


#-----------------------Compare over 10 years-------------------------
yr0 <- 80
yrs <- 90

#Initialize df for saving values
N80_estimates <- NULL

for(i in 1:max(pop_size.df$iteration)){
  #Subset for iteration of interest
  iter <- pop_size.df %>% dplyr::filter(iteration == i)
  
  #Change calculation of lambda
  mean.lam_yr80 <- (iter$Total.adult.pop[yrs]/iter$Total.adult.pop[yr0])^(1/(yrs - yr0))
  
  indv.lambdas <- iter$adult.lambda
  
  #Initialize vectors
  indv.method_yr80.vec = avg.method_yr80.vec <- NULL
  
  #Store N0 and save as first element of each vector
  N0 <- iter$Total.adult.pop[yr0]
  indv.method_yr80.vec[1] = avg.method_yr80.vec[1] <- N0
  
  #Save vector of individual lambda values
  indv.lambdas_yr80 <- iter$adult.lambda
  
  indx <- 1
  
  for(j in yr0:yrs){
    #Calculate mean lambda
    #mean.lam_yr80 <- mean(iter$adult.lambda, na.rm = TRUE)
    
    indx <- indx + 1
    
    #Calculate abundance
    indv.method_yr80.vec[indx] <- round(indv.method_yr80.vec[indx-1] * indv.lambdas_yr80[j+1], 0)
    
    avg.method_yr80.vec[indx] <- round(avg.method_yr80.vec[indx-1] * mean.lam_yr80, 0)
    
  }
  
  combined.df <- tibble(iter$Total.adult.pop[yrs], indv.method_yr80.vec[yrs-yr0+1], avg.method_yr80.vec[yrs-yr0+1]) %>% mutate(iteration = iter$iteration[j],
                                                                                                                               yr0 = yr0,
                                                                                                                               lambda.yrs = yrs-yr0)
  
  N80_estimates <- bind_rows(N80_estimates, combined.df)
  
}

colnames(N80_estimates)[1:3] <- c("true.adult.pop", "indv.method.pop", "avg.method.pop")

N80_estimates


#-----------------------Compare over 5 years-------------------------
yr0 <- 85
yrs <- 90

#Initialize df for saving values
N85_estimates <- NULL

for(i in 1:max(pop_size.df$iteration)){
  #Subset for iteration of interest
  iter <- pop_size.df %>% dplyr::filter(iteration == i)
  
  #Change calculation of lambda
  mean.lam_yr85 <- (iter$Total.adult.pop[yrs]/iter$Total.adult.pop[yr0])^(1/(yrs - yr0))
  
  indv.lambdas <- iter$adult.lambda
  
  #Initialize vectors
  indv.method_yr85.vec = avg.method_yr85.vec <- NULL
  
  #Store N0 and save as first element of each vector
  N0 <- iter$Total.adult.pop[yr0]
  indv.method_yr85.vec[1] = avg.method_yr85.vec[1] <- N0
  
  #Save vector of individual lambda values
  indv.lambdas_yr85 <- iter$adult.lambda
  
  indx <- 1
  
  for(j in yr0:yrs){
    #Calculate mean lambda
    #mean.lam_yr85 <- mean(iter$adult.lambda, na.rm = TRUE)
    
    indx <- indx + 1
    
    #Calculate abundance
    indv.method_yr85.vec[indx] <- round(indv.method_yr85.vec[indx-1] * indv.lambdas_yr85[j+1], 0)
    
    avg.method_yr85.vec[indx] <- round(avg.method_yr85.vec[indx-1] * mean.lam_yr85, 0)
    
  }
  
  combined.df <- tibble(iter$Total.adult.pop[yrs], indv.method_yr85.vec[yrs-yr0+1], avg.method_yr85.vec[yrs-yr0+1]) %>% mutate(iteration = iter$iteration[j],
                                                                                                                               yr0 = yr0,
                                                                                                                               lambda.yrs = yrs-yr0)
  
  N85_estimates <- bind_rows(N85_estimates, combined.df)
  
}

colnames(N85_estimates)[1:3] <- c("true.adult.pop", "indv.method.pop", "avg.method.pop")

N85_estimates


#Combine positive and negative lambda growth
lambda.df.decrease <- bind_rows(N90_estimates, N50_estimates, N70_estimates, N80_estimates, N85_estimates)
lambda.df.increase <- bind_rows(N90_estimates, N50_estimates, N70_estimates, N80_estimates, N85_estimates)

lambda.df.decrease <- lambda.df.decrease %>% mutate(pop.growth = "slight negative")
lambda.df.increase <- lambda.df.increase %>% mutate(pop.growth = "slight positive")

lambda.df <- bind_rows(lambda.df.increase, lambda.df.decrease) 

lambda.df <- lambda.df %>% mutate(relative_bias = ((avg.method.pop - true.adult.pop)/true.adult.pop)*100)

lambda.df %>% group_by(pop.growth, lambda.yrs) %>% 
  summarize(median(relative_bias))


lambda.df %>% 
  ggplot(aes(x = factor(lambda.yrs), y = relative_bias, fill = factor(pop.growth))) +
# geom_violin(draw_quantiles = 0.5, position = position_dodge(0.75), width = 4) +
  #geom_violin(draw_quantiles = 0.5, width = 4) +
  geom_jitter(width = 0.2, aes(col = factor(pop.growth)), alpha = 0.5) +
  geom_boxplot() +
#  geom_hline(yintercept=0, col="black", size=1.25) +
  labs(x = "lambda years", y = "relative bias") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=25),
        legend.title = element_blank())

#Scatter plot of exact relative bias
lambda.df %>% 
  ggplot(aes(x = factor(pop.growth), y = relative_bias, fill = factor(pop.growth))) + 
  geom_jitter() +
  facet_wrap(~lambda.yrs)


#-------------Compare total lambda vs adult lambda----------------
yr0 <- 1
yrs <- 90

pop_size.df <- pop_size.df %>% mutate(juv.pop.size = population_size - Total.adult.pop)

#Initialize df for saving values
adult.lam_90 = total.lam_90 = juv.lam_90 <- NULL

for(i in 1:max(pop_size.df$iteration)){
  #Subset for iteration of interest
  iter <- pop_size.df %>% dplyr::filter(iteration == i)
  
  #Change calculation of lambda
  mean.adult.lam_90 <- (iter$Total.adult.pop[yrs]/iter$Total.adult.pop[1])^(1/(yrs - yr0))
  adult.lam_90 <- c(adult.lam_90, mean.adult.lam_90)

  mean.total.lam_90 <- (iter$population_size[yrs]/iter$population_size[1])^(1/(yrs - yr0))
  total.lam_90 <- c(total.lam_90, mean.total.lam_90)
  
  mean.juv.lam_90 <- (iter$juv.pop.size[yrs]/iter$juv.pop.size[1])^(1/(yrs - yr0))
  juv.lam_90 <- c(juv.lam_90, mean.juv.lam_90)
  
}

#Compare lambda for difference age classes
N90_lambda <- data.frame(total.lam_90, adult.lam_90, juv.lam_90) %>% 
  mutate(total_v_adult = total.lam_90 - adult.lam_90,
         total_v_juv = total.lam_90 - juv.lam_90,
         adult_v_juv = adult.lam_90 - juv.lam_90)

#Scatter plot
N90_lambda %>% pivot_longer(cols = c(total_v_adult, total_v_juv, adult_v_juv),
                            names_to = c("comparison"),
                            values_to = "difference") %>% 
  ggplot(aes(x = difference, y = factor(comparison), fill = factor(comparison))) +
  geom_jitter(aes(col = factor(comparison))) +
  geom_boxplot(alpha = 0.3)