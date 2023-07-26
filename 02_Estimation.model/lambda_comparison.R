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

######################### Lambda comparison #################################
iter1 <- pop_size.df %>% dplyr::filter(iteration == 1)
  
iter1$adult.lambda

#---------------Our method ----------------
all.adults <- NULL
all.adults[1] <- iter1$Total.adult.pop[1]

for(i in 2:nrow(iter1)){
  
  all.adults[i] <- round(iter1$Total.adult.pop[i-1] * iter1$adult.lambda[i], 0)
  
}



########################## MCMC & model parameters #########################
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains
