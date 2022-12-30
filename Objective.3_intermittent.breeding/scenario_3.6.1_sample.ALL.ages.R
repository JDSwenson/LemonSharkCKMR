#Load packages
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

rm(list=ls())

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
sample.props <- 1.5 #Either label this with the percent we want to target (if just one) or if wanting to run over all sample proportions, set as "all"
derived.quantities <- "no"


########################## Objective 1 #########################
#------------------------- Set input file locations -------------------------#
# PopSim.lambda <- "lambda.1" # Can be lambda.1 or lambda.variable
# PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding

#------------------------- Set output file locations -------------------------# 
# MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.validation/Model.output/"
# results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.validation/Model.results/"

#------------------------- Objective 1 common settings -------------------------#
# fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
# jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")
# mating.periodicity <- 1 #number of years between mating for females
#Survival prior info
# survival.prior.mean <- adult.survival
# survival.prior.cv <- 0.05
# survival.prior.sd <- survival.prior.mean * survival.prior.cv
# scenario <- "objective_1_model.validation" #For naming output files and calculating truth.
# model <- "annual.model" #For naming output files
# 
# source("./02_Estimation.model/functions/Obj1.functions.R") #Changed name of script that includes pairwise comparison and other functions
#jags_file = paste0(jags.model_location, "HS.PO_noLambda_annual_model_validation.txt") #Specify JAGS model


#------------------------- Target YOY -------------------------#
# sampling.scheme <- "target.YOY" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
# date.of.PopSim <- "17Nov2022"
# HS.only <- "yes" #Do we only want to filter HS relationships?
# PO.only <- "no" #Do we only want to filter PO relationships?

#------------------------- Sample all juveniles -------------------------#
# sampling.scheme <- "sample.all.juvenile.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
# date.of.PopSim <- "17Nov2022"
# HS.only <- "yes" #Do we only want to filter HS relationships?
# PO.only <- "no" #Do we only want to filter PO relationships?

#------------------------- Sample all ages -------------------------#
# sampling.scheme <- "sample.ALL.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
# date.of.PopSim <- "17Nov2022"
# HS.only <- "no" #Do we only want to filter HS relationships?
# PO.only <- "no" #Do we only want to filter PO relationships?




########################## Objective 2 #########################
#------------------------- Set input file locations -------------------------#
# PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding

#------------------------- Set output file locations -------------------------# 
# MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.2_population.growth/Nov2022/Model.output/"
# results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.2_population.growth/Nov2022/Model.results/"
# outSeeds <- "Seeds2022.04.15"

#------------------------- Objective 2 common settings -------------------------#
# estimation.years <- c(n_yrs - 10, n_yrs - 5, n_yrs)
# fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
# model <- "annual.model" #For naming output files
# source("./02_Estimation.model/functions/Obj2.functions.R") #Changed name of script that includes pairwise comparison and other functions


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
#jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")

#------------------------- Scenario 2.2.1: Small population decline; lambda in model w/ tight prior
# PopSim.lambda <- "lambda.slight.decrease"
# scenario<- "scenario_2.2.1" #For naming output files
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")

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
# sampling.scheme <- "sample.all.juvenile.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
# date.of.PopSim <- "17Nov2022"
# HS.only <- "yes" #Do we only want to filter HS relationships?
# PO.only <- "no" #Do we only want to filter PO relationships?

#------------------------- Sample all ages -------------------------#
# sampling.scheme <- "sample.ALL.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
# date.of.PopSim <- "17Nov2022"
# HS.only <- "no" #Do we only want to filter HS relationships?
# PO.only <- "no" #Do we only want to filter PO relationships?



########################## Objective 3 #########################
#------------------------- Set output file locations -------------------------#
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.3_intermittent.breeding/Nov2022/Model.output/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.3_intermittent.breeding/Nov2022/Model.results/"
outSeeds <- "Seeds2022.04.15"
date.of.PopSim <- "22Nov2022"

#------------------------- Objective 3 common settings -------------------------#
estimation.years <- c(n_yrs - 5, n_yrs)
fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
PopSim.lambda <- "lambda.1"
source("./02_Estimation.model/functions/Obj3.functions.R") #Changed name of script that includes pairwise comparison and other functions
mating.periodicity <- 2 #Overwrite below if triennial breeding


#========================= Scenario 3.1 =========================
#------------------------- Scenario 3.1.1: Biennial breeding; psi = 1; annual model
# PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_3.1.1" #For naming output files
# model <- "annual.model"
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")
# jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")

#------------------------- Scenario 3.1.2: Biennial breeding; psi = 1; multiennial model
# PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_3.1.2" #For naming output files
# model <- "multiennial.model"
# jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")
# jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt") #Specify the HS only model as the model for now. If the likelihood should include POPs, then we will override the model when we select the sampling scheme several lines down.



#========================= Scenario 3.2 =========================
#------------------------- Scenario 3.2.1: Biennial breeding; psi = 0.9; annual model
# PopSim.breeding.schedule <- "biennial.breeding_psi0.90" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_3.2.1" #For naming output files
# model <- "annual.model"
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")
# jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")

#------------------------- Scenario 3.2.2: Biennial breeding; psi = 0.9; multiennial model
# PopSim.breeding.schedule <- "biennial.breeding_psi0.90" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_3.2.2" #For naming output files
# model <- "multiennial.model"
# jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")
# jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt") #Specify the HS only model as the model for now. If the likelihood should include POPs, then we will override the model when we select the sampling scheme several lines down


#========================= Scenario 3.3 =========================
#------------------------- Scenario 3.3.1: Biennial breeding; psi = 0.75; annual model
# PopSim.breeding.schedule <- "biennial.breeding_psi0.75" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_3.3.1" #For naming output files
# model <- "annual.model"
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")
#jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")

#------------------------- Scenario 3.3.2: Biennial breeding; psi = 0.75; multiennial model
# PopSim.breeding.schedule <- "biennial.breeding_psi0.75" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_3.3.2" #For naming output files
# model <- "multiennial.model"
# jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")
# jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt") #Specify the HS only model as the model for now. If the likelihood should include POPs, then we will override the model when we select the sampling scheme several lines down



#========================= Scenario 3.4 =========================
#------------------------- Scenario 3.4.1: Biennial breeding; psi = 0.50; annual model
# PopSim.breeding.schedule <- "biennial.breeding_psi0.50" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_3.4.1" #For naming output files
# model <- "annual.model"
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")
# jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")

#------------------------- Scenario 3.4.2: Biennial breeding; psi = 0.50; multiennial model
# PopSim.breeding.schedule <- "biennial.breeding_psi0.50" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_3.4.2" #For naming output files
# model <- "multiennial.model"
# jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")
# jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt") #Specify the HS only model as the model for now. If the likelihood should include POPs, then we will override the model when we select the sampling scheme several lines down


#========================= Scenario 3.5: Triennial breeding =========================
#------------------------- Scenario 3.5.1: Triennial breeding; psi = 1; annual model
# PopSim.breeding.schedule <- "triennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_3.5.1" #For naming output files
# model <- "annual.model"
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")
#jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")
#mating.periodicity <- 3

#------------------------- Scenario 3.5.2: Triennial breeding; psi = 1; multiennial model
# PopSim.breeding.schedule <- "triennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_3.5.2" #For naming output files
# model <- "multiennial.model"
# jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")
# jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt") #Specify the HS only model as the model for now. If the likelihood should include POPs, then we will override the model when we select the sampling scheme several lines down
# mating.periodicity <- 3

#========================= Scenario 3.6: Biennial breeding w/ stochasticity
#------------------------- Scenario 3.6.1: Biennial breeding; psi = 1; annual model; 10% on-cycle breeders fail to breed; 10% off-cycle breeders do breed.
date.of.PopSim <- "27Dec2022"
PopSim.breeding.schedule <- "biennial.breeding_psi1_stochasticCycles" #Can be annual.breeding or biennial.breeding
scenario <- "scenario_3.6.1" #For naming output files
model <- "annual.model"
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")
jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")

#------------------------- Scenario 3.6.2: Biennial breeding; psi = 1; multiennial model; 10% on-cycle breeders fail to breed; 10% off-cycle breeders do breed.
# date.of.PopSim <- "27Dec2022"
# PopSim.breeding.schedule <- "biennial.breeding_psi1_stochasticCycles" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_3.6.2" #For naming output files
# model <- "multiennial.model"
# jags_params = c("Nf", "psi", "Nm", "survival", "lambda", "Nfb1", "Nfb2") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")
# jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt") #Specify the HS only model as the model for now. If the likelihood should include POPs, then we will override the model when we select the sampling scheme several lines down.
# derived.quantities = "yes"

#========================== Sampling ==========================
#------------------------- Target YOY -------------------------#
#sampling.scheme <- "target.YOY" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
#HS.only <- "yes" #Do we only want to filter HS relationships?

#------------------------- Sample all juveniles -------------------------#
#sampling.scheme <- "sample.all.juvenile.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
#HS.only <- "yes" #Do we only want to filter HS relationships?

#------------------------- Sample all ages -------------------------#
sampling.scheme <- "sample.ALL.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
HS.only <- "no" #Do we only want to filter HS relationships?
PO.only <- "no" #Do we only want to filter PO relationships?
#jags_file <- paste0(jags.model_location, "HS.PO_narrowLambda_Skip_model.txt") #Need to use a different model if sampling all age classes AND using the biennial model



########################## Objective 4 #########################
#------------------------- Set output file locations ------------------------- 
# MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.4_aging.error/Nov2022/Model.output/"
# results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.4_aging.error/Nov2022/Model.results/"
# outSeeds <- "Seeds2022.04.15"
# date.of.PopSim <- "22Nov2022"

#------------------------- Objective 4 common settings -------------------------
# estimation.years <- c(n_yrs - 5, n_yrs)
# fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
# PopSim.lambda <- "lambda.1"
# source("./02_Estimation.model/functions/Obj4.functions.R") #Changed name of script that includes pairwise comparison and other functions

#------------------------- Scenario 4.1: Biennial breeding; psi = 1; aging error 5% CV
# age.cv <- 0.05
#mating.periodicity <- 2
# PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_4.1_5CV" #For naming output files
# model <- "multiennial.model"
# jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")


#------------------------- Scenario 4.2: Biennial breeding; psi = 1; aging error 10% CV
# age.cv <- 0.10
#mating.periodicity <- 2
# PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_4.2_10CV" #For naming output files
# model <- "multiennial.model"
# jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")


#------------------------- Scenario 4.3: Biennial breeding; psi = 1; aging error 20% CV
# age.cv <- 0.20
#mating.periodicity <- 2
# PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_4.3_20CV" #For naming output files
# model <- "multiennial.model"
# jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")

#========================== Sampling - biennial breeding
#------------------------- Target YOY -------------------------#
# sampling.scheme <- "target.YOY" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
# HS.only <- "yes" #Do we only want to filter HS relationships?

#------------------------- Sample all juveniles -------------------------#
# sampling.scheme <- "sample.all.juvenile.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
# HS.only <- "yes" #Do we only want to filter HS relationships?

#------------------------- Sample all ages -------------------------#
# sampling.scheme <- "sample.ALL.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
# HS.only <- "no" #Do we only want to filter HS relationships?
# PO.only <- "no" #Do we only want to filter PO relationships?



#------------------------- Scenario 4.4: Annual breeding; aging error 5% CV
# age.cv <- 0.05
#mating.periodicity <- 1
# PopSim.breeding.schedule <- "annual_breeding" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_4.4_5CV" #For naming output files
# model <- "annual.model"
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")


#------------------------- Scenario 4.5: Biennial breeding; psi = 1; aging error 10% CV
# age.cv <- 0.10
#mating.periodicity <- 1
# PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_4.5_10CV" #For naming output files
# model <- "annual.model"
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")


#------------------------- Scenario 4.6: Biennial breeding; psi = 1; aging error 20% CV
# age.cv <- 0.20
# mating.periodicity <- 1
# PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding
# scenario <- "scenario_4.6_20CV" #For naming output files
# model <- "annual.model"
# jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
# estimated.parameters <- paste0(jags_params, collapse = ",")


#========================== Sampling - annual breeding
#------------------------- Target YOY -------------------------#
# sampling.scheme <- "target.YOY" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
# date.of.PopSim <- "17Nov2022" #Need to use different PopSim date for annual breeding
# HS.only <- "yes" #Do we only want to filter HS relationships?

#------------------------- Sample all juveniles -------------------------#
# sampling.scheme <- "sample.all.juvenile.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
# date.of.PopSim <- "17Nov2022" #Need to use different PopSim date for annual breeding
# HS.only <- "yes" #Do we only want to filter HS relationships?

#------------------------- Sample all ages -------------------------#
# sampling.scheme <- "sample.ALL.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
# date.of.PopSim <- "17Nov2022" #Need to use different PopSim date for annual breeding
# HS.only <- "no" #Do we only want to filter HS relationships?
# PO.only <- "no" #Do we only want to filter PO relationships?





########################## Read in sampling and other dataframes #########################
samples.df <- readRDS(file = paste0(PopSim.location, "sample.info_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))

samples.df %>% group_by(age.x) %>% summarize(n())

pop_size.df <- readRDS(file = paste0(PopSim.location, "pop.size_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))

truth.df <- readRDS(file = paste0(PopSim.location, "truth_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))

n_yrs <- max(pop_size.df$year)

########################## MCMC & model parameters #########################
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains


########################## Start simulation loop #########################
(iterations <- max(samples.df$iteration))

if(sample.props == "all"){
#Run with all sample sizes
(sample.prop_vec <- samples.df %>% distinct(sample.prop) %>% pull(sample.prop))
  (sim.samples.1 <- paste0(sample.prop_vec[1], "_prop.sampled"))
  (sim.samples.2 <- paste0(sample.prop_vec[2], "_prop.sampled"))
  (sim.samples.3 <- paste0(sample.prop_vec[3], "_prop.sampled"))
  (sim.samples.4 <- paste0(sample.prop_vec[4], "_prop.sampled"))
  
  } else if(is.numeric(sample.props) == TRUE){
#Run with just one sample size
    (sample.prop_vec <- samples.df %>%
       distinct(sample.prop) %>%
       dplyr::filter(sample.prop == sample.props) %>%
       pull(sample.prop)) 
    
    (sim.samples.1 <- paste0(sample.prop_vec[1], "_prop.sampled"))
    sim.samples.2 <- NULL
    sim.samples.3 <- NULL
    sim.samples.4 <- NULL
    }

# Initialize arrays for saving results
 results <- NULL
 sims.list.1 <- NULL
 sims.list.2 <- NULL
 sims.list.3 <- NULL
 sims.list.4 <- NULL
 mom.comps.tibble <- NULL
 dad.comps.tibble <- NULL

 for(est in 1:length(estimation.years)){
   #Set estimation year to 10 years in the past, or 0
   estimation.year <- estimation.years[est]
   
   iterations <- 2
   
 for(iter in 1:iterations) {
   #  set.seed(rseeds[iter])
   sim.start <- Sys.time()
   #rseed <- sample(1:1000000,1)
   set.seed(rseeds[iter])
   rseed <- rseeds[iter]

   for(s in 1:length(sample.prop_vec)){
     
     sample.proportion <- sample.prop_vec[s]
     sample.df_all.info <- samples.df %>% dplyr::filter(iteration == iter, sample.prop == sample.proportion)
     
     #Save total sample size
     sample.size.iter <- sample.df_all.info %>% 
       distinct(sample.size.total) %>%
       pull(sample.size.total)

        noDups.list <- split.dups(sample.df_all.info)
        first.capture <- noDups.list[[1]]
        later.capture <- noDups.list[[2]]
    
        #If misassigning ages
        if(exists("age.cv") == TRUE){
        set.seed(rseeds[iter])
        samples.miss <- misassign.ages(later.capture)
        filter1.out <- filter.samples(samples.miss) #Filter for full sibs
        } else {
          filter1.out <- filter.samples(later.capture) #Filter for full sibs
        }
        
 
    #Remove full sibs
    PO.samps.list <- filter1.out[[1]] #Output is a list where each list element corresponds to the offspring birth year and contains the potential parents and offspring for that year.
    HS.samps.df <- filter1.out[[2]] #Output is just the dataframe of samples but filtered for full siblings
    NoFullSibs.df <- filter1.out[[3]] #Save NoFullSibs.df so can downsample PO samples if needed
    
    #-------------Construct pairwise comparison matrix--------------
    #Input is 1) a list of potential parents and offspring for each offspring birth year, and 2) a dataframe with all samples, filtered to remove full siblings
    pairwise.out <- build.pairwise(filtered.samples.PO.list = PO.samps.list, filtered.samples.HS.df = HS.samps.df)
    
    #Save output as different dataframes; includes both HS and PO relationships (but can filter below)
    mom_comps.all <- pairwise.out[[1]]  
    dad_comps.all <- pairwise.out[[2]]
    positives.HS <- pairwise.out[[3]]

    head(mom_comps.all)
    head(dad_comps.all)
    
    mom_comps.all %>% group_by(type) %>% 
      summarize(sum(yes))
    
    dad_comps.all %>% group_by(type) %>% 
      summarize(sum(yes))

#If we only want to use one type of relationship, we'll filter here
    if(HS.only == "yes"){
    #Uncomment below to only run the HS model
    mom_comps.all <- mom_comps.all %>% filter(type == "HS")
    dad_comps.all <- dad_comps.all %>% filter(type == "HS")
    } else if(PO.only == "yes"){
    #Uncomment below to only run the PO model
    mom_comps.all <- mom_comps.all %>% filter(type == "PO")
    dad_comps.all <- dad_comps.all %>% filter(type == "PO")
    }
    
    if(sum(mom_comps.all$yes) == 0 | sum(dad_comps.all$yes) == 0){
      next
    } else {
      
      #HSPs detected
      mom.HSPs <- mom_comps.all %>% dplyr::filter(type == "HS") %>% 
        summarize(HSPs = sum(yes)) %>% 
        pull(HSPs)
      
      dad.HSPs <- dad_comps.all %>% dplyr::filter(type == "HS") %>% 
        summarize(HSPs = sum(yes)) %>% 
        pull(HSPs)
      
      if(HS.only != "yes"){
      #POPs detected
      mom.POPs <- mom_comps.all %>% dplyr::filter(type == "PO") %>% 
        summarize(POPs = sum(yes)) %>% 
        pull(POPs)
      
      dad.POPs <- dad_comps.all %>% dplyr::filter(type == "PO") %>% 
        summarize(POPs = sum(yes)) %>% 
        pull(POPs)      
      }
      
    ####------------------------ Fit CKMR model ----------------####
      
      #Define JAGS data and model, and run the MCMC engine
      set.seed(rseed)
      source("./02_Estimation.model/functions/RunJAGS_collated.R")

      #Calculate truth
      Nf.truth <- pop_size.df %>% dplyr::filter(iteration == iter,                   
                                                year >= ref.year.mom) %>% 
        summarize(females = round(mean(Female.adult.pop), 0)) %>% 
        pull(females)
      
      Nm.truth <- pop_size.df %>% dplyr::filter(iteration == iter,
                                                year >= ref.year.dad) %>% 
        summarize(males = round(mean(Male.adult.pop), 0)) %>% 
        pull(males)
      
    #Calculate expectations
    pop.size.tibble <- pop_size.df %>% dplyr::filter(iteration == iter)
    Exp <- calc.Exp(mom_comps.all, dad_comps.all)
    mom.Exp.HS <- Exp[[1]]
    mom.Exp.PO <- Exp[[2]]
    dad.Exp.HS <- Exp[[3]]
    dad.Exp.PO <- Exp[[4]]
    
    sampled.mothers <- unique(sample.df_all.info$mother.x)
    sampled.fathers <- unique(sample.df_all.info$father.x)
    
    #If there's no lambda in the model, then the true value of abundance is the mean over the estimated years. Here, we make this correction for the scenarios that do not include lambda in the model; for those that do, the values in truth.df are correct.
    if(scenario == "objective_1_model.validation" | 
       scenario == "scenario_2.1.1" | 
       scenario == "scenario_2.1.2" | 
       scenario == "scenario_2.1.3"){
    #Make truth equal to mean over years
    truth.iter <- truth.df %>% dplyr::filter(iteration == iter) %>% 
      mutate(all.truth = ifelse(parameter == "Nf", Nf.truth,
                                ifelse(parameter == "Nm", Nm.truth, all.truth)))
    } else {
      truth.iter <- truth.df %>% dplyr::filter(iteration == iter)
    }

    #If running with derived quantities    
    if(derived.quantities == "yes"){
      Nfb1.temp <- truth.iter %>% dplyr::filter(parameter == "Nf") %>% 
        mutate(parameter = "Nfb1", all.truth = breed.truth)
      
      Nfb2.temp <- truth.iter %>% dplyr::filter(parameter == "Nf") %>% 
        mutate(parameter = "Nfb2", all.truth = breed.truth)
      
      truth.iter <- truth.iter %>% bind_rows(Nfb1.temp, Nfb2.temp)
    }

    #Save some values related to samples so they can be added to the dataframe of results
      samples.iter <- samples.df %>% dplyr::filter(iteration == iter, sample.prop == sample.proportion) %>% 
      distinct(seed, iteration, sample.prop = sample.proportion, sample.size.total)

          
    results.temp <- model.summary2 %>% left_join(truth.iter, by = c("parameter", "iteration", "seed")) %>% 
      left_join(samples.iter, by = c("iteration", "seed"))
    
    if(HS.only == "yes"){
      if(length(jags_params) == 3){
      results.temp <- results.temp %>% 
        mutate(HSPs_detected = c(mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs),
               HSPs_expected = c(mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS),
               scenario = scenario,
               est.yr = estimation.year)
      } else if(length(jags_params) == 4){
        results.temp <- results.temp %>% 
          mutate(HSPs_detected = c(mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs + dad.HSPs),
                 HSPs_expected = c(mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS + dad.Exp.HS),
                 scenario = scenario,
                 est.yr = estimation.year)
      } else if(length(jags_params) == 5){
        results.temp <- results.temp %>% 
          mutate(HSPs_detected = c(mom.HSPs, mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs + dad.HSPs),
                 HSPs_expected = c(mom.Exp.HS, mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS + dad.Exp.HS),
                 scenario = scenario,
                 est.yr = estimation.year)
      } else if(length(jags_params) == 7){
        results.temp <- results.temp %>% 
          mutate(HSPs_detected = c(mom.HSPs, mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs, mom.HSPs),
                 HSPs_expected = c(mom.Exp.HS, mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS, mom.Exp.HS),
                 scenario = scenario,
                 est.yr = estimation.year)
      }
      
    } else if(HS.only != "yes"){
      if(length(jags_params) == 3){
      results.temp <- results.temp %>% 
        mutate(HSPs_detected = c(mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs),
               HSPs_expected = c(mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS),
               POPs_detected = c(mom.POPs, dad.POPs, mom.POPs + dad.POPs),
               POPs_expected = c(mom.Exp.PO, dad.Exp.PO, mom.Exp.PO + dad.Exp.PO),
               scenario = scenario,
               est.yr = estimation.year)
      } else if(length(jags_params) == 4){
        results.temp <- results.temp %>% 
          mutate(HSPs_detected = c(mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs + dad.HSPs),
                 HSPs_expected = c(mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS + dad.Exp.HS),
                 POPs_detected = c(mom.POPs, dad.POPs, mom.POPs + dad.POPs, mom.POPs + dad.POPs),
                 POPs_expected = c(mom.Exp.PO, dad.Exp.PO, mom.Exp.PO + dad.Exp.PO, mom.Exp.PO + dad.Exp.PO),
                 scenario = scenario,
                 est.yr = estimation.year)
      } else if(length(jags_params) == 5){
        results.temp <- results.temp %>% 
          mutate(HSPs_detected = c(mom.HSPs, mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs + dad.HSPs),
                 HSPs_expected = c(mom.Exp.HS, mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS + dad.Exp.HS),
                 POPs_detected = c(mom.POPs, mom.POPs, dad.POPs, mom.POPs + dad.POPs, mom.POPs + dad.POPs),
                 POPs_expected = c(mom.Exp.PO, mom.Exp.PO, dad.Exp.PO, mom.Exp.PO + dad.Exp.PO, mom.Exp.PO + dad.Exp.PO),
                 scenario = scenario,
                 est.yr = estimation.year)
      } else if(length(jags_params) == 7){
        results.temp <- results.temp %>% 
          mutate(HSPs_detected = c(mom.HSPs, mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs, mom.HSPs),
                 HSPs_expected = c(mom.Exp.HS, mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS, mom.Exp.HS),
                 POPs_detected = c(mom.POPs, mom.POPs, dad.POPs, mom.POPs + dad.POPs, mom.POPs + dad.POPs, mom.POPs, mom.POPs),
                 POPs_expected = c(mom.Exp.PO, mom.Exp.PO, dad.Exp.PO, mom.Exp.PO + dad.Exp.PO, mom.Exp.PO + dad.Exp.PO, mom.Exp.PO, mom.Exp.PO),
                 scenario = scenario,
                 est.yr = estimation.year)
      }
    }
    
    
    results <- rbind(results, results.temp)
    
    #Save mom and dad pairwise comparison dataframes
    mom_comps.all <- mom_comps.all %>% mutate(iteration = iter,
                                              sample.prop = sample.proportion,
                                              sample.size = sample.size.iter,
                                              seed = rseed)
    mom.comps.tibble <- rbind(mom.comps.tibble, mom_comps.all)
    
    dad_comps.all <- dad_comps.all %>% mutate(iteration = iter,
                                              sample.prop = sample.proportion,
                                              sample.size = sample.size.iter,
                                              seed = rseed)
    dad.comps.tibble <- rbind(dad.comps.tibble, dad_comps.all)
    
  } # End if/else statement
  } # end loop over sample sizes
    
  #-----------------Save output files iteratively--------------------
   if(iter %% 100 == 0){
     
#Results
    write.table(results, file = paste0(temp_location, results_prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", sampling.scheme, ".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
# 
#    #Model output for diagnostics
     saveRDS(sims.list.1, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.1, "_", MCMC.settings, "_", scenario, "_", model, "_", sampling.scheme))
# 
     saveRDS(sims.list.2, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.2, "_", MCMC.settings, "_", scenario, "_", model, "_", sampling.scheme))
# # 
     saveRDS(sims.list.3, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.3, "_", MCMC.settings, "_", scenario, "_", model, "_", sampling.scheme))
# #    
     saveRDS(sims.list.4, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.4, "_", MCMC.settings, "_", scenario, "_", model, "_", sampling.scheme))
# 
#    #Save pairwise comparisons matrices
    saveRDS(mom.comps.tibble, file = paste0(temp_location, mom.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", sampling.scheme))
#    
    saveRDS(dad.comps.tibble, file = paste0(temp_location, dad.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", sampling.scheme))

   }
    
      sim.end <- Sys.time()
   
   iter.time <- round(as.numeric(difftime(sim.end, sim.start, units = "mins")), 1)
   cat(paste0("\n Finished iteration ", iter, ". \n Took ", iter.time, " minutes \n\n"))
   } # end loop over iterations
 }# end loop over estimation years
  

    #If using all individuals for Nf truth, instead of breeders
   results2 <- results %>%
     mutate(relative_bias = round(((Q50 - all.truth)/all.truth)*100, 1)) %>% #Can change truth to breed.truth if looking for number of active breeders
     mutate(in_interval = ifelse(HPD2.5 < all.truth & all.truth < HPD97.5, "Y", "N")) %>%
     as_tibble()
   
#Within HPD interval?
results2 %>% group_by(sample.prop, parameter) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100)

#Median relative bias by sample size
 results2 %>% group_by(sample.prop, parameter, est.yr) %>% 
   dplyr::summarize(median.bias = median(relative_bias), n = n()) %>% 
   dplyr::arrange(desc(median.bias))

 #-----------------------------Save major output files---------------------------------------------
 #Home computer: Dell Precision
 
 #Save model estimates
write.table(results2, file = paste0(results_location, results_prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", sampling.scheme, ".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
 
 #Save draws from posterior for model diagnostics 
 saveRDS(sims.list.1, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.1, "_", MCMC.settings, "_", scenario, "_", model, "_", sampling.scheme)) #Sample size 1
 
  saveRDS(sims.list.2, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.2, "_", MCMC.settings, "_", scenario, "_", model, "_", sampling.scheme)) #Sample size 2
 # 
  saveRDS(sims.list.3, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.3, "_", MCMC.settings, "_", scenario, "_", model, "_", sampling.scheme)) #Sample size 3
 # 
  saveRDS(sims.list.4, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.4, "_", MCMC.settings, "_", scenario, "_", model, "_", sampling.scheme)) #Sample size 4
 
 #Save final pairwise comparison matrices
 saveRDS(mom.comps.tibble, file = paste0(results_location, mom.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", sampling.scheme))
 
 saveRDS(dad.comps.tibble, file = paste0(results_location, dad.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", sampling.scheme))

 if(exists("age.cv") == TRUE){
   #Save samples dataframe with misassigned ages
   saveRDS(samples.miss, file = paste0(results_location, samples.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", sampling.scheme))
}
#-------------Quick viz of results--------------#
#Box plot of relative bias
# ggplot(data=results2, aes(x=factor(sample.proportion))) +
#   geom_boxplot(aes(y=relative_bias, fill=parameter)) +
#   ylim(-100, 100) +
#   geom_hline(yintercept=0, col="black", size=1.25) +
#   annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
#   labs(x="Sample size", y="Relative bias", title="Relative Bias by sample size") +
#   scale_fill_brewer(palette="Set2") +
#   font("title", size = 10, face = "bold")



####################### End ##############################