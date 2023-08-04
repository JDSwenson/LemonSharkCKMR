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
PopSim.location <- "G://My Drive/Personal_Drive/R/CKMR/Population.simulations/Peer_review/" #Location of population simulation output files
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

############Specify common input prefixes####################
inSeeds <- "Seeds2022.04.15" #Seeds used for population simulation
date.of.PopSim <- "28Jul2023"

###########Specify which simulations to focus on########################
#s.scheme <- "target.YOY" #can be "target.YOY", "sample.all.juvenile.ages", or "sample.ALL.ages"
sample.props <- "all" #Either label this with the percent we want to target if just one (e.g., 1.5)) or if wanting to run over all sample proportions, set as "all"
objective <- 1 #Can be any number 1-5
scenario <- "scenario_1" #See Excel sheet with simulation scenarios: Simulation_log_key_UPDATED.xlsx on Google Drive
sample.scheme.vec <- c("target.YOY", "sample.all.juvenile.ages", "sample.ALL.ages")
est.yr.tests <- 1 #Can be 1 or 3. If 1, that means we will only estimate abundance for the birth year of the second oldest individual in the dataset; if 3, then we will estimate abundance for 10 years before that and the present as well.


if(objective == 1){
########################## Objective 1 #########################
#------------------------- Set input file locations -------------------------#
 PopSim.lambda <- "lambda.1" # Can be lambda.1 or lambda.variable
 PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding

#------------------------- Set output file locations -------------------------#
 MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.validation/Model.output/"
 results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.validation/Model.results/"

#------------------------- Objective 1 common settings -------------------------#
 fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
 jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 mating.periodicity <- 1 #number of years between mating for females
 
 survival.prior.mean <- adult.survival
 survival.prior.cv <- 0.05
 survival.prior.sd <- survival.prior.mean * survival.prior.cv
 scenario <- "objective_1_model.validation" #For naming output files and calculating truth.
 model <- "annual.model" #For naming output files
 
source("./02_Estimation.model/functions/Obj1.functions.R") #Changed name of script that includes pairwise comparison and other functions
jags_file = paste0(jags.model_location, "HS.PO_noLambda_annual_model_validation.txt") 

HS.only <- "no" #Do we only want to filter HS relationships?
PO.only <- "no" #Do we only want to filter PO relationships?

cat(paste0("Testing model validation"))

} else if(objective ==2){
########################## Objective 2 #########################
#------------------------- Set input file locations -------------------------#
PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding

#------------------------- Set output file locations -------------------------# 
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.2_population.growth/Nov2022/Model.output/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.2_population.growth/Nov2022/Model.results/"

#------------------------- Objective 2 common settings -------------------------#
estimation.years <- c(n_yrs - 10, n_yrs - 5, n_yrs)
fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
model <- "annual.model" #For naming output files
source("./02_Estimation.model/functions/Obj2.functions.R") #Changed name of script that includes pairwise comparison and other functions

cat(paste0("Testing population growth"))

if(scenario %in% c("scenario_2.1.1", "scenario_2.1.2", "scenario_2.1.3")){
#========================= Scenario 2.1 =========================
jags_file = paste0(jags.model_location, "HS.PO_noLambda_annual_model.txt") #Annual JAGS model w/o lambda

cat(paste0("with a naive model"))

if(scenario == "scenario_2.1.1"){
#------------------------- Scenario 2.1.1: Small population decline; no lambda in model
 PopSim.lambda <- "lambda.slight.decrease"
 scenario<- "scenario_2.1.1" #For naming output files
 jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")

 cat(paste0("and a slightly decreasing population."))
 
 } else if(scenario == "scenario_2.1.2"){

#------------------------- Scenario 2.1.2: Small population growth; no lambda in model
 PopSim.lambda <- "lambda.slight.increase"
 scenario<- "scenario_2.1.2" #For naming output files
 jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")

 cat(paste0("and a slightly increasing population."))
 
 } else if(scenario == "scenario_2.1.3"){

#------------------------- Scenario 2.1.3: Substantial population decline; no lambda in model
 PopSim.lambda <- "lambda.extreme" 
 scenario<- "scenario_2.1.3" #For naming output files
 jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 
 cat(paste0("and a severely decreasing population."))
 
 }
} else if(scenario %in% c("scenario_2.2.1", "scenario_2.2.2", "scenario_2.2.3")){

#========================= Scenario 2.2 =========================
jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")

cat(paste0("with an adapted model and narrow prior"))

if(scenario == "scenario_2.2.1"){
#------------------------- Scenario 2.2.1: Small population decline; lambda in model w/ tight prior
 PopSim.lambda <- "lambda.slight.decrease"
 scenario<- "scenario_2.2.1" #For naming output files
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")

 cat(paste0("and a slightly decreasing population."))
 
} else if(scenario == "scenario_2.2.2"){
#------------------------- Scenario 2.2.2: Small population growth; lambda in model w/ tight prior
 PopSim.lambda <- "lambda.slight.increase"
 scenario<- "scenario_2.2.2" #For naming output files
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")

 cat(paste0("and a slightly increasing population."))
 
} else if(scenario == "scenario_2.2.3"){
 
#------------------------- Scenario 2.2.3: Substantial population decline; lambda in model w/ tight prior
 PopSim.lambda <- "lambda.extreme" 
 scenario<- "scenario_2.2.3" #For naming output files
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 
 cat(paste0("and a severely decreasing population."))
 
 }
} else if(scenario %in% c("scenario_2.3.1", "scenario_2.3.2", "scenario_2.3.3")){

#========================= Scenario 2.3 =========================
jags_file = paste0(jags.model_location, "HS.PO_wideLambda_annual_model.txt")

cat(paste0("with an adapted model and diffuse prior"))

if(scenario == "scenario_2.3.1"){
#------------------------- Scenario 2.3.1: Small population decline; lambda in model w/ wide prior
 PopSim.lambda <- "lambda.slight.decrease"
 scenario<- "scenario_2.3.1" #For naming output files
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")

 cat(paste0("and a slightly decreasing population."))
 
} else if(scenario == "scenario_2.3.2"){
 
#------------------------- Scenario 2.3.2: Small population growth; lambda in model w/ wide prior
 PopSim.lambda <- "lambda.slight.increase"
 scenario<- "scenario_2.3.2" #For naming output files
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 
 cat(paste0("and a slightly increasing population."))

} else if(scenario == "scenario_2.3.3"){
 
#------------------------- Scenario 2.3.3: Substantial population decline; lambda in model w/ wide prior
 PopSim.lambda <- "lambda.extreme" 
 scenario<- "scenario_2.3.3" #For naming output files
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 
 cat(paste0("and a severely decreasing population."))
 }
}

} else if(objective == 3){
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

cat(paste0("Testing multiennial breeding"))

if(scenario == "scenario_3.1.1"){
#========================= Scenario 3.1 =========================
#------------------------- Scenario 3.1.1: Biennial breeding; psi = 1; annual model
 PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_3.1.1" #For naming output files
 model <- "annual.model"
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")
 
 cat(paste0("with 100% biennial breeders \nand a naive model."))

} else if(scenario == "scenario_3.1.2"){
#------------------------- Scenario 3.1.2: Biennial breeding; psi = 1; multiennial model
 PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_3.1.2" #For naming output files
 model <- "multiennial.model"
 jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt") #Specify the HS only model as the model for now. If the likelihood should include POPs, then we will override the model when we select the sampling scheme several lines down.
 
 cat(paste0("with 100% biennial breeders \nand a multiennial model."))

} else if(scenario == "scenario_3.2.1"){
#========================= Scenario 3.2 =========================
#------------------------- Scenario 3.2.1: Biennial breeding; psi = 0.9; annual model
 PopSim.breeding.schedule <- "biennial.breeding_psi0.90" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_3.2.1" #For naming output files
 model <- "annual.model"
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")
 
 cat(paste0("with 90% biennial breeders \nand a naive model."))

} else if(scenario == "scenario_3.2.2"){
 
#------------------------- Scenario 3.2.2: Biennial breeding; psi = 0.9; multiennial model
 PopSim.breeding.schedule <- "biennial.breeding_psi0.90" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_3.2.2" #For naming output files
 model <- "multiennial.model"
 jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt") #Specify the HS only model as the model for now. If the likelihood should include POPs, then we will override the model when we select the sampling scheme several lines down
 
 cat(paste0("with 90% biennial breeders \nand a multiennial model."))

} else if(scenario == "scenario_3.3.1"){
#========================= Scenario 3.3 =========================
#------------------------- Scenario 3.3.1: Biennial breeding; psi = 0.75; annual model
 PopSim.breeding.schedule <- "biennial.breeding_psi0.75" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_3.3.1" #For naming output files
 model <- "annual.model"
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")

cat(paste0("with 75% biennial breeders \nand a naive model."))

} else if(scenario == "scenario_3.3.2"){
#------------------------- Scenario 3.3.2: Biennial breeding; psi = 0.75; multiennial model
 PopSim.breeding.schedule <- "biennial.breeding_psi0.75" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_3.3.2" #For naming output files
 model <- "multiennial.model"
 jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt") #Specify the HS only model as the model for now. If the likelihood should include POPs, then we will override the model when we select the sampling scheme several lines down
 
 cat(paste0("with 75% biennial breeders \nand a multiennial model."))

} else if(scenario == "scenario_3.4.1"){

#========================= Scenario 3.4 =========================
#------------------------- Scenario 3.4.1: Biennial breeding; psi = 0.50; annual model
 PopSim.breeding.schedule <- "biennial.breeding_psi0.50" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_3.4.1" #For naming output files
 model <- "annual.model"
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")
 
 cat(paste0("with 50% biennial breeders \nand a naive model."))

} else if(scenario == "scenario_3.4.2"){
 
#------------------------- Scenario 3.4.2: Biennial breeding; psi = 0.50; multiennial model
 PopSim.breeding.schedule <- "biennial.breeding_psi0.50" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_3.4.2" #For naming output files
 model <- "multiennial.model"
 jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt") #Specify the HS only model as the model for now. If the likelihood should include POPs, then we will override the model when we select the sampling scheme several lines down
 
 cat(paste0("with 50% biennial breeders \nand a multiennial model."))

} else if(scenario == "scenario_3.5.1"){
#========================= Scenario 3.5: Triennial breeding =========================
#------------------------- Scenario 3.5.1: Triennial breeding; psi = 1; annual model
 PopSim.breeding.schedule <- "triennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_3.5.1" #For naming output files
 model <- "annual.model"
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")
mating.periodicity <- 3

cat(paste0("with 100% triennial breeders \nand a naive model."))

} else if(scenario == "scenario_3.5.2"){

#------------------------- Scenario 3.5.2: Triennial breeding; psi = 1; multiennial model
 PopSim.breeding.schedule <- "triennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_3.5.2" #For naming output files
 model <- "multiennial.model"
 jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt") #Specify the HS only model as the model for now. If the likelihood should include POPs, then we will override the model when we select the sampling scheme several lines down
 mating.periodicity <- 3
 
 cat(paste0("with 100% triennial breeders \nand a multiennial model."))

} else if(scenario == "scenario_3.6.1"){
 
#========================= Scenario 3.6: Biennial breeding w/ stochasticity
#------------------------- Scenario 3.6.1: Biennial breeding; psi = 1; annual model; 10% on-cycle breeders fail to breed; 10% off-cycle breeders do breed.
  date.of.PopSim <- "27Dec2022"
  PopSim.breeding.schedule <- "biennial.breeding_psi1_stochasticCycles" #Can be annual.breeding or biennial.breeding
  scenario <- "scenario_3.6.1" #For naming output files
  jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
  estimated.parameters <- paste0(jags_params, collapse = ",")
  jags_file = paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt")
  
  cat(paste0("with stochastic biennial breeders \nand a naive model."))

} else if(scenario == "scenario_3.6.2"){
#------------------------- Scenario 3.6.2: Biennial breeding; psi = 1; multiennial model; 10% on-cycle breeders fail to breed; 10% off-cycle breeders do breed.
date.of.PopSim <- "27Dec2022"
PopSim.breeding.schedule <- "biennial.breeding_psi1_stochasticCycles" #Can be annual.breeding or biennial.breeding
scenario <- "scenario_3.6.2" #For naming output files
model <- "multiennial.model"
jags_params = c("Nf", "psi", "Nm", "survival", "lambda", "Nfb1", "Nfb2") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")
jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt") #Specify the HS only model as the model for now. If the likelihood should include POPs, then we will override the model when we select the sampling scheme several lines down.
derived.quantities = "yes"

cat(paste0("with stochastic biennial breeders \nand a multiennial model."))
}

} else if(objective ==4){
########################## Objective 4 #########################
#------------------------- Set output file locations ------------------------- 
 MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.4_aging.error/Nov2022/Model.output/"
 results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.4_aging.error/Nov2022/Model.results/"
 outSeeds <- "Seeds2022.04.15"
 date.of.PopSim <- "22Nov2022"

#------------------------- Objective 4 common settings -------------------------
 estimation.years <- c(n_yrs - 5, n_yrs)
 fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
 PopSim.lambda <- "lambda.1"
 source("./02_Estimation.model/functions/Obj4.functions.R") #Changed name of script that includes pairwise comparison and other functions

 cat(paste0("Testing aging error"))
 
if(scenario == "scenario_4.1"){
 
#------------------------- Scenario 4.1: Biennial breeding; psi = 1; aging error 5% CV
 age.cv <- 0.05
mating.periodicity <- 2
 PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_4.1_5CV" #For naming output files
 model <- "multiennial.model"
 jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")

 cat(paste0("with biennial breeders and 5% CV."))
 
}else if(scenario == "scenario_4.2"){
 
#------------------------- Scenario 4.2: Biennial breeding; psi = 1; aging error 10% CV
 age.cv <- 0.10
mating.periodicity <- 2
 PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_4.2_10CV" #For naming output files
 model <- "multiennial.model"
 jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 
 cat(paste0("with biennial breeders and 10% CV."))

}else if(scenario == "scenario_4.3"){
  
#------------------------- Scenario 4.3: Biennial breeding; psi = 1; aging error 20% CV
 age.cv <- 0.20
mating.periodicity <- 2
 PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_4.3_20CV" #For naming output files
 model <- "multiennial.model"
 jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 
 cat(paste0("with biennial breeders and 20% CV."))

}else if(scenario == "scenario_4.4"){
 
#------------------------- Scenario 4.4: Annual breeding; aging error 5% CV
 age.cv <- 0.05
mating.periodicity <- 1
 PopSim.breeding.schedule <- "annual_breeding" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_4.4_5CV" #For naming output files
 model <- "annual.model"
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 
 cat(paste0("with annual breeders and 5% CV."))

}else if(scenario == "scenario_4.5"){

#------------------------- Scenario 4.5: Biennial breeding; psi = 1; aging error 10% CV
 age.cv <- 0.10
mating.periodicity <- 1
 PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_4.5_10CV" #For naming output files
 model <- "annual.model"
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 
 cat(paste0("with annual breeders and 10% CV."))

}else if(scenario == "scenario_4.6"){

#------------------------- Scenario 4.6: Biennial breeding; psi = 1; aging error 20% CV
 age.cv <- 0.20
 mating.periodicity <- 1
 PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding
 scenario <- "scenario_4.6_20CV" #For naming output files
 model <- "annual.model"
 jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 
 cat(paste0("with annual breeders and 20% CV."))
}
}


########################## Read in sampling and other dataframes #########################
#Read in sample dataframe and filter for the relevant sampling scheme
samples.df_all <- readRDS(file = paste0(PopSim.location, "sample.info_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule))

#Check that the sample dataframe was appropriately filtered
samples.df_all %>% group_by(sampling.scheme, age.x) %>% summarize(n())

#Read in dataframe with population size
pop_size.df <- readRDS(file = paste0(PopSim.location, "pop.size_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule))

#Should be 500
pop_size.df %>% dplyr::filter(year == 85) %>% nrow()
n_yrs <- max(pop_size.df$year)

#Read in dataframe with true values (most of these will be recalculated later)
truth.df <- readRDS(file = paste0(PopSim.location, "truth_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule))



#### TO DO: ADJUST CALCULATIONS FOR TRUTH DATAFRAME ####
# New truth dataframes are just survival and psi; want to: 
# 1) make tibble with reference years and iterations (and sampling scheme?),
# 2) subset pop_size.df for estimation years and pivot to long format, making the values in column "parameter" Nf, Nm, Nfb and Nmb.
# 3) add column for estimation year (ref.yr + 1 - call it "year"),
# 4) left_join with pop_size.df based on columns "iteration" and "year" and pull out true M and F abundance and breed abundance and store in columns all.truth and breed.truth.
#THEN
# 5) loop over each estimation year and calculate the true value of lambda
# 6) make any necessary df adjustments and add these values to the tibble




########################## MCMC & model parameters #########################
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains

(iterations <- max(samples.df_all$iteration))

if(sample.props == "all"){
#Run with all sample sizes
(sample.prop_vec <- samples.df_all %>% distinct(sample.prop) %>% pull(sample.prop))
  sim.samples.1 <- paste0(sample.prop_vec[1], "_prop.sampled")
  sim.samples.2 <- paste0(sample.prop_vec[2], "_prop.sampled")
  sim.samples.3 <- paste0(sample.prop_vec[3], "_prop.sampled")
  sim.samples.4 <- paste0(sample.prop_vec[4], "_prop.sampled")
  
  cat(paste0("Will run model with all sampling intensities"))
  
  } else if(is.numeric(sample.props) == TRUE){
#Run with just one sample size
    (sample.prop_vec <- samples.df_all %>%
       distinct(sample.prop) %>%
       dplyr::filter(sample.prop == sample.props) %>%
       pull(sample.prop)) 
    
    (sim.samples.1 <- paste0(sample.prop_vec[1], "_prop.sampled"))
    sim.samples.2 <- NULL
    sim.samples.3 <- NULL
    sim.samples.4 <- NULL
    
    cat(paste0("Will run model with sampling intensity of ", sample.props, "%", " per year."))
    }

########################## Start simulation loop #########################

# Initialize arrays for saving results
 results <- NULL
 sims.list.1 <- NULL
 sims.list.2 <- NULL
 sims.list.3 <- NULL
 sims.list.4 <- NULL
 mom.comps.tibble <- NULL
 dad.comps.tibble <- NULL

 # for(est in 1:length(estimation.years)){
 #   #Set estimation year to 10 years in the past, or 0
 #   estimation.year <- estimation.years[est]
 #   
 
 #iterations <- 2
   
 for(s.scheme in sample.scheme.vec){
   
   #Subset for samples from focal sampling scheme
   samples.df <- samples.df_all %>% dplyr::filter(sampling.scheme == s.scheme)
   
   #First iteration gives NA for reference year, but we can infer what it should be.
   if(s.scheme %in% c("sample.all.juvenile.ages", "sample.ALL.ages")){
     samples.df %>% replace_na(list(ref.yr = 76))

   } else {
     samples.df %>% replace_na(list(ref.yr = 87))
     
   }
 
   #Should just be one sampling scheme represented
   samples.df %>% dplyr::count(sampling.scheme)  
   
 for(iter in 1:iterations) {
   #  set.seed(rseeds[iter])
   sim.start <- Sys.time()
   #rseed <- sample(1:1000000,1)
   set.seed(rseeds[iter])
   rseed <- rseeds[iter]
   
   for(s in 1:length(sample.prop_vec)){
     
     sample.proportion <- sample.prop_vec[s]
     sample.df_all.info <- samples.df %>% dplyr::filter(iteration == iter, sample.prop == sample.proportion)
     
     #Save estimation year
     est.year.calibrate <- sample.df_all.info %>% distinct(ref.yr) %>% 
       mutate(estimation.yr = ref.yr + 1) %>% 
       pull(estimation.yr)
     
     #Make vector of estimation years
     estimation.years <- c(est.year.calibrate - 10, est.year.calibrate, n_yrs)
     
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
    full.sibs <- filter1.out[[4]] #Full siblings
    diff.cohort.sibs <- filter1.out[[5]] #Full siblings from different cohorts
    
    
    #Loop over estimation years
    for(est in 1:est.yr.tests){
      estimation.year <- estimation.years[est]
    
    
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

########### Calculate truth ###################
      #Calculate truth for lambda
      lambda.pop.df <- pop_size.df %>% dplyr::filter(iteration == iter)
      
      lambda.truth <- (lambda.pop.df$Total.adult.pop[n_yrs]/lambda.pop.df$Total.adult.pop[estimation.year])^(1/(n_yrs - estimation.year))
      
      lambda.truth.df <- tibble(estimation.year = estimation.year,
                                iteration = iter,
                                seed = rseed,
                                parameter = "lambda",
                                all.truth = lambda.truth)
      
      
      #Make dataframe of truth
      truth.iter <- pop_size.df %>% dplyr::filter(iteration == iter,
                                                year == estimation.year) %>% 
        mutate(Nfb1 = Num.mothers,
               Nmb1 = Num.fathers) %>%
        dplyr::select(Nf = Female.adult.pop,
                      Nfb1,
                      Nfb2 = Num.mothers,
                      Nm = Male.adult.pop,
                      Nmb1,
                      Nmb2 = Num.fathers,
                      estimation.year = year,
                      iteration = iteration,
                      seed = seed) %>% 
        pivot_longer(cols = starts_with("N"),
                     names_to = "parameter",
                     values_to = "all.truth") %>% 
        bind_rows(lambda.truth.df) %>% 
        mutate(sampling.scheme = s.scheme, 
               sample.prop = sample.proportion)
      
    #Calculate expectations
      #Define Nf and Nm as objects for calculations
      Nf.truth <- truth.iter %>% dplyr::filter(parameter == "Nf") %>% pull(all.truth)
      
      Nm.truth <- truth.iter %>% dplyr::filter(parameter == "Nm") %>% pull(all.truth)
      
      #Subset for just the present iteration
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
      Nf.truth <- pop.size.tibble %>% dplyr::filter(year >= estimation.year) %>% 
        summarize(mean.female.pop = mean(Female.adult.pop)) %>% 
        pull(mean.female.pop)
      
      Nfb.truth <- pop.size.tibble %>% dplyr::filter(year >= estimation.year) %>% 
        summarize(mean.mothers = mean(Num.mothers)) %>% 
        pull(mean.mothers)
      
      Nm.truth <- pop.size.tibble %>% dplyr::filter(year >= estimation.year) %>% 
        summarize(mean.male.pop = mean(Male.adult.pop)) %>% 
        pull(mean.male.pop)
      
      Nmb.truth <- pop.size.tibble %>% dplyr::filter(year >= estimation.year) %>% 
        summarize(mean.fathers = mean(Num.fathers)) %>% 
        pull(mean.fathers)
      
      truth.iter <- truth.iter %>% mutate(all.truth = ifelse(parameter == "Nf", Nf.truth, 
                                               ifelse(parameter == "Nfb1" | parameter == "Nfb2", Nfb.truth, 
                                                      ifelse(parameter == "Nm", Nm.truth, 
                                                             ifelse(parameter == "Nmb1" | parameter == "Nmb2", Nmb.truth, all.truth)))))
    }
    
    
    
    (truth.iter_all.params <- truth.df %>% dplyr::filter(sampling.scheme == s.scheme,
                               iteration == iter,
                               sample.prop == sample.proportion) %>% 
      mutate(estimation.year = ref.yr + 1) %>% 
      dplyr::select(-c(ref.yr, surv_min, surv_max, surv_AM, population.growth)) %>% 
      bind_rows(truth.iter) %>% 
      dplyr::select(parameter, all.truth, estimation.year, iteration, sampling.scheme, sample.prop, seed) %>% #Change order of columns
      dplyr::arrange(parameter))


      ##### PICK UP HERE ON 8/4/2023
    #Need to merge model.summary2 and truth.iter_all.params
    #Also export some of the bigger chunks of code to another script and source
    
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
    
  } # End if/else statement to account for no yes comparisons
  } #end loop over estimation years
   } #end loop over sample sizes
    
  #-----------------Save output files iteratively--------------------
   if(iter %% 100 == 0){
     
#Results
    write.table(results, file = paste0(temp_location, results_prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", s.scheme, ".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
# 
#    #Model output for diagnostics
     saveRDS(sims.list.1, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.1, "_", MCMC.settings, "_", scenario, "_", model, "_", s.scheme))
# 
     saveRDS(sims.list.2, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.2, "_", MCMC.settings, "_", scenario, "_", model, "_", s.scheme))
# # 
     saveRDS(sims.list.3, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.3, "_", MCMC.settings, "_", scenario, "_", model, "_", s.scheme))
# #    
     saveRDS(sims.list.4, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.4, "_", MCMC.settings, "_", scenario, "_", model, "_", s.scheme))

#    #Save pairwise comparisons matrices
    saveRDS(mom.comps.tibble, file = paste0(temp_location, mom.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", s.scheme))
#    
    saveRDS(dad.comps.tibble, file = paste0(temp_location, dad.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", s.scheme))

   }
    
      sim.end <- Sys.time()
   
   iter.time <- round(as.numeric(difftime(sim.end, sim.start, units = "mins")), 1)
   cat(paste0("\n Finished iteration ", iter, ". \n Took ", iter.time, " minutes \n\n"))
   } # end loop over sampling schemes
 }# end loop over iterations
  

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
write.table(results2, file = paste0(results_location, results_prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, ".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
 
 #Save draws from posterior for model diagnostics 
 saveRDS(sims.list.1, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.1, "_", MCMC.settings, "_", scenario, "_", model)) #Sample size 1
 
  saveRDS(sims.list.2, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.2, "_", MCMC.settings, "_", scenario, "_", model)) #Sample size 2
 # 
  saveRDS(sims.list.3, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.3, "_", MCMC.settings, "_", scenario, "_", model)) #Sample size 3
 # 
  saveRDS(sims.list.4, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.4, "_", MCMC.settings, "_", scenario, "_", model)) #Sample size 4
 
 #Save final pairwise comparison matrices
 saveRDS(mom.comps.tibble, file = paste0(results_location, mom.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model))
 
 saveRDS(dad.comps.tibble, file = paste0(results_location, dad.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model))

 if(exists("age.cv") == TRUE){
   #Save samples dataframe with misassigned ages
   saveRDS(samples.miss, file = paste0(results_location, samples.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model))
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