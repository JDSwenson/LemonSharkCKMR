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

rm(list=ls())

######################### Specify common output prefixes #########################
results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
jags.model.prefix <- "CKMR.JAGS_"
mom.comps.prefix <- "comparisons/mom.comps"
dad.comps.prefix <- "comparisons/dad.comps"

######################### Specify other common settings #########################
PopSim.location <- "G://My Drive/Personal_Drive/R/CKMR/Population.simulations/Nov2022/" #Location of population simulation output files
temp_location <- "~/R/working_directory/temp_results/" #Location to save temporary files in case run gets killed
jags.model_location <- "G://My Drive/Personal_Drive/R/CKMR/JAGS_models/" #Location of JAGS models

today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
date.of.simulation <- today

inSeeds <- "Seeds2022.04.15" #Seeds used for population simulation


########################## Objective 1 #########################
#------------------------- Set input file locations ------------------------- 
PopSim.lambda <- "lambda.1" # Can be lambda.1 or lambda.variable
PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding

#------------------------- Set output file locations ------------------------- 
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.validation/Model.output/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.validation/Model.results/"
outSeeds <- "Seeds2022.04.15"

#------------------------- Objective 1 common settings -------------------------
fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")
mating.periodicity <- 1 #number of years between mating for females
#Survival prior info
survival.prior.mean <- adult.survival
survival.prior.cv <- 0.05
survival.prior.sd <- survival.prior.mean * survival.prior.cv
scenario <- "objective_1_model.validation" #For naming output files
model <- "annual.model" #For naming output files

source("./Objective.1_model.validation/functions/Obj1.functions.R") #Changed name of script that includes pairwise comparison and other functions

#------------------------- Target YOY -------------------------#
sampling.scheme <- "target.YOY" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
date.of.PopSim <- "17Nov2022"
HS.only <- "yes" #Do we only want to filter HS relationships?
PO.only <- "no" #Do we only want to filter PO relationships?

#------------------------- Sample all juveniles -------------------------#
sampling.scheme <- "sample.all.juvenile.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
date.of.PopSim <- "17Nov2022"
HS.only <- "yes" #Do we only want to filter HS relationships?
PO.only <- "no" #Do we only want to filter PO relationships?

#------------------------- Sample all ages -------------------------#
sampling.scheme <- "sample.ALL.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
date.of.PopSim <- "17Nov2022"
HS.only <- "no" #Do we only want to filter HS relationships?
PO.only <- "no" #Do we only want to filter PO relationships?




########################## Objective 2 #########################
#------------------------- Set input file locations -------------------------
PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding

#------------------------- Set output file locations ------------------------- 
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.2_population.growth/Nov2022/Model.output/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.2_population.growth/Nov2022/Model.results/"
outSeeds <- "Seeds2022.04.15"

#------------------------- Objective 2 common settings -------------------------
fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
model <- "annual.model" #For naming output files

#========================= Scenario 2.1 =========================
#------------------------- Scenario 2.1.1: Small population decline; no lambda in model
PopSim.lambda <- "lambda.slight.decrease"
scenario<- "scenario_2.1.1" #For naming output files
jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 2.1.2: Small population growth; no lambda in model
PopSim.lambda <- "lambda.slight.increase"
scenario<- "scenario_2.1.2" #For naming output files
jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 2.1.3: Substantial population decline; no lambda in model
PopSim.lambda <- "lambda.extreme" 
scenario<- "scenario_2.1.3" #For naming output files
jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")


#========================= Scenario 2.2 =========================
#------------------------- Scenario 2.2.1: Small population decline; lambda in model w/ tight prior
PopSim.lambda <- "lambda.slight.decrease"
scenario<- "scenario_2.2.1" #For naming output files
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 2.2.2: Small population growth; lambda in model w/ tight prior
PopSim.lambda <- "lambda.slight.increase"
scenario<- "scenario_2.2.2" #For naming output files
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 2.2.3: Substantial population decline; lambda in model w/ tight prior
PopSim.lambda <- "lambda.extreme" 
scenario<- "scenario_2.2.3" #For naming output files
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")


#========================= Scenario 2.3 =========================
#------------------------- Scenario 2.3.1: Small population decline; lambda in model w/ wide prior
PopSim.lambda <- "lambda.slight.decrease"
scenario<- "scenario_2.3.1" #For naming output files
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 2.3.2: Small population growth; lambda in model w/ wide prior
PopSim.lambda <- "lambda.slight.increase"
scenario<- "scenario_2.3.2" #For naming output files
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 2.3.3: Substantial population decline; lambda in model w/ wide prior
PopSim.lambda <- "lambda.extreme" 
scenario<- "scenario_2.3.3" #For naming output files
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")

#========================== Sampling
#------------------------- Target YOY -------------------------#
sampling.scheme <- "target.YOY" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
date.of.PopSim <- "17Nov2022"
HS.only <- "yes" #Do we only want to filter HS relationships?
PO.only <- "no" #Do we only want to filter PO relationships?

#------------------------- Sample all juveniles -------------------------#
sampling.scheme <- "sample.all.juvenile.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
date.of.PopSim <- "17Nov2022"
HS.only <- "yes" #Do we only want to filter HS relationships?
PO.only <- "no" #Do we only want to filter PO relationships?

#------------------------- Sample all ages -------------------------#
sampling.scheme <- "sample.ALL.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
date.of.PopSim <- "17Nov2022"
HS.only <- "no" #Do we only want to filter HS relationships?
PO.only <- "no" #Do we only want to filter PO relationships?



########################## Objective 3 #########################
#------------------------- Set output file locations ------------------------- 
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.3_intermittent.breeding/Nov2022/Model.output/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.3_intermittent.breeding/Nov2022/Model.results/"
outSeeds <- "Seeds2022.04.15"

#------------------------- Objective 3 common settings -------------------------
fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
PopSim.lambda <- "lambda.1"

#========================= Scenario 3.1 =========================
#------------------------- Scenario 3.1.1: Biennial breeding; psi = 1; annual model
PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
scenario <- "scenario_3.1.1" #For naming output files
model <- "annual.model"
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 3.1.2: Biennial breeding; psi = 1; multiennial model
PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
scenario <- "scenario_3.1.2" #For naming output files
model <- "multiennial.model"
jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")


#========================= Scenario 3.2 =========================
#------------------------- Scenario 3.2.1: Biennial breeding; psi = 0.9; annual model
PopSim.breeding.schedule <- "biennial.breeding_psi0.90" #Can be annual.breeding or biennial.breeding
scenario <- "scenario_3.2.1" #For naming output files
model <- "annual.model"
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 3.2.2: Biennial breeding; psi = 0.9; multiennial model
PopSim.breeding.schedule <- "biennial.breeding_psi0.90" #Can be annual.breeding or biennial.breeding
scenario <- "scenario_3.2.2" #For naming output files
model <- "multiennial.model"
jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")


#========================= Scenario 3.3 =========================
#------------------------- Scenario 3.3.1: Biennial breeding; psi = 0.9; annual model
PopSim.breeding.schedule <- "biennial.breeding_psi0.75" #Can be annual.breeding or biennial.breeding
scenario <- "scenario_3.3.1" #For naming output files
model <- "annual.model"
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 3.3.2: Biennial breeding; psi = 0.9; multiennial model
PopSim.breeding.schedule <- "biennial.breeding_psi0.75" #Can be annual.breeding or biennial.breeding
scenario <- "scenario_3.3.2" #For naming output files
model <- "multiennial.model"
jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")



#========================= Scenario 3.4 =========================
#------------------------- Scenario 3.4.1: Biennial breeding; psi = 0.9; annual model
PopSim.breeding.schedule <- "biennial.breeding_psi0.50" #Can be annual.breeding or biennial.breeding
scenario <- "scenario_3.4.1" #For naming output files
model <- "annual.model"
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 3.4.2: Biennial breeding; psi = 0.9; multiennial model
PopSim.breeding.schedule <- "biennial.breeding_psi0.50" #Can be annual.breeding or biennial.breeding
scenario <- "scenario_3.4.2" #For naming output files
model <- "multiennial.model"
jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")


#========================= Scenario 3.5 =========================
#------------------------- Scenario 3.5.1: Biennial breeding; psi = 0.9; annual model
PopSim.breeding.schedule <- "triennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
scenario <- "scenario_3.5.1" #For naming output files
model <- "annual.model"
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")

#------------------------- Scenario 3.5.2: Biennial breeding; psi = 0.9; multiennial model
PopSim.breeding.schedule <- "triennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
scenario <- "scenario_3.5.2" #For naming output files
model <- "multiennial.model"
jags_params = c("Nf", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")


#========================== Sampling
#------------------------- Target YOY -------------------------#
sampling.scheme <- "target.YOY" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
date.of.PopSim <- "22Nov2022"
HS.only <- "yes" #Do we only want to filter HS relationships?
PO.only <- "no" #Do we only want to filter PO relationships?

#------------------------- Sample all juveniles -------------------------#
sampling.scheme <- "sample.all.juvenile.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
date.of.PopSim <- "22Nov2022"
HS.only <- "yes" #Do we only want to filter HS relationships?
PO.only <- "no" #Do we only want to filter PO relationships?

#------------------------- Sample all ages -------------------------#
sampling.scheme <- "sample.ALL.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
date.of.PopSim <- "22Nov2022"
HS.only <- "no" #Do we only want to filter HS relationships?
PO.only <- "no" #Do we only want to filter PO relationships?




########################## Objective 4 #########################
#------------------------- Set output file locations ------------------------- 
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.4_aging.error/Nov2022/Model.output/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.4_aging.error/Nov2022/Model.results/"
outSeeds <- "Seeds2022.04.15"

#------------------------- Objective 3 common settings -------------------------
fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
PopSim.lambda <- "lambda.1"

#========================= Scenario 3.1 =========================
#------------------------- Scenario 3.1.1: Biennial breeding; psi = 1; annual model
PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
scenario <- "scenario_3.1.1" #For naming output files
model <- "annual.model"
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")




#--------------------Leslie matrix parameters----------------------------
#These may not be relevant for model validation
rseeds <- readRDS("rseeds_2022.04.15.rda")

adult.survival <- 0.825 # CHANGED FROM 0.825; Adult survival
repro.age <- 12 # set age of reproductive maturity
max.age <- 50 #set the maximum age allowed in the simulation


#---------------------- Read in sampling and other dataframes --------------------------
samples.df <- readRDS(file = paste0(PopSim.location, "sample.info_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))

samples.df %>% group_by(age.x) %>% summarize(n())

pop_size.df <- readRDS(file = paste0(PopSim.location, "pop.size_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))

truth.df <- readRDS(file = paste0(PopSim.location, "truth_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))

n_yrs <- max(pop_size.df$year)
estimation.year <- n_yrs - 5

#----------------------- MCMC & model parameters ----------------------#
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains


####-------------- Start simulation loop ----------------------
(iterations <- max(samples.df$iteration))
(sample.sizes <- samples.df %>% distinct(sample.prop) %>% pull(sample.prop)) #Subset for sample size of 1%


# Initialize arrays for saving results
 results <- NULL
 sims.list.1 <- NULL
 sims.list.2 <- NULL
 sims.list.3 <- NULL
 sims.list.4 <- NULL
 mom.comps.tibble <- NULL
 dad.comps.tibble <- NULL

 (sim.samples.1 <- paste0(sample.sizes[1], "prop.sampled"))
 (sim.samples.2 <- paste0(sample.sizes[2], "prop.sampled"))
 (sim.samples.3 <- paste0(sample.sizes[3], "prop.sampled"))
 (sim.samples.4 <- paste0(sample.sizes[4], "prop.sampled"))


 for(iter in 1:iterations) {
   #  set.seed(rseeds[iter])
   sim.start <- Sys.time()
   #rseed <- sample(1:1000000,1)
   set.seed(rseeds[iter])
   rseed <- rseeds[iter]

   for(s in 1:length(sample.sizes)){
     
     sample.size <- sample.sizes[s]
     sample.df_all.info <- samples.df %>% dplyr::filter(iteration == iter, sample.prop == sample.size)

        noDups.list <- split.dups(sample.df_all.info)
        first.capture <- noDups.list[[1]]
        later.capture <- noDups.list[[2]]
    
 
    #Remove full sibs
    filter1.out <- filter.samples(later.capture) #Filter for full sibs
    PO.samps.list <- filter1.out[[1]] #Output is a list where each list element corresponds to the offspring birth year and contains the potential parents and offspring for that year.
    HS.samps.df <- filter1.out[[2]] #Output is just the dataframe of samples but filtered for full siblings
    NoFullSibs.df <- filter1.out[[3]] #Save NoFullSibs.df so can downsample PO samples if needed
    
    #-------------Construct pairwise comparison matrix--------------
    #Input is 1) a list of potential parents and offspring for each offspring birth year, and 2) a dataframe with all samples, filtered to remove full siblings
    pairwise.out <- build.pairwise(filtered.samples.PO.list = PO.samps.list, filtered.samples.HS.df = HS.samps.df)
    
    #Save output as different dataframes; includes both HS and PO relationships (but can filter below)
    #Can uncomment to include/exclude different comparisons
    mom_comps.all <- pairwise.out[[1]] #%>% 
      #dplyr::filter(mort.yrs < repro.age & yes >= 1) %>% 
      # dplyr::select(ref.year, all, yes, mort.yrs, type) %>% 
      # group_by(mort.yrs) %>% 
      # summarize(all = sum(all), yes = sum(yes)) %>% 
      # mutate(pop.growth.yrs = 0, ref.year = 90, type = "HS|PO")
    
    dad_comps.all <- pairwise.out[[2]] #%>% 
      # dplyr::filter(mort.yrs < repro.age & yes >=1) %>% 
      # dplyr::select(ref.year, all, yes, mort.yrs, type) %>% 
      # group_by(mort.yrs) %>% 
      # summarize(all = sum(all), yes = sum(yes)) %>% 
      # mutate(pop.growth.yrs = 0, ref.year = 90, type = "HS|PO")
    positives.HS <- pairwise.out[[3]]

    head(mom_comps.all)
    head(dad_comps.all)
    
    mom_comps.all %>% group_by(type) %>% 
      summarize(sum(yes))
    
    dad_comps.all %>% group_by(type) %>% 
      summarize(sum(yes))
    
    #Calculate psi truth
    psi.df <- samples.df %>% dplyr::filter(iteration == iter) %>% #Use all samples from this iteration to calculate the truth
      distinct(indv.name, .keep_all = TRUE) %>% #Make sure we aren't using duplicated individuals
      group_by(repro.strategy) %>%
      summarize(number = n())
    
    (psi.truth <- round(1 - psi.df[psi.df$repro.strategy == "non-conformist",2]/sum(psi.df$number), 3) %>% 
      pull(number)) #Calculate number of non-conformists over number of conformists in samples for each iteration dataset
    

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
      
      mom.HSPs <- sum(mom_comps.all$yes)
      dad.HSPs <- sum(dad_comps.all$yes)
      
      
    # ####------------------------ Fit CKMR model ----------------####
    #Define JAGS data and model, and run the MCMC engine
      set.seed(rseed)
    source("Objective.1_model.validation/functions/Obj1.1_run.JAGS_HS.only.R")

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
    
    
    #Make truth equal to mean over years
    truth.iter <- truth.df %>% dplyr::filter(iteration == iter) %>% 
      mutate(all.truth = ifelse(parameter == "Nf", Nf.truth,
                                ifelse(parameter == "Nm", Nm.truth, all.truth)))

      samples.iter <- samples.df %>% dplyr::filter(iteration == iter, sample.prop == sample.size) %>% 
      distinct(seed, iteration, sample.prop.juvs = sample.prop, sample.size.juvs)
    
    results.temp <- model.summary2 %>% left_join(truth.iter, by = c("parameter", "iteration", "seed")) %>% 
      left_join(samples.iter, by = c("iteration", "seed")) %>% 
      mutate(HSPs_detected = c(mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs),
             HSPs_expected = c(mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS),
             purpose = purpose)
    
    results <- rbind(results, results.temp)
    
    
    #Save mom and dad pairwise comparison dataframes
    mom_comps.all <- mom_comps.all %>% mutate(iteration = iter,
                                              sample.prop.juvs = sample.size,
                                              sample.size.juvs = nrow(HS.samps.df),
                                              seed = rseed)
    mom.comps.tibble <- rbind(mom.comps.tibble, mom_comps.all)
    
    dad_comps.all <- dad_comps.all %>% mutate(iteration = iter,
                                              sample.prop.juvs = sample.size,
                                              sample.size.juvs = nrow(HS.samps.df),
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
  

    #If using all individuals for Nf truth, instead of breeders
   results2 <- results %>%
     mutate(relative_bias = round(((Q50 - all.truth)/all.truth)*100, 1)) %>% #Can change truth to breed.truth if looking for number of active breeders
     mutate(in_interval = ifelse(HPD2.5 < all.truth & all.truth < HPD97.5, "Y", "N")) %>%
     as_tibble()
   
#Within HPD interval?
results2 %>% group_by(sample.prop.juvs, parameter) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100)

#Median relative bias by sample size
 results2 %>% group_by(sample.prop.juvs, parameter) %>% 
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


#-------------Quick viz of results--------------#
#Box plot of relative bias
# ggplot(data=results2, aes(x=factor(sample.prop.juvs))) +
#   geom_boxplot(aes(y=relative_bias, fill=parameter)) +
#   ylim(-100, 100) +
#   geom_hline(yintercept=0, col="black", size=1.25) +
#   annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
#   labs(x="Sample size", y="Relative bias", title="Relative Bias by sample size") +
#   scale_fill_brewer(palette="Set2") +
#   font("title", size = 10, face = "bold")



####################### End ##############################