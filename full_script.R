#Load packages
library(tidyverse) # safe to ignore conflicts with filter() and lag()
library(Rlab)
library(runjags)
library(postpack)
library(coda)
library(tidyverse) # safe to ignore conflicts with filter() and lag()
library(MASS)
library(popbio)
library(mpmtools)
library(ggpubr)
library(rjags)
library(R2jags)
library(jagsUI)
library(FSA)


rm(list=ls())

#source("functions/Dovi_IBS_SB_test.assign.conformity_mat12.R") #All biennial breeders reproduce for the first time at age 12
source("./01_Data.generating.model/functions/Dovi_IBS_SB_test.assign.conformity_mat12OR13_liz.R") #Half of biennial breeders reproduce for the first time at age 12; the other half at age 13.

#################### Set output file locations and labels ####################
temp_location <- "output/"
output.location <- "output/"
parents_prefix <- "parents.breakdown"
sample_prefix <- "sample.info"
pop.size.prefix <- "pop.size"
truth.prefix <- "truth"
PopSim.lambda <- "lambda.1" # Can be lambda.1 or lambda.variable
PopSim.breeding.schedule <- "annual.breeding" #Can be annual.breeding or biennial.breeding

#------------------------- Set output file locations -------------------------#
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.validation/Model.output/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.1_model.validation/Model.results/"

#################### Simulation parameters ####################
init.adult.pop.size <- 1000 # CHANGED FROM 3000; Initial adult population size
init.prop.female <- 0.5 # proportion of the initial population size that is female
birth.sex.ratio <- c(0.5,0.5) # The probability that each baby is F:M - has to add up to 1
YOY.survival <- 0.7 # CHANGED FROM 0.8; young of year survival
juvenile.survival <- 0.8 # CHANGED FROM 0.9; juvenile survival
Adult.survival <- 0.825 # CHANGED FROM 0.825; Adult survival
repro.age <- 12 # set age of reproductive maturity
max.age <- maxAge <- 50 #set the maximum age allowed in the simulation
num.mates <- c(1:3) #1  #c(1:3) #vector of potential number of mates per mating
f <- (1-Adult.survival)/(YOY.survival * juvenile.survival^11) # adult fecundity at equilibrium if no age truncation



#################### Breeding schedule ######################
#------------------------------ Biennial ------------------------------
mating.periodicity <- 2 #number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.

#============================== psi 1 ==============================
  breeding.schedule <- "biennial.breeding_psi1"
  non.conformists <- 0
  percent.breed_off.cycle <- 0 #Percent of off-cycle mothers that will breed each year
  percent.skip_on.cycle <- 0 #Percent of on-cycle mothers that will skip breeding each year

# Adjust fecundity ==============================================================
## for effective number of breeders each year, mating cycle, number of mates ====
#(from liz) ====
psi <- 1-non.conformists
ff <- mating.periodicity/(mating.periodicity-psi*mating.periodicity+psi)*f/init.prop.female/mean(num.mates)
ff
# ====
#################### Population growth ####################
#------------------------------Stable------------------------------
population.growth <- "lambda.1"

#################### Sampling scheme ######################
#============================== target YOY ==============================
sampling.scheme <- "target.YOY"

#============================== sample all juveniles ==============================
#sampling.scheme <- "sample.all.juvenile.ages"

#============================== sample all ages ==============================
#sampling.scheme <- "sample.ALL.ages"
 
#-------------------Set date and load seeds----------------------------
today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
date.of.simulation <- today

rseeds <- readRDS("rseeds_2022.04.15.rda")
seeds <- "Seeds2022.04.15"
# create some seeds: 
#rseeds <- trunc(1529375630*runif(1000, min=0, max=1 ) )
#rseeds[1]  <- 746160703

#----------------------- DATA-GENERATING MODEL --------------------
# Note on sequencing: Birthdays happen at beginning of each year, followed by mating, then death (i.e. a female who dies in year 10 can still give birth and have surviving pups in year 10)

#Stable age distribution
props <- rep(NA, max.age+1)
props[1] <- f
props[2] <- f * YOY.survival
for (y in 3:(repro.age+1)) props[y] <- props[y-1] * juvenile.survival #+1 because of age 0 individuals

for (y in (repro.age+2):(max.age+1)) props[y] <- props[y-1] * Adult.survival #+2 because of age 0 individuals
prop.Adult <- sum(props[(repro.age+1):(max.age+1)])/sum(props)
Nages <- round(props[-1] * init.adult.pop.size) 
init.pop.size <- sum(Nages) # all ages except YOYs

#Set length of simulation and estimation year
burn.in <- 40 # number of years to use as simulation burn in period
Num.years <- 50 # The number of years to run in the simulation beyond the burn in
n_yrs <- burn.in + Num.years #Total number of simulation years
estimation.year <- n_yrs - 5 # Set year of estimation for truth calculations


#--------------------- Sampling parameters ---------------------
sample.years <- c(n_yrs - c(3:0)) #For four years of sampling
sample.vec.prop <- c(.5, 1, 1.5, 2)
sample.vec.prop <- 1.5

####-------------- Prep simulation ----------------------
# Moved sampling below so extract different sample sizes from same population
iterations <- 5  # 1 just to look at output     500 #Number of iterations to loop over


# Initialize arrays for saving results
 results <- NULL
 sims.list.1 <- NULL
 sims.list.2 <- NULL
 sims.list.3 <- NULL
 sims.list.4 <- NULL
 sample.info <- NULL
 parents.tibble_all <- NULL
 pop.size.tibble_all <- NULL
 mom.comps.tibble <- NULL
 dad.comps.tibble <- NULL
 truth.all <- NULL

 sim.samples.1 <- paste0(sample.vec.prop[1], "prop.sampled")
 sim.samples.2 <- paste0(sample.vec.prop[2], "prop.sampled")
 sim.samples.3 <- paste0(sample.vec.prop[3], "prop.sampled")
 sim.samples.4 <- paste0(sample.vec.prop[4], "prop.sampled")

 #
 #---------------------Initialize array from previous checkpoint--------------------------
#Results
#  results <- read_csv(paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose, "_iter_", iter, ".csv"))
# # 
# # #Model output for diagnostics
#  sims.list.1 <- readRDS(file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose))
# # 
#  sims.list.2 <- readRDS(file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose))
# # 
#  sims.list.3 <- readRDS(file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose))
# 
# # Detailed info on samples and parents to examine in more detail
#  sample.info <- readRDS(file = paste0(temp_location, sample_prefix, "_", date.of.simulation, "_", seeds, "_", purpose))

 ####-------------- Start simulation loop ----------------------
 for(iter in 1:iterations) {
   
   if(population.growth == "lambda.extreme"){
     juvenile.survival <- 0.8 # Reset at the beginning of each simulation
     Adult.survival <- 0.825 # Reset at the beginning of each simulation
     
     (Ma <- -log(Adult.survival)) #Mortality of adults
     (Mj <- -log(juvenile.survival)) #Mortality of juveniles
     (Sa <- exp(-Ma)) #Survival of adults
     (Sj <- exp(-Mj)) #Survival of juveniles
     
     set.seed(rseeds[iter])
     Mx <- runif(1, min = 0.05, max = 0.1) #Extra mortality
     (Sa <- exp(-Ma - Mx)) #Survival of adults
     (Sj <- exp(-Mj - Mx)) #Survival of juveniles
   }
   
   sim.start <- Sys.time()
   rseed <- rseeds[iter]
   set.seed(rseed)


  #Run individual based simulation
  out <- simulate.pop(init.pop.size = init.pop.size, 
               init.prop.female = init.prop.female,
               Nages = Nages,
               mating.periodicity = mating.periodicity,
               repro.age = repro.age,
               YOY.survival = YOY.survival,
               juvenile.survival = juvenile.survival,
               Adult.survival = Adult.survival,
               max.age = max.age,
               num.mates = num.mates,
               ff = ff,
               burn.in = burn.in,
               Num.years = Num.years)

  #Save simulation output as objects
  loopy.list <- out[[1]] #List of dataframes for each year of simulation
  pop.size.tibble <- out[[2]] %>%  #population parameters for each year of simulation
    as_tibble() %>% 
    mutate(seed = rseed, iteration = iter)

  parents.tibble <- out[[3]] %>% 
    dplyr::filter(year >= 50) %>% #Tibble for each parent for each year to check the distribution later
    mutate(seed = rseed, iteration = iter)
  
  #organize results and calculate summary statistics from the simulation
  source("./01_Data.generating.model/functions/query_results_PopSim.R")
  
  #-----------------------Collect samples-------------------------
  #Loop over sample sizes stored in sample.vec  
  for(samps in 1:length(sample.vec.prop)){
    
    sample.prop <- sample.vec.prop[samps]
    
    #Initialize sample dataframes
    sample.df_all.info <- NULL
    sample.df_temp <- NULL
    
    #Sample population each year in sample.years and make dataframe of samples with all metadata
    set.seed(rseed)
    
        for(i in sample.years){
          
          #Extract the relevant row from the pop size dataframe
          pop.size.yr <- pop.size.tibble %>% dplyr::filter(year == i)
          
          sample.size <- pop.size.yr %>%
            mutate(sample.size = round(population_size*(sample.prop/100), 0)) %>%
            pull(sample.size)
          
          if(sampling.scheme == "target.YOY"){ #If targeting YOY for juvenile samples
            #Set number of samples to a specific proportion of the population
            
            #Sample YOY only for half-sib analysis
            sample.df_temp <- loopy.list[[i]] %>% mutate(capture.year = i) %>%
              dplyr::filter(age.x == 0) %>%
              dplyr::slice_sample(n = sample.size) # #Sample each year WITHOUT replacement
            
            
            } else if(sampling.scheme == "sample.all.juvenile.ages"){ #If sampling juveniles
              sample.df_temp <- loopy.list[[i]] %>% dplyr::filter(age.x < repro.age & age.x > 0) %>% 
                mutate(capture.year = i) %>%
                dplyr::slice_sample(n = sample.size) #Sample each year WITHOUT replacement (doesn't affect cross-year sampling since it's in a loop)
              
              
              } else if(sampling.scheme == "sample.ALL.ages"){
              
              sample.df_temp <- loopy.list[[i]] %>% mutate(capture.year = i) %>%
                dplyr::slice_sample(n = sample.size) #Sample each year WITHOUT replacement (doesn't affect cross-year sampling since it's in a loop)
        
      }      

      #Combine samples from all years
        sample.df_all.info <- rbind(sample.df_all.info, sample.df_temp)
    }
    

    sampled.mothers <- unique(sample.df_all.info$mother.x)
    sampled.fathers <- unique(sample.df_all.info$father.x)
    
    #Compile results and summary statistics from simulation to compare estimates
    source("./01_Data.generating.model/functions/PopSim_truth.R")
    
    #Save info for samples to examine in more detail
    sample.df_all.info <- sample.df_all.info %>% 
      mutate(sample.size.yr = sample.size,
             sampling.scheme = sampling.scheme,
             iteration = iter,
             seed = rseed,
             sample.prop = sample.prop) %>% 
      mutate(sample.size.total = sample.size.yr * length(sample.years))
    
    sample.info <- rbind(sample.info, sample.df_all.info) %>% 
      as_tibble()
  
  } # end loop over sample sizes

  #-----------------Save output files iteratively--------------------
  
  sim.samples.1 <- paste0(sample.vec.prop[1], "prop.sampled")
  sim.samples.2 <- paste0(sample.vec.prop[2], "prop.sampled")
  sim.samples.3 <- paste0(sample.vec.prop[3], "prop.sampled")
  sim.samples.4 <- paste0(sample.vec.prop[4], "prop.sampled")
  
  #Save parents tibble
  parents.tibble_all <- bind_rows(parents.tibble_all, parents.tibble)
  
#  saveRDS(parents.tibble_all, file = paste0(temp_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule, "_iter_", iter))
  
  # Detailed info on population size
  pop.size.tibble_all <- bind_rows(pop.size.tibble_all, pop.size.tibble)
  
#  saveRDS(pop.size.tibble_all, file = paste0(temp_location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule, "_iter_", iter))

  #True values
  truth.all <- bind_rows(truth.all, true.values)
  
#  saveRDS(truth.all, file = paste0(temp_location, truth.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule, "_iter_", iter))
  
# Detailed info on samples and parents to examine in more detail
#  saveRDS(sample.info, file = paste0(temp_location, sample_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule, "_iter_", iter, "_", sampling.scheme))


      sim.end <- Sys.time()
   
   iter.time <- round(as.numeric(difftime(sim.end, sim.start, units = "mins")), 1)
   cat(paste0("Finished iteration ", iter, ". \n Took ", iter.time, " minutes"))
   } # end loop over iterations
  
 #------------------------- Objective 1 common settings -------------------------#
 fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
 jags_params = c("Nf", "Nm", "survival") #List the parameters to be estimated
 estimated.parameters <- paste0(jags_params, collapse = ",")
 mating.periodicity <- 1 #number of years between mating for females
 
 survival.prior.mean <- Adult.survival
 survival.prior.cv <- 0.05
 survival.prior.sd <- survival.prior.mean * survival.prior.cv
 model <- "annual.model" #For naming output files
 # 
 source("./02_Estimation.model/functions/Obj1.functions.R") #Changed name of script that includes pairwise comparison and other functions
 jags_file = paste0(jags.model_location, "HS.PO_noLambda_annual_model_validation.txt") 
 
 sampling.scheme <- "sample.ALL.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
 date.of.PopSim <- "17Nov2022"
 HS.only <- "no" #Do we only want to filter HS relationships?
 PO.only <- "no" #Do we only want to filter PO relationships?
 
 estimation.years <- c(n_yrs - 5, n_yrs)
 fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
 PopSim.lambda <- "lambda.1"
 source("./02_Estimation.model/functions/Obj3.functions.R") #Changed name of script that includes pairwise comparison and other functions
 mating.periodicity <- 2 #Overwrite below if triennial breeding
 
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
 