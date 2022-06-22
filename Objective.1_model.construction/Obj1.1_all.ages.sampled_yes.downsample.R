#### LEMON SHARKS: DOVI'S IBS MODEL

#Load packages
library(tidyverse) # safe to ignore conflicts with filter() and lag()
library(MASS)
library(popbio)
library(mpmtools, lib.loc = ".")
library(ggpubr)
library(rjags)
library(R2jags)
library(jagsUI)
library(Rlab)
library(runjags)
library(postpack)
library(coda)

rm(list=ls())

source("./01_MAIN_scripts/functions/Dovi_IBS_SB_test.assign.conformity.R")
source("./01_MAIN_scripts/functions/pairwise_comparisons_HS.PO_SB.R")

#----------------Set output file locations ------------------------------
temp_location <- "~/R/working_directory/temp_results/"
MCMC_location <- "./output/Model.output/"
jags.model_location <- "./output/models/"
results_location <- "./output/Model.results/"

results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
jags.model.prefix <- "CKMR.JAGS_"
parents_prefix <- "parents_breakdown/CKMR_parents.breakdown"
sample.prefix <- "sample_info/CKMR_sample.info"
pop.size.prefix <- "pop_size/CKMR_pop.size"
mom.comps.prefix <- "comparisons/mom.comps"
dad.comps.prefix <- "comparisons/dad.comps"


#-------------------Set simulation settings----------------------------
script_name <- "psi1_0.05non.conform_all.ages.sampled_no.downsample_UpdatedEquation.R" #Copy name of script here
primary_goal <- "Test new equation from Liz" #Why am I running this simulation? Provide details

question1 <- "Does the generalized updated equation from Liz improve estimates of abundance and psi?"
question2 <- "Are survival estimates improved by using the same values as the Leslie matrix?"
question3 <- "Better to use targeted sampling of YOY or sample all age classes"
purpose <- "psi1_0.05non.conform_all.ages.sampled_no.downsample_UpdatedEquation" #For naming output files
today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
date.of.simulation <- today

target.YOY <- "no" #For juvenile samples, do we only want to target YOY for each year of sampling?
down_sample <- "yes" #Do we want to downsample to achieve close to max.HSPs?
max.HSPs <- 150
max.POPs <- 150
HS.only <- "yes" #Do we only want to filter HS relationships?
PO.only <- "no" #Do we only want to filter PO relationships? These two are mutually exclusive; cannot have "yes" for both
fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none" and the full model will run, estimating all parameters. If fixing specific parameters, then list them here, and manually change in the run.JAGS_HS.PO_SB.R script
jags_params = c("Nfb", "psi", "Nm", "surv", "lam")
estimated.parameters <- paste0(jags_params, collapse = ",")

#rseeds <- sample(1:1000000,iterations)
#save(rseeds, file = "rseeds_2022.04.15.rda")

#Save paths and file labels as objects
load("rseeds_2022.04.15.rda")
seeds <- "Seeds2022.04.15"


#----------------------- DATA-GENERATING MODEL --------------------
# Note on sequencing: Births happen at beginning of each year, followed by deaths 
# (i.e. a female who dies in year 10 can still give birth and have surviving pups in year 10)

#--------------------Simulation parameters----------------------------
init.adult.pop.size <- 1000 # CHANGED FROM 3000; Initial adult population size
init.prop.female <- 0.5 # proportion of the initial population size that is female
birth.sex.ratio <- c(0.5,0.5) # The probability that each baby is F:M - has to add up to 1
YOY.survival <- 0.7 # CHANGED FROM 0.8; young of year survival
juvenile.survival <- 0.8 # CHANGED FROM 0.9; juvenile survival
Adult.survival <- 0.825 # CHANGED FROM 0.825; Adult survival
repro.age <- 12 # set age of reproductive maturity
max.age <- maxAge <- 50 #set the maximum age allowed in the simulation
mating.periodicity <- 2 #number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
non.conformists <- 0.05 #proportion of off-year breeders to randomly include off their breeding cycle - want to change this to non.conformists
num.mates <- c(1:3) #vector of potential number of mates per mating
#avg.num.offspring <- 3 # NOT USED? CHANGED FROM 3; set the average number of offspring per mating (from a poisson distribution)

f <- (1-Adult.survival)/(YOY.survival * juvenile.survival^11) # adult fecundity at equilibrium if no age truncation
ff <- f/init.prop.female * mating.periodicity/mean(num.mates) # female fecundity per breeding cycle
ff
ff <- ff*(1-non.conformists) #Change female fecundity per breeding cycle to account for non-conformists
ff

#Stable age distribution
props <- rep(NA, max.age+1)
props[1] <- f
props[2] <- f * YOY.survival
for (y in 3:(repro.age+1)) props[y] <- props[y-1] * juvenile.survival
#props[repro.age+1] <- props[repro.age] * juvenile.survival + Adult.survival

for (y in (repro.age+2):(max.age+1)) props[y] <- props[y-1] * Adult.survival
prop.Adult <- sum(props[(repro.age+1):(max.age+1)])/sum(props)
Nages <- round(props[-1] * init.adult.pop.size) 
init.pop.size <- sum(Nages) # all ages except YOYs

#Set length of simulation and estimation year
burn.in <- 40 # number of years to use as simulation burn in period
Num.years <- 50 # The number of years to run in the simulation beyond the burn in
n_yrs <- burn.in + Num.years #Total number of simulation years
estimation.year <- n_yrs - 5 # Set year of estimation

#rseeds <- sample(1:1000000,iterations)
#load("rseeds_2022.03.23.rda")


#-----------------Leslie Matrix parameters--------------------
#To set a prior on lambda, we will run a Leslie matrix, assuming the following information for survival and fecundity
leslie.survival <- Adult.survival
leslie.fecundity <- ff/mating.periodicity
surv.cv <- 0.1 #What is the CV on survival?
fec.cv <- 0.1 #What is the CV on fecundity?
corr.vec <- c(0, -0.25, -0.5)
n.draws <- 100 #Number of draws from a multivariate normal distribution
#Run leslie matrix to generate priors for lambda and survival
source("./01_MAIN_scripts/functions/Leslie_matrix_source.R")

#Check values from Leslie matrix
lambda.sd
lambda.df
mean(lambda.df$mean.survival)
min(lambda.df$mean.survival)
max(lambda.df$mean.survival)
mean(lambda.df$mean.fecundity)
min(lambda.df$mean.fecundity)
max(lambda.df$mean.fecundity)


#--------------------- Sampling parameters ---------------------
sample.years <- c(n_yrs - c(3:0)) #For two years of sampling
#sample.years <- n_yrs #One year of sampling
#sample.size <- 300 #sample size per year
#sample.vec.juvs <- c(50, 100, 150, 200) #vector to sample over per year
#sample.vec.adults <- c(sample.vec.juvs/5)
#sample.vec.total <- sample.vec.juvs + sample.vec.adults
sample.vec.prop <- c(.5, 1, 1.5, 2)

#----------------------- MCMC parameters ----------------------#
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains


####---------------Update and save simulation log-------------------####
simulation.df <- tibble(script_name = script_name,
                        primary_goal = primary_goal,
                        question1 = question1,
                        question2 = question2,
                        question3 = question3,
                        purpose = purpose,
                        date.of.simulation = date.of.simulation,
                        target.YOY = target.YOY,
                        down_sample = down_sample,
                        max.HSPs = max.HSPs,
                        max.POPs = max.POPs,
                        HS.only = HS.only,
                        PO.only = PO.only,
                        fixed.parameters = fixed.parameters,
                        estimated.parameters = estimated.parameters,
                        seeds = seeds,
                        thinning_rate = nt,
                        posterior_samples = ni,
                        burn_in = nb,
                        years_sampled = length(sample.years),
                        breeding_periodicity = mating.periodicity,
                        non_conformists = non.conformists,
                        survival_cv = surv.cv,
                        fecundity_cv = fec.cv
)

#Save simulation settings in Simulation_log
   # simulation.log <- read_csv("Simulation_log.csv")
   #  tail(simulation.log)
   #  simulation.log_updated <- bind_rows(simulation.log, simulation.df) #Combine old simulation settings with these
   #  write_csv(simulation.log_updated, file = "Simulation_log.csv") #Save the updated simulation log

####-------------- Start simulation loop ----------------------
# Moved sampling below so extract different sample sizes from same population
iterations <- 100 #Number of iterations to loop over


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

 sim.samples.1 <- paste0(sample.vec.prop[1], "prop.sampled")
 sim.samples.2 <- paste0(sample.vec.prop[2], "prop.sampled")
 sim.samples.3 <- paste0(sample.vec.prop[3], "prop.sampled")
 sim.samples.4 <- paste0(sample.vec.prop[4], "prop.sampled")

 ####Initialize array from previous checkpoint
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
#  sample.info <- readRDS(file = paste0(temp_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))


#rseed.pop <- sample(1:100000000, size = 1)
#rseed.pop <- 87625053 #Want to use the same population for all simulations
#set.seed(rseed.pop)

 for(iter in 1:iterations) {
   #  set.seed(rseeds[iter])
   sim.start <- Sys.time()
   #rseed <- sample(1:1000000,1)
   set.seed(rseeds[iter])
   rseed <- rseeds[iter]

   source("./01_MAIN_scripts/functions/Dovi_IBS_SB_test.assign.conformity.R")
   
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

  pop.size.tibble_all <- bind_rows(pop.size.tibble_all, pop.size.tibble)
  
  parents.tibble <- out[[3]] %>% 
    dplyr::filter(year >= 50)#Tibble for each parent for each year to check the distribution later
  
  parents.tibble_all <- bind_rows(parents.tibble_all, parents.tibble)
  
  #Save parents tibble
  saveRDS(parents.tibble_all, file = paste0(temp_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose, "_iter_", iter))
  
  # Detailed info on population size
  saveRDS(pop.size.tibble_all, file = paste0(temp_location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", purpose, "_iter_", iter))

  
  # saveRDS(loopy.list, file = "../loopy.list")
  # saveRDS(parents.tibble, file = "../parents.tibble")
  # saveRDS(pop.size.tibble, file = "../pop.size.tibble")
  # loopy.list <- readRDS("loopy.list")
  # parents.tibble <- readRDS(parents.tibble, file = paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose))
  #pop.size.tibble <- readRDS(pop.size.tibble, file = paste0(results_location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))
  # 
  #organize results and calculate summary statistics from the simulation
  source("./01_MAIN_scripts/functions/query_results_SB.R")
  
  #-----------------------Collect samples-------------------------
  #Loop over sample sizes stored in sample.vec  
  for(samps in 1:length(sample.vec.prop)){
    
    juv.sample.prop <- sample.vec.prop[samps]
    rent.sample.prop <- juv.sample.prop/5
    
    #Initialize sample dataframes
    sample.df_all.info <- NULL
    sample.df_temp.off = sample.df_temp.rents <- NULL
    
    #Sample population each year in sample.years and make dataframe of samples with all metadata
    set.seed(rseed)
    
        for(i in sample.years){
      
      pop.size.yr <- pop.size.tibble %>% dplyr::filter(year == i) %>% 
        mutate(Total.juvenile.pop = population_size - Total.adult.pop)
      
      #Set number of juvenile samples to a specific proportion of the population  
      sample.size.juvs <- pop.size.yr %>% 
        mutate(sample.size = round(Total.juvenile.pop*(juv.sample.prop/100)), 0) %>% 
        pull(sample.size)
      
      #Supplement with parents
      sample.size.rents <- pop.size.yr %>% 
        mutate(sample.size = round(Total.adult.pop*(rent.sample.prop/100)), 0) %>% 
        pull(sample.size)
      
      total.prop.sampled <- round((sample.size.juvs + sample.size.rents)/(pop.size.yr$population_size), 2)
      
      if(target.YOY == "yes"){ #If targeting YOY for juvenile samples
      #Sample YOY only for half-sib analysis
      sample.df_temp.off <- loopy.list[[i]] %>% mutate(capture.year = i) %>% 
        dplyr::filter(age.x == 0) %>% 
        dplyr::slice_sample(n = sample.size.juvs) # #Sample each year WITHOUT replacement (doesn't affect cross-year sampling since it's in a loop)

      #Sample reproductively mature adults only for parent-offspring analysis
      sample.df_temp.rents <- loopy.list[[i]] %>% mutate(capture.year = i) %>% 
        dplyr::filter(age.x >= repro.age) %>% 
        dplyr::slice_sample(n = sample.size.rents)
      } else { #If indiscriminately sampling juveniles
        
        sample.df_temp.off <- loopy.list[[i]] %>% mutate(capture.year = i) %>% 
          dplyr::filter(age.x < repro.age) %>% 
          dplyr::slice_sample(n = sample.size.juvs) # #Sample each year WITHOUT replacement (doesn't affect cross-year sampling since it's in a loop)
        
        #Sample reproductively mature adults only for parent-offspring analysis
        sample.df_temp.rents <- loopy.list[[i]] %>% mutate(capture.year = i) %>% 
          dplyr::filter(age.x >= repro.age) %>% 
          dplyr::slice_sample(n = sample.size.rents)
      }      

      #Combine all
        sample.df_all.info <- rbind(sample.df_all.info, sample.df_temp.off, sample.df_temp.rents)
    }
    

        noDups.list <- split.dups(sample.df_all.info)
        first.capture <- noDups.list[[1]]
        later.capture <- noDups.list[[2]]
    
 
    source("./01_MAIN_scripts/functions/pairwise_comparisons_HS.PO_SB.R") #Already loaded above; here for troubleshooting
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
    psi.truth <- calc.psi(loopy.list, mom_comps.all)

    
    #source("./01_MAIN_scripts/functions/pairwise_comparisons_HS.PO_SB.R") #Already loaded above; here for troubleshooting
    ####-----------------------------Downsample if more than max.HSPs------------------------------------####
    if(down_sample == "yes"){
      set.seed(rseed)
      
    down.samples <- downsample(mom_comps.all, dad_comps.all, HS.samps.df, PO.samps.list) #calculates approximate number of samples needed to achieve max.HSPs for filtering below
    HS.samps.df.down <- down.samples[[1]]
    PO.samps.list.down <- down.samples[[2]]
    
    # (HS_target.samples <- down.samples[[1]])
    # HS_props <- down.samples[[2]]
    # HS_pos.mom.comps <- down.samples[[3]] #Number of positives for m
    # HS_pos.dad.comps <- down.samples[[4]]
    # PO_target.samples.per.yr <- down.samples[[5]]
    # PO_pos.mom.comps <- down.samples[[6]]
    # PO_pos.dad.comps <- down.samples[[7]]
    # PO.props.vec <- down.samples[[8]]

    #down-sample to achieve a more reasonable number of kin (if necessary)
      
      #Make a new pairwise comparison matrix
    pairwise.out2 <- build.pairwise(filtered.samples.PO.list = PO.samps.list.down, filtered.samples.HS.df = HS.samps.df.down)

    #Save output as different dataframes; includes both HS and PO relationships (but can filter below)
    mom_comps.all <- pairwise.out2[[1]]
    dad_comps.all <- pairwise.out2[[2]]

    head(mom_comps.all)
    head(dad_comps.all)

    mom_comps.all %>% group_by(type) %>%
      summarize(sum(yes))

    dad_comps.all %>% group_by(type) %>%
      summarize(sum(yes))
    } #End down sample
    
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
    
      
    # ####------------------------ Fit CKMR model ----------------####
    #Define JAGS data and model, and run the MCMC engine
      set.seed(rseed)
    source("01_MAIN_scripts/functions/run.JAGS_HS.only_SB_UPDATED.R")

    #Calculate expectations
    Exp <- calc.Exp(mom_comps.all, dad_comps.all)
    mom.Exp.HS <- Exp[[1]]
    mom.Exp.PO <- Exp[[2]]
    dad.Exp.HS <- Exp[[3]]
    dad.Exp.PO <- Exp[[4]]
    
    sampled.mothers <- unique(sample.df_all.info$mother.x)
    sampled.fathers <- unique(sample.df_all.info$father.x)
    
    #Compile results and summary statistics from simulation to compare estimates
    source("01_MAIN_scripts/functions/compile.results_HS.only_SB.R")
    
    #-----------------Loop end-----------------------------
    #Bind results from previous iterations with current iteration
    (results.temp <- cbind(estimates, metrics) %>% 
       mutate(purpose = purpose))
    results <- rbind(results, results.temp)
    
    
    #Save info for samples to examine in more detail
    sample.df_all.info <- sample.df_all.info %>% mutate(iteration = iter, 
                                                        sample.size.juvs = Total.juv.samples, 
                                                        sample.size.rents = Total.adult.samples,
                                                        seed = rseed)
    sample.info <- rbind(sample.info, sample.df_all.info)
  
    #Save mom and dad pairwise comparison dataframes
    mom_comps.all <- mom_comps.all %>% mutate(iteration = iter,
                                              sample.prop.juvs = Juv_sample_prop,
                                              sample.prop.rents = Adult_sample_prop,
                                              sample.size.juvs = Total.juv.samples, 
                                              sample.size.rents = Total.adult.samples,
                                              seed = rseed)
    mom.comps.tibble <- rbind(mom.comps.tibble, mom_comps.all)
    
    dad_comps.all <- dad_comps.all %>% mutate(iteration = iter,
                                              sample.prop.juvs = Juv_sample_prop,
                                              sample.prop.rents = Adult_sample_prop,
                                              sample.size.juvs = Total.juv.samples, 
                                              sample.size.rents = Total.adult.samples,
                                              seed = rseed)
    dad.comps.tibble <- rbind(dad.comps.tibble, dad_comps.all)
    
  } # end loop over sample sizes
  }
    
  #-----------------Save output files iteratively--------------------
  
  sim.samples.1 <- paste0(sample.vec.prop[1], "prop.sampled")
  sim.samples.2 <- paste0(sample.vec.prop[2], "prop.sampled")
  sim.samples.3 <- paste0(sample.vec.prop[3], "prop.sampled")
  sim.samples.4 <- paste0(sample.vec.prop[4], "prop.sampled")
  
#Results
    write.table(results, file = paste0(temp_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose, "_iter_", iter, ".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
# 
#    #Model output for diagnostics
     saveRDS(sims.list.1, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose))
# 
    saveRDS(sims.list.2, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose))
# 
    saveRDS(sims.list.3, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose))
#    
    saveRDS(sims.list.4, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.4, "_", MCMC.settings, "_", purpose))
# 
# # Detailed info on samples and parents to examine in more detail
    saveRDS(sample.info, file = paste0(temp_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))
# 
#    #Save pairwise comparisons matrices
    saveRDS(mom.comps.tibble, file = paste0(temp_location, mom.comps.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))
#    
    saveRDS(dad.comps.tibble, file = paste0(temp_location, dad.comps.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))

      sim.end <- Sys.time()
   
   iter.time <- round(as.numeric(difftime(sim.end, sim.start, units = "mins")), 1)
   cat(paste0("Finished iteration ", iter, ". \n Took ", iter.time, " minutes"))
   } # end loop over iterations
  
  
  
########## Save and check results ##########
#Calculate relative bias for all estimates
#If using breeding individuals for Nf truth
 # results2 <- results %>%
 #   mutate(relative_bias = ifelse(parameter == "Nfb", round(((Q50 - breed.truth)/breed.truth)*100, 1),
 #                                 round(((Q50 - all.truth)/all.truth)*100, 1))) %>% 
 #   mutate(in_interval = ifelse(parameter == "Nfb", 
 #                               ifelse(HPD2.5 < breed.truth & breed.truth < HPD97.5, "Y", "N"),
 #                               ifelse(HPD2.5 < all.truth & all.truth < HPD97.5, "Y", "N"))) %>%
 #   mutate(total_samples = total_juvenile_samples + total_adult_samples) %>% 
 #   as_tibble()

 
   #If using all individuals for Nf truth, instead of breeders
   results2 <- results %>%
     mutate(relative_bias = round(((Q50 - all.truth)/all.truth)*100, 1)) %>% #Can change truth to breed.truth if looking for number of active breeders
     mutate(in_interval = ifelse(HPD2.5 < all.truth & all.truth < HPD97.5, "Y", "N")) %>%
     mutate(total_samples = tota_juvenile_samples + total_adult_samples) %>%
     as_tibble()
   
#Within HPD interval?
results2 %>% group_by(prop_sampled_juvs, parameter, purpose) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100)

#Median relative bias by sample size
 results2 %>% group_by(prop_sampled_juvs, parameter, purpose) %>% 
   dplyr::summarize(median.bias = median(relative_bias), n = n()) %>% 
   dplyr::arrange(desc(median.bias))

 # #View results
 # results2 %>% dplyr::select(parameter, 
 #                            Q2.5,
 #                            Q50,
 #                            Q97.5,
 #                            mean,
 #                            sd,
 #                            prop_sampled_juvs,
 #                            prop_sampled_adults,
 #                            total_samples,
 #                            purpose,
 #                            HPD2.5,
 #                            HPD97.5,
 #                            truth,
 #                            relative_bias,
 #                            in_interval,
 #                            Exp_POPs,
 #                            POPs_detected,
 #                            Exp_HSPs,
 #                            HSPs_detected,
 #                            iteration) %>% 
 #   #dplyr::arrange(parameter, iteration, total_samples) %>% 
 #   dplyr::arrange(desc(relative_bias)) %>% 
 #   View()
   
 
 #Mean number of parents detected
 #Median relative bias by sample size
 # results2 %>% group_by(total_samples, parameter) %>% 
 #   dplyr::summarize(mean = mean(parents_detected), n = n())
 
 
 #-----------------------------Save major output files---------------------------------------------
 #Home computer: Dell Precision
 
 #Save model estimates
write.table(results2, file = paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose, ".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
 
 #Save draws from posterior for model diagnostics 
 saveRDS(sims.list.1, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose)) #Sample size 1
 
 saveRDS(sims.list.2, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose)) #Sample size 2
 
 saveRDS(sims.list.3, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose)) #Sample size 3
 
 saveRDS(sims.list.4, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.4, "_", MCMC.settings, "_", purpose)) #Sample size 4
 
 #Save detailed info about samples from population
 saveRDS(sample.info, file = paste0(results_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))
 
 #Save final pairwise comparison matrices
 saveRDS(mom.comps.tibble, file = paste0(results_location, mom.comps.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))
 
 saveRDS(dad.comps.tibble, file = paste0(results_location, dad.comps.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))
 
 #Save parents tibble
 saveRDS(parents.tibble_all, file = paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose))
 
 # Detailed info on population size
 saveRDS(pop.size.tibble_all, file = paste0(results_location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))
 
#To read in RDS file
#pp <- readRDS("~/R/working_directory/temp_results/neutralGrowth_estSurv_iteration_5_samplesize_800")


#MGHPCC
#write.table(results2, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/Dovi_lambdaModel_06_22.2021_neutralPopGrowth.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)


#-------------Quick viz of results--------------#
#Box plot of relative bias
ggplot(data=results2, aes(x=factor(prop_sampled_juvs))) +
  geom_boxplot(aes(y=relative_bias, fill=parameter)) +
  ylim(-100, 100) +
  geom_hline(yintercept=0, col="black", size=1.25) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  labs(x="Sample size", y="Relative bias", title="Relative Bias by sample size") +
  scale_fill_brewer(palette="Set2") +
  font("title", size = 10, face = "bold")



####################### End ##############################