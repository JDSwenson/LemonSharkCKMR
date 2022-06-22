#### LEMON SHARKS: DOVI'S IBS MODEL

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

source("./01_MAIN_scripts/functions/Dovi_IBS_SB_test.assign.conformity.R")
source("./01_MAIN_scripts/functions/pairwise_comparisons_HS.PO_SB.R")

#----------------Set output file locations and labels ------------------------------
temp_location <- "~/R/working_directory/temp_results/"
output.location <- "G://My Drive/Personal_Drive/R/CKMR/Population.simulations/"
parents_prefix <- "parents.breakdown"
sample_prefix <- "sample.info"
pop.size.prefix <- "pop.size"
truth.prefix <- "truth"
population.growth <- "lambda.variable"
sampling.scheme <- "target.YOY"

#-------------------Set and save population simulation settings----------------------------
script_name <- "population_simulation_lambda.variable_target.YOY.R" #Copy name of script here

today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
date.of.simulation <- today

target.YOY <- "yes" #For juvenile samples, do we only want to target YOY for each year of sampling?
HS.only <- "yes"

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
ref.year <- 75 #First year for calculations of survival and lambda

f <- (1-Adult.survival)/(YOY.survival * juvenile.survival^11) # adult fecundity at equilibrium if no age truncation
ff1 <- f/init.prop.female * mating.periodicity/mean(num.mates) # female fecundity per breeding cycle
ff1
ff1 <- ff1*(1-non.conformists) #Change female fecundity per breeding cycle to account for non-conformists
ff1

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
pop.sim.df <- tibble(script_name = script_name,
                     population.growth = population.growth,
                     sampling.scheme = sampling.scheme,
                     date.of.simulation = date.of.simulation,
                     target.YOY = target.YOY,
                     seeds = seeds,
                     years_sampled = length(sample.years),
                     YOY.survival = YOY.survival,
                     juvenile.survival = juvenile.survival,
                     adult.survival = Adult.survival,
                     repro.age = repro.age,
                     max.age = max.age,
                     breeding_periodicity = mating.periodicity,
                     non_conformists = non.conformists,
                     female.fecundity = paste0(ff1, "+/- 0.5"),
                     n.sim.yrs = n_yrs
)

#Save simulation settings in Simulation_log
  # popSim.log <- read_csv("Population_simulation_log.csv") %>%
  #   mutate(female.fecundity = as.character(female.fecundity)) #Read in simulation log
  # tail(popSim.log)
  # simulation.log_updated <- bind_rows(popSim.log, pop.sim.df) #Combine old simulation settings with these
  # write_csv(simulation.log_updated, file = "Population_simulation_log.csv") #Save the updated simulation log

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
 truth.all <- NULL

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
#  sample.info <- readRDS(file = paste0(temp_location, sample_prefix, "_", date.of.simulation, "_", seeds, "_", purpose))


#rseed.pop <- sample(1:100000000, size = 1)
#rseed.pop <- 87625053 #Want to use the same population for all simulations
#set.seed(rseed.pop)

 for(iter in 1:iterations) {
   #  set.seed(rseeds[iter])
   sim.start <- Sys.time()
   #rseed <- sample(1:1000000,1)
   set.seed(rseeds[iter])
   rseed <- rseeds[iter]
   
   ff <- sample(c(ff1-0.5, ff1, ff1+0.5), size = 1)

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

  parents.tibble <- out[[3]] %>% 
    dplyr::filter(year >= 50) %>% #Tibble for each parent for each year to check the distribution later
    mutate(seed = rseed, iteration = iter)
  
  #organize results and calculate summary statistics from the simulation
  source("./01_MAIN_scripts/functions/query_results_PopSim.R")
  
  #-----------------------Collect samples-------------------------
  #Loop over sample sizes stored in sample.vec  
  for(samps in 1:length(sample.vec.prop)){
    
    juv.sample.prop <- sample.vec.prop[samps]
    rent.sample.prop <- juv.sample.prop/5
    
    #Initialize sample dataframes
    sample.df_all.info <- NULL
    sample.df_temp.off = sample.df_temp.rents <- NULL
    
    Juv.samples <- NULL
    Adult.samples <- NULL
    
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
          dplyr::slice_sample(n = sample.size.juvs)  #Sample each year WITHOUT replacement (doesn't affect cross-year sampling since it's in a loop)
        
        #Sample reproductively mature adults only for parent-offspring analysis
        sample.df_temp.rents <- loopy.list[[i]] %>% mutate(capture.year = i) %>% 
          dplyr::filter(age.x >= repro.age) %>% 
          dplyr::slice_sample(n = sample.size.rents)
      }      

      #Combine all
        sample.df_all.info <- rbind(sample.df_all.info, sample.df_temp.off, sample.df_temp.rents)
        Juv.samples <- c(Juv.samples, sample.size.juvs)
        Adult.samples <- c(Adult.samples, sample.size.rents)
    }
    

    sampled.mothers <- unique(sample.df_all.info$mother.x)
    sampled.fathers <- unique(sample.df_all.info$father.x)
    Total.juv.samples <- sum(Juv.samples)
    Total.adult.samples <- if(HS.only == "yes") 0 else sum(Adult.samples)
    
    #Compile results and summary statistics from simulation to compare estimates
    source("01_MAIN_scripts/functions/PopSim_truth.R")
    
    #Save info for samples to examine in more detail
    sample.df_all.info <- sample.df_all.info %>% 
      mutate(sample.size.juvs = Total.juv.samples,
             sample.size.rents = Total.adult.samples,
             sampling.scheme = sampling.scheme,
             iteration = iter,
             seed = rseed,
             sample.prop = juv.sample.prop)
    
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
  
  saveRDS(parents.tibble_all, file = paste0(temp_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_iter_", iter))
  
  # Detailed info on population size
  pop.size.tibble_all <- bind_rows(pop.size.tibble_all, pop.size.tibble)
  
  saveRDS(pop.size.tibble_all, file = paste0(temp_location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_iter_", iter))

  #True values
  truth.all <- bind_rows(truth.all, estimates)
  
  saveRDS(truth.all, file = paste0(temp_location, truth.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_iter_", iter))
  
# Detailed info on samples and parents to examine in more detail
  saveRDS(sample.info, file = paste0(temp_location, sample_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_iter_", iter, "_", sampling.scheme))


      sim.end <- Sys.time()
   
   iter.time <- round(as.numeric(difftime(sim.end, sim.start, units = "mins")), 1)
   cat(paste0("Finished iteration ", iter, ". \n Took ", iter.time, " minutes"))
   } # end loop over iterations
  
   #-----------------------------Save major output files---------------------------------------------
 #Save detailed info about samples from population
 saveRDS(sample.info, file = paste0(output.location, sample_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", sampling.scheme))

 #Save parents tibble
 saveRDS(parents.tibble_all, file = paste0(output.location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", sampling.scheme))
 
 # Detailed info on population size
 saveRDS(pop.size.tibble_all, file = paste0(output.location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", sampling.scheme))
 
 # Truth
 saveRDS(truth.all, file = paste0(output.location, truth.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", sampling.scheme))
 