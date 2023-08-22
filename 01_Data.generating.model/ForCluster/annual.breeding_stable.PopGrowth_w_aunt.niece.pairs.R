#Load packages
library(tidyverse) # safe to ignore conflicts with filter() and lag()
library(MASS)
library(popbio)
#library(mpmtools)  # not at CRAN;
# install.packages("devtools")
#devtools::install_github("BruceKendall/mpmtools")
library(ggpubr)
# The rjags package is just an interface to the JAGS library
# Make sure you have installed JAGS-4.x.y.exe (for any x >=0, y>=0) from
# http://www.sourceforge.net/projects/mcmc-jags/files
# library(rjags)
# library(R2jags)
# library(jagsUI)
library(Rlab)
library(runjags)
library(postpack)
library(coda)

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

#################### Toggles for  population simulation options ####################
b.schedule <- "annual" #annual or biennial or triennial
input.psi <- 1 #Can use any number here; 1 if wanting every individual to have the same schedule
offcycle.breeding <- "no" #Options are "yes" or "no"
input.popGrowth <- "stable" #Options are "stable", "slight increase", "slight decline", or "severe decline"
#sampling.scheme <- "target.YOY"
#sampling.scheme <- "sample.all.juvenile.ages"
#sampling.scheme <- "sample.ALL.ages"


#Specify simulation details based on above input
source("./01_Data.generating.model/functions/specify.simulation.R")

#-------------------Set date and load seeds----------------------------
today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
date.of.simulation <- today

rseeds <- readRDS("rseeds_2022.04.15.rda")
seeds <- "Seeds2022.04.15"
# create some seeds: 
#rseeds <- trunc(1529375630*runif(1000, min=0, max=1 ) )
#rseeds[1]  <- 746160703

#----------------------- DATA-GENERATING MODEL --------------------
# Note on sequencing: Birthdays happen at beginning of each year, followed by mating, then death (i.e. a female who dies in year 10 can still give birth and have surviving pups born in year 10)

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
sample.scheme.vec <- c("target.YOY", "sample.all.juvenile.ages", "sample.ALL.ages")
sample.vec.prop <- c(.5, 1, 1.5, 2)


####-------------- Prep simulation ----------------------
# Moved sampling below so extract different sample sizes from same population
iterations <- 2  # 1 just to look at output     500 #Number of iterations to loop over


# Initialize arrays for saving results
 results <- NULL
 sims.list.1 <- NULL
 sims.list.2 <- NULL
 sims.list.3 <- NULL
 sims.list.4 <- NULL
 sample.info <- NULL
 aunt.unc_nephew.niece_info <- NULL
 parents.tibble_all <- NULL
 pop.size.tibble_all <- NULL
 mom.comps.tibble <- NULL
 dad.comps.tibble <- NULL
 truth.all <- NULL
 decoy_HSPs <- NULL

 sim.samples.1 <- paste0(sample.vec.prop[1], "prop.sampled")
 sim.samples.2 <- paste0(sample.vec.prop[2], "prop.sampled")
 sim.samples.3 <- paste0(sample.vec.prop[3], "prop.sampled")
 sim.samples.4 <- paste0(sample.vec.prop[4], "prop.sampled")

 #
 #---------------------Initialize array from previous checkpoint--------------------------
# Results
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
     
   } #End lambda.extreme conditional statement
   
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

  #-------------------Save simulation output as objects------------------------#
  #List of dataframes for each year of simulation
  loopy.list <- out[[1]] 

  #Population parameters for each year of simulation
  pop.size.tibble <- out[[2]] %>% 
    as_tibble() %>% 
    mutate(seed = rseed, iteration = iter)
  
  #Save dataframe of parents' output
  parents.tibble <- out[[3]] %>% 
    dplyr::filter(year >= 50) %>% #Tibble for each parent for each year to check the distribution later
    mutate(seed = rseed, iteration = iter)
  
  #organize results and calculate summary statistics from the simulation
  source("./01_Data.generating.model/functions/query_results_PopSim.R")
  
  ref.tibble <- NULL #Initialize dataframe for reference years i.e. second earliest birth date of juveniles in the dataset. This changes based on the sampling scheme and iteration (though it's usually the same across iterations), and is used to calculate the truth.
  
  #-----------------------Collect samples-------------------------
  #Loop over sample sizes stored in sample.vec  
  for(samps in 1:length(sample.vec.prop)){

    #Save proportion of population to sample    
    sample.prop <- sample.vec.prop[samps]
    
    #Initialize dataframes that will store output from EACH sampling scheme and ALL sampling schemes
    sample.df_all.info <- NULL #Final df with info from ALL sampling schemes, for ALL years
    sample.df_all.info_temp <- NULL #Df with 

    aunt.unc_niece.nephew_pw.comps.all <- NULL #Final df with info from ALL sampling schemes, for ALL years
    aunt.unc_niece.nephew_pw.comps.all_temp <- NULL
    
    
    #Sample population each year in sample.years and make dataframe of samples with all metadata
    set.seed(rseed)
    
          for(j in 1:length(sample.scheme.vec)){
            
            #Set sampling scheme to current "target.samples"
            target.samples <- sample.scheme.vec[j]
            
            #Sample over every sampling year with the chosen sampling scheme. Output is the combined samples from all years.
            samples.out <- draw.samples(target.samples = target.samples) #Saves a list of output files
            
            sample.df_all.info_temp1 <- samples.out[[1]] %>% #Saves all samples from this sampling scheme from ALL sampling years
              as_tibble()
            
            aunt.unc_niece.nephew_pw.comps.all_temp1 <- samples.out[[2]] %>% #Saves all aunt|uncle / niece|nephew pairs from this sampling scheme from ALL sampling years
              as_tibble()
            
            sample.size <- samples.out[[3]] #Sames object sample.size for later use
            
            
            #Calculate reference year
            ref.tibble <- calculate.ref.year()
            
            sample.df_all.info_temp <- bind_rows(sample.df_all.info_temp, sample.df_all.info_temp1) #Iteratively stores sample info from EACH sampling scheme. Output object contains sample info for ALL samples from ALL sampling schemes for this iteration.
            
            aunt.unc_niece.nephew_pw.comps.all_temp <- bind_rows(aunt.unc_niece.nephew_pw.comps.all_temp, aunt.unc_niece.nephew_pw.comps.all_temp1) #Iteratively stores aunt/niece info from EACH sampling scheme. Output object contains sample info for ALL samples from ALL sampling schemes for this iteration.
    } #End loop over sampling schemes
    

    sampled.mothers <- unique(sample.df_all.info$mother.x)
    sampled.fathers <- unique(sample.df_all.info$father.x)
    
    #Compile results and summary statistics from simulation to compare estimates
    source("./01_Data.generating.model/functions/PopSim_truth.R")
   
    #Contains sample info for ALL samples from ALL sampling schemes for this iteration
    #Rename columns for row bind with sample.info.
    sample.df_all.info <- sample.df_all.info_temp %>% 
      mutate(sample.size.yr = sample.size,
             #sampling.scheme = sampling.scheme,
             iteration = iter,
             seed = rseed,
             sample.prop = sample.prop) %>% 
      mutate(sample.size.total = sample.size.yr * length(sample.years))
    
    #Save info from ALL sampling schemes over ALL sampling years from every iteration. Output object contains sample info from ALL sampling schemes over ALL sampling years over ALL iterations.
    sample.info <- rbind(sample.info, sample.df_all.info) %>% 
      as_tibble()

    #Contains sample info for ALL samples from ALL sampling schemes for this iteration
    #Add columns for row bind with sample.info.
    aunt.unc_niece.nephew_pw.comps.all  <-  aunt.unc_niece.nephew_pw.comps.all %>% 
      bind_rows(aunt.unc_niece.nephew_pw.comps.all_temp) %>% 
      mutate(iteration = iter,
             seed = rseed,
             sample.prop = sample.prop)
    
    #Save info from ALL sampling schemes over ALL sampling years from every iteration. Output object contains aunt|uncle / niece|nephew comparisons and info from ALL sampling schemes over ALL sampling years over ALL iterations.
    aunt.unc_nephew.niece_info <- rbind(aunt.unc_nephew.niece_info, aunt.unc_niece.nephew_pw.comps.all) %>% 
      as_tibble()
  
  } # end loop over sample proportions

  ref.tibble #just view the file
  
  #Adds ref year to each sample depending on how it was collected
  sample.info2 <- sample.info %>% left_join(ref.tibble, by = c("sampling.scheme", "sample.prop", "iteration")) 
  
  #-----------------Save output files iteratively--------------------
  
  sim.samples.1 <- paste0(sample.vec.prop[1], "prop.sampled")
  sim.samples.2 <- paste0(sample.vec.prop[2], "prop.sampled")
  sim.samples.3 <- paste0(sample.vec.prop[3], "prop.sampled")
  sim.samples.4 <- paste0(sample.vec.prop[4], "prop.sampled")
  
  #Save parents tibble from ALL sample proportions
  parents.tibble_all <- bind_rows(parents.tibble_all, parents.tibble)
  
#  saveRDS(parents.tibble_all, file = paste0(temp_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule, "_iter_", iter))
  
  # Detailed info on population size from ALL sample proportions
  pop.size.tibble_all <- bind_rows(pop.size.tibble_all, pop.size.tibble)
  
#  saveRDS(pop.size.tibble_all, file = paste0(temp_location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule, "_iter_", iter))

  #True values from ALL sample proportions
  truth.all <- bind_rows(truth.all, true.values)

  #Aunt|uncle / nice|nephew comparisons and info from ALL sample proportions  
  decoy_HSPs <-  bind_rows(decoy_HSPs, aunt.unc_nephew.niece_info)
  
#  saveRDS(truth.all, file = paste0(temp_location, truth.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule, "_iter_", iter))
  
# Detailed info on samples and parents to examine in more detail
#  saveRDS(sample.info, file = paste0(temp_location, sample_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule, "_iter_", iter, "_", sampling.scheme))


      sim.end <- Sys.time()
   
   iter.time <- round(as.numeric(difftime(sim.end, sim.start, units = "mins")), 1)
   cat(paste0("Finished iteration ", iter, ". \n Took ", iter.time, " minutes"))
   } # end loop over iterations
  
   #-----------------------------Save major output files---------------------------------------------
 #Save detailed info about samples from population
 saveRDS(sample.info2, file = paste0(output.location, sample_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule))

 #Save parents tibble
 saveRDS(parents.tibble_all, file = paste0(output.location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule))
 
 # Detailed info on population size
 saveRDS(pop.size.tibble_all, file = paste0(output.location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule))
 
 # Truth
 saveRDS(truth.all, file = paste0(output.location, truth.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule))
 
 #Charlatan HSPs
 saveRDS(decoy_HSPs, file = paste0(output.location, "aunt.uncle_niece.nephews_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule))