#Load packages
#library(tidyverse) # safe to ignore conflicts with filter() and lag()
library(dplyr)
library(tidyr)
library(stringr)
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
offcycle.breeding <- "yes" #Options are "yes" or "no"
input.popGrowth <- "slight decline" #Options are "stable", "slight increase", "slight decline", or "severe decline"
#sampling.scheme <- "target.YOY"
#sampling.scheme <- "sample.all.juvenile.ages"
#sampling.scheme <- "sample.ALL.ages"



#################### Breeding schedule ######################
if(b.schedule == "annual"){
#------------------------------ Annual ------------------------------
 breeding.schedule <- "annual.breeding"
 mating.periodicity <- 1 #number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
 non.conformists <- 0
 ff <- f/init.prop.female * mating.periodicity/mean(num.mates) # female fecundity per breeding cycle
 ff
 print("breeding schedule is annual")
 
} else if(b.schedule == "biennial"){
  
#------------------------------ Biennial ------------------------------
mating.periodicity <- 2 #number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
print("breeding schedule is biennial")

if(input.psi == 1){
#============================== psi 1 ==============================
  breeding.schedule <- "biennial.breeding_psi1"
  non.conformists <- 0
  
  print("psi equals 1")
  
  if(offcycle.breeding == "yes"){
    
    percent.breed_off.cycle <- 0.1 #Percent of off-cycle mothers that will breed each year
    percent.skip_on.cycle <- 0.1 #Percent of on-cycle mothers that will skip breeding each year
    
    print(paste0("There are ", percent.breed_off.cycle*100, "% of females that will breed off-cycle and ", percent.skip_on.cycle*100, "% of females that will fail when on-cycle."))
    
  }
} else if(input.psi != 1){
  
  breeding.schedule <- paste0("biennial.breeding_psi", input.psi)
  non.conformists <- 1 - input.psi #proportion of off-year breeders to randomly include off their breeding cycle
  
  print(paste0("psi equals ", 1 - input.psi))
}

} else if(b.schedule == "triennial"){
  
  mating.periodicity <- 3 #number of years between mating; assigned to an individual and sticks with them through their life.
  breeding.schedule <- "triennial.breeding_psi1"
  non.conformists <- 0
  
  print("breeding schedule is triennial")
}

if(b.schedule == "biennial" | b.schedule == "triennial"){
  # Adjust fecundity ==============================================================
  ## for effective number of breeders each year, mating cycle, number of mates ====
  #(from liz) ====
  psi <- 1-non.conformists
  ff <- mating.periodicity/(mating.periodicity-psi*mating.periodicity+psi)*f/init.prop.female/mean(num.mates)
  ff
}

#################### Population growth ####################
if(input.popGrowth == "stable"){
#------------------------------Stable------------------------------
population.growth <- "lambda.1"

print("population growth is stable")

} else if(input.popGrowth == "slight increase"){

#------------------------------Slight increase------------------------------
 population.growth <- "lambda.slight.increase"
# ff.shift <- ff+0.5 #Increase fecundity to slightly increase population growth - only works for annual breeding
 
 (Ma <- -log(Adult.survival)) #Mortality of adults
 (Mj <- -log(juvenile.survival)) #Mortality of juveniles
 (Sa <- exp(-Ma)) #Survival of adults
 (Sj <- exp(-Mj)) #Survival of juveniles
 
 Mx <- -0.01 #Extra mortality (make negative so it reduces mortality)
 (Sa <- exp(-Ma - Mx)) #Survival of adults
 (Sj <- exp(-Mj - Mx)) #Survival of juveniles
 
 Adult.survival <- Sa
 juvenile.survival <- Sj
 
 cat(paste0("Adult survival is ", round(Adult.survival, 0), ";\nJuvenile survival is ", round(juvenile.survival,0), ";\nPopulation will increase in size"))
 
} else if(input.popGrowth == "slight decline"){

#------------------------------Slight decline------------------------------
 population.growth <- "lambda.slight.decrease"
# ff.shift <- ff-0.5 #Decrease fecundity to slightly decrease population growth - only works for annual breeding
 (Ma <- -log(Adult.survival)) #Mortality of adults
 (Mj <- -log(juvenile.survival)) #Mortality of juveniles
 (Sa <- exp(-Ma)) #Survival of adults
 (Sj <- exp(-Mj)) #Survival of juveniles
 
 Mx <- 0.01 #Extra mortality
 (Sa <- exp(-Ma - Mx)) #Survival of adults
 (Sj <- exp(-Mj - Mx)) #Survival of juveniles
 
 Adult.survival <- Sa
 juvenile.survival <- Sj
 
 cat(paste0("Adult survival is ", round(Adult.survival, 0), ";\nJuvenile survival is ", round(juvenile.survival,0),";\nPopulation will decrease in size"))
 
} else if(input.popGrowth == "severe decline"){

#------------------------------Substantial decrease------------------------------
population.growth <- "lambda.extreme" #Mortality will change later inside the loop

cat(paste0("Adult survival is ", round(Adult.survival,0), ";\nJuvenile survival is ", round(juvenile.survival,0),";\nPopulation will decline rapidly"))
}


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
iterations <- 500  # 1 just to look at output     500 #Number of iterations to loop over


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
     
     Adult.survival <- Sa
     juvenile.survival <- Sj
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
          
          for(j in 1:length(sample.scheme.vec)){
            
            #Set sampling scheme to 
            target.samples <- sample.scheme.vec[j]
            
            if(target.samples == "target.YOY"){ #If targeting YOY for juvenile samples
            
              #Sample YOY only for half-sib analysis
              sample.df_temp <- loopy.list[[i]] %>% mutate(capture.year = i,
                                                           sampling.scheme = target.samples) %>%
                dplyr::filter(age.x == 0) %>%
                dplyr::slice_sample(n = sample.size) # #Sample each year WITHOUT replacement
            
            
            } else if(target.samples == "sample.all.juvenile.ages"){ #If sampling juveniles
              sample.df_temp <- loopy.list[[i]] %>% dplyr::filter(age.x < repro.age & age.x > 0) %>% 
                mutate(capture.year = i,
                       sampling.scheme = target.samples) %>%
                dplyr::slice_sample(n = sample.size) #Sample each year WITHOUT replacement (doesn't affect cross-year sampling since it's in a loop)
              
              
              } else if(target.samples == "sample.ALL.ages"){
              
              sample.df_temp <- loopy.list[[i]] %>% mutate(capture.year = i,
                                                           sampling.scheme = target.samples) %>%
                dplyr::slice_sample(n = sample.size) #Sample each year WITHOUT replacement (doesn't affect cross-year sampling since it's in a loop)
        
      }      

      #Combine samples from all years
        sample.df_all.info <- rbind(sample.df_all.info, sample.df_temp)
          }
    }
    

    sampled.mothers <- unique(sample.df_all.info$mother.x)
    sampled.fathers <- unique(sample.df_all.info$father.x)
    
    #Compile results and summary statistics from simulation to compare estimates
    source("./01_Data.generating.model/functions/PopSim_truth.R")
    
    #Save info for samples to examine in more detail
    sample.df_all.info <- sample.df_all.info %>% 
      mutate(sample.size.yr = sample.size,
             #sampling.scheme = sampling.scheme,
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
  
   #-----------------------------Save major output files---------------------------------------------
 #Save detailed info about samples from population
 saveRDS(sample.info, file = paste0(output.location, sample_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule))
 
 #Save parents tibble
 saveRDS(parents.tibble_all, file = paste0(output.location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule))
 
 # Detailed info on population size
 saveRDS(pop.size.tibble_all, file = paste0(output.location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule))
 
 # Truth
 saveRDS(truth.all, file = paste0(output.location, truth.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule))