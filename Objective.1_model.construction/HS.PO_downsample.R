#### LEMON SHARKS: DOVI'S IBS MODEL

#Load packages
library(tidyverse) # safe to ignore conflicts with filter() and lag()
#library(popbio) #not needed
#library(mpmtools) #not needed
library(ggpubr)
library(rjags)
library(R2jags)
library(jagsUI)
library(Rlab)
library(runjags)
library(postpack)
library(coda)

rm(list=ls())

source("./01_MAIN_scripts/functions/Dovi_IBS.R")
source("./01_MAIN_scripts/functions/pairwise_comparisons_HS.PO_downsample.R")
source("./01_MAIN_scripts/functions/remove_dups.R")

#----------------Set output file locations ------------------------------
#Check results from model diagnostics
today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
date.of.simulation <- today

#rseeds <- sample(1:1000000,iterations)
#save(rseeds, file = "rseeds_2022.03.23.rda")

#Save paths and file labels as objects
load("rseeds_2022.03.23.rda")
seeds <- "Seeds2022.03.23"
purpose <- "HS.PO_downsample"
temp_location <- "~/R/working_directory/temp_results/"
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Model.validation/Model.output/"
jags.model_location <- "G://My Drive/Personal_Drive/R/CKMR/Model.validation/models/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Model.validation/Model.results/"

results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
jags.model.prefix <- "HS.PO_"
parents_prefix <- "CKMR_parents.breakdown"
sample.prefix <- "CKMR_sample.info"
pop.size.prefix <- "CKMR_pop.size"

#----------------------- DATA-GENERATING MODEL --------------------
# Note on sequencing: Births happen at beginning of each year, followed by deaths 
# (i.e. a female who dies in year 10 can still give birth and have surviving pups in year 10)

#--------------------Simulation parameters----------------------------

init.adult.pop.size <- 1000 # CHANGED FROM 3000; Initial adult population size
init.prop.female <- 0.5 # proportion of the initial population size that is female
birth.sex.ratio <- c(0.5,0.5) # The probability that each baby is F:M - has to add up to 1
YOY.survival <- 0.7 # CHANGED FROM 0.8; young of year survival
juvenile.survival <- 0.8 # CHANGED FROM 0.9; juvenile survival
Adult.survival <- 0.825 # CHANGED FROM 0.9; Adult survival
repro.age <- 12 # set age of reproductive maturity
max.age <- maxAge <- 50 #set the maximum age allowed in the simulation
mating.periodicity <- 1 #number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
num.mates <- c(1:3) #vector of potential number of mates per mating
#avg.num.offspring <- 3 # NOT USED? CHANGED FROM 3; set the average number of offspring per mating (from a poisson distribution)

f <- (1-Adult.survival)/(YOY.survival * juvenile.survival^11) # adult fecundity at equilibrium if no age truncation
ff <- f/init.prop.female * mating.periodicity/mean(num.mates) # female fecundity per breeding cycle
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

#----------------------- MCMC parameters ----------------------
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains

#--------------------- Sampling parameters ---------------------
sample.years <- c(n_yrs - c(3:0)) #For two years of sampling
#sample.years <- n_yrs #One year of sampling
#sample.size <- 300 #sample size per year
sample.vec <- c(50, 150, 250) #vector to sample over per year


####-------------- Start simulation loop ----------------------
# Moved sampling below so extract different sample sizes from same population
iterations <- 100#Number of iterations to loop over


# Initialize arrays for saving results
 results <- NULL
 sample.info <- NULL
 rents.info <- NULL
 sims.list.1 <- NULL
 sims.list.2 <- NULL
 sims.list.3 <- NULL
 pop.size.tibble <- NULL

sim.samples.1 <- paste0(sample.vec[1]*length(sample.years), ".samples")
sim.samples.2 <- paste0(sample.vec[2]*length(sample.years), ".samples")
sim.samples.3 <- paste0(sample.vec[3]*length(sample.years), ".samples")

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


for(iter in 1:iterations) {
#  set.seed(rseeds[iter])
  sim.start <- Sys.time()
  #rseed <- sample(1:1000000,1)
  set.seed(rseeds[iter])
  rseed <- rseeds[iter]

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
  pop.size.temp <- out[[2]] %>%  #population parameters for each year of simulation
    as_tibble() %>% 
    mutate(iter = iter)

  parents.tibble <- out[[3]] %>% #Tibble for each parent for each year to check the distribution later
    mutate(iteration = iter) %>% 
    dplyr::filter(year > 50) #Just make the dataframe smaller
  
  #organize results and calculate summary statistics from the simulation
  source("./01_MAIN_scripts/functions/query_results.R")

  #-----------------------Collect samples-------------------------
  #Loop over sample sizes stored in sample.vec  
  for(samps in 1:length(sample.vec)){
        sample.size <- sample.vec[samps] #Specify sample size
    
    #Initialize sample dataframes
    sample.df_all.info <- NULL
    sample.df_temp <- NULL
    
    #Sample population each year in sample.years and make dataframe of samples with all metadata
    for(i in sample.years){
      sample.df_temp <- loopy.list[[i]] %>% mutate(capture.year = i) %>% 
        dplyr::slice_sample(n = sample.size) #%>% #Sample each year WITHOUT replacement (doesn't affect cross-year sampling since it's in a loop)
      #In case wanting to filter for one offspring per parent per cohort  
      #dplyr::distinct(mother.x, birth.year, .keep_all = TRUE) %>% 
      #dplyr::distinct(father.x, birth.year, .keep_all = TRUE)
      sample.df_all.info <- rbind(sample.df_all.info, sample.df_temp)
    }
    
    
    #Keep just one instance of each individual (to avoid self-recapture) and sort by birth year so when we make the pairwise comparison matrix, Ind_1 is always older than Ind_2 (bc it comes first)
    # sample.df_all.info <- sample.df_all.info %>% distinct(indv.name, .keep_all = TRUE) %>% 
    #   dplyr::arrange(birth.year) %>% 
    #   as_tibble()

    source("./01_MAIN_scripts/functions/remove_dups.R")

        noDups.list <- split.dups(sample.df_all.info)
        first.capture <- noDups.list[[1]]
        later.capture <- noDups.list[[2]]
    
     sampled.mothers <- unique(sample.df_all.info$mother.x)
     sampled.fathers <- unique(sample.df_all.info$father.x)
     
    source("./01_MAIN_scripts/functions/pairwise_comparisons_HS.PO_downsample.R")
    #Remove full sibs
    filter1.out <- filter.samples(later.capture) #Filter for full sibs
    PO.samps.list <- filter1.out[[1]] #Output is a list where each list element corresponds to the offspring birth year and contains the potential parents and offspring for that year.
    HS.samps.df <- filter1.out[[2]] #Output is just the dataframe of samples but filtered for full siblings
    
    
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
    
    #Downsample if more than max.HSPs
    max.HSPs <- 150
    down.samples <- downsample(mom_comps.all, dad_comps.all, HS.samps.df) #calculates approximate number of samples needed to achieve max.HSPs for filtering below
    (target.samples <- down.samples[[1]])
    props <- down.samples[[2]]
    pos.mom.comps <- down.samples[[3]] #Number of positives for m
    pos.dad.comps <- down.samples[[4]]
    
    #downs-ample to achieve a more reasonable number of kin (if necessary)
    if(pos.mom.comps + pos.dad.comps > max.HSPs){
      HS.samps.df.down <- HS.samps.df %>% slice_sample(n = target.samples, weight_by = props)
        #Make a new pairwise comparison matrix
    pairwise.out2 <- build.pairwise(filtered.samples.PO.list = PO.samps.list, filtered.samples.HS.df = HS.samps.df.down)
    
    #Save output as different dataframes; includes both HS and PO relationships (but can filter below)
    mom_comps.all <- pairwise.out2[[1]]
    dad_comps.all <- pairwise.out2[[2]]
    
    head(mom_comps.all)
    head(dad_comps.all)
    
    mom_comps.all %>% group_by(type) %>% 
      summarize(sum(yes))
    
    dad_comps.all %>% group_by(type) %>% 
      summarize(sum(yes))
    }
    ####------------------------ Fit CKMR model ----------------####
    #Uncomment below to only run the HS model
    #mom_comps.all <- mom_comps.all %>% filter(type == "HS")
    #dad_comps.all <- dad_comps.all %>% filter(type == "HS")

    #Uncomment below to only run the PO model
    #mom_comps.all <- mom_comps.all %>% filter(type == "PO")
    #dad_comps.all <- dad_comps.all %>% filter(type == "PO")

    #Define JAGS data and model, and run the MCMC engine   
    source("01_MAIN_scripts/functions/run.JAGS_HS.PO.R")
    
    #Compile results and summary statistics from simulation to compare estimates
    source("01_MAIN_scripts/functions/compile.results_HS.PO.R")
    
    #-----------------Loop end-----------------------------
    #Bind results from previous iterations with current iteration
    (results.temp <- cbind(estimates, metrics))
    results <- rbind(results, results.temp)
    
    
    #Save info for samples to examine in more detail
    sample.df_all.info <- sample.df_all.info %>% mutate(iteration = iter, sample.size = sample.size, seed = rseed)
    sample.info <- rbind(sample.info, sample.df_all.info)
    
    #Save parents info
    rents.info <- rbind(rents.info, parents.tibble)
    
    #Save survival info
    pop.size.temp <- pop.size.temp  %>% 
      dplyr::filter(year >= 50)

  } # end loop over sample sizes
  
  
  #-----------------Save output files iteratively--------------------
  #in case R crashes or computer shuts down
  pop.size.tibble <- rbind(pop.size.tibble, pop.size.temp)
  
  sim.samples.1 <- paste0(sample.vec[1]*length(sample.years), ".samples")
  sim.samples.2 <- paste0(sample.vec[2]*length(sample.years), ".samples")
  sim.samples.3 <- paste0(sample.vec[3]*length(sample.years), ".samples")
  
#Results
   write.table(results, file = paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose, "_iter_", iter, ".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

   #Model output for diagnostics
    saveRDS(sims.list.1, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose))

   saveRDS(sims.list.2, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose))

   saveRDS(sims.list.3, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose))

# Detailed info on samples and parents to examine in more detail
   saveRDS(sample.info, file = paste0(temp_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))

   saveRDS(rents.info, file = paste0(temp_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose))
   
   # Detailed info on survival
   saveRDS(pop.size.tibble, file = paste0(temp_location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))
   sim.end <- Sys.time() 
   
   iter.time <- round(as.numeric(difftime(sim.end, sim.start, units = "mins")), 1)
   cat(paste0("Finished iteration ", iter, ". \n Took ", iter.time, " minutes"))
} # end loop over iterations

########## Save and check results ##########
#Calculate relative bias for all estimates
results2 <- results %>% 
  mutate(relative_bias = round(((Q50 - truth)/truth)*100,1)) %>%
  mutate(in_interval = ifelse(HPD2.5 < truth & truth < HPD97.5, "Y", "N"))

#Within HPD interval?
results2 %>% group_by(total_samples, parameter) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100)

#Median relative bias by sample size
 results2 %>% group_by(total_samples, parameter) %>% 
   dplyr::summarize(median = median(relative_bias), n = n())

 #Mean number of parents detected
 #Median relative bias by sample size
 results2 %>% group_by(total_samples, parameter) %>% 
   dplyr::summarize(mean = mean(parents_detected), n = n())
 
 
 #-----------------------------Save major output files---------------------------------------------
 #Home computer: Dell Precision
 
 #Save model estimates
write.table(results2, file = paste0(results_location, results_prefix, "_", date.of.simulation, "_", seeds, "_", purpose, ".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
 
 #Save draws from posterior for model diagnostics 
 saveRDS(sims.list.1, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose)) #Sample size 1
 
 saveRDS(sims.list.2, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose)) #Sample size 2
 
 saveRDS(sims.list.3, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", seeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose)) #Sample size 3
 
 #Save detailed info about samples from population
 saveRDS(sample.info, file = paste0(results_location, sample.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))
 
 #Save detailed info about parents
 saveRDS(rents.info, file = paste0(results_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", purpose))
 
 #Save detailed info about survival
 saveRDS(survival.df, file = paste0(results_location, survival.prefix, "_", date.of.simulation, "_", seeds, "_", purpose))

#To read in RDS file
#pp <- readRDS("~/R/working_directory/temp_results/neutralGrowth_estSurv_iteration_5_samplesize_800")


#MGHPCC
#write.table(results2, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/Dovi_lambdaModel_06_22.2021_neutralPopGrowth.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)


#-------------Quick viz of results--------------#
#Box plot of relative bias
ggplot(data=results2, aes(x=factor(total_samples))) +
  geom_boxplot(aes(y=relative_bias, fill=parameter)) +
  ylim(-100, 100) +
  geom_hline(yintercept=0, col="black", size=1.25) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  labs(x="Sample size", y="Relative bias", title="Relative Bias by sample size") +
  scale_fill_brewer(palette="Set2") +
  font("title", size = 10, face = "bold")



####################### End ##############################