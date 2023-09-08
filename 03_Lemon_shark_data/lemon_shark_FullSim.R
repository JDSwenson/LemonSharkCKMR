###NEXT TO DO: 
# 1) Dig into the code in detail. Think the five year code is almost there.
# 2) Make code that uses entire dataset to estimate abundance in each year
# 2b) Might need to include a derived quantity that estimates Nt using N0 * lambda ^ t
# 3) Run five year code on a loop


#Load packages
library(tidyverse) # safe to ignore conflicts with filter() and lag()
library(MASS)
library(popbio)
library(mpmtools)  # not at CRAN;
# install.packages("devtools")
#devtools::install_github("BruceKendall/mpmtools")
library(ggpubr)
# The rjags package is just an interface to the JAGS library
# Make sure you have installed JAGS-4.x.y.exe (for any x >=0, y>=0) from
# http://www.sourceforge.net/projects/mcmc-jags/files
library(rjags)
library(R2jags)
library(jagsUI)
library(Rlab)
library(runjags)
library(postpack)
library(coda)

rm(list = ls())
options(dplyr.summarise.inform = FALSE)


#source("functions/Dovi_IBS_SB_test.assign.conformity_mat12.R") #All biennial breeders reproduce for the first time at age 12
source("./03_Lemon_shark_data/functions/Obj5.functions.R") #Half of biennial breeders reproduce for the first time at age 12; the other half at age 13.

#################### Set output file locations and labels ####################
temp_location <- "output/"
output.location <- "output/"
parents_prefix <- "parents.breakdown"
sample_prefix <- "sample.info"
pop.size.prefix <- "pop.size"
truth.prefix <- "truth"
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.output/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/"
jags.model_location <- "G://My Drive/Personal_Drive/R/CKMR/JAGS_models/" #Location of JAGS models

results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
jags.model.prefix <- "CKMR.JAGS_"
mom.comps.prefix <- "mom.comps"
dad.comps.prefix <- "dad.comps"
samples.prefix <- "samples/samples.missassigned"

PopSim.breeding.schedule <- "biennial.breeding_psi1" #Can be annual.breeding or biennial.breeding
model <- "multiennial.model"
jags_params = c("Nf", "Nfb1", "Nfb2", "psi", "Nm", "survival", "lambda") #List the parameters to be estimated
HS.only <- "yes"

########################## MCMC & model parameters #########################
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains

jags_file <- paste0(jags.model_location, "HS.only_narrowLambda_Skip_model.txt")


#################### Simulation parameters ####################
init.adult.pop.size <- 100 # CHANGED FROM 3000; Initial adult population size
init.prop.female <- 0.5 # proportion of the initial population size that is female
birth.sex.ratio <- c(0.5, 0.5) # The probability that each baby is F:M - has to add up to 1
YOY.survival <- 0.7 # CHANGED FROM 0.8; young of year survival
juvenile.survival <- 0.8 # CHANGED FROM 0.9; juvenile survival
Adult.survival = adult.survival <- 0.825 # CHANGED FROM 0.825; Adult survival
max.age <- maxAge <- 50 #set the maximum age allowed in the simulation

repro.age <- 12 # set age of reproductive maturity
mating.periodicity <- 2
num.mates <- c(1:3) #1  #c(1:3) #vector of potential number of mates per mating
f <- (1 - Adult.survival) / (YOY.survival * juvenile.survival ^ 11) # adult fecundity at equilibrium if no age truncation
non.conformists <- 0
psi <- 1 - non.conformists

ff <- mating.periodicity / (mating.periodicity - psi * mating.periodicity + psi) * f / init.prop.female / mean(num.mates)

ff

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
props <- rep(NA, max.age + 1)
props[1] <- f
props[2] <- f * YOY.survival
for (y in 3:(repro.age + 1)) props[y] <- props[y - 1] * juvenile.survival #+1 because of age 0 individuals

for (y in (repro.age + 2):(max.age + 1)) props[y] <- props[y - 1] * Adult.survival #+2 because of age 0 individuals
prop.Adult <- sum(props[(repro.age + 1):(max.age + 1)]) / sum(props)
Nages <- round(props[-1] * init.adult.pop.size)
init.pop.size <- sum(Nages) # all ages except YOYs

#Set length of simulation and estimation year
burn.in <- 40 # number of years to use as simulation burn in period
Num.years <- 50 # The number of years to run in the simulation beyond the burn in
n_yrs <- burn.in + Num.years #Total number of simulation years

#--------------------- Sampling parameters ---------------------
sample.years <- c(n_yrs - c(25:0)) #Twenty years of sampling


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

sim.samples.1 <- "NinetyPercent_sampled"


####-------------- Start simulation loop ----------------------
for (iter in 1:iterations) {
  sim.start <- Sys.time()
  rseed <- rseeds[iter]
  set.seed(rseed)
  
  
  #Run individual based simulation
  out <- simulate.pop(
    init.pop.size = init.pop.size,
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
    Num.years = Num.years
  )
  
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
  source("./03_Lemon_shark_data/functions/query_results_PopSim.R")
  
  ref.tibble <- NULL #Initialize dataframe for reference years i.e. second earliest birth date of juveniles in the dataset. This changes based on the sampling scheme and iteration (though it's usually the same across iterations), and is used to calculate the truth.
  
  #-----------------------Collect samples-------------------------
  sample.prop <- 90 #Will divide by 100 before setting proportion, so want this to be the percent being targeted (i.e. 90 instead of 0.9).
  
  #Initialize dataframes that will store output from EACH sampling scheme and ALL sampling schemes
  sample.df_all.info <- NULL #Final df with info from ALL sampling schemes, for ALL years
  sample.df_all.info_temp <- NULL #Df with
  
  aunt.unc_niece.nephew_pw.comps.all <- NULL #Final df with info from ALL sampling schemes, for ALL years
  aunt.unc_niece.nephew_pw.comps.all_temp <- NULL
  
  
  #Sample population each year in sample.years and make dataframe of samples with all metadata
  set.seed(rseed)
  
  #Set sampling scheme to current "target.samples"
  target.samples <- "target.YOY"
  
  #Reset dataframes that will be used within sample functions to save info from EACH sampling scheme for ALL years
  sample.df_all.info_temp1 <- NULL
  aunt.unc_niece.nephew_pw.comps.all_temp1 <- NULL
  
  #Sample over every sampling year with the chosen sampling scheme. Output is the combined samples from all years.
  samples.out <- draw.samples(target.samples = target.samples) #Saves a list of output files
  sample.size <- samples.out[[3]]
  
  sample.info <- samples.out[[1]] %>%
    as_tibble() %>% 
    mutate(
      sample.size.yr = sample.size,
      iteration = iter,
      seed = rseed,
      sample.prop = sample.prop
    ) %>%
    mutate(sample.size.total = sample.size.yr * length(sample.years))
  
  aunt.unc_niece.nephew_pw.comps.all <- samples.out[[2]] %>% 
    as_tibble() %>% 
    mutate(iteration = iter,
           seed = rseed,
           sample.prop = sample.prop)


sampled.mothers <- unique(sample.df_all.info$mother.x)
sampled.fathers <- unique(sample.df_all.info$father.x)


#UPDATE 09/08/2023 Need to include a model fit to the full dataset either before or after the loop, then integrate the results.

#------------ Run loop of sampling over five year interval ------------------#
for(block in 70:n_yrs){
  
  samples.df <- sample.info %>% dplyr::filter(capture.year >= block - 4 & capture.year <= block)
  
  #Calculate reference year
  ref.year <- calculate.ref.year(sample.info = samples.df)
  est.year.calibrate <- ref.year + 1
  
  #Compile results and summary statistics from simulation to compare estimates
  source("./03_Lemon_shark_data/functions/PopSim_truth.R")
  
  #Adds ref year to each sample depending on how it was collected
  samples.df <- samples.df %>% mutate(ref.yr = ref.year)
  
  estimation.year <- block
  
  #-------------Construct pairwise comparison matrix--------------
  #Input is 1) a list of potential parents and offspring for each offspring birth year, and 2) a dataframe with all samples, filtered to remove full siblings
  pairwise.out <- build.pairwise(filtered.samples.HS.df = samples.df)
  
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
  
  mom_comps.all <- mom_comps.all %>% filter(type == "HS")
  dad_comps.all <- dad_comps.all %>% filter(type == "HS")
  
  #HSPs detected
  mom.HSPs <- mom_comps.all %>% summarize(HSPs = sum(yes)) %>% 
    pull(HSPs)
  
  dad.HSPs <- dad_comps.all %>% summarize(HSPs = sum(yes)) %>% 
    pull(HSPs)
  

  ####------------------------ Fit CKMR model ----------------####
  
  #Define JAGS data and model, and run the MCMC engine
  set.seed(rseed)
  source("./02_Estimation.model/functions/RunJAGS_collated.R")

  ########### Calculate truth ###################
  #Calculate truth for lambda
  lambda.truth <- (pop.size.tibble$Total.adult.pop[n_yrs]/pop.size.tibble$Total.adult.pop[est.year.calibrate])^(1/(n_yrs - est.year.calibrate))
  
  lambda.truth.df <- tibble(estimation.year = estimation.year,
                            iteration = iter,
                            seed = rseed,
                            parameter = "lambda",
                            all.truth = lambda.truth)
  
  #Make dataframe of truth
  truth.iter <- pop.size.tibble %>% dplyr::filter(iteration == iter,
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
    mutate(T0 = est.year.calibrate)
  
  #Define Nf and Nm as objects for calculations
  Nf.truth <- truth.iter %>% dplyr::filter(parameter == "Nf") %>% pull(all.truth)
  
  Nm.truth <- truth.iter %>% dplyr::filter(parameter == "Nm") %>% pull(all.truth)
  
  #Subset for just the present iteration
  Exp <- calc.Exp(mom_comps.all, dad_comps.all)
  mom.Exp.HS <- Exp[[1]]
  mom.Exp.PO <- Exp[[2]]
  dad.Exp.HS <- Exp[[3]]
  dad.Exp.PO <- Exp[[4]]
  
  sampled.mothers <- unique(sample.df_all.info$mother.x)
  sampled.fathers <- unique(sample.df_all.info$father.x)
  
  
  (truth.iter_all.params <- truth.iter %>% 
      mutate(estimation.year = estimation.year) %>% 
      dplyr::mutate(T0 = ref.year + 1) %>% 
      dplyr::select(parameter, all.truth, estimation.year, T0, iteration, seed) %>% #Change order of columns
      dplyr::arrange(parameter))
  
  
  results.temp <- model.summary2 %>% left_join(truth.iter_all.params, by = c("parameter", "iteration", "seed")) %>% 
    dplyr::select(parameter,
                  Q50,
                  all.truth,
                  HPD2.5,
                  HPD97.5,
                  mean,
                  sd,
                  Q2.5,
                  Q97.5,
                  Rhat,
                  estimation.year,
                  T0,
                  iteration,
                  neff,
                  seed) %>% 
    mutate(estimation.yr = block) %>% 
    mutate(HSPs_detected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.HSPs, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.HSPs, mom.HSPs + dad.HSPs)),
             HSPs_expected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.Exp.HS, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.Exp.HS, mom.Exp.HS + dad.Exp.HS)))
    
  results <- rbind(results, results.temp)
  
  
  #Save mom and dad pairwise comparison dataframes
  mom_comps.all <- mom_comps.all %>% mutate(iteration = iter,
                                            seed = rseed)
  mom.comps.tibble <- rbind(mom.comps.tibble, mom_comps.all)
  
  dad_comps.all <- dad_comps.all %>% mutate(iteration = iter,
                                            seed = rseed)
  dad.comps.tibble <- rbind(dad.comps.tibble, dad_comps.all)
  
  
  
  
  
#-----------------Store output files iteratively--------------------

#Save parents tibble from ALL sample proportions for EACH iteration
parents.tibble_all <- bind_rows(parents.tibble_all, parents.tibble)

#  saveRDS(parents.tibble_all, file = paste0(temp_location, parents_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule, "_iter_", iter))

# Detailed info on population size for each iteration
pop.size.tibble_all <- bind_rows(pop.size.tibble_all, pop.size.tibble)

#  saveRDS(pop.size.tibble_all, file = paste0(temp_location, pop.size.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule, "_iter_", iter))

#True values for each iteration
truth.all <- bind_rows(truth.all, true.values)

#Save Aunt|uncle / nice|nephew comparisons with different name
decoy_HSPs <- aunt.unc_nephew.niece_info

#  saveRDS(truth.all, file = paste0(temp_location, truth.prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule, "_iter_", iter))

# Detailed info on samples and parents to examine in more detail
#  saveRDS(sample.info, file = paste0(temp_location, sample_prefix, "_", date.of.simulation, "_", seeds, "_", population.growth, "_", breeding.schedule, "_iter_", iter, "_", sampling.scheme))


sim.end <- Sys.time()

iter.time <-
  round(as.numeric(difftime(sim.end, sim.start, units = "mins")), 1)
cat(paste0("Finished iteration ", iter, ". \n Took ", iter.time, " minutes"))
} # end loop over iterations

#-----------------------------Save major output files---------------------------------------------
#Save detailed info about samples from population
saveRDS(
  sample.info2,
  file = paste0(
    output.location,
    sample_prefix,
    "_",
    date.of.simulation,
    "_",
    seeds,
    "_",
    population.growth,
    "_",
    breeding.schedule
  )
)

#Save parents tibble
saveRDS(
  parents.tibble_all,
  file = paste0(
    output.location,
    parents_prefix,
    "_",
    date.of.simulation,
    "_",
    seeds,
    "_",
    population.growth,
    "_",
    breeding.schedule
  )
)

# Detailed info on population size
saveRDS(
  pop.size.tibble_all,
  file = paste0(
    output.location,
    pop.size.prefix,
    "_",
    date.of.simulation,
    "_",
    seeds,
    "_",
    population.growth,
    "_",
    breeding.schedule
  )
)

# Truth
saveRDS(
  truth.all,
  file = paste0(
    output.location,
    truth.prefix,
    "_",
    date.of.simulation,
    "_",
    seeds,
    "_",
    population.growth,
    "_",
    breeding.schedule
  )
)

#Charlatan HSPs
saveRDS(
  decoy_HSPs,
  file = paste0(
    output.location,
    "aunt.uncle_niece.nephews_",
    date.of.simulation,
    "_",
    seeds,
    "_",
    population.growth,
    "_",
    breeding.schedule
  )
)
}