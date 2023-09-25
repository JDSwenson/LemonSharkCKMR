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


source("./03_Lemon_shark_data/functions/Obj5.functions.R") #Added population simulation to this script.

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

PopSim.breeding.schedule <- "biennial.breeding" #Can be annual.breeding or biennial.breeding
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

percent.skip_on.cycle <- 0.1
percent.breed_off.cycle <- 0.1

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
sample.years <- c(n_yrs - c(20:0)) #Twenty years of sampling


####-------------- Prep simulation ----------------------
# Moved sampling below so extract different sample sizes from same population
iterations <- 100


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
for (iter in 1:iterations) { #Restart morning of 09/20/2023 - made it through iteration 19
  sim.start <- Sys.time()
  rseed <- rseeds[iter]
  set.seed(rseed)
  
  
  #Run individual based simulation - population will increase for 10 years starting at year 70, stay approximately stable for 10 years, then decrease for the last five years
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
  
  rm(out)
  
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

  #Check on samples
  sample.info %>% dplyr::count(age.x) #Should all be 0
  sample.info %>% dplyr::count(birth.year) #Should be relatively even numbers across years
    
  #Don't need this ...
  # aunt.unc_niece.nephew_pw.comps.all <- samples.out[[2]] %>% 
  #   as_tibble() %>% 
  #   mutate(iteration = iter,
  #          seed = rseed,
  #          sample.prop = sample.prop)


sampled.mothers <- unique(sample.info$mother.x)
sampled.fathers <- unique(sample.info$father.x)

#First year of estimation will be four years after the first year of sampling
first.est.year <- sample.years + 4

#------------ Run loop of sampling over five year interval ------------------#
for(block in first.est.year:n_yrs){
  
  #Estimate in the present
  estimation.year <- block

  #-------------Split samples--------------
  #Blocked  samples - five year
  samples.df.block5 <- sample.info %>% dplyr::filter(capture.year >= block - 4 & capture.year <= block)
  
  #Blocked  samples - three year
  samples.df.block3 <- sample.info %>% dplyr::filter(capture.year >= block - 2 & capture.year <= block)
  
  #All samples
  samples.df.all <- sample.info %>% dplyr::filter(capture.year <= block)
  
#--------------Calculate reference year--------------
  #Five year block
  ref.year.block5 <- calculate.ref.year(sample.info = samples.df.block5)
  est.year.calibrate.block5 <- ref.year.block5 + 1
  
  #Three year block
  ref.year.block3 <- calculate.ref.year(sample.info = samples.df.block3)
  est.year.calibrate.block3 <- ref.year.block3 + 1
  
  #All samples
  ref.year.all <- calculate.ref.year(sample.info = samples.df.all)
  est.year.calibrate.all <- ref.year.all + 1
  
  
  #Adds ref year to each sample
  #Five year block
  samples.df.block5 <- samples.df.block5 %>% mutate(ref.yr = ref.year.block5)
  
  #Three year block
  samples.df.block3 <- samples.df.block3 %>% mutate(ref.yr = ref.year.block3)
  
  #All samples
  samples.df.all <- samples.df.all %>% mutate(ref.yr = ref.year.all)
  

  ####------------------------ Split block and all ----------------
  #----------------Start with blocked data ---------------------#
  #---------------Three year block------------------------
  #-------------Construct pairwise comparison matrix--------------
  pairwise.out <- build.pairwise(filtered.samples.HS.df = samples.df.block3)
  mom_comps.all.block3 <- pairwise.out[[1]]  
  dad_comps.all.block3 <- pairwise.out[[2]]
  positives.HS.block3 <- pairwise.out[[3]]
  
  mom_comps.all.block3 <- mom_comps.all.block3 %>% filter(type == "HS")
  dad_comps.all.block3 <- dad_comps.all.block3 %>% filter(type == "HS")
  mom.HSPs.block3 <- mom_comps.all.block3 %>% summarize(HSPs = sum(yes)) %>%
    pull(HSPs)
  dad.HSPs.block3 <- dad_comps.all.block3 %>% summarize(HSPs = sum(yes)) %>%
    pull(HSPs)
  
  
  
  #yrs <- c(estimation.year:n_yrs)
  ref.year.mom <- min(mom_comps.all.block3$ref.year)
  ref.year.dad <- min(dad_comps.all.block3$ref.year)
  est.year.calibrate <- est.year.calibrate.block3
  
  #Create vectors of data for JAGS
  #Mom
  #HS - even years
  mom_comps.HS_on <- mom_comps.all.block3 %>% dplyr::filter(type == "HS", BI == "on")
  mom.mort.yrs_HS.on <- mom_comps.HS_on$mort.yrs
  mom.popGrowth.yrs_HS.on <- mom_comps.HS_on$pop.growth.yrs
  mom.n.comps_HS.on <- mom_comps.HS_on$all
  mom.positives_HS.on <- mom_comps.HS_on$yes
  mom.yrs_HS.on <- nrow(mom_comps.HS_on)
  
  #Mom
  #HS - off years
  mom_comps.HS_off <- mom_comps.all.block3 %>% dplyr::filter(type == "HS", BI == "off")
  mom.mort.yrs_HS.off <- mom_comps.HS_off$mort.yrs
  mom.popGrowth.yrs_HS.off <- mom_comps.HS_off$pop.growth.yrs
  mom.n.comps_HS.off <- mom_comps.HS_off$all
  mom.positives_HS.off <- mom_comps.HS_off$yes
  mom.yrs_HS.off <- nrow(mom_comps.HS_off)
  
  #Dad
  dad.mort.yrs <- dad_comps.all.block3$mort.yrs
  dad.popGrowth.yrs <- dad_comps.all.block3$pop.growth.yrs
  dad.n.comps <- dad_comps.all.block3$all
  dad.positives <- dad_comps.all.block3$yes
  dad.yrs <- nrow(dad_comps.all.block3)
  
  cat("Running model using three year blocks of samples\n")
  
  set.seed(rseed)
  cat(paste0("Using samples from year ", block - 2, " to ", block, "\n"))
  
  source("./03_Lemon_shark_data/functions/Obj5_run.JAGS_HS.only.R")
  
  
  
  ########### Calculate truth ###################
  #Calculate truth for survival and psi
  #results in a file called "true.values" that contains the truth for survival and psi
  est.year.calibrate <- est.year.calibrate.block3
  ref.year <- ref.year.block3
  source("./03_Lemon_shark_data/functions/PopSim_truth.R")
  
  #Pull lambda truth
  lambda.truth <- (pop.size.tibble$adult.lambda[block])
  lambda.truth.df <- tibble(estimation.year = estimation.year,
                            iteration = iter,
                            seed = rseed,
                            parameter = "lambda",
                            all.truth = lambda.truth)
  
  #Make dataframe of truth
  truth.iter <- pop.size.tibble %>% dplyr::filter(year == estimation.year) %>% 
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
    bind_rows(lambda.truth.df, true.values) %>% 
    mutate(T0 = est.year.calibrate)
  
  #Define Nf and Nm as objects for calculations
  Nf.truth <- truth.iter %>% dplyr::filter(parameter == "Nf") %>% pull(all.truth)
  
  Nm.truth <- truth.iter %>% dplyr::filter(parameter == "Nm") %>% pull(all.truth)
  
  #Subset for just the present iteration
  Exp.block3 <- calc.Exp(mom_comps.all.block3, dad_comps.all.block3)
  mom.Exp.HS.block3 <- Exp.block3[[1]]
  mom.Exp.PO.block3 <- Exp.block3[[2]]
  dad.Exp.HS.block3 <- Exp.block3[[3]]
  dad.Exp.PO.block3 <- Exp.block3[[4]]
  
  sampled.mothers <- unique(sample.df_all.info$mother.x)
  sampled.fathers <- unique(sample.df_all.info$father.x)
  
  
  results.temp.block3 <- model.summary2 %>% left_join(truth.iter, by = c("parameter", "iteration", "seed")) %>% 
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
    mutate(HSPs_detected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.HSPs.block3, 
                                  ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.HSPs.block3, mom.HSPs.block3 + dad.HSPs.block3)),
           HSPs_expected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.Exp.HS.block3, 
                                  ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.Exp.HS.block3, mom.Exp.HS.block3 + dad.Exp.HS.block3))) %>% 
    mutate(time_window = "three year block")
  
  gc()
  #---------------Five year block------------------------
  #-------------Construct pairwise comparison matrix--------------
  pairwise.out <- build.pairwise(filtered.samples.HS.df = samples.df.block5)
  mom_comps.all.block5 <- pairwise.out[[1]]  
  dad_comps.all.block5 <- pairwise.out[[2]]
  positives.HS.block5 <- pairwise.out[[3]]
  
  mom_comps.all.block5 <- mom_comps.all.block5 %>% filter(type == "HS")
  dad_comps.all.block5 <- dad_comps.all.block5 %>% filter(type == "HS")
  mom.HSPs.block5 <- mom_comps.all.block5 %>% summarize(HSPs = sum(yes)) %>%
    pull(HSPs)
  dad.HSPs.block5 <- dad_comps.all.block5 %>% summarize(HSPs = sum(yes)) %>%
    pull(HSPs)
  
  #yrs <- c(estimation.year:n_yrs)
  ref.year.mom <- min(mom_comps.all.block5$ref.year)
  ref.year.dad <- min(dad_comps.all.block5$ref.year)
  est.year.calibrate <- est.year.calibrate.block5
  
  #Create vectors of data for JAGS
  #Mom
  #HS - even years
  mom_comps.HS_on <- mom_comps.all.block5 %>% dplyr::filter(type == "HS", BI == "on")
  mom.mort.yrs_HS.on <- mom_comps.HS_on$mort.yrs
  mom.popGrowth.yrs_HS.on <- mom_comps.HS_on$pop.growth.yrs
  mom.n.comps_HS.on <- mom_comps.HS_on$all
  mom.positives_HS.on <- mom_comps.HS_on$yes
  mom.yrs_HS.on <- nrow(mom_comps.HS_on)
  
  #Mom
  #HS - off years
  mom_comps.HS_off <- mom_comps.all.block5 %>% dplyr::filter(type == "HS", BI == "off")
  mom.mort.yrs_HS.off <- mom_comps.HS_off$mort.yrs
  mom.popGrowth.yrs_HS.off <- mom_comps.HS_off$pop.growth.yrs
  mom.n.comps_HS.off <- mom_comps.HS_off$all
  mom.positives_HS.off <- mom_comps.HS_off$yes
  mom.yrs_HS.off <- nrow(mom_comps.HS_off)
  
  #Dad
  dad.mort.yrs <- dad_comps.all.block5$mort.yrs
  dad.popGrowth.yrs <- dad_comps.all.block5$pop.growth.yrs
  dad.n.comps <- dad_comps.all.block5$all
  dad.positives <- dad_comps.all.block5$yes
  dad.yrs <- nrow(dad_comps.all.block5)
  
  cat("Running model using five year blocks of samples\n")
  
  set.seed(rseed)
  cat(paste0("Using samples from year ", block - 4, " to ", block, "\n"))
  
  source("./03_Lemon_shark_data/functions/Obj5_run.JAGS_HS.only.R")

  ########### Calculate truth ###################
  #Calculate truth for survival and psi
  #results in a file called "true.values" that contains the truth for survival and psi
  est.year.calibrate <- est.year.calibrate.block5
  ref.year <- ref.year.block5
  source("./03_Lemon_shark_data/functions/PopSim_truth.R")
  
  #Calculate truth for lambda
  lambda.truth <- (pop.size.tibble$adult.lambda[block])
  
  lambda.truth.df <- tibble(estimation.year = estimation.year,
                            iteration = iter,
                            seed = rseed,
                            parameter = "lambda",
                            all.truth = lambda.truth)
  
  #Make dataframe of truth
  truth.iter <- pop.size.tibble %>% dplyr::filter(year == estimation.year) %>% 
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
    bind_rows(lambda.truth.df, true.values) %>% 
    mutate(T0 = est.year.calibrate)
  
  #Define Nf and Nm as objects for calculations
  Nf.truth <- truth.iter %>% dplyr::filter(parameter == "Nf") %>% pull(all.truth)
  
  Nm.truth <- truth.iter %>% dplyr::filter(parameter == "Nm") %>% pull(all.truth)
  
  #Subset for just the present iteration
  Exp.block5 <- calc.Exp(mom_comps.all.block5, dad_comps.all.block5)
  mom.Exp.HS.block5 <- Exp.block5[[1]]
  mom.Exp.PO.block5 <- Exp.block5[[2]]
  dad.Exp.HS.block5 <- Exp.block5[[3]]
  dad.Exp.PO.block5 <- Exp.block5[[4]]
  
  sampled.mothers <- unique(sample.df_all.info$mother.x)
  sampled.fathers <- unique(sample.df_all.info$father.x)
  
  
  results.temp.block5 <- model.summary2 %>% left_join(truth.iter, by = c("parameter", "iteration", "seed")) %>% 
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
    mutate(HSPs_detected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.HSPs.block5, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.HSPs.block5, mom.HSPs.block5 + dad.HSPs.block5)),
             HSPs_expected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.Exp.HS.block5, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.Exp.HS.block5, mom.Exp.HS.block5 + dad.Exp.HS.block5))) %>% 
    mutate(time_window = "five year block")
  
  gc()
  
  #----------------Use all data ---------------------
  #-------------Construct pairwise comparison matrix--------------
  pairwise.out <- build.pairwise(filtered.samples.HS.df = samples.df.all)
  mom_comps.all.all <- pairwise.out[[1]]  
  dad_comps.all.all <- pairwise.out[[2]]
  positives.HS.all <- pairwise.out[[3]]
  
  mom_comps.all.all <- mom_comps.all.all %>% filter(type == "HS")
  dad_comps.all.all <- dad_comps.all.all %>% filter(type == "HS")
  mom.HSPs.all <- mom_comps.all.all %>% summarize(HSPs = sum(yes)) %>%
    pull(HSPs)
  dad.HSPs.all <- dad_comps.all.all %>% summarize(HSPs = sum(yes)) %>%
    pull(HSPs)
  
  #yrs <- c(estimation.year:n_yrs)
  ref.year.mom <- min(mom_comps.all.all$ref.year)
  ref.year.dad <- min(dad_comps.all.all$ref.year)
  
  #Create vectors of data for JAGS
  #Mom
  #HS - even years
  mom_comps.HS_on <- mom_comps.all.all %>% dplyr::filter(type == "HS", BI == "on")
  mom.mort.yrs_HS.on <- mom_comps.HS_on$mort.yrs
  mom.popGrowth.yrs_HS.on <- mom_comps.HS_on$pop.growth.yrs
  mom.n.comps_HS.on <- mom_comps.HS_on$all
  mom.positives_HS.on <- mom_comps.HS_on$yes
  mom.yrs_HS.on <- nrow(mom_comps.HS_on)
  
  #Mom
  #HS - off years
  mom_comps.HS_off <- mom_comps.all.all %>% dplyr::filter(type == "HS", BI == "off")
  mom.mort.yrs_HS.off <- mom_comps.HS_off$mort.yrs
  mom.popGrowth.yrs_HS.off <- mom_comps.HS_off$pop.growth.yrs
  mom.n.comps_HS.off <- mom_comps.HS_off$all
  mom.positives_HS.off <- mom_comps.HS_off$yes
  mom.yrs_HS.off <- nrow(mom_comps.HS_off)
  
  #Dad
  dad.mort.yrs <- dad_comps.all.all$mort.yrs
  dad.popGrowth.yrs <- dad_comps.all.all$pop.growth.yrs
  dad.n.comps <- dad_comps.all.all$all
  dad.positives <- dad_comps.all.all$yes
  dad.yrs <- nrow(dad_comps.all.all)
  
  cat("Running model using all samples\n")
  
  set.seed(rseed)
  cat(paste0("Using all samples with estimation year ", block, "\n"))
  
  source("./03_Lemon_shark_data/functions/Obj5_run.JAGS_HS.only.R")
  
  ########### Calculate truth ###################
  #Calculate truth for survival and psi
  #results in a file called "true.values" that contains the truth for survival and psi
  est.year.calibrate <- est.year.calibrate.all
  ref.year <- ref.year.all
  source("./03_Lemon_shark_data/functions/PopSim_truth.R")
  
  #Calculate truth for lambda
  lambda.truth <- (pop.size.tibble$adult.lambda[block])
  
  lambda.truth.df <- tibble(estimation.year = estimation.year,
                            iteration = iter,
                            seed = rseed,
                            parameter = "lambda",
                            all.truth = lambda.truth)
  
  #Make dataframe of truth
  truth.iter <- pop.size.tibble %>% dplyr::filter(year == estimation.year) %>% 
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
    bind_rows(lambda.truth.df, true.values) %>% 
    mutate(T0 = est.year.calibrate)
  
  #Define Nf and Nm as objects for calculations
  Nf.truth <- truth.iter %>% dplyr::filter(parameter == "Nf") %>% pull(all.truth)
  
  Nm.truth <- truth.iter %>% dplyr::filter(parameter == "Nm") %>% pull(all.truth)
  
  #Subset for just the present iteration
  Exp.all <- calc.Exp(mom_comps.all.all, dad_comps.all.all)
  mom.Exp.HS.all <- Exp.all[[1]]
  mom.Exp.PO.all <- Exp.all[[2]]
  dad.Exp.HS.all <- Exp.all[[3]]
  dad.Exp.PO.all <- Exp.all[[4]]
  
  sampled.mothers <- unique(sample.df_all.info$mother.x)
  sampled.fathers <- unique(sample.df_all.info$father.x)
  
  
  results.temp.all <- model.summary2 %>% left_join(truth.iter, by = c("parameter", "iteration", "seed")) %>% 
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
    mutate(HSPs_detected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.HSPs.all, 
                                  ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.HSPs.all, mom.HSPs.all + dad.HSPs.all)),
           HSPs_expected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.Exp.HS.all, 
                                  ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.Exp.HS.all, mom.Exp.HS.all + dad.Exp.HS.all))) %>% 
    mutate(time_window = "all samples")
  
  
  ####------------------------ Combine block and all ----------------####
  #Save mom and dad pairwise comparison dataframes
  #Five year block
  mom_comps.all.block5 <- mom_comps.all.block5 %>% mutate(iteration = iter,
                                            seed = rseed,
                                            time_window = "five year block")
  
  dad_comps.all.block5 <- dad_comps.all.block5 %>% mutate(iteration = iter,
                                                        seed = rseed,
                                                        time_window = "five year block")
  
  
#Three year block  
  mom_comps.all.block3 <- mom_comps.all.block3 %>% mutate(iteration = iter,
                                                          seed = rseed,
                                                          time_window = "three year block")
  
  dad_comps.all.block3 <- dad_comps.all.block3 %>% mutate(iteration = iter,
                                                          seed = rseed,
                                                          time_window = "three year block")
  
  
#All samples
  mom_comps.all.all <- mom_comps.all.all %>% mutate(iteration = iter,
                                                        seed = rseed,
                                                        time_window = "all samples")
  
  
  dad_comps.all.all <- dad_comps.all.all %>% mutate(iteration = iter,
                                                        seed = rseed,
                                                        time_window = "five year all")
  
  mom.comps.tibble <- rbind(mom.comps.tibble, mom_comps.all.block3, mom_comps.all.block5, mom_comps.all.all)
  
  dad.comps.tibble <- rbind(dad.comps.tibble, dad_comps.all.block3, dad_comps.all.block5, dad_comps.all.all)
  
  
#-----------------Store output files iteratively--------------------
#Results
results <- rbind(results, results.temp.block3, results.temp.block5, results.temp.all)
  
#Save parents tibble from ALL sample proportions for EACH iteration
#parents.tibble_all <- bind_rows(parents.tibble_all, parents.tibble)

# Detailed info on population size for each iteration
pop.size.tibble_all <- bind_rows(pop.size.tibble_all, pop.size.tibble)

#True values for each iteration
truth.all <- bind_rows(truth.all, true.values)
} #End loop over buh-locks


sim.end <- Sys.time()

iter.time <-
  round(as.numeric(difftime(sim.end, sim.start, units = "mins")), 1)
cat(paste0("Finished iteration ", iter, ". \n Took ", iter.time, " minutes"))

gc()
} # end loop over iterations

#If using all individuals for Nf truth, instead of breeders
results2 <- results %>%
  mutate(relative_bias = round(((Q50 - all.truth)/all.truth)*100, 1)) %>% #Can change truth to breed.truth if looking for number of active breeders
  mutate(in_interval = ifelse(HPD2.5 < all.truth & all.truth < HPD97.5, "Y", "N")) %>%
  as_tibble()

#Within HPD interval?
results2 %>% group_by(parameter, estimation.year, time_window) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100)

#Median relative bias by sample size
results2 %>% group_by(parameter, estimation.year, time_window) %>% 
  dplyr::summarize(median.bias = median(relative_bias), n = n()) %>% 
  dplyr::arrange(desc(median.bias))


#-----------------------------Save major output files---------------------------------------------
#Home computer: Dell Precision
#Save model estimates
write.table(results2, file = paste0(results_location, "lemon_shark_time_series_sims_results.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

#Save pop size estimates
saveRDS(pop.size.tibble_all, file = paste0(results_location, "lemon_shark_time_series_sims_popSize"))

#Save draws from posterior for model diagnostics 
#saveRDS(sims.list.1, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.1, "_", MCMC.settings, "_", scenario, "_", model)) #Sample size 1

#Save final pairwise comparison matrices
saveRDS(mom.comps.tibble, file = paste0(results_location, "lemon_shark_time_series_sims_mom.comps"))

saveRDS(dad.comps.tibble, file = paste0(results_location, "lemon_shark_time_series_sims_dad.comps"))

#-------------Quick viz of results--------------#
results2 %>% ggplot(aes(x = relative_bias,
           y = parameter)) + 
  geom_density_ridges(alpha = 0.6) + 
  labs(x="Relative bias", y="Parameter") +
  #  theme_ridges() + 
  facet_wrap(~time_window) + 
  theme_bw() +
  xlim(-100, 100) +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed", size = 1.1) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "none")
