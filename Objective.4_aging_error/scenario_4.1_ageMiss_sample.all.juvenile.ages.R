#Objective 1.1: Validate model using informed priors

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

source("./Objective.3_intermittent.breeding/functions/Obj3.functions.R") #Changed name of script that includes pairwise comparison and other functions

#----------------Set input file locations ------------------------------
PopSim.location <- "G://My Drive/Personal_Drive/R/CKMR/Population.simulations/"
PopSim.lambda <- "lambda.1" # Can be lambda.1 or lambda.variable
PopSim.breeding.schedule <- "biennial.breeding" #Can be annual.breeding or biennial.breeding
Sampling.scheme <- "sample.all.juvenile.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
date.of.PopSim <- "19Jul2022" # 11Jul2022
inSeeds <- "Seeds2022.04.15"

#----------------Set output file locations ------------------------------
temp_location <- "~/R/working_directory/temp_results/"
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.3_intermittent.breeding/Model.output/"
jags.model_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.3_intermittent.breeding/models/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/Objective.3_intermittent.breeding/Model.results/"

results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
jags.model.prefix <- "CKMR.JAGS_"
mom.comps.prefix <- "comparisons/mom.comps"
dad.comps.prefix <- "comparisons/dad.comps"


#-------------------Set simulation settings and scenario info----------------------------
script_name <- "scenario_3.2.1_SB_sample.all.juvenile.ages.R" #Copy name of script here
primary_goal <- "Test model performance with biennial breeding" #Why am I running this simulation? Provide details

question1 <- "How does a base-case CKMR model perform with biennial breeding?"
question2 <- "Do we need to account for this in a CKMR model for elasmobranchs?"
question3 <- "How does the model perform with and without lambda?"
purpose <- "scenario_3.2.1_SB_sample.all.juvenile.ages" #For naming output files
today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
date.of.simulation <- today

target.YOY <- "no" #For juvenile samples, do we only want to target YOY for each year of sampling?
down_sample <- "no" #Do we want to downsample to achieve close to max.HSPs?
max.HSPs <- NA
max.POPs <- NA
HS.only <- "yes" #Do we only want to filter HS relationships?
PO.only <- "no" #Do we only want to filter PO relationships? These two are mutually exclusive; cannot have "yes" for both
fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none" and the full model will run, estimating all parameters. If fixing specific parameters, then list them here, and manually change in the run.JAGS_HS.PO_SB.R script
jags_params = c("Nf", "Nm", "survival", "psi")
estimated.parameters <- paste0(jags_params, collapse = ",")

#rseeds <- sample(1:1000000,iterations)
#save(rseeds, file = "rseeds_2022.04.15.rda")

#Save paths and file labels as objects
rseeds <- readRDS("rseeds_2022.04.15.rda")
outSeeds <- "Seeds2022.04.15"


#--------------------Leslie matrix parameters----------------------------
#These may not be relevant for model validation

adult.survival <- 0.825 # CHANGED FROM 0.825; Adult survival
repro.age <- 12 # set age of reproductive maturity
max.age <- 50 #set the maximum age allowed in the simulation
mating.periodicity <- 2 #number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.

#---------------------- Read in sampling and other dataframes --------------------------
samples.df <- readRDS(file = paste0(PopSim.location, "sample.info_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", Sampling.scheme)) %>% 
  rename(total.sample.size = sample.size.juvs)

samples.df %>% group_by(age.x) %>% summarize(n())

pop_size.df <- readRDS(file = paste0(PopSim.location, "pop.size_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", Sampling.scheme))

truth.df <- readRDS(file = paste0(PopSim.location, "truth_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", Sampling.scheme))

n_yrs <- max(pop_size.df$year)

estimation.year <- n_yrs - 5

#----------------------- MCMC & model parameters ----------------------#
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains

#Survival prior info
survival.prior.mean <- NA
survival.prior.cv <- NA
survival.prior.sd <- NA
survival.prior.info <- "Uniform: 0.5 - 0.95"

#Lambda prior info
lambda.prior.mean <- NA
lambda.prior.cv <- NA
lambda.prior.sd <- NA
lambda.prior.info <- NA

#psi prior
psi.prior.info <- "Uniform: 0.5 - 0.99"

#abundance prior
abundance.prior.info <- "diffuse Normal w diffuse Uniform hyperprior"


####---------------Update and save simulation log-------------------####
#Will want to change for Objective 2 to include CVs and misspecified parameter values#
model_settings.df <- tibble(script_name = script_name,
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
                            seeds = outSeeds,
                            thinning_rate = nt,
                            posterior_samples = ni,
                            burn_in = nb,
                            survival.prior.info = survival.prior.mean,
                            lambda.prior.info = lambda.prior.mean,
                            psi.prior = psi.prior.info,
                            abundance.prior = abundance.prior.info
)

#Save simulation settings in Simulation_log
 # model.log <- read_csv("model_settings.log.csv")
 # tail(model.log)
 # (model.log_updated <- bind_rows(model.log, model_settings.df)) #Combine old simulation settings with these
 # write_csv(model.log_updated, file = "model_settings.log.csv") #Save the updated simulation log

####-------------- Start simulation loop ----------------------
(iterations <- max(samples.df$iteration))
(sample.sizes <- samples.df %>% 
    distinct(sample.prop) %>% 
    dplyr::filter(sample.prop == 1.5) %>% 
    pull(sample.prop)) #Subset for sample size of 1.5%


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
    

        ####-----------------------------Downsample if more than max.HSPs------------------------------------####
    if(down_sample == "yes"){
      set.seed(rseed)
      
    down.samples <- downsample(mom_comps.all, dad_comps.all, HS.samps.df, PO.samps.list) #calculates approximate number of samples needed to achieve max.HSPs for filtering below
    HS.samps.df.down <- down.samples[[1]]
    PO.samps.list.down <- down.samples[[2]]
    
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
      
      #HSPs detected
      mom.HSPs <- mom_comps.all %>% dplyr::filter(type == "HS") %>% 
        summarize(HSPs = sum(yes)) %>% 
        pull(HSPs)
      
      dad.HSPs <- dad_comps.all %>% dplyr::filter(type == "HS") %>% 
        summarize(HSPs = sum(yes)) %>% 
        pull(HSPs)
      
      #POPs detected
      # mom.POPs <- mom_comps.all %>% dplyr::filter(type == "PO") %>% 
      #   summarize(POPs = sum(yes)) %>% 
      #   pull(POPs)
      # 
      # dad.POPs <- dad_comps.all %>% dplyr::filter(type == "PO") %>% 
      #   summarize(POPs = sum(yes)) %>% 
      #   pull(POPs)
      
    # ####------------------------ Fit CKMR model ----------------####
    #Define JAGS data and model, and run the MCMC engine
      set.seed(rseed)
    source("Objective.3_intermittent.breeding/functions/scenario_3.2.1_run.JAGS_HS.only.R")

      #Calculate truth
      Nf.truth <- pop_size.df %>% dplyr::filter(iteration == iter,
                                                year == estimation.year) %>%
        pull(Female.adult.pop)
      
      Nm.truth <- pop_size.df %>% dplyr::filter(iteration == iter,
                                                year == estimation.year) %>%
        pull(Male.adult.pop)      
      
    #Calculate expectations
    pop.size.tibble <- pop_size.df %>% dplyr::filter(iteration == iter)
    Exp <- calc.Exp(mom_comps.all, dad_comps.all)
    (mom.Exp.HS <- Exp[[1]])
#    (mom.Exp.PO <- Exp[[2]])
    (dad.Exp.HS <- Exp[[3]])
#    (dad.Exp.PO <- Exp[[4]])
    
    sampled.mothers <- unique(sample.df_all.info$mother.x)
    sampled.fathers <- unique(sample.df_all.info$father.x)
    
    
    
    #Make truth equal to mean over years
    truth.iter <- truth.df %>% dplyr::filter(iteration == iter) %>% 
      mutate(all.truth = ifelse(parameter == "Nf", Nf.truth,
                                ifelse(parameter == "Nm", Nm.truth, all.truth)))

      samples.iter <- samples.df %>% dplyr::filter(iteration == iter, sample.prop == sample.size) %>% 
      distinct(seed, iteration, sample.prop.juvs = sample.prop, total.sample.size)
    
    results.temp <- model.summary2 %>% left_join(truth.iter, by = c("parameter", "iteration", "seed")) %>% 
      left_join(samples.iter, by = c("iteration", "seed")) %>% 
      mutate(HSPs_detected = c(mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs),
             HSPs_expected = c(mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS),
             purpose = purpose,
             est.yr = estimation.year)
    
    results <- rbind(results, results.temp)
    
    
    #Save mom and dad pairwise comparison dataframes
    mom_comps.all <- mom_comps.all %>% mutate(iteration = iter,
                                              sample.prop.all = sample.size,
                                              sample.size.juvs = nrow(HS.samps.df),
                                              seed = rseed,
                                              est.yr = estimation.year)
    mom.comps.tibble <- rbind(mom.comps.tibble, mom_comps.all)
    
    dad_comps.all <- dad_comps.all %>% mutate(iteration = iter,
                                              sample.prop.all = sample.size,
                                              sample.size.juvs = nrow(HS.samps.df),
                                              seed = rseed,
                                              est.yr = estimation.year)
    dad.comps.tibble <- rbind(dad.comps.tibble, dad_comps.all)
    
  } # End if/else statement
  } # end loop over sample sizes
   
if(iter %% 100 == 0){
  #-----------------Save output files iteratively--------------------
 

  #Results
    write.table(results, file = paste0(temp_location, results_prefix, "_", date.of.simulation, "_", outSeeds, "_", purpose, "_", estimation.year, "_est.yr_", "_iter_", iter, ".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
# 
#    #Model output for diagnostics
     saveRDS(sims.list.1, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose, "_", estimation.year, "_est.yr"))
# 
#     saveRDS(sims.list.2, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose, "_", estimation.year, "_est.yr"))
# #
#     saveRDS(sims.list.3, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose, "_", estimation.year, "_est.yr"))
# #
#     saveRDS(sims.list.4, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.4, "_", MCMC.settings, "_", purpose, "_", estimation.year, "_est.yr"))
# 
#    #Save pairwise comparisons matrices
    saveRDS(mom.comps.tibble, file = paste0(temp_location, mom.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", purpose, "_", estimation.year, "_est.yr"))
#    
    saveRDS(dad.comps.tibble, file = paste0(temp_location, dad.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", purpose, "_", estimation.year, "_est.yr"))

 }
  
  sim.end <- Sys.time()
  
  iter.time <- round(as.numeric(difftime(sim.end, sim.start, units = "mins")), 1)
  cat(paste0("\n Finished iteration ", iter, ". \n Took ", iter.time, " minutes \n\n"))
    
 }  # end loop over iterations
  
  
  
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
     as_tibble()
   
#Within HPD interval?
results2 %>% group_by(sample.prop.juvs, parameter, purpose) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100)

#Median relative bias by sample size
 results2 %>% group_by(sample.prop.juvs, parameter, purpose) %>% 
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
write.table(results2, file = paste0(results_location, results_prefix, "_", date.of.simulation, "_", outSeeds, "_", purpose, ".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
 
 #Save draws from posterior for model diagnostics 
 saveRDS(sims.list.1, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.1, "_", MCMC.settings, "_", purpose)) #Sample size 1
 
 #  saveRDS(sims.list.2, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.2, "_", MCMC.settings, "_", purpose)) #Sample size 2
 # #
 #  saveRDS(sims.list.3, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.3, "_", MCMC.settings, "_", purpose)) #Sample size 3
 # #
 #  saveRDS(sims.list.4, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.4, "_", MCMC.settings, "_", purpose)) #Sample size 4
 
 #Save final pairwise comparison matrices
 saveRDS(mom.comps.tibble, file = paste0(results_location, mom.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", purpose))
 
 saveRDS(dad.comps.tibble, file = paste0(results_location, dad.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", purpose))
 
#To read in RDS file
#pp <- readRDS("~/R/working_directory/temp_results/neutralGrowth_estSurv_iteration_5_samplesize_800")


#MGHPCC
#write.table(results2, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/Dovi_lambdaModel_06_22.2021_neutralPopGrowth.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)


#-------------Quick viz of results--------------#
#Box plot of relative bias
ggplot(data=results2, aes(x=factor(sample.prop.juvs))) +
  geom_boxplot(aes(y=relative_bias, fill=parameter)) +
  ylim(-100, 100) +
  geom_hline(yintercept=0, col="black", size=1.25) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  labs(x="Sample size", y="Relative bias", title="Relative Bias by sample size") +
  scale_fill_brewer(palette="Set2") +
  font("title", size = 10, face = "bold")



####################### End ##############################