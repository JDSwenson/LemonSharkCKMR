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
library(FSA)

rm(list=ls())

######################### Specify common output prefixes #########################
results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
jags.model.prefix <- "CKMR.JAGS_"
mom.comps.prefix <- "comparisons/mom.comps"
dad.comps.prefix <- "comparisons/dad.comps"
samples.prefix <- "samples/samples.missassigned"

######################### Specify other common settings #########################
temp_location <- "~/R/working_directory/temp_results/" #Location to save temporary files in case run gets killed
jags.model_location <- "G://My Drive/Personal_Drive/R/CKMR/JAGS_models/" #Location of JAGS models

#Focus on objective 2 because this should have the model and functions we want
source("./02_Estimation.model/functions/Obj2.functions.R") #Changed name of script that includes pairwise comparison and other functions
today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
date.of.simulation <- today

#I mean may as well
rseeds <- readRDS("rseeds_2022.04.15.rda")
outSeeds <- "Seeds2022.04.15"
adult.survival <- 0.825 # CHANGED FROM 0.825; Adult survival
repro.age <- 12 # set age of reproductive maturity
max.age <- 50 #set the maximum age allowed in the simulation


#------------------------- Set output file locations -------------------------#
 MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/SER/Model.output/"
 results_location <- "G://My Drive/Personal_Drive/R/CKMR/SER/Model.results/"

 estimation.years <- c(n_yrs - 10, n_yrs - 5, n_yrs)
 fixed.parameters <- "none" #List the fixed parameters here; if none, then leave as "none".
 model <- "annual.model" #For naming output files
 sampling.scheme <- "sample.ALL.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
 HS.only <- "no" #Do we only want to filter HS relationships?
 PO.only <- "no" #Do we only want to filter PO relationships?
 jags_file <- paste0(jags.model_location, "HS.PO_narrowLambda_annual_model.txt") #Need to use a different model if sampling all age classes AND using the biennial model


########################## Read in sampling and other dataframes #########################
# samples.df <- readRDS(file = paste0(PopSim.location, "sample.info_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))
# 
# samples.df %>% group_by(age.x) %>% summarize(n())
# 
# pop_size.df <- readRDS(file = paste0(PopSim.location, "pop.size_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))
# 
# truth.df <- readRDS(file = paste0(PopSim.location, "truth_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule, "_", sampling.scheme))

#Assuming population simulation was just run and is still in working environment
samples.df <- sample.info
truth.df <- truth.all
pop_size.df <- pop.size.tibble_all
parents.df <- parents.tibble_all
n_yrs <- max(pop_size.df$year)

########################## MCMC & model parameters #########################
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains
jags_params = c("Nf", "Nm", "survival", "lambda") #List the parameters to be estimated
estimated.parameters <- paste0(jags_params, collapse = ",")
derived.quantities <- "no"

estimation.years <- c(n_yrs - 5, n_yrs)

########################## Start simulation loop #########################
(iterations <- max(samples.df$iteration))

# Initialize arrays for saving results
 results <- NULL
 mom.comps.tibble <- NULL
 dad.comps.tibble <- NULL

 for(est in 1:length(estimation.years)){
   #Set estimation year to 10 years in the past, or 0
   estimation.year <- estimation.years[est]
   
 for(iter in 1:iterations) {
   #  set.seed(rseeds[iter])
   sim.start <- Sys.time()
   #rseed <- sample(1:1000000,1)
   set.seed(rseeds[iter])
   rseed <- rseeds[iter]

     sample.df_all.info <- samples.df %>% dplyr::filter(iteration == iter)
     
     #Save total sample size
     sample.size.iter <- sample.df_all.info %>% 
       distinct(sample.size.total) %>%
       pull(sample.size.total)

        noDups.list <- split.dups(sample.df_all.info)
        first.capture <- noDups.list[[1]]
        later.capture <- noDups.list[[2]]
  
        filter1.out <- filter.samples(later.capture) #Filter for full sibs
 
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
    
      truth.iter <- truth.df %>% dplyr::filter(iteration == iter)

    #If running with derived quantities    
    if(derived.quantities == "yes"){
      Nfb1.temp <- truth.iter %>% dplyr::filter(parameter == "Nf") %>% 
        mutate(parameter = "Nfb1", all.truth = breed.truth)
      
      Nfb2.temp <- truth.iter %>% dplyr::filter(parameter == "Nf") %>% 
        mutate(parameter = "Nfb2", all.truth = breed.truth)
      
      truth.iter <- truth.iter %>% bind_rows(Nfb1.temp, Nfb2.temp)
    }

    #Save some values related to samples so they can be added to the dataframe of results
      samples.iter <- samples.df %>% dplyr::filter(iteration == iter) %>% 
      distinct(seed, iteration)

          
    results.temp <- model.summary2 %>% left_join(truth.iter, by = c("parameter", "iteration", "seed")) %>% 
      left_join(samples.iter, by = c("iteration", "seed"))
    
    if(HS.only == "yes"){
      if(length(jags_params) == 3){
      results.temp <- results.temp %>% 
        mutate(HSPs_detected = c(mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs),
               HSPs_expected = c(mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS),
               est.yr = estimation.year)
      } else if(length(jags_params) == 4){
        results.temp <- results.temp %>% 
          mutate(HSPs_detected = c(mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs + dad.HSPs),
                 HSPs_expected = c(mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS + dad.Exp.HS),
                 est.yr = estimation.year)
      } else if(length(jags_params) == 5){
        results.temp <- results.temp %>% 
          mutate(HSPs_detected = c(mom.HSPs, mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs + dad.HSPs),
                 HSPs_expected = c(mom.Exp.HS, mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS + dad.Exp.HS),
                 est.yr = estimation.year)
      } else if(length(jags_params) == 7){
        results.temp <- results.temp %>% 
          mutate(HSPs_detected = c(mom.HSPs, mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs, mom.HSPs),
                 HSPs_expected = c(mom.Exp.HS, mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS, mom.Exp.HS),
                 est.yr = estimation.year)
      }
      
    } else if(HS.only != "yes"){
      if(length(jags_params) == 3){
      results.temp <- results.temp %>% 
        mutate(HSPs_detected = c(mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs),
               HSPs_expected = c(mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS),
               POPs_detected = c(mom.POPs, dad.POPs, mom.POPs + dad.POPs),
               POPs_expected = c(mom.Exp.PO, dad.Exp.PO, mom.Exp.PO + dad.Exp.PO),
               est.yr = estimation.year)
      } else if(length(jags_params) == 4){
        results.temp <- results.temp %>% 
          mutate(HSPs_detected = c(mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs + dad.HSPs),
                 HSPs_expected = c(mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS + dad.Exp.HS),
                 POPs_detected = c(mom.POPs, dad.POPs, mom.POPs + dad.POPs, mom.POPs + dad.POPs),
                 POPs_expected = c(mom.Exp.PO, dad.Exp.PO, mom.Exp.PO + dad.Exp.PO, mom.Exp.PO + dad.Exp.PO),
                 est.yr = estimation.year)
      } else if(length(jags_params) == 5){
        results.temp <- results.temp %>% 
          mutate(HSPs_detected = c(mom.HSPs, mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs + dad.HSPs),
                 HSPs_expected = c(mom.Exp.HS, mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS + dad.Exp.HS),
                 POPs_detected = c(mom.POPs, mom.POPs, dad.POPs, mom.POPs + dad.POPs, mom.POPs + dad.POPs),
                 POPs_expected = c(mom.Exp.PO, mom.Exp.PO, dad.Exp.PO, mom.Exp.PO + dad.Exp.PO, mom.Exp.PO + dad.Exp.PO),
                 est.yr = estimation.year)
      } else if(length(jags_params) == 7){
        results.temp <- results.temp %>% 
          mutate(HSPs_detected = c(mom.HSPs, mom.HSPs, dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs + dad.HSPs, mom.HSPs, mom.HSPs),
                 HSPs_expected = c(mom.Exp.HS, mom.Exp.HS, dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS + dad.Exp.HS, mom.Exp.HS, mom.Exp.HS),
                 POPs_detected = c(mom.POPs, mom.POPs, dad.POPs, mom.POPs + dad.POPs, mom.POPs + dad.POPs, mom.POPs, mom.POPs),
                 POPs_expected = c(mom.Exp.PO, mom.Exp.PO, dad.Exp.PO, mom.Exp.PO + dad.Exp.PO, mom.Exp.PO + dad.Exp.PO, mom.Exp.PO, mom.Exp.PO),
                 est.yr = estimation.year)
      }
    }
    
    
    results <- rbind(results, results.temp)
    
    #Save mom and dad pairwise comparison dataframes
    mom_comps.all <- mom_comps.all %>% mutate(iteration = iter,
                                              sample.size = sample.size.iter,
                                              seed = rseed)
    mom.comps.tibble <- rbind(mom.comps.tibble, mom_comps.all)
    
    dad_comps.all <- dad_comps.all %>% mutate(iteration = iter,
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
# #    #Model output for diagnostics
#      saveRDS(sims.list.1, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.1, "_", MCMC.settings, "_", scenario, "_", model, "_", sampling.scheme))
# # 
#      saveRDS(sims.list.2, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.2, "_", MCMC.settings, "_", scenario, "_", model, "_", sampling.scheme))
# # # 
#      saveRDS(sims.list.3, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.3, "_", MCMC.settings, "_", scenario, "_", model, "_", sampling.scheme))
# # #    
#      saveRDS(sims.list.4, file = paste0(temp_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.4, "_", MCMC.settings, "_", scenario, "_", model, "_", sampling.scheme))
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
results2 %>% group_by(sample.prop, parameter) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100)

#Median relative bias by sample size
 results2 %>% group_by(sample.prop, parameter, est.yr) %>% 
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

 if(exists("age.cv") == TRUE){
   #Save samples dataframe with misassigned ages
   saveRDS(samples.miss, file = paste0(results_location, samples.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", sampling.scheme))
}
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



####################### End ##############################