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
MCMC_location <- "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.output/"
results_location <- "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/"
results_prefix <- "CKMR_results"
MCMC_prefix <- "CKMR_modelout"
jags.model.prefix <- "CKMR.JAGS_"
mom.comps.prefix <- "mom.comps"
dad.comps.prefix <- "dad.comps"
samples.prefix <- "samples/samples.missassigned"

######################### Specify other common settings #########################
PopSim.location <- "G://My Drive/Personal_Drive/R/CKMR/Population.simulations/Peer_review/" #Location of population simulation output files
temp_location <- "~/R/working_directory/temp_results/" #Location to save temporary files in case run gets killed
jags.model_location <- "G://My Drive/Personal_Drive/R/CKMR/JAGS_models/" #Location of JAGS models
today <- format(Sys.Date(), "%d%b%Y") # Store date for use in file name
date.of.simulation <- today
rseeds <- readRDS("rseeds_2022.04.15.rda")
outSeeds <- "Seeds2022.04.15"

adult.survival <- 0.825 # CHANGED FROM 0.825; Adult survival
repro.age <- 12 # set age of reproductive maturity
max.age <- 50 #set the maximum age allowed in the simulation
n_yrs <- 90 #Number of years the simulation was run for

############Specify common input prefixes####################
inSeeds <- "Seeds2022.04.15" #Seeds used for population simulation
#date.of.PopSim <- "03Aug2023" #Most common date for population simulations: 03Aug2023
date.of.PopSim <- "06Sep2023" #On 22Aug2023 I re-ran the stable population growth/annual breeding simulation, but identified aunt/niece and uncle/nephew pairs. Re-ran AGAIN on 06Sep2023 to iron out glitches with the code.

###########Specify which simulations to focus on########################
#s.scheme <- "target.YOY" #can be "target.YOY", "sample.all.juvenile.ages", or "sample.ALL.ages"
sample.props <- 1.5 #Either label this with the percent we want to target if just one (e.g., 1.5)) or if wanting to run over all sample proportions, set as "all"
objective <- 2 #Can be any number for the objectives (1-5)
scenario <- "scenario_2.2.4" #See Excel sheet with simulation scenarios: Simulation_log_key_UPDATED.xlsx on Google Drive
sample.scheme.vec <- c("target.YOY", "sample.all.juvenile.ages", "sample.ALL.ages")
#sample.scheme.vec <- c("sample.ALL.ages") #If wanting to just run one
est.yr.tests <- 4 #Can be 1 or 4. If 1, that means we will only estimate abundance for the birth year of the second oldest individual in the dataset; if 4, then we will estimate abundance for 10 years before that, the present, and five years before the present. If running with the base-case CKMR model, then should set to 1

#Assume we're not including aunt/niece pairs, but need to define the object. The specify.simulation code will adjust this setting if we are including aunt/niece pairs and if we want to filter based on age at repro.
test.decoys <- "no"

#Specify simulation details based on inputs above
source("./02_Estimation.model/functions/specify.simulation.R")

########################## Read in sampling and other dataframes #########################
#Read in sample dataframe and filter for the relevant sampling scheme
samples.df_all <- readRDS(file = paste0(PopSim.location, "sample.info_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule))

#Check that the sample dataframe was appropriately filtered
samples.df_all %>% group_by(sampling.scheme, age.x) %>% summarize(n())

#Read in dataframe with population size
pop_size.df <- readRDS(file = paste0(PopSim.location, "pop.size_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule))

pop_size.df %>% dplyr::filter(year == 85) %>% nrow() #Should be 500
n_yrs <- max(pop_size.df$year) #Save number of simulation years

pop_size.df %>% tail(15)
pop_size.df %>% tail(15) %>% dplyr::select(adult.lambda)

#Read in dataframe with true values (some of these will be recalculated later)
truth.df <- readRDS(file = paste0(PopSim.location, "truth_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule))

#Double-check that things are cool - everything should be 500
truth.df %>% group_by(sampling.scheme, sample.prop, parameter) %>% 
  summarize(n())

if(test.decoys == "yes"){

  #Read in dataframe with aunt|uncle/niece|nephew pairs masquerading as half-sibs
  imposters.df <- readRDS(file = paste0(PopSim.location, "aunt.uncle_niece.nephews_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule)) %>% 
    dplyr::distinct(.keep_all = TRUE)
  
}

########################## MCMC & model parameters #########################
ni <- 40000 # number of post-burn-in samples per chain
nb <- 50000 # number of burn-in samples
nt <- 20     # thinning rate
nc <- 2      # number of chains

(iterations <- max(samples.df_all$iteration))

if(sample.props == "all"){
#Run with all sample sizes
(sample.prop_vec <- samples.df_all %>% distinct(sample.prop) %>% pull(sample.prop))
  sim.samples.1 <- paste0(sample.prop_vec[1], "_prop.sampled")
  sim.samples.2 <- paste0(sample.prop_vec[2], "_prop.sampled")
  sim.samples.3 <- paste0(sample.prop_vec[3], "_prop.sampled")
  sim.samples.4 <- paste0(sample.prop_vec[4], "_prop.sampled")
  
  cat(paste0("Will run model with all sampling intensities"))
  
  } else if(is.numeric(sample.props) == TRUE){
#Run with just one sample size
    (sample.prop_vec <- samples.df_all %>%
       distinct(sample.prop) %>%
       dplyr::filter(sample.prop == sample.props) %>%
       pull(sample.prop)) 
    
    (sim.samples.1 <- paste0(sample.prop_vec[1], "_prop.sampled"))
    sim.samples.2 <- NULL
    sim.samples.3 <- NULL
    sim.samples.4 <- NULL
    
    cat(paste0("Will run model with sampling intensity of ", sample.props, "%", " per year."))
    }

########################## Start simulation loop #########################

# Initialize arrays for saving results
 results <- NULL
 sims.list.1 <- NULL
 sims.list.2 <- NULL
 sims.list.3 <- NULL
 sims.list.4 <- NULL
 mom.comps.tibble <- NULL
 dad.comps.tibble <- NULL

 # for(est in 1:length(estimation.years)){
 #   #Set estimation year to 10 years in the past, or 0
 #   estimation.year <- estimation.years[est]
 #   
   
 for(s.scheme in sample.scheme.vec){
   
   #Specify the model we'll be using
   jags_file <- specify.model()
   
   #Subset for samples from focal sampling scheme
   samples.df <- samples.df_all %>% dplyr::filter(sampling.scheme == s.scheme)

   HS.only <- "yes" #Specify that we'll only use half-siblings. Will change below if the sampling scheme includes adults.
   PO.only <- "no"
      
   #First iteration gives NA for reference year, but we can infer what it should be based on the sampling scheme
   if(s.scheme %in% c("sample.all.juvenile.ages", "sample.ALL.ages")){
     samples.df <- samples.df %>% replace_na(list(ref.yr = 76))
   } else {
    samples.df <- samples.df %>% replace_na(list(ref.yr = 87))
   }
   
   #Change a couple of settings if sampling adults
   if(s.scheme == "sample.ALL.ages"){
     HS.only <- "no"
     if(model == "multiennial.model"){
       jags_file <- paste0(jags.model_location, "HS.PO_narrowLambda_Skip_model.txt") #We are not testing the multiennial model with a changing population so we do not need to have a scenario where we use a wide lambda prior.
     }
   }
   
   #Should just be one sampling scheme represented
   samples.df %>% dplyr::count(sampling.scheme)  
   
 for(iter in 1:iterations) {
   #  set.seed(rseeds[iter])
   sim.start <- Sys.time()
   #rseed <- sample(1:1000000,1)
   set.seed(rseeds[iter])
   rseed <- rseeds[iter]
   
   for(s in 1:length(sample.prop_vec)){
     
     sample.proportion <- sample.prop_vec[s]
     sample.df_all.info <- samples.df %>% dplyr::filter(iteration == iter, sample.prop == sample.proportion)
     
     #Save estimation year
     est.year.calibrate <- sample.df_all.info %>% distinct(ref.yr) %>% 
       mutate(estimation.yr = ref.yr + 1) %>% 
       pull(estimation.yr)
     
     #Make vector of estimation years: t0, t0-10, five years into the past, present
     estimation.years <- c(est.year.calibrate, est.year.calibrate - 10, n_yrs - 5, n_yrs)
     
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
        
 
    PO.samps.list <- filter1.out[[1]] #Output is a list where each list element corresponds to the offspring birth year and contains the potential parents and offspring for that year.
    HS.samps.df <- filter1.out[[2]] #Output is just the dataframe of samples but filtered for full siblings i.e. kept one of the two
    NoFullSibs.df <- filter1.out[[3]] #Save NoFullSibs.df so can downsample if desired
    full.sibs <- filter1.out[[4]] #Number of full siblings
    diff.cohort.sibs <- filter1.out[[5]] #Number of full siblings from different cohorts
    
    
    #Loop over estimation years
    for(est in 1:est.yr.tests){
      estimation.year <- estimation.years[est]
      
      cat(paste0("Working on iteration ", iter, " for sampling scheme ", s.scheme, " with sampling intensity ", sample.proportion, " and estimation year ", estimation.year, "\n"))
    
    
    #-------------Construct pairwise comparison matrix--------------
    #Input is 1) a list of potential parents and offspring for each offspring birth year, and 2) a dataframe with all samples, filtered to remove full siblings
      #ADDED inclusion of aunt|uncle and niece|nephew pairs 09/08/2023 to the build.pairwise function. Will only run if test.decoys = "yes"
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

########### Calculate truth ###################
      #Calculate truth for lambda
      lambda.pop.df <- pop_size.df %>% dplyr::filter(iteration == iter)
      
      lambda.truth <- (lambda.pop.df$Total.adult.pop[n_yrs]/lambda.pop.df$Total.adult.pop[est.year.calibrate])^(1/(n_yrs - est.year.calibrate))
      
      lambda.truth.df <- tibble(estimation.year = estimation.year,
                                iteration = iter,
                                seed = rseed,
                                parameter = "lambda",
                                all.truth = lambda.truth)
      
      
      #Make dataframe of truth
      truth.iter <- pop_size.df %>% dplyr::filter(iteration == iter,
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
        mutate(sampling.scheme = s.scheme, 
               sample.prop = sample.proportion,
               T0 = est.year.calibrate)
      
    #Calculate expectations
      #Define Nf and Nm as objects for calculations
      Nf.truth <- truth.iter %>% dplyr::filter(parameter == "Nf") %>% pull(all.truth)
      
      Nm.truth <- truth.iter %>% dplyr::filter(parameter == "Nm") %>% pull(all.truth)
      
      #Subset for just the present iteration
    pop.size.tibble <- pop_size.df %>% dplyr::filter(iteration == iter)
    Exp <- calc.Exp(mom_comps.all, dad_comps.all)
    mom.Exp.HS <- Exp[[1]]
    mom.Exp.PO <- Exp[[2]]
    dad.Exp.HS <- Exp[[3]]
    dad.Exp.PO <- Exp[[4]]
    
    sampled.mothers <- unique(sample.df_all.info$mother.x)
    sampled.fathers <- unique(sample.df_all.info$father.x)
    
    #Calculate truth
    results.temp <- calc.truth(truth.df = truth.df)
    
    
    results <- rbind(results, results.temp)
    
    #Save mom and dad pairwise comparison dataframes
    mom_comps.all <- mom_comps.all %>% mutate(iteration = iter,
                                              sampling.scheme = s.scheme,
                                              sample.prop = sample.proportion,
                                              sample.size = sample.size.iter,
                                              seed = rseed)
    mom.comps.tibble <- rbind(mom.comps.tibble, mom_comps.all)
    
    dad_comps.all <- dad_comps.all %>% mutate(iteration = iter,
                                              sampling.scheme = s.scheme,
                                              sample.prop = sample.proportion,
                                              sample.size = sample.size.iter,
                                              seed = rseed)
    dad.comps.tibble <- rbind(dad.comps.tibble, dad_comps.all)
    
  } # End if/else statement to account for no yes comparisons
  } #end loop over estimation years
   } #end loop over sample sizes
    
  #-----------------Save output files iteratively--------------------
   #if(iter %% 100 == 0){
     
#Results
    write.table(results, file = paste0(temp_location, results_prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", s.scheme, ".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
# 
#    #Model output for diagnostics
     saveRDS(sims.list.1, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.1, "_", MCMC.settings, "_", scenario, "_", model, "_", s.scheme))
# 
     saveRDS(sims.list.2, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.2, "_", MCMC.settings, "_", scenario, "_", model, "_", s.scheme))
# # 
     saveRDS(sims.list.3, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.3, "_", MCMC.settings, "_", scenario, "_", model, "_", s.scheme))
# #    
     saveRDS(sims.list.4, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.4, "_", MCMC.settings, "_", scenario, "_", model, "_", s.scheme))

#    #Save pairwise comparisons matrices
    saveRDS(mom.comps.tibble, file = paste0(temp_location, mom.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", s.scheme))
#    
    saveRDS(dad.comps.tibble, file = paste0(temp_location, dad.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, "_", s.scheme))

#   }
    
      sim.end <- Sys.time()
   
   iter.time <- round(as.numeric(difftime(sim.end, sim.start, units = "mins")), 1)
   cat(paste0("\n Finished iteration ", iter, " with sampling scheme: ", s.scheme, ". \n Took ", iter.time, " minutes \n\n"))
   } # end loop over sampling schemes
 }# end loop over iterations
  

    #If using all individuals for Nf truth, instead of breeders
   results2 <- results %>%
     mutate(relative_bias = round(((Q50 - all.truth)/all.truth)*100, 1)) %>% #Can change truth to breed.truth if looking for number of active breeders
     mutate(in_interval = ifelse(HPD2.5 < all.truth & all.truth < HPD97.5, "Y", "N")) %>%
     as_tibble()
   
#Within HPD interval?
results2 %>% group_by(sample.prop, parameter, estimation.year, sampling.scheme) %>% 
  dplyr::summarize(percent_in_interval = sum(in_interval == "Y")/n() * 100)

#Median relative bias by sample size
 results2 %>% group_by(sampling.scheme, sample.prop, parameter, estimation.year) %>% 
   dplyr::summarize(median.bias = median(relative_bias), n = n()) %>% 
   dplyr::arrange(desc(median.bias))


 #-----------------------------Save major output files---------------------------------------------
 #Home computer: Dell Precision
 
 #Save model estimates
write.table(results2, file = paste0(results_location, results_prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model, ".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
 
 #Save draws from posterior for model diagnostics 
 saveRDS(sims.list.1, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.1, "_", MCMC.settings, "_", scenario, "_", model)) #Sample size 1
 
  saveRDS(sims.list.2, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.2, "_", MCMC.settings, "_", scenario, "_", model)) #Sample size 2
 # 
  saveRDS(sims.list.3, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.3, "_", MCMC.settings, "_", scenario, "_", model)) #Sample size 3
 # 
  saveRDS(sims.list.4, file = paste0(MCMC_location, MCMC_prefix, "_", date.of.simulation, "_", outSeeds, "_", sim.samples.4, "_", MCMC.settings, "_", scenario, "_", model)) #Sample size 4
 
 #Save final pairwise comparison matrices
 saveRDS(mom.comps.tibble, file = paste0(results_location, mom.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model))
 
 saveRDS(dad.comps.tibble, file = paste0(results_location, dad.comps.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model))

 if(exists("age.cv") == TRUE){
   #Save samples dataframe with misassigned ages
   saveRDS(samples.miss, file = paste0(results_location, samples.prefix, "_", date.of.simulation, "_", outSeeds, "_", scenario, "_", model))
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