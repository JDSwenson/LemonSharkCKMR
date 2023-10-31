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
date.of.PopSim <- "03Aug2023" #Most common date for population simulations: 03Aug2023
#date.of.PopSim <- "06Sep2023" #On 22Aug2023 I re-ran the stable population growth/annual breeding simulation, but identified aunt/niece and uncle/nephew pairs. Re-ran AGAIN on 06Sep2023 to iron out glitches with the code.

###########Specify which simulations to focus on########################
#s.scheme <- "target.YOY" #can be "target.YOY", "sample.all.juvenile.ages", or "sample.ALL.ages"
sample.props <- "all" #Either label this with the percent we want to target if just one (e.g., 1.5)) or if wanting to run over all sample proportions, set as "all"
objective <- 2 #Can be any number for the objectives (1-5)
sample.scheme.vec <- c("sample.all.juvenile.ages") #If wanting to just run one
est.yr.tests <- 1 #Can be 1 or 4. If 1, that means we will only estimate abundance for the birth year of the second oldest individual in the dataset; if 4, then we will estimate abundance for 10 years before that, the present, and five years before the present. If running with the base-case CKMR model, then should set to 1

#Assume we're not including aunt/niece pairs, but need to define the object. The specify.simulation code will adjust this setting if we are including aunt/niece pairs.
test.decoys <- "no"


########################## Read in sampling and other dataframes #########################
#Read in sample dataframe and filter for the relevant sampling scheme
#samples.df_all <- readRDS(file = paste0(PopSim.location, "sample.info_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule))

#Check that the sample dataframe was appropriately filtered
#samples.df_all %>% group_by(sampling.scheme, age.x) %>% summarize(n())

#Read in dataframe with population size
scenario <- "scenario_2.2.2"

#Specify simulation details based on inputs above
source("./02_Estimation.model/functions/specify.simulation.R")

pop_size_increase.df <- readRDS(file = paste0(PopSim.location, "pop.size_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule)) %>% 
  mutate(population.growth = "slightly increasing")

pop_size_increase.df %>% dplyr::filter(year == 85) %>% nrow() #Should be 500
n_yrs <- max(pop_size_increase.df$year) #Save number of simulation years

pop_size_increase.df %>% tail(15)
pop_size_increase.df %>% tail(15) %>% dplyr::select(adult.lambda)


#Read in dataframe with population size
scenario <- "scenario_2.2.1" 

#Specify simulation details based on inputs above
source("./02_Estimation.model/functions/specify.simulation.R")

pop_size_decrease.df <- readRDS(file = paste0(PopSim.location, "pop.size_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule)) %>% 
  mutate(population.growth = "slightly decreasing")

pop_size_decrease.df %>% dplyr::filter(year == 85) %>% nrow() #Should be 500
n_yrs <- max(pop_size_decrease.df$year) #Save number of simulation years

pop_size_decrease.df %>% tail(15)
pop_size_decrease.df %>% tail(15) %>% dplyr::select(adult.lambda)

#Read in dataframe with population size
scenario <- "scenario_2.2.3" 

#Specify simulation details based on inputs above
source("./02_Estimation.model/functions/specify.simulation.R")

pop_size_severe.df <- readRDS(file = paste0(PopSim.location, "pop.size_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule)) %>% 
  mutate(population.growth = "severe decline")

pop_size_severe.df %>% dplyr::filter(year == 85) %>% nrow() #Should be 500
n_yrs <- max(pop_size_severe.df$year) #Save number of simulation years


#Read in dataframe with population size
scenario <- "scenario_2.2.4"

#Specify simulation details based on inputs above
source("./02_Estimation.model/functions/specify.simulation.R")

pop_size_stable.df <- readRDS(file = paste0(PopSim.location, "pop.size_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule)) %>% 
  mutate(population.growth = "stable")

pop_size_stable.df %>% dplyr::filter(year == 85) %>% nrow() #Should be 500
n_yrs <- max(pop_size_stable.df$year) #Save number of simulation years



pop.size_all <- bind_rows(pop_size_increase.df, pop_size_decrease.df, pop_size_severe.df, pop_size_stable.df)

pop.size_all %>% write_rds(file = "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/pop.size_all")

#Read in dataframe with true values (some of these will be recalculated later)
truth.df <- readRDS(file = paste0(PopSim.location, "truth_", date.of.PopSim, "_", inSeeds, "_", PopSim.lambda, "_", PopSim.breeding.schedule))

#Double-check that things are cool - everything should be 500
truth.df %>% group_by(sampling.scheme, sample.prop, parameter) %>% 
  summarize(n())




#######Compare this for increasing and declining populations###############
truth.df %>% dplyr::filter(sampling.scheme == "sample.all.juvenile.ages", sample.prop == 1.5) %>% dplyr::group_by(parameter) %>% summarize(sd = sd(all.truth), mean = mean(all.truth), min = min(surv_min), max = max(surv_max))





#Run the below twice: once after loading files from population increase, and once after loading files from population decrease

######################### Lambda comparison #################################
#--------------------------Compare over 90 years-------------------------
yr0 <- 1
yrs <- 90

#Initialize df for saving values
N90_estimates <- NULL

for(i in 1:max(pop_size.df$iteration)){
  #Subset for iteration of interest
  iter <- pop_size.df %>% dplyr::filter(iteration == i)

    #Change calculation of lambda
  mean.lam_yr1 <- (iter$Total.adult.pop[yrs]/iter$Total.adult.pop[1])^(1/(yrs - yr0))
  
  indv.lambdas <- iter$adult.lambda
  
  #Initialize vectors
  indv.method_yr1 = avg.method_yr1 <- NULL
  
  for(j in 2:length(indv.lambdas)){
    
    #Store N0 and save as first element of each vector
    N0 <- iter$Total.adult.pop[1]
    indv.method_yr1[1] = avg.method_yr1[1] <- N0
    
    #Save vector of individual lambda values
    indv.lambdas_yr1 <- iter$adult.lambda
    
    #Calculate mean lambda -- removed after Liz suggested a different approach
#    mean.lam_yr1 <- mean(iter$adult.lambda, na.rm = TRUE)
    
    #Calculate abundance
    indv.method_yr1[j] <- round(indv.method_yr1[j-1] * indv.lambdas_yr1[j], 0)
    avg.method_yr1[j] <- round(avg.method_yr1[j-1] * mean.lam_yr1, 0)
    
  }
  
  combined.df <- tibble(iter$Total.adult.pop[90], indv.method_yr1[90], avg.method_yr1[90]) %>% mutate(iteration = iter$iteration[j],
                                                                                                          yr0 = yr0,
                                                                                                          lambda.yrs = yrs-yr0)
  
  N90_estimates <- bind_rows(N90_estimates, combined.df)

  }
  
colnames(N90_estimates)[1:3] <- c("true.adult.pop", "indv.method.pop", "avg.method.pop")

N90_estimates



#-----------------------Compare over 40 years-------------------------
yr0 <- 50
yrs <- 90

#Initialize df for saving values
N50_estimates <- NULL

for(i in 1:max(pop_size.df$iteration)){
  #Subset for iteration of interest
  iter <- pop_size.df %>% dplyr::filter(iteration == i)
  
  #Change calculation of lambda
  mean.lam_yr50 <- (iter$Total.adult.pop[yrs]/iter$Total.adult.pop[yr0])^(1/(yrs - yr0))
  
  indv.lambdas <- iter$adult.lambda
  
  #Initialize vectors
  indv.method_yr50.vec = avg.method_yr50.vec <- NULL
  
  #Store N0 and save as first element of each vector
  N0 <- iter$Total.adult.pop[yr0]
  indv.method_yr50.vec[1] = avg.method_yr50.vec[1] <- N0
  
  #Save vector of individual lambda values
  indv.lambdas_yr50 <- iter$adult.lambda

  indx <- 1
    
  for(j in yr0:yrs){
    #Calculate mean lambda
    #mean.lam_yr50 <- mean(iter$adult.lambda, na.rm = TRUE)
    
    indx <- indx + 1
    
    #Calculate abundance
    indv.method_yr50.vec[indx] <- round(indv.method_yr50.vec[indx-1] * indv.lambdas_yr50[j+1], 0)

    avg.method_yr50.vec[indx] <- round(avg.method_yr50.vec[indx-1] * mean.lam_yr50, 0)

  }
  
  combined.df <- tibble(iter$Total.adult.pop[yrs], indv.method_yr50.vec[yrs-yr0+1], avg.method_yr50.vec[yrs-yr0+1]) %>% mutate(iteration = iter$iteration[j],
                                                                                                          yr0 = yr0,
                                                                                                          lambda.yrs = yrs-yr0)
  
  N50_estimates <- bind_rows(N50_estimates, combined.df)
  
}

colnames(N50_estimates)[1:3] <- c("true.adult.pop", "indv.method.pop", "avg.method.pop")

N50_estimates

#-----------------------Compare over 20 years-------------------------
yr0 <- 70
yrs <- 90

#Initialize df for saving values
N70_estimates <- NULL

for(i in 1:max(pop_size.df$iteration)){
  #Subset for iteration of interest
  iter <- pop_size.df %>% dplyr::filter(iteration == i)
  
  #Change calculation of lambda
  mean.lam_yr70 <- (iter$Total.adult.pop[yrs]/iter$Total.adult.pop[yr0])^(1/(yrs - yr0))
  
  indv.lambdas <- iter$adult.lambda
  
  #Initialize vectors
  indv.method_yr70.vec = avg.method_yr70.vec <- NULL
  
  #Store N0 and save as first element of each vector
  N0 <- iter$Total.adult.pop[yr0]
  indv.method_yr70.vec[1] = avg.method_yr70.vec[1] <- N0
  
  #Save vector of individual lambda values
  indv.lambdas_yr70 <- iter$adult.lambda
  
  indx <- 1
  
  for(j in yr0:yrs){
    #Calculate mean lambda
    #mean.lam_yr70 <- mean(iter$adult.lambda, na.rm = TRUE)
    
    indx <- indx + 1
    
    #Calculate abundance
    indv.method_yr70.vec[indx] <- round(indv.method_yr70.vec[indx-1] * indv.lambdas_yr70[j+1], 0)
    
    avg.method_yr70.vec[indx] <- round(avg.method_yr70.vec[indx-1] * mean.lam_yr70, 0)
    
  }
  
  combined.df <- tibble(iter$Total.adult.pop[yrs], indv.method_yr70.vec[yrs-yr0+1], avg.method_yr70.vec[yrs-yr0+1]) %>% mutate(iteration = iter$iteration[j],
                                                                                                                               yr0 = yr0,
                                                                                                                               lambda.yrs = yrs-yr0)
  
  N70_estimates <- bind_rows(N70_estimates, combined.df)
  
}

colnames(N70_estimates)[1:3] <- c("true.adult.pop", "indv.method.pop", "avg.method.pop")

N70_estimates


#-----------------------Compare over 10 years-------------------------
yr0 <- 80
yrs <- 90

#Initialize df for saving values
N80_estimates <- NULL

for(i in 1:max(pop_size.df$iteration)){
  #Subset for iteration of interest
  iter <- pop_size.df %>% dplyr::filter(iteration == i)
  
  #Change calculation of lambda
  mean.lam_yr80 <- (iter$Total.adult.pop[yrs]/iter$Total.adult.pop[yr0])^(1/(yrs - yr0))
  
  indv.lambdas <- iter$adult.lambda
  
  #Initialize vectors
  indv.method_yr80.vec = avg.method_yr80.vec <- NULL
  
  #Store N0 and save as first element of each vector
  N0 <- iter$Total.adult.pop[yr0]
  indv.method_yr80.vec[1] = avg.method_yr80.vec[1] <- N0
  
  #Save vector of individual lambda values
  indv.lambdas_yr80 <- iter$adult.lambda
  
  indx <- 1
  
  for(j in yr0:yrs){
    #Calculate mean lambda
    #mean.lam_yr80 <- mean(iter$adult.lambda, na.rm = TRUE)
    
    indx <- indx + 1
    
    #Calculate abundance
    indv.method_yr80.vec[indx] <- round(indv.method_yr80.vec[indx-1] * indv.lambdas_yr80[j+1], 0)
    
    avg.method_yr80.vec[indx] <- round(avg.method_yr80.vec[indx-1] * mean.lam_yr80, 0)
    
  }
  
  combined.df <- tibble(iter$Total.adult.pop[yrs], indv.method_yr80.vec[yrs-yr0+1], avg.method_yr80.vec[yrs-yr0+1]) %>% mutate(iteration = iter$iteration[j],
                                                                                                                               yr0 = yr0,
                                                                                                                               lambda.yrs = yrs-yr0)
  
  N80_estimates <- bind_rows(N80_estimates, combined.df)
  
}

colnames(N80_estimates)[1:3] <- c("true.adult.pop", "indv.method.pop", "avg.method.pop")

N80_estimates


#-----------------------Compare over 5 years-------------------------
yr0 <- 85
yrs <- 90

#Initialize df for saving values
N85_estimates <- NULL

for(i in 1:max(pop_size.df$iteration)){
  #Subset for iteration of interest
  iter <- pop_size.df %>% dplyr::filter(iteration == i)
  
  #Change calculation of lambda
  mean.lam_yr85 <- (iter$Total.adult.pop[yrs]/iter$Total.adult.pop[yr0])^(1/(yrs - yr0))
  
  indv.lambdas <- iter$adult.lambda
  
  #Initialize vectors
  indv.method_yr85.vec = avg.method_yr85.vec <- NULL
  
  #Store N0 and save as first element of each vector
  N0 <- iter$Total.adult.pop[yr0]
  indv.method_yr85.vec[1] = avg.method_yr85.vec[1] <- N0
  
  #Save vector of individual lambda values
  indv.lambdas_yr85 <- iter$adult.lambda
  
  indx <- 1
  
  for(j in yr0:yrs){
    #Calculate mean lambda
    #mean.lam_yr85 <- mean(iter$adult.lambda, na.rm = TRUE)
    
    indx <- indx + 1
    
    #Calculate abundance
    indv.method_yr85.vec[indx] <- round(indv.method_yr85.vec[indx-1] * indv.lambdas_yr85[j+1], 0)
    
    avg.method_yr85.vec[indx] <- round(avg.method_yr85.vec[indx-1] * mean.lam_yr85, 0)
    
  }
  
  combined.df <- tibble(iter$Total.adult.pop[yrs], indv.method_yr85.vec[yrs-yr0+1], avg.method_yr85.vec[yrs-yr0+1]) %>% mutate(iteration = iter$iteration[j],
                                                                                                                               yr0 = yr0,
                                                                                                                               lambda.yrs = yrs-yr0)
  
  N85_estimates <- bind_rows(N85_estimates, combined.df)
  
}

colnames(N85_estimates)[1:3] <- c("true.adult.pop", "indv.method.pop", "avg.method.pop")

N85_estimates


#Combine positive and negative lambda growth
#lambda.df.decrease <- bind_rows(N90_estimates, N50_estimates, N70_estimates, N80_estimates, N85_estimates)
lambda.df.increase <- bind_rows(N90_estimates, N50_estimates, N70_estimates, N80_estimates, N85_estimates)

lambda.df.decrease <- lambda.df.decrease %>% mutate(pop.growth = "slight negative")
lambda.df.increase <- lambda.df.increase %>% mutate(pop.growth = "slight positive")

lambda.df <- bind_rows(lambda.df.increase, lambda.df.decrease) 

lambda.df <- lambda.df %>% mutate(relative_bias = ((avg.method.pop - true.adult.pop)/true.adult.pop)*100)

lambda.df %>% group_by(pop.growth, lambda.yrs) %>% 
  summarize(median(relative_bias))

lambda.df %>% write_rds(file = "G://My Drive/Personal_Drive/R/CKMR/output_peer_review/Model.results/lambda.df")


lambda.df %>% 
  ggplot(aes(x = factor(lambda.yrs), y = relative_bias, fill = factor(pop.growth))) +
# geom_violin(draw_quantiles = 0.5, position = position_dodge(0.75), width = 4) +
  #geom_violin(draw_quantiles = 0.5, width = 4) +
  geom_jitter(width = 0.2, aes(col = factor(pop.growth)), alpha = 0.5) +
  geom_boxplot() +
#  geom_hline(yintercept=0, col="black", size=1.25) +
  labs(x = "lambda years", y = "relative bias") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=25),
        legend.title = element_blank())

#Scatter plot of exact relative bias
lambda.df %>% 
  ggplot(aes(x = factor(pop.growth), y = relative_bias, fill = factor(pop.growth))) + 
  geom_jitter() +
  facet_wrap(~lambda.yrs)


#-------------Compare total lambda vs adult lambda----------------
yr0 <- 1
yrs <- 90

pop_size.df <- pop_size.df %>% mutate(juv.pop.size = population_size - Total.adult.pop)

#Initialize df for saving values
adult.lam_90 = total.lam_90 = juv.lam_90 <- NULL

for(i in 1:max(pop_size.df$iteration)){
  #Subset for iteration of interest
  iter <- pop_size.df %>% dplyr::filter(iteration == i)
  
  #Change calculation of lambda
  mean.adult.lam_90 <- (iter$Total.adult.pop[yrs]/iter$Total.adult.pop[1])^(1/(yrs - yr0))
  adult.lam_90 <- c(adult.lam_90, mean.adult.lam_90)

  mean.total.lam_90 <- (iter$population_size[yrs]/iter$population_size[1])^(1/(yrs - yr0))
  total.lam_90 <- c(total.lam_90, mean.total.lam_90)
  
  mean.juv.lam_90 <- (iter$juv.pop.size[yrs]/iter$juv.pop.size[1])^(1/(yrs - yr0))
  juv.lam_90 <- c(juv.lam_90, mean.juv.lam_90)
  
}

#Compare lambda for difference age classes
N90_lambda <- data.frame(total.lam_90, adult.lam_90, juv.lam_90) %>% 
  mutate(total_v_adult = total.lam_90 - adult.lam_90,
         total_v_juv = total.lam_90 - juv.lam_90,
         adult_v_juv = adult.lam_90 - juv.lam_90)

#Scatter plot
N90_lambda %>% pivot_longer(cols = c(total_v_adult, total_v_juv, adult_v_juv),
                            names_to = c("comparison"),
                            values_to = "difference") %>% 
  ggplot(aes(x = difference, y = factor(comparison), fill = factor(comparison))) +
  geom_jitter(aes(col = factor(comparison))) +
  geom_boxplot(alpha = 0.3)