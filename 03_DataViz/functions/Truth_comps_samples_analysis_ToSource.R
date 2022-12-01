library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(coda)
library(runjags)
library(ggmcmc)
library(viridis)
library(scales)
library(ggridges)

#rm(list=ls())
#today <- format(Sys.Date(), "%d%b%Y")

#----------------Set input file locations - same as truth_calcs script ------------------------------
#Population simulation files
PopSim.location <- "G://My Drive/Personal_Drive/R/CKMR/Population.simulations/"
PopSim.lambda.1 <- "lambda.1" # Can be lambda.1 or lambda.variable
PopSim.lambda.variable <- "lambda.variable"
PopSim.lambda.extreme <- "lambda.extreme"
PopSim.annual.breeding <- "annual.breeding" #Can be annual.breeding or biennial.breeding
PopSim.biennial.breeding <- "biennial.breeding_NoNonConform"
Sampling.scheme.YOY <- "target.YOY" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
Sampling.scheme.juvs <- "sample.all.juvenile.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
Sampling.scheme.all <- "sample.ALL.ages" # Can be sample.all.juvenile.ages, target.YOY, or sample.ALL.ages
lambda_date.1 <- "19Jul2022" #Annual breeding; all ages and target YOY, both lambda.1 and lambda.variable
lambda_date.2 <- "28Jul2022" #Biennial breeding w/o nonconformists; all ages and target YOY
lambda_date.3 <- "08Aug2022" #Annual and biennial breeding w/o nonconformists: juvenile ages - lambda.1, lambda.variable, lambda.extreme
lambda_date.4 <- "24Jul2022" #Annual breeding; all ages and target YOY, lambda.extreme

inSeeds <- "Seeds2022.04.15"

date.of.comps_2.2.2 <- "16Aug2022"

#----------------------Objective 2----------------------------#
#Confirmed that pairwise comparison matrices are the same with scenario 2.2.2 and 2.3.2
#YOY
mom.comps_lambda.extreme_YOY <- readRDS(file = paste0(objective_2_results_location, "scenario_2.2.2/", mom.comps.prefix, "_", date.of.comps_2.2.2, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_target.YOY_change.est.yr")) %>% 
  dplyr::select(-BI)


dad.comps_lambda.extreme_YOY <- readRDS(file = paste0(objective_2_results_location, "scenario_2.2.2/", dad.comps.prefix, "_", date.of.comps_2.2.2, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_target.YOY_change.est.yr"))

all.comps_lambda.extreme_YOY <- rbind(mom.comps_lambda.extreme_YOY, dad.comps_lambda.extreme_YOY) %>% 
  mutate(sampling.scheme = "target YOY") %>% 
  dplyr::filter(type == "HS")

#Juveniles
mom.comps_lambda.extreme_juvs <- readRDS(file = paste0(objective_2_results_location, "scenario_2.2.2/", mom.comps.prefix, "_", date.of.comps_2.2.2, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_sample.all.juvenile.ages_change.est.yr")) %>% 
  dplyr::select(-BI)

dad.comps_lambda.extreme_juvs <- readRDS(file = paste0(objective_2_results_location, "scenario_2.2.2/", dad.comps.prefix, "_", date.of.comps_2.2.2, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_sample.all.juvenile.ages_change.est.yr"))

all.comps_lambda.extreme_juvs <- rbind(mom.comps_lambda.extreme_juvs, dad.comps_lambda.extreme_juvs) %>% 
  mutate(sampling.scheme = "sample all juvenile ages") %>% 
  dplyr::filter(type == "HS")


#All age classes
mom.comps_lambda.extreme_all.ages <- readRDS(file = paste0(objective_2_results_location, "scenario_2.2.2/", mom.comps.prefix, "_", date.of.comps_2.2.2, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_sample.ALL.ages_change.est.yr")) %>% 
  dplyr::select(-BI)

dad.comps_lambda.extreme_all.ages <- readRDS(file = paste0(objective_2_results_location, "scenario_2.2.2/", dad.comps.prefix, "_", date.of.comps_2.2.2, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_sample.ALL.ages_change.est.yr"))

all.comps_lambda.extreme_all.ages <- rbind(mom.comps_lambda.extreme_all.ages, dad.comps_lambda.extreme_all.ages) %>% 
  mutate(sampling.scheme = "sample all age classes")


#All
lambda.extreme_comps.all_temp <- bind_rows(all.comps_lambda.extreme_YOY, all.comps_lambda.extreme_juvs, all.comps_lambda.extreme_all.ages) %>% 
  mutate(population.growth = factor("severe decline"),
         prop.sampled = factor(sample.prop.all),
         sampling.scheme = factor(sampling.scheme))
  
obj2_results4Join <- obj2_results %>% dplyr::filter(population.growth == "severe decline") %>% 
  dplyr::select(parameter, sampling.scheme, population.growth, iteration, est.yr, est.yr.lab, prop.sampled, relative_bias, in_interval)

#Each estimation year has the same seed, so the comparisons should all be the same
lambda.extreme_comps.all <- obj2_results4Join %>% inner_join(lambda.extreme_comps.all_temp, by = c("sampling.scheme", "population.growth", "iteration", "est.yr", "prop.sampled")) %>% 
  mutate(ref.yrs = ref.year - est.yr) #CHANGED FROM abs(ref.year - est.yr)

###Summarize and make dataframes for plotting
lambda.extreme_comps4viz <- lambda.extreme_comps.all %>% dplyr::filter(yes != 0) %>% 
  uncount(yes) %>% 
  mutate(est.yr.lab = factor(est.yr.lab, levels = c("present", "5 years past", "10 years past")))


#See how frequently lambda is estimated to be negative
#Had to do some weird summarize/mutate wizardry to get the dataframe formatted correctly. Doesn't like ifelse statements in summarize.
obj2_lam.results <- obj2_results %>% dplyr::filter(population.growth != "stable", 
                               parameter == "lambda") %>% 
  group_by(lambda.prior, population.growth, sampling.scheme, est.yr) %>% 
  dplyr::summarize(percent.correct = round(sum(mean < 1)/n() * 100, 0)) %>% #This will be accurate for negative lambda but not positive
  mutate(percent.correct = ifelse(population.growth == "slight increase", 100-percent.correct, percent.correct)) %>% #positive lambda is the inverse of the percent.correct column
  mutate(percent.wrong = 100 - percent.correct) %>% 
  dplyr::filter(lambda.prior == "diffuse (0.80 - 1.20)")

obj2.temp1 <- obj2_lam.results %>% dplyr::select(-c(percent.wrong)) %>% 
  dplyr::rename(percent = percent.correct) %>% 
  mutate(status = factor("correct")) %>% 
  mutate(degree = ifelse(population.growth == "slight decline" | population.growth == "slight increase", "slight change", "substantial decline")) %>% 
  group_by(sampling.scheme, est.yr, status, degree) %>% 
  dplyr::summarize(percent.summ = round(mean(percent)))

obj2.temp2 <- obj2_lam.results %>% dplyr::select(-c(percent.correct)) %>% 
  dplyr::rename(percent = percent.wrong) %>% 
  mutate(status = factor("incorrect")) %>% 
  mutate(degree = ifelse(population.growth == "slight decline" | population.growth == "slight positive", "slight change", "substantial decline")) %>% 
  group_by(sampling.scheme, est.yr, status, degree) %>% 
  dplyr::summarize(percent.summ = round(mean(percent)))

obj2_lam.results_4viz <- obj2.temp1 %>% bind_rows(obj2.temp2) %>% 
  mutate(status = factor(status, levels = c("incorrect", "correct"))) #Reverse the order for better viz with the bar plot
