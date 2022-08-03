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

#----------------------Objective 2----------------------------#
#Confirmed that the comparisons are the same for scenarios 1.1 and 1.2
#Need to calculate mean truth over years that comparisons span
#YOY
scenario_2.2.2_mom.comps.YOY <- readRDS(file = paste0(objective_2_results_location, mom.comps.prefix, "_", date.of.simulation.scenario_2.2Y, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_target.YOY_change.est.yr")) %>% 
  dplyr::select(-BI)

scenario_2.2.2_dad.comps.YOY <- readRDS(file = paste0(objective_2_results_location, dad.comps.prefix, "_", date.of.simulation.scenario_2.2Y, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_target.YOY_change.est.yr"))

scenario_2.2.2_all.comps.YOY <- rbind(scenario_2.2.2_mom.comps.YOY, scenario_2.2.2_dad.comps.YOY) %>% 
  mutate(sampling.scheme = "target YOY") %>% 
  dplyr::filter(type == "HS")

#Juveniles
scenario_2.2.2_mom.comps.juvs <- readRDS(file = paste0(objective_2_results_location, mom.comps.prefix, "_", date.of.simulation.scenario_2.2J, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_sample.all.juvenile.ages_change.est.yr")) %>% 
  dplyr::select(-BI)

scenario_2.2.2_dad.comps.juvs <- readRDS(file = paste0(objective_2_results_location, dad.comps.prefix, "_", date.of.simulation.scenario_2.2J, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_sample.all.juvenile.ages_change.est.yr"))

scenario_2.2.2_all.comps.juvs <- rbind(scenario_2.2.2_mom.comps.juvs, scenario_2.2.2_dad.comps.juvs) %>%
  mutate(sampling.scheme = "sample all juvenile ages")
  

#All age classes
scenario_2.2.2_mom.comps.all <- readRDS(file = paste0(objective_2_results_location, mom.comps.prefix, "_", date.of.simulation.scenario_2.2A, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_sample.ALL.ages_change.est.yr")) %>% 
  dplyr::select(-BI)

scenario_2.2.2_dad.comps.all <- readRDS(file = paste0(objective_2_results_location, dad.comps.prefix, "_", date.of.simulation.scenario_2.2A, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_sample.ALL.ages_change.est.yr"))

scenario_2.2.2_all.comps.all <- rbind(scenario_2.2.2_mom.comps.all, scenario_2.2.2_dad.comps.all) %>% 
  mutate(sampling.scheme = "sample all age classes")

#All
scenario_2.2.2_comps.all_temp <- bind_rows(scenario_2.2.2_all.comps.YOY, scenario_2.2.2_all.comps.juvs, scenario_2.2.2_all.comps.all) %>% 
  mutate(population.growth = factor("extreme negative"),
         prop.sampled = factor(sample.prop.all),
         sampling.scheme = factor(sampling.scheme))
  
obj2_results4Join <- obj2_results %>% dplyr::filter(population.growth == "extreme negative") %>% 
  dplyr::select(parameter, sampling.scheme, population.growth, iteration, est.yr, est.yr.lab, prop.sampled, relative_bias, in_interval)

#Each estimation year has the same seed, so the comparisons should all be the same
scenario_2.2.2_comps.all <- obj2_results4Join %>% inner_join(scenario_2.2.2_comps.all_temp, by = c("sampling.scheme", "population.growth", "iteration", "est.yr", "prop.sampled")) %>% 
  mutate(ref.yrs = abs(ref.year - est.yr))

###Summarize and make dataframes for plotting
scenario_2.2.2_comps4viz <- scenario_2.2.2_comps.all %>% dplyr::filter(yes != 0) %>% 
  uncount(yes)
