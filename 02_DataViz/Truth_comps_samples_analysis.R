library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(coda)
library(runjags)
library(ggmcmc)
library(viridis)
library(scales)
library(ggridges)

rm(list=ls())
today <- format(Sys.Date(), "%d%b%Y")

#----------------------Objective 2----------------------------#
#Confirmed that the comparisons are the same for scenarios 1.1 and 1.2
#Need to calculate mean truth over years that comparisons span
#YOY
scenario_2.2.2_mom.comps.YOY <- readRDS(file = paste0(objective_2_results_location, mom.comps.prefix, "_", date.of.simulation.scenario_2.2Y, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_target.YOY_change.est.yr")) %>% 
  dplyr::select(-BI)

scenario_2.2.2_dad.comps.YOY <- readRDS(file = paste0(objective_2_results_location, dad.comps.prefix, "_", date.of.simulation.scenario_2.2Y, "_", inSeeds, "_scenario_2.2.2_lambda.extreme_target.YOY_change.est.yr"))

scenario_2.2.2_all.comps.YOY <- rbind(scenario_2.2.2_mom.comps.YOY, scenario_2.2.2_dad.comps.YOY) %>% 
  mutate(sampling.scheme = "target YOY")

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

#Make example pairwise comparison matrix
scenario_2.2.2_mom.comps.all %>% dplyr::select(ref.year, all, yes, mort.yrs, type, parent, pop.growth.yrs, BI) %>% head(20) %>% write_csv(file = "example.pairwise.csv")


#See how frequently lambda is estimated to be negative
obj2_results %>% dplyr::filter(population.growth == "slight negative" | population.growth == "slight positive") %>% 
  dplyr::filter(lambda.prior %in% c("diffuse"),
                               parameter == "lambda") %>% 
  group_by(lambda.prior, sampling.scheme) %>% 
  dplyr::summarize(percent.correct = sum(mean < 1)/n() * 100)
