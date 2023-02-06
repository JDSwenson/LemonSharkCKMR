#Calculate truth for Objectives 1 and 2
#----------------Set input file locations ------------------------------
#Population simulation files
PopSim.location <- "G://My Drive/Personal_Drive/R/CKMR/Population.simulations/Pre_Nov2022/"
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

#---------------------------Objective 1 truth calculations----------------------#
#Confirmed that the population size is the same for all sampling schemes
obj1.popsize <- readRDS(file = paste0(PopSim.location, PopSim.lambda.1,"/annual_breeding/pop.size_", lambda_date.1, "_", inSeeds, "_", PopSim.lambda.1, "_", PopSim.annual.breeding, "_", Sampling.scheme.YOY))

#Confirmed that the comparisons are the same for scenarios 1.1 and 1.2
#Need to calculate mean truth over years that comparisons span
#YOY
scenario_1.1_mom.comps.YOY <- readRDS(file = paste0(objective_1_results_location, mom.comps.prefix, "_", date.of.simulation.scenario_1.1Y, "_", inSeeds, "_", file.scenario_1.1Y)) %>% 
  dplyr::select(-BI)

scenario_1.1_dad.comps.YOY <- readRDS(file = paste0(objective_1_results_location, dad.comps.prefix, "_", date.of.simulation.scenario_1.1Y, "_", inSeeds, "_", file.scenario_1.1Y))

obj1_all.comps.YOY <- rbind(scenario_1.1_mom.comps.YOY, scenario_1.1_dad.comps.YOY)

#Juveniles
scenario_1.1_mom.comps.juvs <- readRDS(file = paste0(objective_1_results_location, mom.comps.prefix, "_", date.of.simulation.scenario_1.1J, "_", inSeeds, "_", file.scenario_1.1J)) %>% 
  dplyr::select(-BI)

scenario_1.1_dad.comps.juvs <- readRDS(file = paste0(objective_1_results_location, dad.comps.prefix, "_", date.of.simulation.scenario_1.1J, "_", inSeeds, "_", file.scenario_1.1J))

obj1_all.comps.juvs <- rbind(scenario_1.1_mom.comps.juvs, scenario_1.1_dad.comps.juvs)

#All age classes
scenario_1.1_mom.comps.all <- readRDS(file = paste0(objective_1_results_location, mom.comps.prefix, "_", date.of.simulation.scenario_1.1A, "_", inSeeds, "_", file.scenario_1.1A)) %>% 
  dplyr::select(-BI)

scenario_1.1_dad.comps.all <- readRDS(file = paste0(objective_1_results_location, dad.comps.prefix, "_", date.of.simulation.scenario_1.1A, "_", inSeeds, "_", file.scenario_1.1A))

obj1_all.comps.all <- rbind(scenario_1.1_mom.comps.all, scenario_1.1_dad.comps.all)


(YOY.ref <- obj1_all.comps.YOY %>% group_by(iteration) %>% summarize(min.ref = min(ref.year)))
(juvs.ref <- obj1_all.comps.juvs %>% group_by(iteration) %>% summarize(min.ref = min(ref.year)))
(all.ref <- obj1_all.comps.all %>% group_by(iteration) %>% summarize(min.ref = min(ref.year)))


obj1.truth.YOY.temp = obj1.truth.juvs.temp = obj1.truth.all.temp = obj1.truth.df <- NULL

#Assume all sampling schemes have the same number of iterations
for(i in 1:nrow(YOY.ref)){
#YOY
  reference.yr.YOY <- YOY.ref %>% dplyr::filter(iteration == i) %>% 
    pull(min.ref)
  
  obj1.truth.YOY.temp <- obj1.popsize %>% dplyr::filter(iteration == i,
                                             year >= reference.yr.YOY) %>% 
    summarize(male.truth = mean(Male.adult.pop),
              female.truth = mean(Female.adult.pop)) %>% 
    mutate(iteration = i,
           sampling.scheme = "target YOY")

  #YOY
  reference.yr.juvs <- juvs.ref %>% dplyr::filter(iteration == i) %>% 
    pull(min.ref)
  
  obj1.truth.juvs.temp <- obj1.popsize %>% dplyr::filter(iteration == i,
                                              year >= reference.yr.juvs) %>% 
    summarize(male.truth = mean(Male.adult.pop),
              female.truth = mean(Female.adult.pop)) %>% 
    mutate(iteration = i,
           sampling.scheme = "sample all juvenile ages")

  #all
  reference.yr.all <- all.ref %>% dplyr::filter(iteration == i) %>% 
    pull(min.ref)
  
  obj1.truth.all.temp <- obj1.popsize %>% dplyr::filter(iteration == i,
                                              year >= reference.yr.all) %>% 
    summarize(male.truth = mean(Male.adult.pop),
              female.truth = mean(Female.adult.pop)) %>% 
    mutate(iteration = i,
           sampling.scheme = "sample all age classes")
  
  obj1.truth.df <- rbind(obj1.truth.df, obj1.truth.YOY.temp, obj1.truth.juvs.temp, obj1.truth.all.temp)
}


#---------------------------Objective 2 truth calculations----------------------#
#obj2.1 is variable lambda but no lambda parameter in the model; combine with obj1.2 with stable lambda
#obj2.2.1 is variable lambda w/ lambda parameter in the model and estimating in year 85; combine with obj2.2.2 for estimates in year 80 and 90 as well.

#Confirmed that the population size is the same for all sampling schemes
obj2.popsize.neutral <- readRDS(file = paste0(PopSim.location, PopSim.lambda.1, "/annual_breeding/pop.size_", lambda_date.1, "_", inSeeds, "_", PopSim.lambda.1, "_", PopSim.annual.breeding, "_", Sampling.scheme.YOY)) %>% 
  dplyr::filter(year %in% c(80, 85, 90)) %>% 
  dplyr::mutate(population.growth = "stable")

obj2.popsize.variable <- readRDS(file = paste0(PopSim.location, PopSim.lambda.variable, "/pop.size_", lambda_date.1, "_", inSeeds, "_", PopSim.lambda.variable, "_", PopSim.annual.breeding, "_", Sampling.scheme.YOY)) %>% 
  dplyr::filter(year %in% c(80, 85, 90)) %>% 
  dplyr::mutate(population.growth = ifelse(iteration <= 500, "slight negative", "slight positive"))

obj2.popsize.extreme <- readRDS(file = paste0(PopSim.location, PopSim.lambda.variable, "/pop.size_", lambda_date.4, "_", inSeeds, "_", PopSim.lambda.extreme, "_", PopSim.annual.breeding, "_", Sampling.scheme.YOY)) %>% 
  dplyr::filter(year %in% c(80, 85, 90)) %>% 
  dplyr::mutate(population.growth = "extreme negative")


obj2.truth.df <- rbind(obj2.popsize.neutral, obj2.popsize.variable, obj2.popsize.extreme) %>% 
  dplyr::rename(est.yr = year) %>% 
  dplyr::select(est.yr, iteration, seed, Male.adult.pop, Female.adult.pop, population.growth)


#---------------------------Objective 3 truth calculations----------------------#
#Confirmed that the population size is the same for all sampling schemes
#No non-conformists
obj3.popsize.df <- readRDS(file = paste0(PopSim.location, PopSim.lambda.1, "/biennial_breeding/pop.size_", lambda_date.2, "_", inSeeds, "_", PopSim.lambda.1, "_", PopSim.biennial.breeding, "_", Sampling.scheme.YOY)) %>% 
  dplyr::filter(year %in% c(80, 85, 90)) %>% 
  dplyr::mutate(population.growth = "neutral") 

#Add this when I revist calculating psi with Liz
#%>% mutate(psi.truth = (-(Num.mothers/Female.adult.pop - 0.5)*2) - 1)

obj3.truth.df_psi1 <- obj3.popsize.df %>%
  dplyr::rename(est.yr = year) %>% 
  dplyr::select(est.yr, iteration, seed, Male.adult.pop, Female.adult.pop, Num.mothers, Num.fathers, population.growth)


#----------------------------Calculate truth of psi-----------------------------#
options(dplyr.summarise.inform = FALSE)
# 
# #----------psi 0.5--------------------
# (psi_0.5_samplefiles <- list.files(path = paste0(PopSim.location, "lambda.1/biennial_breeding/plus_NonConformists/sample.info/"),
#                             recursive = FALSE,
#                             pattern = "*psi0.5_*",
#                             full.names = TRUE))
# 
# psi_0.5 <- map_dfr(psi_0.5_samplefiles, readRDS) %>% 
#   mutate(sim.psi = 0.5) %>% 
#   dplyr::filter(sex == "F")
# 
# #-----------psi 0.75----------------
# (psi_0.75_samplefiles <- list.files(path = paste0(PopSim.location, "lambda.1/biennial_breeding/plus_NonConformists/sample.info/"),
#                                    recursive = FALSE,
#                                    pattern = "*psi0.75_*",
#                                    full.names = TRUE))
# 
# psi_0.75 <- map_dfr(psi_0.75_samplefiles, readRDS) %>% 
#   mutate(sim.psi = 0.75) %>% 
#   dplyr::filter(sex == "F")
# 
# 
# #-----------psi 0.9----------------
# (psi_0.9_samplefiles <- list.files(path = paste0(PopSim.location, "lambda.1/biennial_breeding/plus_NonConformists/sample.info/"),
#                                     recursive = FALSE,
#                                     pattern = "*psi0.9_*",
#                                     full.names = TRUE))
# 
# psi_0.9 <- map_dfr(psi_0.9_samplefiles, readRDS) %>% 
#   mutate(sim.psi = 0.9) %>% 
#   dplyr::filter(sex == "F")
# 
# psi.truth.temp =psi.truth.df_0.5 = psi.truth.df_0.75 = psi.truth.df_0.9 <- NULL
# 
# for(iter in 1:max(psi_0.5$iteration)){
# #Calculate psi truth
# psi.df <- psi_0.5 %>% dplyr::filter(iteration == iter) %>% #Use all samples from this iteration to calculate the truth
#   distinct(indv.name, .keep_all = TRUE) %>% #Make sure we aren't using duplicated individuals
#   group_by(repro.strategy, sim.psi) %>%
#   summarize(number = n())
# 
# (psi.truth.temp <- round(1 - psi.df[psi.df$repro.strategy == "non-conformist",3]/sum(psi.df$number), 3) %>% 
#     mutate(iteration = iter,
#            sim.psi = 0.5) %>% 
#   dplyr::rename(psi.truth = number))
#   
# psi.truth.df_0.5 <- bind_rows(psi.truth.df_0.5, psi.truth.temp)
# print(paste0("Finished iteration ", iter))
# }
# 
# for(iter in 1:max(psi_0.75$iteration)){
#   #Calculate psi truth
#   psi.df <- psi_0.75 %>% dplyr::filter(iteration == iter) %>% #Use all samples from this iteration to calculate the truth
#     distinct(indv.name, .keep_all = TRUE) %>% #Make sure we aren't using duplicated individuals
#     group_by(repro.strategy, sim.psi) %>%
#     summarize(number = n())
#   
#   (psi.truth.temp <- round(1 - psi.df[psi.df$repro.strategy == "non-conformist",3]/sum(psi.df$number), 3) %>% 
#       mutate(iteration = iter,
#              sim.psi = 0.75) %>% 
#     dplyr::rename(psi.truth = number))
#   
#   psi.truth.df_0.75 <- bind_rows(psi.truth.df_0.75, psi.truth.temp)
#   print(paste0("Finished iteration ", iter))
# }
# 
# for(iter in 1:max(psi_0.9$iteration)){
#   #Calculate psi truth
#   psi.df <- psi_0.9 %>% dplyr::filter(iteration == iter) %>% #Use all samples from this iteration to calculate the truth
#     distinct(indv.name, .keep_all = TRUE) %>% #Make sure we aren't using duplicated individuals
#     group_by(repro.strategy, sim.psi) %>%
#     summarize(number = n())
#   
#   (psi.truth.temp <- round(1 - psi.df[psi.df$repro.strategy == "non-conformist",3]/sum(psi.df$number), 3) %>%
#       mutate(iteration = iter,
#              sim.psi = 0.9) %>%
#     dplyr::rename(psi.truth = number))
#   
#   psi.truth.df_0.9 <- bind_rows(psi.truth.df_0.9, psi.truth.temp)
#   print(paste0("Finished iteration ", iter))
#   }
# 
# #Add psi = 1 to make join easier
# psi.truth.df_1 <- tibble(psi.truth = 1, iteration = c(1:500), sim.psi = 1)
#   
# psi_all.truth <- bind_rows(psi.truth.df_0.5, psi.truth.df_0.75, psi.truth.df_0.9, psi.truth.df_1) %>% 
#   as_tibble()
# 
#write_rds(psi_all.truth, file = paste0(objective_3_results_location, "psi_all.truth"))

#psi_all.truth <- read_rds(file = paste0(objective_3_results_location, "psi_all.truth"))
