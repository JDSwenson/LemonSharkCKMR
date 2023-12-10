mom.comps1 <- readRDS(file = paste0(results_location, "mom.comps_04Oct2023_Seeds2022.04.15_scenario_3.1.2_multiennial.model")) %>% 
  mutate(scenario = "scenario_3.1.2")

mom.comps2 <- readRDS(file = paste0(results_location, "mom.comps_04Oct2023_Seeds2022.04.15_scenario_3.2.2_multiennial.model")) %>% 
  mutate(scenario = "scenario_3.2.2")

mom.comps3 <- readRDS(file = paste0(results_location, "mom.comps_04Oct2023_Seeds2022.04.15_scenario_3.3.2_multiennial.model")) %>% 
  mutate(scenario = "scenario_3.3.2")

mom.comps4 <- readRDS(file = paste0(results_location, "mom.comps_04Oct2023_Seeds2022.04.15_scenario_3.4.2_multiennial.model")) %>% 
  mutate(scenario = "scenario_3.4.2")

mom.comps5 <- readRDS(file = paste0(results_location, "mom.comps_04Oct2023_Seeds2022.04.15_scenario_3.5.2_multiennial.model")) %>% 
  mutate(scenario = "scenario_3.5.2")

mom.comps6 <- readRDS(file = paste0(results_location, "mom.comps_04Oct2023_Seeds2022.04.15_scenario_3.6.2_multiennial.model")) %>% 
  mutate(scenario = "scenario_3.6.2")

mom.comps7 <- readRDS(file = paste0(results_location, "mom.comps_04Oct2023_Seeds2022.04.15_scenario_3.7.2_multiennial.model")) %>% 
  mutate(scenario = "scenario_3.7.2")

mom.comps.all <- bind_rows(mom.comps1, mom.comps2, mom.comps3, mom.comps4, mom.comps5, mom.comps6, mom.comps7)

mom.comps.all %>% dplyr::count(scenario)

mom.comps.all %>% saveRDS(file = paste0(results_location, "objective3_mom.comps.all"))

ls_dataSim_allSibs %>% dplyr::filter(parameter == "lambda", Q50 >= 1.29)

dad.comps <- readRDS(file = paste0(results_location, "lemon_shark_sims/peer_review_rd2/lemon_shark_time_series_sims_dad.comps_Conn5_allSibs_temp"))

mom.comps <- readRDS(file = paste0(results_location, "lemon_shark_sims/peer_review_rd2/lemon_shark_time_series_sims_mom.comps_Conn5_allSibs_temp"))

mom.comps %>% dplyr::filter(iteration == 1) %>% mutate(estimation.yr = ref.year + pop.growth.yrs) %>% 
  group_by(time_window, BI, pop.growth.yrs) %>% 
  summarize(sum(yes))

mom.comps %>% dplyr::filter(pop.growth.yrs == 19) %>% 
  group_by(BI) %>% 
  summarize(sum(yes))

obj2_results %>% dplyr::filter(scenario %in% c("scenario_2.3.3", "scenario_2.4.3")) %>% dplyr::filter(parameter == "lambda") %>% 
  dplyr::select(parameter, Q50, all.truth, scenario, estimation.sim.label, model.label)

obj2_results %>% dplyr::filter(scenario %in% c("scenario_2.3.3", "scenario_2.4.3"))%>% 
  dplyr::filter(parameter == "lambda") %>% 
  group_by(scenario) %>% 
  summarize(min(Q50), max(Q50))



ls_data %>% dplyr::filter(years_sampled == 3)

obj2_results %>% dplyr::filter(parameter %in% c("lambda", "survival"),
                               lambda.prior == "diffuse (0.80 - 1.20)") %>% 
  group_by(parameter, est.yr, prop.sampled, sampling.scheme, population.growth) %>% 
  summarize(n())

#Data Viz - check why sample sizes don't match
test.all <- s4.1.1_all.ages %>% group_by(iteration) %>% summarize(num = n())
test.juvs <- s4.1.1_all.juvs %>% group_by(iteration) %>% summarize(num = n())
test.YOY <- s4.1.1_YOY %>% group_by(iteration) %>% summarize(num = n())

#Population simulations
born90 <- loopy.list[[90]] %>% dplyr::filter(birth.year == 90)
born89 <- loopy.list[[89]] %>% dplyr::filter(birth.year == 89)
born88 <- loopy.list[[89]] %>% dplyr::filter(birth.year == 88)
born87 <- loopy.list[[88]] %>% dplyr::filter(birth.year == 87)

sum(born90$mother.x %in% born89$mother.x == T)
sum(born90$mother.x %in% born88$mother.x == T)
sum(born90$mother.x %in% born87$mother.x == T)




#Results
obj1_results %>% dplyr::filter(sampling.scheme == "sample all age classes") %>% 
  group_by(parameter, prop.sampled.lab) %>% 
  summarize(mean(POPs_detected))

obj1_results %>% dplyr::filter(sampling.scheme == "sample all age classes") %>% 
  group_by(parameter, prop.sampled.lab) %>% 
  summarize(mean(mean.adult.lambda),
            sd(mean.adult.lambda))

obj2_results %>% dplyr::filter(population.growth == "extreme negative", est.yr == 85) %>% 
  distinct(iteration, .keep_all = TRUE) %>%
  summarize(mean(mean.adult.lambda.last10),
            sd(mean.adult.lambda.last10))


obj3_results %>% dplyr::filter(sampling.scheme == "sample all age classes") %>% 
  group_by(parameter, prop.sampled.lab) %>% 
  summarize(mean(mean.adult.lambda),
            sd(mean.adult.lambda))


#---------------------Check mean number of offspring by reproductive strategy----------------
mean.off.conf = mean.off.temp <- NULL

for(i in 1:length(loopy.list)){
  
  moms.off <- loopy.list[[i]] %>% dplyr::filter(age.x == 0) %>% 
    group_by(mother.x) %>% 
    summarize(offspring = n()) %>% 
    dplyr::rename(indv.name = mother.x)
  
  mom.info <- loopy.list[[i]] %>% inner_join(moms.off, by = "indv.name") %>% 
    dplyr::select(indv.name, repro.strategy, offspring)
  
  mean.off.temp <- mom.info %>% group_by(repro.strategy) %>% 
    summarize(mean.off = mean(offspring)) %>% 
    mutate(year = i)

  mean.off.conf <- bind_rows(mean.off.conf, mean.off.temp)    
}

#------------------Check the years in which conformists and non-conformists bred-----------------#
test.pop <- loopy.list[50:90] %>% map_dfr(bind_rows)

parents.tibble %>% dplyr::rename(indv.name = parent) %>% 
  inner_join(test.pop, by = "indv.name") %>% 
  dplyr::select(indv.name, num.off, year, repro.strategy) %>% 
  dplyr::distinct(indv.name, year, .keep_all = TRUE) %>% 
  dplyr::filter(indv.name == "aekruksyblrgpenfazbc")

#-----------------Check proportion of conformists vs non-conformists----------------#

table(loopy.list[[90]]$repro.strategy)
