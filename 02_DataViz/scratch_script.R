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
