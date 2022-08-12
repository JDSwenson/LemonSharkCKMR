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


