pop.size.1.summary <- pop.size.1 %>% dplyr::filter(year == 85) %>% 
  dplyr::distinct(year, iter, .keep_all = TRUE) %>%
  dplyr::select(!Total.adult.pop) %>% 
 # dplyr::rename(Nm = Male.adult.pop, Nf = Female.adult.pop) %>% 
  pivot_longer(cols = ends_with("adult.pop"), names_to = "col", values_to = "truth.N") %>% 
  mutate(parameter = ifelse(col == "Male.adult.pop", "Nm", "Nf")) %>% 
  dplyr::select(parameter, truth.N, iteration = iter)

pop.size.1.lam <- pop.size.1 %>% dplyr::filter(year >= 85) %>% 
  group_by(iter) %>% 
  summarize(true.lam = mean(adult.lambda)) %>% 
  dplyr::rename(iteration = iter)


results.1 <- results.1 %>% left_join(pop.size.1.summary, by = c("iteration", "parameter")) %>% 
  left_join(pop.size.1.lam, by = "iteration") %>% 
  mutate(truth = ifelse(parameter == "Nf" | parameter == "Nm", truth.N, 
                        ifelse(parameter == "lam", true.lam, truth))) 


pop.size.2.summary <- pop.size.2 %>% dplyr::filter(year == 85) %>% 
  dplyr::distinct(year, iter, .keep_all = TRUE) %>%
  dplyr::select(!Total.adult.pop) %>% 
  # dplyr::rename(Nm = Male.adult.pop, Nf = Female.adult.pop) %>% 
  pivot_longer(cols = ends_with("adult.pop"), names_to = "col", values_to = "truth.N") %>% 
  mutate(parameter = ifelse(col == "Male.adult.pop", "Nm", "Nf")) %>% 
  dplyr::select(parameter, truth.N, iteration = iter)


results.2 <- results.2 %>% left_join(pop.size.2.summary, by = c("iteration", "parameter")) %>% 
  left_join(pop.size.2.lam, by = "iteration") %>% 
  mutate(truth = ifelse(parameter == "Nf" | parameter == "Nm", truth.N, 
                        ifelse(parameter == "lam", true.lam, truth))) 

pop.size.2.lam <- pop.size.2 %>% dplyr::filter(year >= 85) %>% 
  group_by(iter) %>% 
  summarize(true.lam = mean(adult.lambda)) %>% 
  dplyr::rename(iteration = iter)
