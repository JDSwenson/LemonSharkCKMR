########## Compile and report results #########
# Combine above to make dataframe with truth and estimates side-by-side
# store years from youngest sibling in comparisons to end of study
yrs <- c(estimation.year:n_yrs)
ref.year <- min(mom_comps.all$ref.year, dad_comps.all$ref.year)

mean.adult.lambda <- mean(adult.lambda[ref.year:n_yrs], na.rm=T) # mean Lambda over years of estimation for adults ### HOW DO WE KNOW THIS?
#Extract true values from year of estimation (ie estimation.year)
Mom_truth <- round(pop.size.temp$Female.adult.pop[estimation.year],0) # True Nf
Dad_truth <- round(pop.size.temp$Male.adult.pop[estimation.year], 0) # True Nm
surv_mean <- round(mean(sVec[ref.year:n_yrs]), 4) # True adult survival over estimation period
#Adult_truth <- round(pop.size.temp$Total.adult.pop[estimation.year], 0) # Used for sex-aggregated model
lam_truth <- round(mean.adult.lambda, 4)
Mom_min <- min(pop.size.temp$Female.adult.pop[estimation.year:n_yrs]) #Minimum Nf over estimation period
Mom_max <- max(pop.size.temp$Female.adult.pop[estimation.year:n_yrs]) #Maximum Nf over estimation period
Dad_min <- min(pop.size.temp$Male.adult.pop[estimation.year:n_yrs]) #Minimum Nm over estimation period
Dad_max <- max(pop.size.temp$Male.adult.pop[estimation.year:n_yrs]) #Maximum Nm over estimation period
#Adult_min <- min(pop.size.temp$Total.adult.pop[estimation.year:n_yrs]) # Used for sex-aggregated model
#Adult_max <- max(pop.size.temp$Total.adult.pop[estimation.year:n_yrs]) # Used for sex-aggregated model
surv_min <- min(sVec[estimation.year:n_yrs]) #Minimum survival over estimation period
surv_max <- max(sVec[estimation.year:n_yrs]) #Maximum survival over estimation period
lam_min <- min(adult.lambda[ref.year:n_yrs]) #Minimum lambda over estimation period
lam_max <- max(adult.lambda[ref.year:n_yrs]) #Maximum lambda over estimation period

#Create dataframe of estimates and truth
estimates <- model.summary2 %>% 
  mutate(min = c(Mom_min, Dad_min, surv_min, lam_min), max = c(Mom_max, Dad_max, surv_max, lam_max), truth = c(Mom_truth, Dad_truth, surv_mean, lam_truth)) %>%
  as_tibble()

#Extract more metrics that can help with troubleshooting and visualization
total_samples <- sample.size * length(sample.years) # total samples
pop_size_mean <- round(mean(pop.size.temp$population_size[estimation.year:n_yrs]),0) #Mean TOTAL population size over estimation period

mom.PO.matches <- mom_comps.all %>% filter(type == "PO") %>%
  summarize(matches = sum(yes)) %>% 
  pull(matches)

mom.HS.matches <- mom_comps.all %>% filter(type == "HS") %>%
  summarize(matches = sum(yes)) %>% 
  pull(matches)

dad.PO.matches <- dad_comps.all %>% filter(type == "PO") %>%
  summarize(matches = sum(yes)) %>% 
  pull(matches)

dad.HS.matches <- dad_comps.all %>% filter(type == "HS") %>%
  summarize(matches = sum(yes)) %>% 
  pull(matches)

#Bind metrics together
metrics <- cbind(c(mom.PO.matches, 
                   dad.PO.matches, 
                   rep(mom.PO.matches + dad.PO.matches,
                       times = n_params-2)), # number of positive IDs i.e. half-sibs; subtract 2 for sex-specific abundance parameters
                 c(mom.HS.matches,
                   dad.HS.matches,
                   rep(mom.HS.matches + dad.HS.matches,
                       times = n_params-2)),
                 c(length(sampled.mothers), 
                   length(sampled.fathers), 
                   rep(length(sampled.mothers) + length(sampled.fathers), 
                       times = n_params-2)), #number of unique sampled parents
                 c(rep(mean.adult.lambda, times = n_params)), # mean lambda over estimation period
                 c(rep(total_samples, times = n_params)), # total samples
                 c(rep(iter, times = n_params))) #iteration
colnames(metrics) <- c("POPs_detected", "HSPs_detected", "unique_parents_in_sample", "mean_adult_lambda", "total_samples", "iteration")