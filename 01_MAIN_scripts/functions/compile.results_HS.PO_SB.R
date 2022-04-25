########## Compile and report results #########
# Combine above to make dataframe with truth and estimates side-by-side
# store years from youngest sibling in comparisons to end of study

####-------------------Compile results when estimating all parameters----------------
#Extract true values from year of estimation (ie estimation.year)
if(abundance.only != "yes"){
Mom.all_truth <- round(pop.size.tibble$Female.adult.pop[estimation.year],0) # True Nf
Dad.all_truth <- round(pop.size.tibble$Male.adult.pop[estimation.year], 0) # True Nm
Mom.breed_truth <- round(pop.size.tibble$Num.mothers[estimation.year],0) # True Nf
Dad.breed_truth <- round(pop.size.tibble$Num.fathers[estimation.year], 0) # True Nm
pb_truth <- Mom.breed_truth/Mom.all_truth #True percentage of breeders
surv_mean <- round(mean(sVec[ref.year:n_yrs]), 4) # True adult survival over estimation period
#Adult_truth <- round(pop.size.tibble$Total.adult.pop[estimation.year], 0) # Used for sex-aggregated model
lam_truth <- round(mean.adult.lambda, 4)
psi_truth <- psi.truth
Mom_min <- min(pop.size.tibble$Female.adult.pop[estimation.year:n_yrs]) #Minimum Nf over estimation period
Mom_max <- max(pop.size.tibble$Female.adult.pop[estimation.year:n_yrs]) #Maximum Nf over estimation period
Dad_min <- min(pop.size.tibble$Male.adult.pop[estimation.year:n_yrs]) #Minimum Nm over estimation period
Dad_max <- max(pop.size.tibble$Male.adult.pop[estimation.year:n_yrs]) #Maximum Nm over estimation period
#Adult_min <- min(pop.size.tibble$Total.adult.pop[estimation.year:n_yrs]) # Used for sex-aggregated model
#Adult_max <- max(pop.size.tibble$Total.adult.pop[estimation.year:n_yrs]) # Used for sex-aggregated model
surv_min <- min(sVec[estimation.year:n_yrs]) #Minimum survival over estimation period
surv_max <- max(sVec[estimation.year:n_yrs]) #Maximum survival over estimation period
lam_min <- min(adult.lambda[ref.year:n_yrs]) #Minimum lambda over estimation period
lam_max <- max(adult.lambda[ref.year:n_yrs]) #Maximum lambda over estimation period

#Create dataframe of estimates and truth
estimates <- model.summary2 %>%
  mutate(truth = c(Mom.all_truth, Mom.breed_truth, Dad.all_truth, surv_mean, lam_truth, psi_truth, pb_truth)) %>%
  as_tibble()

#Extract more metrics that can help with troubleshooting and visualization
Juv_total_samples <- sample.size.juvs * length(sample.years) # total samples
Adult_total_samples <- sample.size.rents * length(sample.years)
pop_size_mean <- round(mean(pop.size.tibble$population_size[estimation.year:n_yrs]),0) #Mean TOTAL population size over estimation period

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
metrics <- cbind(c(rep(mom.Exp.PO, times = 2), #for Nfa, Nfb
                   dad.Exp.PO, #For Nm
                   rep(mom.Exp.PO + dad.Exp.PO, #For surv, lam, psi, and pb
                       times = n_params-3)),
                 c(rep(mom.PO.matches, times = 2),
                   dad.PO.matches,
                   rep(mom.PO.matches + dad.PO.matches,
                       times = n_params-3)), # number of positive IDs i.e. half-sibs; subtract 2 for sex-specific abundance parameters
                 c(rep(mom.Exp.HS, times = 2),
                   dad.Exp.HS,
                   rep(mom.Exp.HS + dad.Exp.HS,
                       times = n_params-3)),
                 c(rep(mom.HS.matches, times = 2),
                   dad.HS.matches,
                   rep(mom.HS.matches + dad.HS.matches,
                       times = n_params-3)),
                 c(rep(length(sampled.mothers), times = 2),
                   length(sampled.fathers),
                   rep(length(sampled.mothers) + length(sampled.fathers),
                       times = n_params-3)), #number of unique sampled parents
                 c(rep(mean.adult.lambda, times = n_params)), # mean lambda over estimation period
                 c(rep(Juv_total_samples, times = n_params)), # Juv total samples
                 c(rep(Adult_total_samples, times = n_params)), # Adult total samples
                 c(rep(iter, times = n_params)), #iteration
                 c(rep(rseed, times = n_params)))
colnames(metrics) <- c("Exp_POPs", "POPs_detected", "Exp_HSPs", "HSPs_detected", "unique_parents_in_sample", "mean_adult_lambda", "Juvenile_samples", "Adult_samples", "iteration", "seed")
}

####---------------------Compile results when fixing parameters------------------
#Extract true values from year of estimation (ie estimation.year)
if(abundance.only == "yes"){
Mom.all_truth <- round(pop.size.tibble$Female.adult.pop[estimation.year],0) # True Nf
Dad.all_truth <- round(pop.size.tibble$Male.adult.pop[estimation.year], 0) # True Nm
Mom.breed_truth <- round(pop.size.tibble$Num.mothers[estimation.year],0) # True Nf
Dad.breed_truth <- round(pop.size.tibble$Num.fathers[estimation.year], 0) # True Nm
surv_mean <- round(mean(sVec[ref.year:n_yrs]), 4) # True adult survival over estimation period
#Adult_truth <- round(pop.size.tibble$Total.adult.pop[estimation.year], 0) # Used for sex-aggregated model
lam_truth <- round(mean.adult.lambda, 4)
psi_truth <- psi.truth
Mom_min <- min(pop.size.tibble$Female.adult.pop[estimation.year:n_yrs]) #Minimum Nf over estimation period
Mom_max <- max(pop.size.tibble$Female.adult.pop[estimation.year:n_yrs]) #Maximum Nf over estimation period
Dad_min <- min(pop.size.tibble$Male.adult.pop[estimation.year:n_yrs]) #Minimum Nm over estimation period
Dad_max <- max(pop.size.tibble$Male.adult.pop[estimation.year:n_yrs]) #Maximum Nm over estimation period
#Adult_min <- min(pop.size.tibble$Total.adult.pop[estimation.year:n_yrs]) # Used for sex-aggregated model
#Adult_max <- max(pop.size.tibble$Total.adult.pop[estimation.year:n_yrs]) # Used for sex-aggregated model
surv_min <- min(sVec[estimation.year:n_yrs]) #Minimum survival over estimation period
surv_max <- max(sVec[estimation.year:n_yrs]) #Maximum survival over estimation period
lam_min <- min(adult.lambda[ref.year:n_yrs]) #Minimum lambda over estimation period
lam_max <- max(adult.lambda[ref.year:n_yrs]) #Maximum lambda over estimation period
#Create dataframe of estimates and truth
estimates <- model.summary2 %>% 
  mutate(min = c(Mom_min, Dad_min), 
         max = c(Mom_max, Dad_max), 
         truth = c(Mom.all_truth, Dad.all_truth),
         breed.truth = c(Mom.breed_truth, Dad.breed_truth)) %>%
  as_tibble()

#Extract more metrics that can help with troubleshooting and visualization
Juv_total_samples <- sample.size.juvs * length(sample.years) # total samples
Adult_total_samples <- sample.size.rents * length(sample.years)
pop_size_mean <- round(mean(pop.size.tibble$population_size[estimation.year:n_yrs]),0) #Mean TOTAL population size over estimation period

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
metrics <- cbind(c(mom.Exp.PO,
                   dad.Exp.PO,
                   rep(mom.Exp.PO + dad.Exp.PO,
                       times = n_params-2)),
                 c(mom.PO.matches, 
                   dad.PO.matches, 
                   rep(mom.PO.matches + dad.PO.matches,
                       times = n_params-2)), # number of positive IDs i.e. half-sibs; subtract 2 for sex-specific abundance parameters
                 c(mom.Exp.HS,
                   dad.Exp.HS,
                   rep(mom.Exp.HS + dad.Exp.HS,
                       times = n_params-2)),
                 c(mom.HS.matches,
                   dad.HS.matches,
                   rep(mom.HS.matches + dad.HS.matches,
                       times = n_params-2)),
                 c(length(sampled.mothers), 
                   length(sampled.fathers), 
                   rep(length(sampled.mothers) + length(sampled.fathers), 
                       times = n_params-2)), #number of unique sampled parents
                 c(rep(mean.adult.lambda, times = n_params)), # mean lambda over estimation period
                 c(rep(Juv_total_samples, times = n_params)), # Juv total samples
                 c(rep(Adult_total_samples, times = n_params)), # Adult total samples
                 c(rep(iter, times = n_params)), #iteration
                 c(rep(rseed, times = n_params))) 
colnames(metrics) <- c("Exp_POPs", "POPs_detected", "Exp_HSPs", "HSPs_detected", "unique_parents_in_sample", "mean_adult_lambda", "Juvenile_samples", "Adult_samples", "iteration", "seed")
}
