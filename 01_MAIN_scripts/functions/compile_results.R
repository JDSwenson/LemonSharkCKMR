# Combine above to make dataframe with truth and estimates side-by-side
# store years from youngest sibling in comparisons to end of study
yrs <- c(est.year:n_yrs)

#Extract true values from year of estimation (ie est.year)
Mom_truth <- round(pop.size$Female.adult.pop[est.year],0) # True Nf
Dad_truth <- round(pop.size$Male.adult.pop[est.year], 0) # True Nm
surv_truth <- round(mean(sVec[est.year:n_yrs]), 4) # True adult survival over estimation period
#Adult_truth <- round(pop.size$Total.adult.pop[est.year], 0) # Used for sex-aggregated model
lam_truth <- round(mean.adult.lambda, 4)
Mom_min <- min(pop.size$Female.adult.pop[est.year:n_yrs]) #Minimum Nf over estimation period
Mom_max <- max(pop.size$Female.adult.pop[est.year:n_yrs]) #Maximum Nf over estimation period
Dad_min <- min(pop.size$Male.adult.pop[est.year:n_yrs]) #Minimum Nm over estimation period
Dad_max <- max(pop.size$Male.adult.pop[est.year:n_yrs]) #Maximum Nm over estimation period
#Adult_min <- min(pop.size$Total.adult.pop[est.year:n_yrs]) # Used for sex-aggregated model
#Adult_max <- max(pop.size$Total.adult.pop[est.year:n_yrs]) # Used for sex-aggregated model
surv_min <- min(sVec[est.year:n_yrs]) #Minimum survival over estimation period
surv_max <- max(sVec[est.year:n_yrs]) #Maximum survival over estimation period
lam_min <- min(adult.lambda[est.year:n_yrs]) #Minimum lambda over estimation period
lam_max <- max(adult.lambda[est.year:n_yrs]) #Maximum lambda over estimation period
mean.num.mothers.total <- round(mean(pop.size$Num.mothers[est.year:n_yrs]), 0) #Mean number of mothers in population over estimation period
mean.num.fathers.total <- round(mean(pop.size$Num.fathers[est.year:n_yrs]), 0) #Mean number of fathers in population over estimation period

#Create dataframe of estimates and truth
estimates <- model.summary2 %>% 
  mutate(min = c(Mom_min, Dad_min, surv_min, lam_min), max = c(Mom_max, Dad_max, surv_max, lam_max), truth = c(Mom_truth, Dad_truth, surv_truth, lam_truth)) %>%
  as_tibble()

#Extract more metrics that can help with troubleshooting and visualization
total_samples <- sample.size * length(sample.years) # total samples
pop_size_mean <- round(mean(pop.size$population_size[est.year:n_yrs]),0) #Mean TOTAL population size over estimation period

#Bind metrics together
metrics <- cbind(c(sum(mom_comps[,3]), sum(dad_comps[,3]), rep(sum(mom_comps[,3]) + sum(dad_comps[3]), times = n_params-2)), # number of positive IDs i.e. half-sibs; subtract 2 for sex-specific abundance parameters
                 c(mean.num.mothers.total, mean.num.fathers.total, rep(mean.num.mothers.total + mean.num.fathers.total, times = n_params - 2)), #Mean number of parents in population over estimation period
                 c(length(sampled.mothers), length(sampled.fathers), rep(length(sampled.mothers) + length(sampled.fathers), times = n_params-2)), #number of unique sampled parents
                 c(rep(mean.adult.lambda, times = n_params)), # mean lambda over estimation period
                 c(rep(total_samples, times = n_params)), # total samples
                 c(rep(pop_size_mean, times = n_params)), # mean population size over estimation period
                 c(rep(iter, times = n_params)),
                 c(rep(N.prior.max, times = n_params))) #iteration
colnames(metrics) <- c("parents_detected", "mean_unique_parents_in_pop", "unique_parents_in_sample", "mean_adult_lambda", "total_samples", "pop_size_mean", "iteration", "N.prior_max")