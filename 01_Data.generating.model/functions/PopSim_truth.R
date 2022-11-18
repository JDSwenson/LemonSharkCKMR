########## Compile and report results #########
#Extract reference year as earliest cohort
ref.year <- sample.df_all.info %>% dplyr::filter(age.x < repro.age) %>% 
  arrange(birth.year) %>% 
  slice_min(birth.year) %>% 
  distinct(birth.year) %>% 
  pull(birth.year)

####-------------------Compile results when estimating all parameters----------------
#### TRUTH ####
#Abundance
Mom.all_truth <- round(pop.size.tibble$Female.adult.pop[estimation.year],0) # True Nf
Dad.all_truth <- round(pop.size.tibble$Male.adult.pop[estimation.year], 0) # True Nm
Mom.breed_truth <- round(pop.size.tibble$Num.mothers[estimation.year],0) # True Nf
Dad.breed_truth <- round(pop.size.tibble$Num.fathers[estimation.year], 0) # True Nm
pb_truth <- Mom.breed_truth/Mom.all_truth #True percentage of female breeders
Mom_min <- min(pop.size.tibble$Female.adult.pop[estimation.year:n_yrs]) #Minimum Nf over estimation period
Mom_max <- max(pop.size.tibble$Female.adult.pop[estimation.year:n_yrs]) #Maximum Nf over estimation period
Dad_min <- min(pop.size.tibble$Male.adult.pop[estimation.year:n_yrs]) #Minimum Nm over estimation period
Dad_max <- max(pop.size.tibble$Male.adult.pop[estimation.year:n_yrs]) #Maximum Nm over estimation period

#Survival
surv_mean <- round(mean(sVec[ref.year:n_yrs]), 4) # True adult survival over sample period
surv_min <- min(sVec[ref.year:n_yrs]) #Minimum survival over sample period
surv_max <- max(sVec[ref.year:n_yrs]) #Maximum survival over sample period

#Lambda
(mean.total.lam <- mean(pop.size.tibble$total.lambda[ref.year:n_yrs], na.rm=T)) # mean population growth for whole population
sd(pop.size.tibble$total.lambda[ref.year:n_yrs], na.rm=T) # sd Lambda
(mean.adult.lambda <- mean(pop.size.tibble$adult.lambda[ref.year:n_yrs], na.rm=T)) # mean population growth for whole population
(mean.adult.lambda_last.10 <- mean(pop.size.tibble$adult.lambda[80:n_yrs], na.rm=T)) # mean population growth for whole population)

lam_truth <- round(mean.adult.lambda, 4)
lam_min <- min(adult.lambda[ref.year:n_yrs]) #Minimum lambda over estimation period
lam_max <- max(adult.lambda[ref.year:n_yrs]) #Maximum lambda over estimation period

#psi
psi_truth <- 1-non.conformists

#Create dataframe of estimates and truth
true.values <- tibble(parameter = c("Nf", "psi", "Nm", "survival", "lambda"),
                    all.truth = c(Mom.all_truth, psi_truth, Dad.all_truth, surv_mean, lam_truth),
                    breed.truth = c(Mom.breed_truth, psi_truth, Dad.breed_truth, surv_mean, lam_truth),
                    mean.adult.lambda = mean.adult.lambda,
                    mean.adult.lambda_last.10 = mean.adult.lambda_last.10,
                    female.fecundity = ff,
                    iteration = iter,
                    seed = rseed,
                    population.growth = population.growth)


if(exists("ff.shift") == TRUE){
  true.values$female.fecundity = ff.shift
}

if(population.growth == "lambda.extreme"){
  true.values <- true.values %>% mutate(adult.survival.adj = Sa,
                                    juvenile.survival.adj = Sj)
  }