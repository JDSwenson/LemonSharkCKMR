########## Compile and report results #########
#Extract reference year as earliest cohort
# ref.year <- sample.df_all.info %>% dplyr::filter(age.x < repro.age) %>%
#   arrange(birth.year) %>%
#   slice_min(birth.year) %>%
#   distinct(birth.year) %>%
#   pull(birth.year)



####-------------------Compile results when estimating all parameters----------------
#### TRUTH ####
# #Abundance
# Mom.all_truth <- round(pop.size.tibble$Female.adult.pop[estimation.year],0) # True Nf
# Dad.all_truth <- round(pop.size.tibble$Male.adult.pop[estimation.year], 0) # True Nm
# Mom.breed_truth <- round(pop.size.tibble$Num.mothers[estimation.year],0) # True Nf
# Dad.breed_truth <- round(pop.size.tibble$Num.fathers[estimation.year], 0) # True Nm
# pb_truth <- Mom.breed_truth/Mom.all_truth #True percentage of female breeders
# Mom_min <- min(pop.size.tibble$Female.adult.pop[estimation.year:n_yrs]) #Minimum Nf over estimation period
# Mom_max <- max(pop.size.tibble$Female.adult.pop[estimation.year:n_yrs]) #Maximum Nf over estimation period
# Dad_min <- min(pop.size.tibble$Male.adult.pop[estimation.year:n_yrs]) #Minimum Nm over estimation period
# Dad_max <- max(pop.size.tibble$Male.adult.pop[estimation.year:n_yrs]) #Maximum Nm over estimation period

#Calculate geometric mean of survival for each reference year
#See: https://www.r-bloggers.com/2021/08/calculate-geometric-mean-in-r/
#ref.vec <- ref.tibble %>% distinct(ref.yr) %>% pull(ref.yr)

ref.tibble2 <- ref.tibble

ref.tibble2$survival <- NA #name "survival" for pivot later
ref.tibble2$surv_min <- NA
ref.tibble2$surv_max <- NA
ref.tibble2$surv_AM <- NA
ref.tibble2$psi <- NA #name "psi" for pivot later

for(r in 1:nrow(ref.tibble2)){
  
  ref.value <- ref.tibble2$ref.yr[r]
  #Geometric mean of survival
  (surv_mean <- exp(mean(log(sVec[ref.value:(n_yrs - 1)]))))
  
  #Arithemetic mean of survival
  (surv_mean_AM <- round(mean(sVec[ref.value:(n_yrs-1)]), 4))
  


##Equivalent way of calculating geometric mean
# sVec2 <- NULL
# 
# #Loop over the survival values for each year; note that survival here represents survival FROM the year it's recorded e..g sVec[2] represents survival from year 2 to year 3
# for(s in ref.value:(n_yrs - 1)){
#   sTemp <- sVec[s]
#   sVec2 <- c(sVec2, sTemp)
# }
# 
# delta <- n_yrs - ref.value
# 
# surv_mean <- (prod(sVec2))^(1/delta)

#surv_mean <- (sVec[n_yrs]/sVec[ref.value])^(1/(n_yrs - ref.value))
surv_min <- min(sVec[ref.value:(n_yrs-1)]) #Minimum survival over sample period
surv_max <- max(sVec[ref.value:(n_yrs-1)]) #Maximum survival over sample period

#Lambda
# (mean.total.lam <- mean(pop.size.tibble$total.lambda[ref.value:n_yrs], na.rm=T)) # mean population growth for whole population
# sd(pop.size.tibble$total.lambda[ref.value:n_yrs], na.rm=T) # sd Lambda
# (mean.adult.lambda <- mean(pop.size.tibble$adult.lambda[ref.value:n_yrs], na.rm=T)) # mean population growth for whole population
# (mean.adult.lambda_last.10 <- mean(pop.size.tibble$adult.lambda[80:n_yrs], na.rm=T)) # mean population growth for whole population)
# 
# lam_truth <- round(mean.adult.lambda, 4)
# lam_min <- min(adult.lambda[ref.value:n_yrs]) #Minimum lambda over estimation period
# lam_max <- max(adult.lambda[ref.value:n_yrs]) #Maximum lambda over estimation period

#psi
psi_truth.vec <- NULL
#Check proportion of conformists and non-conformists that are breeding each year
for(i in ref.value:n_yrs){
  Newmoms.vec <- loopy.list[[i]] %>% dplyr::filter(age.x == 0) %>% pull(mother.x)
  
  # Pot.conf.moms <- loopy.list[[i]] %>% dplyr::filter(sex == "F", 
  #                                                    age.x == repro.age,
  #                                                    repro.strategy == "conformist")
  # 
  # 
  # Newmoms.conf.df <- loopy.list[[i]] %>% dplyr::filter(indv.name %in% Newmoms.vec, 
  #                                                      age.x == repro.age,
  #                                                      repro.strategy == "conformist")
  # 
  # Pot.Nonconf.moms <- loopy.list[[i]] %>% dplyr::filter(sex == "F", 
  #                                                       age.x == repro.age,
  #                                                       repro.strategy == "non-conformist")
  # 
  # 
  # Newmoms.Nonconf.df <- loopy.list[[i]] %>% dplyr::filter(indv.name %in% Newmoms.vec, 
  #                                                         age.x == repro.age,
  #                                                         repro.strategy == "non-conformist")
  
  moms.all.df <- loopy.list[[i]] %>% dplyr::filter(indv.name %in% Newmoms.vec)
  
  psi_truth.vec[i] <- moms.all.df %>% group_by(repro.strategy) %>% 
    summarize(n = n()) %>% 
    mutate(perc.conf = round(n/sum(n), 2)) %>% #percent conformists
    dplyr::filter(repro.strategy == "conformist") %>% 
    pull(perc.conf)
}

psi_truth <- round(mean(psi_truth.vec[ref.value:n_yrs]), 2)

#Add the calculated values to ref.tibble2
ref.tibble2$survival[r] <- surv_mean
ref.tibble2$surv_AM[r] <- surv_mean_AM
ref.tibble2$surv_min[r] <- surv_min
ref.tibble2$surv_max[r] <- surv_max
ref.tibble2$psi[r] <- psi_truth
}

true.values <- ref.tibble2 %>% pivot_longer(cols = c(survival, psi), names_to = "parameter", values_to = "all.truth") %>% 
  mutate(iteration = iter,
         seed = rseed)

#Create dataframe of estimates and truth
# true.values <- tibble(parameter = c("Nf", "psi", "Nm", "survival", "lambda"),
#                     all.truth = c(Mom.all_truth, psi_truth, Dad.all_truth, surv_mean, lam_truth),
#                     breed.truth = c(Mom.breed_truth, psi_truth, Dad.breed_truth, surv_mean, lam_truth),
#                     mean.adult.lambda = mean.adult.lambda,
#                     mean.adult.lambda_last.10 = mean.adult.lambda_last.10,
#                     female.fecundity = ff,
#                     iteration = iter,
#                     seed = rseed,
#                     population.growth = population.growth)


# if(exists("ff.shift") == TRUE){
#   true.values$female.fecundity = ff.shift
# }
# 
# if(population.growth == "lambda.extreme"){
#   true.values <- true.values %>% mutate(adult.survival.adj = Sa,
#                                     juvenile.survival.adj = Sj)
#   }