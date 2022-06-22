########## Compile and report results #########
#---------------Calculate psi truth----------------
  #Calculate breeding interval
  #Initialize sample dataframes
  # BI.df <- NULL
  # BI.df_temp <- NULL
  # 
  # ref.years <- c(75:n_yrs)
  # 
  # for(i in ref.years){ #Extract all YOY for each year sampled
  #   BI.df_temp <- loopy.list[[i]] %>% mutate(capture.year = i) %>% 
  #     dplyr::filter(age.x == 0)
  #   
  #   #Combine all
  #   BI.df <- bind_rows(BI.df, BI.df_temp)
  # }
  # 
  # #Extract unique mothers for each birth year
  # BI.df <- BI.df %>% as_tibble() %>% 
  #   arrange(birth.year) %>% 
  #   distinct(mother.x, birth.year, .keep_all = TRUE)
  # 
  # #-----------------Identify all relatives in population since estimation year and find breeding interval--------------------#
  # pairwise.BI.df <- data.frame(t(combn(BI.df$indv.name, m=2)))  # generate all combinations of the elements of x, taken m at a time.
  # colnames(pairwise.BI.df) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined
  # 
  # #Create dataframe of pairwise comparisons with just individual IDs
  # head(BI.df)
  # #head(pairwise.df)
  # 
  # #Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
  # OlderSib_birthyears.BI <- BI.df %>%
  #   dplyr::select(older.sib = indv.name, 
  #                 older.sib.birth = birth.year, 
  #                 older.sib.age = age.x, 
  #                 older.sib.mom = mother.x,
  #                 older.sib.dad = father.x)
  # 
  # #Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
  # #This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
  # pairwise.BI.all.info <- pairwise.BI.df %>% left_join(OlderSib_birthyears.BI, by = "older.sib") %>% 
  #   as_tibble() %>%  
  #   left_join(BI.df, by = "indv.name") %>% 
  #   dplyr::rename("younger.sib" = indv.name, 
  #                 "younger.sib.birth" = birth.year, 
  #                 "younger.sib.age" = age.x, 
  #                 "younger.sib.mom" = mother.x, 
  #                 "younger.sib.dad" = father.x) %>%
  #   dplyr::select(older.sib, 
  #                 older.sib.birth, 
  #                 older.sib.age, 
  #                 younger.sib, 
  #                 younger.sib.birth, 
  #                 younger.sib.age, 
  #                 older.sib.mom, 
  #                 younger.sib.mom, 
  #                 older.sib.dad, 
  #                 younger.sib.dad) %>% 
  #   dplyr::filter(older.sib.birth != younger.sib.birth)
  # 
  # 
  # #Extract all positive maternal half-sib comparisons from sampled years
  # positives.mom.BI <- pairwise.BI.all.info %>% filter(older.sib.mom == younger.sib.mom) %>% 
  #   mutate(yr.gap = younger.sib.birth - older.sib.birth) %>% 
  #   mutate(BI = ifelse(yr.gap %% 2 == 0, "even", "odd")) %>% #If the year gap between ha;f-sibs is divisible by 2, then call the breeding interval (BI) "even"
  #   dplyr::distinct(older.sib.birth, younger.sib.birth, younger.sib.mom, .keep_all = TRUE) #Duplicative of earlier; making sure cohort size isn't biasing things here
  # 
  # #Everything matches up re: positives and 
  # positives.mom.BI %>% group_by(BI) %>% summarize(n()) #How many positive comparisons for 
  # 
  # skipped.moms <- positives.mom.BI %>% dplyr::filter(BI == "even") %>% 
  #   dplyr::distinct(younger.sib.mom) %>% 
  #   pull(younger.sib.mom) #Should be the same as the older sib mom
  # 
  # annual.moms <- positives.mom.BI %>% dplyr::filter(BI == "odd") %>% 
  #   dplyr::distinct(younger.sib.mom) %>% 
  #   pull(younger.sib.mom) #Should be the same as the older sib mom
  # 
  # #Which moms ONLY reproduced in even years?
  # skipped.only.moms <- skipped.moms[which(!skipped.moms %in% annual.moms)]
  # all.moms <- unique(positives.mom.BI$younger.sib.mom)
  # 
  # #Percent of individuals that only breed when years are evenly spaced
  # (psi.truth <- length(skipped.only.moms)/length(all.moms))
  
  
  

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
surv_mean <- round(mean(sVec[ref.year:n_yrs]), 4) # True adult survival over estimation period
surv_min <- min(sVec[estimation.year:n_yrs]) #Minimum survival over estimation period
surv_max <- max(sVec[estimation.year:n_yrs]) #Maximum survival over estimation period

#Lambda
lam_truth <- round(mean.adult.lambda, 4)
lam_min <- min(adult.lambda[ref.year:n_yrs]) #Minimum lambda over estimation period
lam_max <- max(adult.lambda[ref.year:n_yrs]) #Maximum lambda over estimation period

#psi
psi_truth <- 1-non.conformists

#Create dataframe of estimates and truth
estimates <- tibble(parameter = c("Nf", "psi", "Nm", "survival", "lambda"),
                    all.truth = c(Mom.all_truth, psi_truth, Dad.all_truth, surv_mean, lam_truth),
                    breed.truth = c(Mom.breed_truth, psi_truth, Dad.breed_truth, surv_mean, lam_truth),
                    mean.adult.lambda = mean.adult.lambda,
                    female.fecundity = ff,
                    iteration = iter,
                    seed = rseed,
                    population.growth = population.growth)
