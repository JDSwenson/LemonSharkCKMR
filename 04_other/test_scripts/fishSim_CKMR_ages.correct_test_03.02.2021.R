#rm(list=ls())
#Sex-specific estimates of N;
#Population growth fixed to observed population growth of the whole population.
#set so R doesn't use scientific notation
options("scipen"=100, "digits"=4)

library(fishSim)
library(foreach)
library(parallel)
library(doParallel)
library(optimx)
#Load individual packages because tidyverse won't load on cluster
library(plyr)
library(tidyverse)
library(popbio)
library(mpmtools)

set.seed(47)
#----------Simulation parameters-----------
t_start = 61 #First year we take samples
t_end = 65 #Last year we take samples
n_samp_yr = (t_end-t_start)+1 #Number of years being sampled
sim_yrs = c(1:60) #Number of pre-sample simulation years
n_yrs <- t_end #Total years of simulation + sampling

#----------------Set parameters for altMate-------------------------
firstBreed <- 12 #Using same knife-edge maturity for males and females; want age 12 to be the first year of breeding
maxAge <- 30 #Age after which animals die

#Set male and female curves. They're the same at present.
mat_m = mat_f <- rep(0,maxAge) #creates an empty vector for each age that will be used for maturity
mat_m[firstBreed+1:maxAge] = mat_f[firstBreed+1:maxAge] <- 1  #Knife-edge maturity for males and females, beginning at age 12
maleCurve <- c(mat_m, rep(max(mat_m), 1000)) #Set male age at maturity so breeders are not exhausted
femaleCurve <- c(mat_f, rep(max(mat_m), 1000)) #Set female age at maturity
maxClutch <- 18 #Maximum number of offspring one individual can possibly produce
mat_type <- "ageSex" #Set mat_type for altMate call. Could just do age but shouldn't change anything.

#-----------------Leslie Matrix--------------------
#Make a Leslie Matrix to find parameters that produce zero population growth and extract the stable age distribution

batchSize <- 4 #Average size of clutch (8) divided by two for skipped breeding
fb <- batchSize/2 #Litter size for Leslie Matrix - females only

Leslie_input <- data.frame(
  x = c(0:maxAge), #age
  sx = c(0.6, 0.7, 0.81, rep(0.85, times = 27), 0), #survival
  mx = c(rep(0, times = 12), rep(fb, times = 19)) #age-specific birth rates (female proportion of the population only)
)

A1_pre <- make_Leslie_matrix(Leslie_input) #Make pre-census Leslie Matrix
#View(A1)
A1_post <- pre_to_post(Amat = A1_pre, S0 = .6) #Convert to post-censue Leslie Matrix, using age 0 survival of 0.6; ignore warning - top row does contain fecundity coefficients.
#View(A1_post)

#Calculate dominant eigenvalue (i.e. population growth rate) from transition matrix
lambda1(A1_pre)
lambda1(A1_post) #Make sure the post-census and pre-census values match

#stable_age <- mpmtools::stable_stage(A1_post)
stable_age <- mpmtools::stable_stage(A1_pre) #to set a stable age distribution in fishSim with makeFounders, need to use a pre-breeding census i.e. start with age 1

survCurv <- stable_age #Sets probability of founder cohort belonging to each age class
osr <- c(0.5, 0.5) #Operating sex ratio of 50/50
stocks <- c(1) #stocks (just one)


#-------------------Set parameters for mort-------------------------
Surv_age <- Leslie_input$sx #Survival values from Leslie Matrix
ageMort <- (1 - Surv_age) #Subtract survival values from 1 to get mortality for mort
death_type <- "age"

pop <- 10000 #Set starting population size for makeFounders
perc_adult <- sum(stable_age[12:30]) #Percent of population that is mature i.e. that will be estimated by CKMR
N_a <- round(pop*perc_adult, 0) #Set starting parameter for number of adults (N_a)
N_f <- round(N_a/2, 0) #Set starting parameter value for number of adult females (N_f)
N_m <- round(N_a-N_f, 0) #Set starting parameter value for number of adult females (N_f)
Pars1 <- c(log(N_f),log(N_m)) #Parameters for model 1 (sex-specific)
Pars2 <- log(N_a) #Parameters for model 2 (sex-aggregated)


iterations <- 20

for(samps in 1:4){
  set.seed(100)
  results <- NULL
  #results <- data.frame(matrix(0, nrow = iterations*2, ncol = 7))
  n_samples <- c(60, 90, 120, 150)[samps]
  age_dist <- NULL #Reset the age_dist dataframe
  survival_at_age <- NULL
  
  for(iter in 1:iterations) {
    #set.seed(iter) #set the seed to the same ten numbers so the results are repeatable.
    index <- (iter*2)-1 #Index for filling in the matrix at the end of each loop over n_samples
    #set.seed(47)
    #-------------Loop start-------------------------    
    #makeFounders creates a matrix of specified size (pop) where each row is an individual in the founder popualtion.
    indiv <- makeFounders(pop = pop, osr = osr, stocks = c(1), maxAge = maxAge, survCurv=survCurv)
    
    #PoNG(mateType = "ageSex", femaleCurve = femaleCurve, mortType = "age", ageMort = ageMort, batchSize = batchSize, firstBreed = 0, maxClutch = maxClutch)
    
    #Simulate 40 years of mating
    ##This loop can take a minute depending on population size
    
    Dad_truth <- c()
    Mom_truth <- c()
    All_truth <- c()
    
    #-----------Proportion-at-age-----------------
    age_p <- indiv %>% filter(is.na(DeathY) == TRUE) %>% 
      group_by(AgeLast) %>% 
      dplyr::summarize(yr1 = n()/pop)
    
    #-------------Survival-at-age--------------
    #Numbers alive-at-age for age-specific survival
    alive_at_age <- indiv %>% filter(is.na(DeathY) == TRUE) %>% 
      group_by(AgeLast) %>% 
      dplyr::summarize(N_yr1 = n())
    
    #-----------------Start fishSim loop--------------------  
    for (y in 1:length(sim_yrs)) {
      year <- y
      #Store true values for each year
      All_truth[year] <- nrow(indiv[is.na(indiv[,6]),]) ## the currently-alive population size.
      
      Dad_truth[year] <- indiv %>% filter(Sex == "M" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>%
        dplyr::summarize(n())
      
      Mom_truth[year] <- indiv %>% filter(Sex == "F" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
        dplyr::summarize(n())
      
      ###Proportion-at-age###
      #Save age distribution at the beginning of each year. Column corresponds to age; iteration is added later
      ap <- indiv %>% filter(is.na(DeathY) == TRUE) %>% 
        group_by(AgeLast) %>% 
        dplyr::summarize(age_prop = n()/All_truth[year]) 
      
      names(ap)[2] <- paste0("age_prop_yr", year)
      
      age_p <- full_join(age_p, ap, by = "AgeLast") #Do this separate so columns are added in the right order
      
      ###Survival-at-age### 
      #Numbers alive at age after mating but before dying - to calculate survival
      aa <- indiv %>% filter(is.na(DeathY) == TRUE) %>% 
        group_by(AgeLast) %>% 
        dplyr::summarize(N = n()) 
      
      names(aa)[2] <- paste0("N_yr", year)
      
      alive_at_age <- full_join(alive_at_age, aa, by = "AgeLast") #Do this separate so columns are added in the right order
      
      indiv <- altMate(indiv, batchSize = batchSize, fecundityDist = "poisson", year = y, type = mat_type, maxClutch = maxClutch, singlePaternity = FALSE, maleCurve = maleCurve, femaleCurve = femaleCurve)
      
      indiv <- mort(indiv, type = death_type, year=y, ageMort = ageMort, maxAge = maxAge)
      
      indiv <- birthdays(indiv)
    }
    
    age_p <- age_p[,-2] #Remove duplicate value from first year
    alive_at_age <- alive_at_age[,-3] #Remove duplicate value from first year
    colnames(alive_at_age)[2] <- "N_yr1"
    
    
    # #Repeats the above loop, assigning offspring, deaths, and birthdays, but this time also draws samples each year and records which animals were sampled in indiv
    for(y in c(t_start:t_end)) {
      indiv <- capture(indiv, n=n_samples, year = y, fatal = FALSE) #Capture individuals
      
      year <- y
      #Store true values for each year
      All_truth[year] <- nrow(indiv[is.na(indiv[,6]),]) ## the currently-alive population size.
      
      Dad_truth[year] <- indiv %>% filter(Sex == "M" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>%
        dplyr::summarize(n())
      
      Mom_truth[year] <- indiv %>% filter(Sex == "F" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
        dplyr::summarize(n())
      
      ###Proportion-at-age###
      # #Save age distribution at the beginning of each year. Column corresponds to age; iteration is added later
      ap <- indiv %>% filter(is.na(DeathY) == TRUE) %>% 
        group_by(AgeLast) %>% 
        dplyr::summarize(age_prop = n()/All_truth[year]) 
      
      names(ap)[2] <- paste0("age_prop_yr", year)
      
      age_p <- full_join(age_p, ap, by = "AgeLast") %>% 
        mutate_all(round, 4) #Do this separate so columns are added in the right order
      
      ###Survival-at-age### 
      #Numbers alive at age after mating but before dying - to calculate survival
      aa <- indiv %>% filter(is.na(DeathY) == TRUE) %>% 
        group_by(AgeLast) %>% 
        dplyr::summarize(N = n()) 
      
      names(aa)[2] <- paste0("N_yr", year)
      
      alive_at_age <- full_join(alive_at_age, aa, by = "AgeLast") #Do this separate so columns are added in the right order
      
      indiv <- altMate(indiv, batchSize = batchSize, fecundityDist = "poisson", year = y, type = mat_type, maxClutch = maxClutch, singlePaternity = FALSE, maleCurve = maleCurve, femaleCurve = femaleCurve)
      
      indiv <- mort(indiv, type = death_type, year=y, ageMort = ageMort, maxAge = maxAge)
      
      indiv <- birthdays(indiv)
    }
    #indiv <- capture(indiv, n=100, year = 59, fatal = FALSE) #Capture individuals
    #Store truth values as vectors (dplyr makes them a list)
    Mom_truth <- unlist(Mom_truth) 
    Dad_truth <- unlist(Dad_truth)
    Adult_truth <- Mom_truth + Dad_truth
    All_truth <- unlist(All_truth)
    
    pop_growth_all <- All_truth[2:65]/All_truth[1:64]
    pop_growth_F <- Mom_truth[2:65]/Mom_truth[1:64]
    pop_growth_M <- Dad_truth[2:65]/Dad_truth[1:64]
    
    #Mom_truth
    #Dad_truth
    
    pairs <- findRelativesPar(indiv = indiv, sampled = TRUE) #verified that the number of comparisons is what it should be
    POPs <- pairs[pairs$OneTwo == 1,1:2] ##Parent-Offspring pairs
    HSPs <- pairs[pairs$TwoTwo == 1,1:2] ##Half-Sibling pairs - verified from fishSim vignette
    non_POPs <- pairs[pairs$OneTwo == 0,1:2] ##pairs that are not POPs
    non_HSPs <- pairs[pairs$TwoTwo == 0,1:2] ##pairs that are not half-sibs
    
    # think of this dataframe as storing only information about the younger fish. Stores ID and birth year only (at present) for every individual in indiv, but renames the ID column to younger so it can be joined with HSPs_tbl below)
    youngerbirthyears <- indiv %>%
      select(Me, BirthY, Mum, Dad) %>% 
      rename("younger" = Me, "Young_sib_birth" = BirthY, "Young_sib_mom" = Mum, "Young_sib_dad" = Dad)
    
    #Create dataframe with IDs of sampled half-sibs and columns named appropriately for joins
    HSPs_tbl <- HSPs %>% 
      rename(
        Me = Var1,
        younger = Var2) %>%
      mutate_all(as.character)
    
    #Extract the values of indiv for the older individual in each comparison.
    #Join all the information from indiv with the IDs of the sampled individuals in HSPs_2_tbl. 
    HSPs_2_tbl <- inner_join(indiv, HSPs_tbl, by = "Me") %>%
      left_join(youngerbirthyears, by = "younger")  %>%
      rename("Old_sib_birth" = BirthY, "Old_sib_mom" = Mum, "Old_sib_dad" = Dad) #%>%   select(c(Old_sib_birth, Young_sib_birth, Old_sib_mom, Young_sib_mom, Old_sib_dad, Young_sib_dad))
    
   
    ####Try removing individuals from same cohort comparisons####
    same_cohort_HSPs <- HSPs_2_tbl %>% filter(Old_sib_birth == Young_sib_birth)
    same_cohort_HSPs <- data.frame(Me = c(same_cohort_HSPs$Me, same_cohort_HSPs$younger))
    
    indiv2 <- indiv %>% anti_join(same_cohort_HSPs, by = "Me")
    
    pairs <- findRelativesPar(indiv = indiv2, sampled = TRUE) #verified that the number of comparisons is what it should be
    POPs <- pairs[pairs$OneTwo == 1,1:2] ##Parent-Offspring pairs
    HSPs <- pairs[pairs$TwoTwo == 1,1:2] ##Half-Sibling pairs - verified from fishSim vignette
    non_POPs <- pairs[pairs$OneTwo == 0,1:2] ##pairs that are not POPs
    non_HSPs <- pairs[pairs$TwoTwo == 0,1:2] ##pairs that are not half-sibs
    
    # think of this dataframe as storing only information about the younger fish. Stores ID and birth year only (at present) for every individual in indiv, but renames the ID column to younger so it can be joined with HSPs_tbl below)
    youngerbirthyears <- indiv %>%
      select(Me, BirthY, Mum, Dad) %>% 
      rename("younger" = Me, "Young_sib_birth" = BirthY, "Young_sib_mom" = Mum, "Young_sib_dad" = Dad)
    
    #Create dataframe with IDs of sampled half-sibs and columns named appropriately for joins
    HSPs_tbl <- HSPs %>% 
      rename(
        Me = Var1,
        younger = Var2) %>%
      mutate_all(as.character)
    
    #Extract the values of indiv for the older individual in each comparison.
    #Join all the information from indiv with the IDs of the sampled individuals in HSPs_2_tbl. 
    HSPs_2_tbl <- inner_join(indiv, HSPs_tbl, by = "Me") %>%
      left_join(youngerbirthyears, by = "younger")  %>%
      rename("Old_sib_birth" = BirthY, "Old_sib_mom" = Mum, "Old_sib_dad" = Dad) %>%   select(c(Old_sib_birth, Young_sib_birth, Old_sib_mom, Young_sib_mom, Old_sib_dad, Young_sib_dad))
    
    #Split HSPs dataframe into MHS and PHS pairs
    #Filter out intra-cohort comparisons
    mom_positives <- HSPs_2_tbl[HSPs_2_tbl$Old_sib_mom == HSPs_2_tbl$Young_sib_mom,] %>% 
      select(Old_sib_birth, Young_sib_birth) %>% 
      filter(Young_sib_birth != Old_sib_birth) %>% 
      plyr::count() %>% 
      select(Old_sib_birth, Young_sib_birth, freq)
    
    
    dad_positives <- HSPs_2_tbl[HSPs_2_tbl$Old_sib_dad == HSPs_2_tbl$Young_sib_dad,] %>% 
      select(Old_sib_birth, Young_sib_birth) %>% 
      filter(Young_sib_birth != Old_sib_birth) %>% 
      plyr::count() %>% 
      select(Old_sib_birth, Young_sib_birth, freq)
    
    #Checked that below matches above 02/03/2021
    parent_positives <- HSPs_2_tbl[HSPs_2_tbl$Old_sib_dad == HSPs_2_tbl$Young_sib_dad | HSPs_2_tbl$Old_sib_mom == HSPs_2_tbl$Young_sib_mom,] %>% 
      select(Old_sib_birth, Young_sib_birth) %>% 
      filter(Young_sib_birth != Old_sib_birth) %>% 
      plyr::count() %>% 
      select(Old_sib_birth, Young_sib_birth, freq)
    
    # Rename the columns so that the join functions know what to work with. The column `Me` is still the actual id column. However, we also need the column `younger` so that we can look up the 
    #birth years stored in `youngerbirthyears`
    non_HSPs_tbl <- non_HSPs %>% 
      rename(
        Me = Var1,
        younger = Var2) %>%
      mutate_all(as.character) # convert columns to characters since not all levels of indiv$Me are present in non_POPs_tbl
    #head(non_HSPs_tbl)
    
    ##Extract the values of indiv for the older individual in each comparison.
    #Join all the information from indiv with the IDs of the sampled individuals in non_HSPs_2_tbl. 
    non_HSPs_2_tbl <- inner_join(indiv, non_HSPs_tbl, by = "Me") %>%
      left_join(youngerbirthyears, by = "younger")  %>%
      rename(Old_sib_birth = BirthY) %>% 
      select(c(Old_sib_birth, Young_sib_birth))
    
    #Create table of all negative comparisons, count occurences, and filter intracohort comparisons
    mom_negatives = dad_negatives = parent_negatives <- non_HSPs_2_tbl %>% 
      plyr::count() %>% 
      filter(Young_sib_birth != Old_sib_birth)
    
    #-------------Kinship probabilities - Half-sib-------------------
    #Set the first cohort to the first year for which we have data
    min_est_cohort <- min(non_HSPs_2_tbl$Old_sib_birth, HSPs_2_tbl$Old_sib_birth)
    
    #Not assessing survival - keeping constant at 0.85
    #surv <- 0.85
    #max_age = 30 #max age of lemon sharks
    m_adult_age <- c(12:maxAge) #Set ages at which males and females are mature.
    f_adult_age <- c(12:maxAge)
    #pop_growth_F_mean <- mean(pop_growth_F[min_est_cohort_F:length(pop_growth_F)])
    #pop_growth_M_mean <- mean(pop_growth_M[min_est_cohort_M:length(pop_growth_M)])
    pop_growth_all_mean <- mean(pop_growth_all[min_est_cohort:length(pop_growth_all)])
    
    
    #Calculate mean survival-at-age
    mean_alive_at_age <- alive_at_age[1:30,] %>% rowMeans() %>% 
      round(digits = 0)
    #store mean adult survival
    mean_adult_surv <- mean(mean_alive_at_age[13:30]/mean_alive_at_age[12:29], na.rm = TRUE)
    
    ###Used in model###
    #surv <- 0.85
    surv <- mean_adult_surv
    #lam <- lambda1(A1_pre)
    lam <- pop_growth_all_mean
    
    
    #Wrong way of doing it ...
    #What is the first year for which we have a positive comparison for CKMR
    #CKMR estimates abundance for the younger cohort of each comparison, so we use the younger sibling birth year
    #min_est_cohort_F <- min(mom_positives$Young_sib_birth)
    #min_est_cohort_M <- min(dad_positives$Young_sib_birth)
    #min_est_cohort <- min(min_est_cohort_F, min_est_cohort_M)
    
    #Right way of doing it
    
    #For PC
    source("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/functions/get_P_lemon_HS_sex-specific.R")
    source("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/functions/lemon_neg_log_lik_HS_sex-specific.R")
    source("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/functions/get_P_lemon_HS_sex-aggregated.R")
    source("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/functions/lemon_neg_log_lik_HS_sex-aggregated.R")
    
    # #For cluster
    # source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/get_P_lemon_HS.R")
    # source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/lemon_neg_log_lik_HS_Sex_specific.R")
    # source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/get_P_lemon_HS_TotalA.R")
    # source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/lemon_neg_log_lik_HS_TotalA.R")
    
    #Fit model
    
    #Fit model - optimx version
    CK_fit1 <- optimx(par=Pars1,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Mother=mom_negatives, Negatives_Father=dad_negatives, Pairs_Mother=mom_positives, Pairs_Father=dad_positives, P_Mother=P_Mother, P_Father=P_Father, t_start=t_start, t_end=t_end)
    
    CK_fit2 <- optimx(par=Pars2,fn=lemon_neg_log_lik_TotalA,hessian=TRUE, method="BFGS", Negatives_Parent=parent_negatives, Pairs_Parent=parent_positives, P_Parent=P_Parent, t_start=t_start, t_end=t_end)
    
    #summary(CK_fit1)
    #exp(CK_fit1[1:2])
    
    #compute variance covariance matrix - optimx
    D1=diag(length(Pars1))*c(exp(CK_fit1$p1[1]),exp(CK_fit1$p2[1])) #derivatives of transformations
    VC_trans1 = solve(attr(CK_fit1, "details")["BFGS" ,"nhatend"][[1]])
    VC1 = (t(D1)%*%VC_trans1%*%D1) #delta method
    SE1=round(sqrt(diag(VC1)),0)
    
    
    #compute variance covariance matrix - optimx
    D2=diag(length(Pars2))*exp(CK_fit2$p1[1]) #derivatives of transformations
    VC_trans2 = solve(attr(CK_fit2, "details")["BFGS" ,"nhatend"][[1]])
    VC2 = (t(D2)%*%VC_trans2%*%D2) #delta method
    SE2 = round(sqrt(diag(VC2)),0)
    
    
    #Combine above to make dataframe with truth and estimates side-by-side
    #store years from youngest sibling in comparisons to end of study
    yrs <- c(min(mom_positives$Young_sib_birth, dad_positives$Young_sib_birth, parent_positives$Young_sib_birth):t_end)
    
    estimates <- data.frame(cbind(round(exp(c(CK_fit1$p1[1], CK_fit1$p2[1], CK_fit2$p1[1])),0)), c(SE1, SE2), c("F", "M", "All"))
    estimates <- cbind(estimates, c(Mom_truth[min_est_cohort], Dad_truth[min_est_cohort], Adult_truth[min_est_cohort]))
    colnames(estimates) <- c("CKMR_estimate", "SE", "Sex", "Truth")
    estimates

    total_samples <- n_samples*n_samp_yr
        
    metrics <- cbind(c(sum(mom_positives[,3]), sum(dad_positives[,3]), sum(parent_positives[,3])), c(rep(lam, times = 3)), c(rep(total_samples, times=3)))
    colnames(metrics) <- c("Parents_detected", "Pop_growth", "Samples")
    
    #-----------------Loop end-----------------------------    
    results <- rbind(results, cbind(estimates, metrics))
    
    #     #Remove any random age 31 individuals
    age_p <- age_p %>%
      filter(AgeLast <= 30)
    
    ##Transpose proportion-at-age so ages are different columns
    age_pT <- data.frame(t(age_p[,-1]))
    colnames(age_pT) <- paste0("Age_", unlist(age_p[,1]))
    age_pT$iter <- iter
    
    #Store proportion-at-age for each year of the fishSim simulation for each iteration of the loop.
    age_dist <- rbind(age_dist, age_pT)
    
    
    ##calculate survival for each age
    alive_at_age2 <- data.frame(alive_at_age[-31,-1])
    
    obs_survival <- data.frame()
    for(i in 1:(nrow(alive_at_age2)-1)){
      for(j in 1:(ncol(alive_at_age2)-1)){
        obs_survival[i, j] <- round(alive_at_age2[i+1, j+1]/alive_at_age2[i,j], digits=2)
      }
    }
    
    obs_survivalT <- data.frame(t(obs_survival))
    colnames(obs_survivalT) <- paste0("Age_", c(1:29))
    obs_survivalT$iter <- iter
    survival_at_age <- rbind(survival_at_age, obs_survivalT)
    #rownames(survival_at_age) <- paste0("sim_yr_", c(1:(n_yrs-1)))
    
    print(paste0("finished iteration", iter, " at: ", Sys.time()))
  }
  results <- results %>% 
    mutate(Relative_bias = round(((CKMR_estimate - Truth)/Truth)*100,1))
  
   results %>% group_by(Sex) %>% 
     dplyr::summarize(median = median(Relative_bias), n = n())
  
  write.table(results, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/fishSim_null", total_samples, ".samples_02.10.2021_ages.correct_fishSim.N.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

    write.table(age_dist, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/fishSim_age.distributions_", total_samples, ".samples_02.10.2021_ages.correct_age.dist.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

    write.table(survival_at_age, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/fishSim_survival.at.age_", total_samples, ".samples_02.10.2021_ages.correct_surv.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
}



