#rm(list=ls())
#Sex-specific estimates of N;
#Population growth fixed to observed population growth of the whole population.
#set so R doesn't use scientific notation
options("scipen"=100, "digits"=4)

library(fishSim)
library(foreach)
library(parallel)
library(doParallel)
#library(optimx)
#Load individual packages because tidyverse won't load on cluster
library(plyr)
library(dplyr)
library(tidyr)
library(popbio)
library(mpmtools)

set.seed(47)
#----------Simulation parameters-----------
t_start = 61 #First year we take samples
t_end = 65 #Last year we take samples
#min_est_cohort <- 42 #First year we are estimating abundance
#max_est_cohort <- t_end #Last year we are estimating abundance
#est_yrs <- min_est_cohort:max_est_cohort #All years for which we are estimating abundance
n_yrs <- t_end #Number of years of study

sim_yrs = c(1:60)
n_samp_yr = t_end-t_start +1 #Number of years being sampled
#n_samples_per_yr = 75
#n_samples_total <- n_samples_per_yr * n_samp_yr + n_samples_per_yr #because sampling starts in year t_start, we need to add one more year

#----------------Leslie Matrix and makeFounder parameters---------------
# s_YOY = 0.6    # YOY survival
# s_adult = 0.85  # Adult (mature) survival 
# 
# #Change length.out to ramp up from s_YOY to s_adult, and then edit times so the total number of survival rates is 11
# s_juv = c(round(seq(s_YOY, s_adult, length.out = 3), 2), rep(s_adult, times = 8))    # Juvenile survival - ramp up from YOY survival to adult survival; with 13 columns corresponding to years 0 (YOY) - 12 (mature adults), we set 11 different survival rates, starting with the YOY survival and moving to adult. The for loop in the function below starts at index 2, so doesn't repeat YOY survival
# a_m = 12         # Age at sexual maturity
# b_max = 1.5     #Females recruited into the population. Assume equal sex ratio. Calculated from: average fecundity (6) divided by two (for skipped-breeding) divided by two (for females vs males)

# age_x <- max_age
# 
# 
# init_matrix = function(s_j, s_a, a_mat, age_xx, b_maxx){
#   A_matrix = matrix(data = 0, nrow = age_xx, ncol = age_xx) # dimension matrix and fill with zeros
#   A_matrix[1, (a_mat:age_xx)] = b_maxx # assign birth rates to mature ages in first row of matrix
#   A_matrix[2,1] <- s_YOY
#   for(ii in 2:(a_mat-1)) A_matrix[ii+1, ii] = s_j[ii] # assign juvenile survival rates
#   #A_matrix[age_xx,a_mat] = s_a # adult survival assumed for maturing animals transitioning into plus-group
#   for(jj in a_mat:(age_xx-1)) A_matrix[jj+1, jj] = s_adult # assign juvenile
#   A_matrix[age_xx, age_xx] = s_a # adult survival assumed for plus-group
#   return(A_matrix)
# }
# 
# A = init_matrix(s_j = s_juv, s_a = s_adult, a_mat = a_m, 
#                 age_xx = age_x, b_maxx = b_max) # Call function to initialize matrix A
#View(A)


########### Old way of calculating "stable age distribution"#############
#Surv_age <- rep(0.85,maxAge) #Set survival for each age to 85%
#Surv_age[1] <- 0.6 #Set survival for first year to 60%
#Surv_age[2] <- 0.72 #Set survival for second year to 70%
#Prop_age <- rep(1:maxAge)
#for(iage in 2:maxAge)Prop_age[iage]=Prop_age[iage-1]*Surv_age[iage-1]
#Prop_age <- Prop_age/sum(Prop_age)  #stable age distribution ... supposedly
#survCurv <- Prop_age #Sets probability of founder cohort belonging to each age class


#----------------Set parameters for altMate-------------------------
firstBreed <- 12 #If using same knife-edge maturity for males and females
maxAge <- 30

mat_m <- rep(0,maxAge) #creates an empty vector that has the length of n_ages
mat_m[firstBreed:maxAge]=1  #Knife-edge maturity for males, beginning at age 12

mat_f <- rep(0,maxAge) #creates an empty vector that has the length of n_ages
mat_f[firstBreed:maxAge]=1  #Knife-edge maturity for females, beginning at age 12

maleCurve <- c(mat_m, rep(max(mat_m), 1000)) #Set male age at maturity
femaleCurve <- c(mat_f, rep(max(mat_m), 1000)) #Set female age at maturity
maxClutch <- 18 #Maximum number of offspring one individual can possibly produce
mat_type <- "ageSex"

#-----------------Leslie Matrix--------------------
batchSize <- 4
fb <- batchSize/2 #Batch size for Leslie Matrix - females only

Leslie_input <- data.frame(
  x = c(0:maxAge), #age
  sx = c(0.6, 0.7, 0.81, rep(0.85, times = 27), 0), #survival
  mx = c(rep(0, times = 12), rep(fb, times = 19)) #age-specific birth rates (female proportion of the population only)
)

A1_pre <- make_Leslie_matrix(Leslie_input)
#View(A1)
A1_post <- pre_to_post(Amat = A1_pre, S0 = .6)
#View(A1_post)

#Calculate dominant eigenvalue (i.e. population growth rate) from transition matrix
(lam <- lambda1(A1_pre))
lambda1(A1_post)
#stable_age <- mpmtools::stable_stage(A1_post)

stable_age <- mpmtools::stable_stage(A1_pre) #fishSim seems to use a pre-breeding census for the age distribution

#Calculate stable age structure of the population
# stable_age = stable.stage(A)
# stable_age
survCurv <- stable_age #Sets probability of founder cohort belonging to each age class

Surv_age <- Leslie_input$sx
#rep(0.85,maxAge+1) #Set survival for each age to 85%
#Surv_age[1] <- .6 #Set survival for first year to 60%
#Surv_age[2] <- .72 #Set survival for second year to 72%
osr <- c(0.5, 0.5) #Operating sex ratio
stocks <- c(1) #stocks (empty for now)


#-------------------Set parameters for mort-------------------------
ageMort <- (1 - Surv_age)
#ageMort[maxAge+1] <- ageMort[maxAge]
death_type <- "age"

pop <- 5000
perc_adult <- sum(stable_age[12:30])
N_a <- round(pop*perc_adult, 0)
N_f <- round(N_a/2, 0)
N_m <- round(N_a-N_f, 0)
Pars1 <- c(log(N_f),log(N_m))
Pars2 <- log(N_a)

#Check growth rate to find first year survival rate that produces 0 population growth; store in a variable that will be fixed in the model to estimate N.
#pop_growth <- as.numeric(check_growthrate(mateType = "ageSex", femaleCurve = femaleCurve, mortType = "age", ageMort = ageMort, batchSize = batchSize, maxClutch = maxClutch, maxAge = maxAge))

iterations <- 100

for(samps in 1:4){
  set.seed(47)
  results <- NULL
  #results <- data.frame(matrix(0, nrow = iterations*2, ncol = 7))
  n_samples <- c(30, 60, 90, 120)[samps]
  age_dist <- NULL #Reset the age_dist dataframe
  survival_at_age <- NULL
  
  for(iter in 1:iterations) {
    #set.seed(iter) #set the seed to the same ten numbers so the results are repeatable.
        index <- (iter*2)-1 #Index for filling in the matrix at the end of each loop over n_samples
    #set.seed(47)
#-------------Loop start-------------------------    
#makeFounders creates a matrix of specified size (pop) where each row is an individual in the founder popualtion.
    #makeFounders creates a matrix of specified size (pop) where each row is an individual in the founder popualtion.
    indiv <- makeFounders(pop = pop, osr = osr, stocks = c(1), maxAge = maxAge, survCurv=survCurv)

    #PoNG(mateType = "ageSex", femaleCurve = femaleCurve, mortType = "age", ageMort = ageMort, batchSize = batchSize, firstBreed = 0, maxClutch = maxClutch)
    
    #Simulate 40 years of mating
    ##This loop can take a minute depending on population size
    
    Dad_truth <- c()
    Mom_truth <- c()
    All_truth <- c()

    #-----------Proportion-at-age-----------------
    #Proportion-at-age row-wise
    #age_p <- data.frame(matrix(0, nrow = n_yrs, ncol = max_age + 1))
    #colnames(age_p) <- paste0("age_", c(0:max_age))
    
    age_p <- indiv %>% filter(is.na(DeathY) == TRUE) %>% 
      group_by(AgeLast) %>% 
      summarize(yr1 = n()/N_a)
    
    #Save age distribution at the beginning of each year. Column corresponds to age; iteration is added later
    #age_p[1,2:31] <- indiv %>% filter(is.na(DeathY) == TRUE) %>% 
    #  group_by(AgeLast) %>% 
    #  summarize(age_prop = n()/N_a) %>% 
    #  select(age_prop) %>% 
    #  t()
    
    #-------------Survival-at-age--------------
    #Numbers alive-at-age for age-specific survival
    alive_at_age <- indiv %>% filter(is.na(DeathY) == TRUE) %>% 
      group_by(AgeLast) %>% 
      summarize(N_yr1 = n())
    
    #-----------------Start fishSim loop--------------------  
    for (y in 1:length(sim_yrs)) {
        year <- y
        #Store true values for each year
        All_truth[year] <- nrow(indiv[is.na(indiv[,6]),]) ## the currently-alive population size.
        
        Dad_truth[year] <- indiv %>% filter(Sex == "M" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>%
          summarize(n())
        
        Mom_truth[year] <- indiv %>% filter(Sex == "F" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
          summarize(n())

      ###Proportion-at-age###
      #Save age distribution at the beginning of each year. Column corresponds to age; iteration is added later
      ap <- indiv %>% filter(is.na(DeathY) == TRUE) %>% 
        group_by(AgeLast) %>% 
        summarize(age_prop = n()/All_truth[year]) 
      
      names(ap)[2] <- paste0("age_prop_yr", year)
      
      age_p <- full_join(age_p, ap, by = "AgeLast") #Do this separate so columns are added in the right order
      
      ###Survival-at-age### 
      #Numbers alive at age after mating but before dying - to calculate survival
      aa <- indiv %>% filter(is.na(DeathY) == TRUE) %>% 
        group_by(AgeLast) %>% 
        summarize(N = n()) 
    
      names(aa)[2] <- paste0("N_yr", year)
      
        alive_at_age <- full_join(alive_at_age, aa, by = "AgeLast") #Do this separate so columns are added in the right order
      
      indiv <- altMate(indiv, batchSize = batchSize, fecundityDist = "poisson", year = y, type = mat_type, maxClutch = maxClutch, singlePaternity = FALSE, maleCurve = maleCurve, femaleCurve = femaleCurve)
      
      indiv <- mort(indiv, type = death_type, year=y, ageMort = ageMort, maxAge = maxAge)

      indiv <- birthdays(indiv)
    }

    age_p <- age_p[,-3] #Remove duplicate value from first year
    colnames(age_p)[2] <- "age_prop_yr1"
    alive_at_age <- alive_at_age[,-3] #Remove duplicate value from first year
    colnames(alive_at_age)[2] <- "N_yr1"

    
# #Repeats the above loop, assigning offspring, deaths, and birthdays, but this time also draws samples each year and records which animals were sampled in indiv
     for(y in c(t_start:t_end)) {
    indiv <- capture(indiv, n=n_samples, year = y, fatal = FALSE) #Capture individuals

        year <- y
    #Store true values for each year
    All_truth[year] <- nrow(indiv[is.na(indiv[,6]),]) ## the currently-alive population size.
    
    Dad_truth[year] <- indiv %>% filter(Sex == "M" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>%
      summarize(n())
    
    Mom_truth[year] <- indiv %>% filter(Sex == "F" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
      summarize(n())
    
    ###Proportion-at-age###
    #Save age distribution at the beginning of each year. Column corresponds to age; iteration is added later
    ap <- indiv %>% filter(is.na(DeathY) == TRUE) %>% 
      group_by(AgeLast) %>% 
      summarize(age_prop = n()/All_truth[year]) 
    
    names(ap)[2] <- paste0("age_prop_yr", year)
    
    age_p <- full_join(age_p, ap, by = "AgeLast") %>% 
      mutate_all(round, 4) #Do this separate so columns are added in the right order
    
    ###Survival-at-age### 
    #Numbers alive at age after mating but before dying - to calculate survival
    aa <- indiv %>% filter(is.na(DeathY) == TRUE) %>% 
      group_by(AgeLast) %>% 
      summarize(N = n()) 
    
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
    mis_vec <- c(0, 1, 2)

    indiv_mis <- indiv %>% rename(trueAgeLast = AgeLast) %>% 
      mutate(AgeLast_temp = trueAgeLast - sample(mis_vec, nrow(indiv), replace = TRUE)) %>% 
      mutate(AgeLast = ifelse(AgeLast_temp < 0, 0, AgeLast_temp)) %>% 
      select(Me, Sex, Dad, Mum, BirthY, DeathY, Stock, AgeLast, SampY, trueAgeLast) 
    
    pairs <- findRelativesPar(indiv = indiv_mis, sampled = TRUE)
    POPs <- pairs[pairs$OneTwo == 1,1:2] ##Parent-Offspring pairs
    HSPs <- pairs[pairs$TwoTwo == 1,1:2] ##Half-Sibling pairs - verified from fishSim vignette
    non_POPs <- pairs[pairs$OneTwo == 0,1:2] ##pairs that are not POPs
    non_HSPs <- pairs[pairs$TwoTwo == 0,1:2] ##pairs that are not half-sibs

    # think of this dataframe as storing only information about the younger fish. Stores ID and birth year only (at present) for every indiv2idual in indiv_mis, but renames the ID column to younger so it can be joined with HSPs_tbl below)
    youngerbirthyears <- indiv_mis %>%
      select(Me, BirthY, Mum, Dad) %>% 
      rename("younger" = Me, "Young_sib_birth" = BirthY, "Young_sib_mom" = Mum, "Young_sib_dad" = Dad)
    
    #Create dataframe with IDs of sampled half-sibs and columns named appropriately for joins
    HSPs_tbl <- HSPs %>% 
      rename(
        Me = Var1,
        younger = Var2) %>%
      mutate_all(as.character)
    
    #Extract the values of indiv_mis for the older indiv_misidual in each comparison.
    #Join all the information from indiv with the IDs of the sampled individuals in HSPs_2_tbl. 
    HSPs_2_tbl <- inner_join(indiv_mis, HSPs_tbl, by = "Me") %>%
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
    
    min_est_cohort_F <- min(mom_positives$Young_sib_birth)
    min_est_cohort_M <- min(dad_positives$Young_sib_birth)
    min_est_cohort <- min(min_est_cohort_F, min_est_cohort_M)

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
    non_HSPs_2_tbl <- inner_join(indiv_mis, non_HSPs_tbl, by = "Me") %>%
      left_join(youngerbirthyears, by = "younger")  %>%
      rename(Old_sib_birth = BirthY) %>% 
      select(c(Old_sib_birth, Young_sib_birth))
    
    #Create table of all negative comparisons, count occurences, and filter intracohort comparisons
    mom_negatives = dad_negatives = parent_negatives <- non_HSPs_2_tbl %>% 
        plyr::count() %>% 
        filter(Young_sib_birth != Old_sib_birth)
    
    #-------------Kinship probabilities - Half-sib-------------------
    #Not assessing survival - keeping constant at 0.85
    Surv <- 0.85
    max_age = 30 #max age of lemon sharks
    m_adult_age <- c(12:max_age) #Set ages at which males and females are mature.
    f_adult_age <- c(12:max_age)
    pop_growth_F_mean <- mean(pop_growth_F[min_est_cohort_F:length(pop_growth_F)])
    pop_growth_M_mean <- mean(pop_growth_M[min_est_cohort_M:length(pop_growth_M)])
    pop_growth_all_mean <- mean(pop_growth_all[min_est_cohort:length(pop_growth_all)])
    
    #Calculate mean survival-at-age
    mean_alive_at_age <- alive_at_age[1:30,] %>% rowMeans() %>% 
      round(digits = 0)
    #store mean adult survival
    mean_adult_surv <- mean(mean_alive_at_age[13:30]/mean_alive_at_age[12:29])
    
    #For PC
    # source("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/null/functions/get_P_lemon_HS.R")
    # source("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/null/functions/lemon_neg_log_lik_HS_Sex_specific.R")
    # source("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/null/functions/get_P_lemon_HS_TotalA.R")
    # source("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/null/functions/lemon_neg_log_lik_HS_TotalA.R")
    
    # #For cluster
    source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/get_P_lemon_HS.R")
    source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/lemon_neg_log_lik_HS_Sex_specific.R")
    source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/get_P_lemon_HS_TotalA.R")
    source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/lemon_neg_log_lik_HS_TotalA.R")

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
    
    total_samples <- n_samples * n_samp_yr
    
    metrics <- cbind(c(sum(mom_positives[,3]), sum(dad_positives[,3]), sum(parent_positives[,3])), c(rep(lam, times = 3)), c(rep(total_samples, times=3)))
    colnames(metrics) <- c("Parents_detected", "Pop_growth", "Samples")

    #-----------------Loop end-----------------------------    
    results <- rbind(results, cbind(estimates, metrics))

# results[index,] <- c(round(exp(CK_fit[1]),0), round(SE[1],0), round(Mom_truth[min_est_cohort_F],0), "F", sum(mom_positives[,3]), pop_growth_all_mean, n_samples)
# results[index+1,] <- c(round(exp(CK_fit[2]),0), round(SE[2],0), round(Dad_truth[min_est_cohort_M],0), "M", sum(dad_positives[,3]), pop_growth_all_mean, n_samples)
# colnames(results) <- c("N_est", "SE", "Truth", "Sex", "Parents_detected", "Pop_growth", "Samples")

    #Remove any random age 31 individuals
    age_p <- age_p %>% 
      filter(AgeLast <= 30)

    ##Transpose proportion-at-age so ages are different columns
age_pT <- data.frame(t(age_p[,-1]))
colnames(age_pT) <- paste0("Age_", age_p[,1])
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
   summarize(median = median(Relative_bias), n = n())
  
  write.table(results, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/fishSim_null", total_samples, ".samples_02.03.2021_ages.under_fishSim.N.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

    write.table(age_dist, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/fishSim_age.distributions_", total_samples, ".samples_02.03.2021_ages.under_age.dist.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

    write.table(survival_at_age, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/fishSim_survival.at.age_", total_samples, ".samples_02.03.2021_ages.under_surv.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
}