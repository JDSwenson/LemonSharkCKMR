#fishSim CKMR loop
#Constant abundance, one estimate
#Presently, estimating one abundance for males and one for females, and calling the "truth" the mean for the range of years covered by the estimate (aka earliest year for which we have an estimate)

#6/29/2020 update: made new indiv dataframe to select only cohorts from the years we're estimating. However, when we subset for only the years of interest, then we can't get an estimate for the first year ... need to think about this more and troubleshoot ... 

rm(list=ls())

library(fishSim)
library(foreach)
library(parallel)
library(doParallel)
library(optimx)
library(tidyverse)


##############Set up experiment parameters##############
#Used in prior probability matrix:
#min_est_cohort: sets the columns (younger birth year) that will be populated in the matrix
#est_ages: sets the index to refer to the appropriate value in the vector of parameter values
#max_age: helps determine which cells to populate in the prior probability matrix. 
#f_adult_age: minimum value is used in conjunction with max_age to determine which comparisons should have a non-zero probability.

####Time series - 3 years of estimates, 2 years of samples####
t_start = 39 #First year we take samples
t_end = 41 #Last year we take samples
min_est_cohort <- 40 #First year we are estimating abundance
max_est_cohort <- t_end #Last year we are estimating abundance
est_yrs <- min_est_cohort:max_est_cohort #All years for which we are estimating abundance
est_ages <- length(est_yrs)
n_yrs <- t_end #Number of years of study - used to determine dimensions of prior probability matrix
Surv <- 0.85 #adult survival
max_age = 30 #max age of lemon sharks
m_adult_age <- c(12:30) #Set ages at which males and females are mature.
f_adult_age <- c(13:30)
m_mat <- min(m_adult_age) #Age at maturity for males
f_mat <- min(f_adult_age) #Age at maturity for females


####set up the simulation parameters####
sim_yrs = c(1:(t_start-1))
#makeFounder parameters
osr <- c(0.5, 0.5) #Operating sex ratio
stocks <- c() #stocks (empty for now)
maxAge <- 30 #max age of lemon sharks
Surv_age <- rep(0.85,maxAge) #Set survival for each age to 85%
Surv_age[1] <- 0.55 #Set survival for first year to 55%
Surv_age[2] <- 0.65 #Set survival for second year to 65%
Prop_age <- rep(1:maxAge)
pop <- 5000

for(iage in 2:maxAge)Prop_age[iage]=Prop_age[iage-1]*Surv_age[iage-1]
Prop_age <- Prop_age/sum(Prop_age)  #stable age distribution
survCurv <- Prop_age #Sets probability of founder cohort belonging to each age class

#Set parameters for altMate
batchSize = 6 #average size of brood -- assumes Poisson distribution for fecundity (which is default for altMate)
mat_m <- rep(0,maxAge) #creates an empty vector that has the length of n_ages
mat_m[12:maxAge]=1  #Knife-edge maturity for males, beginning at age 12

mat_f <- rep(0,maxAge) #creates an empty vector that has the length of n_ages
mat_f[13:maxAge] <- 1  #Knife-edge maturity for females, beginning at age 13

maleCurve <- c(mat_m, rep(max(mat_m), 5000)) #Set male age at maturity
femaleCurve <- c(mat_f, rep(max(mat_f), 5000)) #Set female age at maturity
maxClutch <- 18 #Maximum number of offspring one individual can possibly produce
mat_type <- "ageSex"

#Set parameters for mort()
ageMort <- (1 - Surv_age)
ageMort[maxAge+1] <- ageMort[maxAge]
death_type <- "age"

total_abundance <- rep(600, times = length(est_yrs)) #Set a priori guess for prior total adult abundance
N_f <- c(total_abundance/2)
N_m <- c(total_abundance-N_f)
Pars=c(log(N_f),log(N_m))

####Start loop####
sim_start_time <- Sys.time()
print(paste0("Simulation started at ", Sys.time()))
iterations <- 100
n_samp_yr = t_end-t_start #Number of years being sampled

for(samps in 1:2){
  all_est <- data.frame() #Initialize dataframe for storing all iterations
  n_samples_per_yr = c(300, 500)[samps]
  n_samples <- n_samples_per_yr * n_samp_yr + n_samples_per_yr #because sampling starts in year t_start, we need to add one more year
  for(iter in 1:iterations) {

#makeFounders creates a matrix of specified size (pop) where each row is an individual in the founder popualtion.
indiv <- makeFounders(pop = pop, osr = osr, stocks = c(1), maxAge = maxAge, survCurv=survCurv)

#head(indiv)
#hist(as.numeric(indiv[,8])) #histogram of population ages

#Simulate 40 years of mating
#altMate derives the number of offspring produced by drawing from a sampling distribution for each female based on the parameters specified in the function. Fathers are randomly drawn from all mature males in the mother's stock.
##This loop can take a minute depending on population size

Dad_truth <- c()
Mom_truth <- c()

for (y in 1:length(sim_yrs)) {
  indiv <- altMate(indiv, batchSize = batchSize, fecundityDist = "poisson", year = y, type = mat_type, maxClutch = maxClutch, singlePaternity = FALSE, maleCurve = maleCurve, femaleCurve = femaleCurve, firstBreed = 0)
  
  #Store truth values
  Dad_truth[y] <- indiv %>% 
    filter(Sex == "M" &  is.na(DeathY)==TRUE & AgeLast >= m_mat) %>%
    summarize(n())
  
  Mom_truth[y] <- indiv %>% 
    filter(Sex == "F" & is.na(DeathY)==TRUE & AgeLast >= f_mat) %>%
    summarize(n())
  
  indiv <- mort(indiv, type = death_type, year=y, maxAge = maxAge, ageMort = ageMort)
  
  indiv <- birthdays(indiv)
}

#Repeats the above loop, in the following order: 1) breed, 2) extract true abundance of adults that bred, 3) sample, 4) kill individuals, 5) birthday - adds one to AgeLast column
for (y in c(t_start:t_end)) {
  indiv <- altMate(indiv, batchSize = batchSize, fecundityDist = "poisson", year = y, type = mat_type, maxClutch = maxClutch, singlePaternity = FALSE, maleCurve = maleCurve, femaleCurve = femaleCurve, firstBreed = 0)
 
  #Store truth values
  Dad_truth[y] <- indiv %>% 
    filter(Sex == "M" &  is.na(DeathY)==TRUE & AgeLast >= m_mat) %>%
    summarize(n())
  
  Mom_truth[y] <- indiv %>% 
    filter(Sex == "F" & is.na(DeathY)==TRUE & AgeLast >= f_mat) %>%
    summarize(n())
  
  indiv <- capture(indiv, n=n_samples_per_yr, year=y, fatal = FALSE) #Capture individuals
  
  indiv <- indiv %>% 
    mutate(SampY=replace(SampY, BirthY < (min_est_cohort-1), NA))
  
  indiv <- mort(indiv, type = death_type, year=y, maxAge = maxAge, ageMort = ageMort)
  
  indiv <- birthdays(indiv) #Age each individual by one year
}

#Store truth values as vectors (dplyr makes them a list)
Mom_truth <- unlist(Mom_truth) 
Dad_truth <- unlist(Dad_truth)

pairs <- findRelativesPar(indiv = indiv, sampled = TRUE)
POPs <- pairs[pairs$OneTwo == 1,1:2] ##Parent-Offspring pairs
HSPs <- pairs[pairs$TwoTwo == 1,1:2] ##Half-Sibling pairs - verified from fishSim vignette
non_POPs <- pairs[pairs$OneTwo == 0,1:2] ##pairs that are not POPs
non_HSPs <- pairs[pairs$TwoTwo == 0,1:2] ##pairs that are not half-sibs
#relatives <- namedRelatives(pairs)

# think of this dataframe as a renaming of indiv - specifically, the ID column is renamed "younger" so it can be joined with HSPs_tbl below)
youngerbirthyears <- indiv %>%
  select(Me, BirthY, Mum, Dad) %>% 
  rename("younger" = Me, "Young_sib_birth" = BirthY, "Young_sib_mom" = Mum, "Young_sib_dad" = Dad)

#Create dataframe with IDs of sampled half-sibs and columns named appropriately for joins. The second column is the younger individual in each comparison. We name this column "younger" so we can retrieve the appropriate information from the youngerbirthyears dataframe above. We name the first column Me so we can retrieve the appropriate information from the indiv dataframe.
HSPs_tbl <- HSPs %>% 
  rename(
    Me = Var1,
    younger = Var2) %>%
  mutate_all(as.character)

#Extract the values of indiv for the older individual in each comparison.
#Inner join returns all rows from x where there are matching values in y, so returns all rows of indiv that correspond with the older sib IDs stored in HSPs_tbl. Also returns the ID for all the younger individuals in the comparison.

#Join the information from above with the IDs of the sampled individuals in HSPs_2_tbl.
#left_join returns all rows from x and all columns from x and y, so grabs all the information from the renamed indiv dataframe (above) for all the younger sampled individuals. 
HSPs_2_tbl <- inner_join(indiv, HSPs_tbl, by = "Me") %>%
  left_join(youngerbirthyears, by = "younger")  %>%
  rename("Old_sib_birth" = BirthY, "Old_sib_mom" = Mum, "Old_sib_dad" = Dad) %>% 
  select(c(Old_sib_birth, Young_sib_birth, Old_sib_mom, Young_sib_mom, Old_sib_dad, Young_sib_dad))

#Split HSPs dataframe into MHS and PHS pairs
#Filter out 1) intra-cohort comparisons, 2) comparisons where the younger individual is outside the range of years we're estimating, & 3)
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


nrow(HSPs_2_tbl[HSPs_2_tbl$Old_sib_dad == HSPs_2_tbl$Young_sib_dad & HSPs_2_tbl$Old_sib_mom == HSPs_2_tbl$Young_sib_mom,]) #Should be 0 or there are full sibs

# Rename the columns so that the join functions know what to work with. The column `Me` is still the actual id column. However, we also need the column `younger` so that we can look up the 
#birth years stored in `youngerbirthyears`
non_HSPs_tbl <- non_HSPs %>% 
  rename(
    Me = Var1,
    younger = Var2) %>%
  mutate_all(as.character) # convert columns to characters since not all levels of indiv$Me are present in non_POPs_tbl
#head(non_HSPs_tbl)

#Extract the values of indiv for the older individual in each comparison.
#Inner join returns all rows from x where there are matching values in y, so returns all rows of indiv that correspond with the older sib IDs stored in non_HSPs_tbl. Also returns the ID for all the younger individuals in the comparison.

#Join the information from above with the IDs of the sampled individuals in non_HSPs_2_tbl.
#left_join returns all rows from x and all columns from x and y, so grabs all the information from the renamed indiv dataframe (above) for all the younger sampled individuals. 
non_HSPs_2_tbl <- inner_join(indiv, non_HSPs_tbl, by = "Me") %>%
  left_join(youngerbirthyears, by = "younger")  %>%
  rename(Old_sib_birth = BirthY) %>% 
  select(c(Old_sib_birth, Young_sib_birth))

#Create table of all negative comparisons, count occurences, and filter intracohort comparisons
mom_negatives <- non_HSPs_2_tbl %>%   
  filter(Young_sib_birth != Old_sib_birth) %>% 
  plyr::count()

dad_negatives <- non_HSPs_2_tbl %>%   
  filter(Young_sib_birth != Old_sib_birth) %>% 
  plyr::count()

##Define variables and functions
#source("~/R/R_working_dir/CKMR/LemonSharkCKMR/time_series/models/get_P_lemon_HS_time_series_4yrs.R") #For desktop
#source("~/R/R_working_dir/CKMR/LemonSharkCKMR/likelihood_functions/lemon_neg_log_like_HS_Sex_specific.R") #For desktop

source("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/time_series/models/get_P_lemon_HS_time_series_4yrs.R") #for laptop
source("~/R/R_working_dir/CKMR/LemonSharkCKMR_GitHub/likelihood_functions/lemon_neg_log_like_HS_Sex_specific.R") #For laptop

P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)

#P$P_Mother

#Check P function
#P$P_Mother[19,40,45] #Dimensions are parent birth year, parent capture year, offspring birth year (all of which are specified by n_yrs)

#Fit model - optimx
CK_fit <- optimx(par=Pars,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Mother=mom_negatives, Negatives_Father=dad_negatives, Pairs_Mother=mom_positives, Pairs_Father=dad_positives, P_Mother=P_Mother, P_Father=P_Father, n_yrs=n_yrs, t_start=t_start, t_end=t_end)

#summary(CK_fit)
exp(CK_fit[1:4])

#compute variance covariance matrix
#Estimating for two years
#For optimx
D=diag(length(Pars))*c(exp(CK_fit$p1[1]),exp(CK_fit$p2[1]), exp(CK_fit$p3[1]), exp(CK_fit$p4[1])) #derivatives of transformations

#Estimating for three years
#D=diag(length(Pars))*c(exp(CK_fit$p1[1]),exp(CK_fit$p2[1]), exp(CK_fit$p3[1]), exp(CK_fit$p4[1]), exp(CK_fit$p5[1]), exp(CK_fit$p6[1])) #derivatives of transformations

VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
VC = (t(D)%*%VC_trans%*%D) #delta method
SE=round(sqrt(diag(VC)),0)

#Rename columns - change according to years estimating
colnames(CK_fit)[1:4] <- rep(c("N_est_F", "N_est_M"), each = 2)

#Combine above to make dataframe with truth and estimates side-by-side
#store years from youngest sibling in comparisons to end of study
yrs <- c(min(mom_positives$Young_sib_birth, dad_positives$Young_sib_birth):t_end)

estimates <- data.frame(cbind(t(round(exp(CK_fit[1:4]),0)), SE, rep(c("F", "M"), each = 2)), yr = yrs)
colnames(estimates) <- c("CKMR_estimate", "SE", "sex", "yr")

moms_detected <- mom_positives %>% 
  group_by(Young_sib_birth) %>% 
  summarize(sum(freq)) %>% 
  pull(2)

dads_detected <- dad_positives %>% 
  group_by(Young_sib_birth) %>% 
  summarize(sum(freq)) %>% 
  pull(2)

(estimates <- cbind(estimates, truth = c(Mom_truth[yrs], Dad_truth[yrs]), total_samples=rep(n_samples, 4), parents_detected = c(moms_detected, dads_detected)))

####Fit model - nlminb####
CK_fit <- nlminb(start=Pars, objective=lemon_neg_log_lik, Negatives_Mother=mom_negatives, Negatives_Father=dad_negatives, Pairs_Mother=mom_positives, Pairs_Father=dad_positives, P_Mother=P_Mother, P_Father=P_Father, n_yrs=n_yrs, t_start=t_start, t_end=t_end)
D=diag(length(Pars))*c(exp(CK_fit$par[1]),exp(CK_fit$par[2]), exp(CK_fit$par[3]), exp(CK_fit$par[4])) #derivatives of transformations

#Estimating two years - need to edit to match nlminb function requirements
#VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
#VC = (t(D)%*%VC_trans%*%D) #delta method
#SE=round(sqrt(diag(VC)),0)

#Rename columns - change according to years estimating
#colnames(CK_fit)[1:4] <- rep(c("N_est_F", "N_est_M"), each = 2)


#Combine above to make dataframe with truth and estimates side-by-side
#store years from youngest sibling in comparisons to end of study
yrs <- c(min(mom_positives$Young_sib_birth, dad_positives$Young_sib_birth):t_end)

#Add SE to estimates once I can compute it
estimates <- data.frame(cbind(round(exp(CK_fit$par[1:4]),0), rep(c("F", "M"), each = 2)), yr = yrs)
colnames(estimates) <- c("CKMR_estimate", "sex", "yr")

moms_detected <- mom_positives %>% 
  group_by(Young_sib_birth) %>% 
  summarize(sum(freq)) %>% 
  pull(2)

dads_detected <- dad_positives %>% 
  group_by(Young_sib_birth) %>% 
  summarize(sum(freq)) %>% 
  pull(2)

(estimates <- cbind(estimates, truth = c(Mom_truth[yrs], Dad_truth[yrs]), total_samples=rep(n_samples, 4), parents_detected = c(moms_detected, dads_detected)))

all_est <- rbind(all_est, estimates)
row.names(all_est) <- NULL

save(CK_fit, file=paste0("models/Lemon_CKModel_HS_time_series_fishSim_6.29.20_", iter))
print(paste0("finished iteration", iter, " at: ", Sys.time()))
}

write.table(all_est, file = paste0("fishSim_HS_time_series_null_", n_samples, "samps_06.29.2020.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
}

sim_end_time <- Sys.time()
sim_end_time-sim_start_time
#Need to figure out which units for below
#print(paste0("Run time of simulation: ", round(sim_end_time - sim_start_time, 2), ""))
print(paste0("Simulation finished at ", Sys.time()))

###Troubleshooting
tail(indiv)
Moms <- indiv %>% 
  select(Mum, Dad) %>% 
  rename(Me = Mum) %>%
  mutate_all(as.character)

Grouped_Moms <- indiv %>% 
  inner_join(Moms, by = "Me") %>% 
  select(AgeLast, BirthY, DeathY)

Dads <- indiv %>% 
  select(Mum, Dad) %>% 
  rename(Me = Dad) %>%
  mutate_all(as.character)

Grouped_Dads <- indiv %>% 
  inner_join(Dads, by = "Me") %>% 
  select(AgeLast, BirthY, DeathY)
