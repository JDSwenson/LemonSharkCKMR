rm(list=ls())

library(fishSim)
library(foreach)
library(parallel)
library(doParallel)
library(optimx)
library(plyr)
library(dplyr)
library(tidyr)

set.seed(4747)
#---------------Set up experiment parameters--------------------
#Constant Abundance - one estimate
t_start = 41 #First year we take samples
t_end = 45 #Last year we take samples
#min_est_cohort <- 42 #First year we are estimating abundance
#max_est_cohort <- t_end #Last year we are estimating abundance
#est_yrs <- min_est_cohort:max_est_cohort #All years for which we are estimating abundance
n_yrs <- t_end #Number of years of study

#Set up the simulation parameters
sim_yrs = c(1:(t_start-1))
n_samp_yr = t_end-t_start #Number of years being sampled
n_samples_per_yr = 40
n_samples_total <- n_samples_per_yr * n_samp_yr + n_samples_per_yr #because sampling starts in year t_start, we need to add one more year

#makeFounder parameters
osr <- c(0.5, 0.5) #Operating sex ratio
stocks <- c() #stocks (empty for now)
maxAge <- 30 #max age of lemon sharks
Surv_age <- rep(0.85,maxAge) #Set survival for each age to 85%
Surv_age[1] <- 0.55 #Set survival for first year to 60%
Surv_age[2] <- 0.8 #Set survival for second year to 70%
Prop_age <- rep(1:maxAge)
pop <- 5000

for(iage in 2:maxAge)Prop_age[iage]=Prop_age[iage-1]*Surv_age[iage-1]
Prop_age <- Prop_age/sum(Prop_age)  #stable age distribution
survCurv <- Prop_age #Sets probability of founder cohort belonging to each age class

#Set parameters for altMate
batchSize = 3.5 #average size of brood -- assumes Poisson distribution for fecundity (which is default for altMate). Average fecundity (6) divided by 2 for skipped-breeding
firstBreed <- 12 #If using same knife-edge maturity for males and females

mat_m <- rep(0,maxAge) #creates an empty vector that has the length of n_ages
mat_m[12:maxAge]=1  #Knife-edge maturity for males, beginning at age 12

mat_f <- rep(0,maxAge) #creates an empty vector that has the length of n_ages
mat_f[12:maxAge]=1  #Knife-edge maturity for females, beginning at age 12

maleCurve <- c(mat_m, rep(max(mat_m), 1000)) #Set male age at maturity
femaleCurve <- c(mat_f, rep(max(mat_m), 1000)) #Set female age at maturity
maxClutch <- 18 #Maximum number of offspring one individual can possibly produce
mat_type <- "ageSex"

#Set parameters for mort()
ageMort <- (1 - Surv_age)
ageMort[maxAge+1] <- ageMort[maxAge]
death_type <- "age"

total_abundance <- pop
N_f <- c(total_abundance/2)
N_m <- c(total_abundance-N_f)
Pars1=c(log(N_f),log(N_m))
Pars2=c(log(total_abundance))

#-------------Kinship probabilities - Half-sib-------------------
#Not assessing survival - keeping constant at 0.85
surv <- 0.85
m_adult_age <- c(12:maxAge) #Set ages at which males and females are mature.
f_adult_age <- c(12:maxAge)

iterations <- 50

source("~/R/R_working_dir/LemonSharkCKMR_GitHub/00_functions/lemon_funcs_05.10.2021.R")

for(samps in 1:3){
  results <- NULL
  n_samples <- c(50, 75, 100)[samps]
  sampled_ages <- data.frame(matrix(0, nrow = n_samples, ncol = iterations))
  for(iter in 1:iterations) {
    index <- (iter*2)-1

#-------------Loop start-------------------------    
#makeFounders creates a matrix of specified size (pop) where each row is an individual in the founder popualtion.
indiv <- makeFounders(pop = pop, osr = osr, stocks = c(1), maxAge = maxAge, survCurv=survCurv)

#Check growth rate to find first year survival rate that produces 0 population growth
#check_growthrate(mateType = "ageSex", femaleCurve = femaleCurve, mortType = "age", ageMort = ageMort, batchSize = batchSize, firstBreed = 12, maxClutch = maxClutch)

#PoNG(mateType = "ageSex", femaleCurve = femaleCurve, mortType = "age", ageMort = ageMort, batchSize = batchSize, firstBreed = 0, maxClutch = maxClutch)

#Simulate 40 years of mating
#altMate derives the number of offspring produced by drawing from a sampling distribution for each female based on the parameters specified in the function. Fathers are randomly drawn from all mature males in the mother's stock.
##This loop can take a minute depending on population size

Dad_truth <- c()
Mom_truth <- c()
All_truth <- c()

for (y in 1:length(sim_yrs)) {
  #Store true values for each year
  Dad_truth[y] <- indiv %>% 
    filter(Sex == "M" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
    dplyr::summarize(n())
  
  Mom_truth[y] <- indiv %>% 
    filter(Sex == "F" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
    dplyr::summarize(n())
  indiv <- altMate(indiv, batchSize = batchSize, fecundityDist = "poisson", year = y, type = mat_type, maxClutch = maxClutch, singlePaternity = FALSE, maleCurve = maleCurve, femaleCurve = femaleCurve, firstBreed = firstBreed)
  indiv <- mort(indiv, type = death_type, year=y, maxAge = maxAge, ageMort = ageMort) 
  indiv <- birthdays(indiv)
}

nrow(indiv[is.na(indiv[,6]),]) ## the currently-alive population size.

#Repeats the above loop, assigning offspring, deaths, and birthdays, but this time also draws samples each year and records which animals were sampled in indiv
for (y in c(t_start:t_end)) {
  #Store true values for each year
  Dad_truth[y] <- indiv %>% 
    filter(Sex == "M" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
    dplyr::summarize(n())
  
  Mom_truth[y] <- indiv %>% 
    filter(Sex == "F" & AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
    dplyr::summarize(n())
  
  indiv <- altMate(indiv, batchSize = batchSize, fecundityDist = "poisson", year = y, type = mat_type, maxClutch = maxClutch, singlePaternity = FALSE, maleCurve = maleCurve, femaleCurve = femaleCurve, firstBreed = firstBreed)
  indiv <- mort(indiv, type = death_type, year=y, maxAge = maxAge, ageMort = ageMort)
  
  indiv <- birthdays(indiv) #Age each individual by one year
  indiv <- capture(indiv, n=n_samples, year = y, fatal = FALSE) #Capture individuals
}

#Store truth values as vectors (dplyr makes them a list)
Mom_truth <- unlist(Mom_truth) 
Dad_truth <- unlist(Dad_truth)
All_truth <- Mom_truth + Dad_truth

#Mom_truth
#Dad_truth

pairs <- findRelativesPar(indiv = indiv, sampled = TRUE)
POPs <- pairs[pairs$OneTwo == 1,1:2] ##Parent-Offspring pairs
HSPs <- pairs[pairs$TwoTwo == 1,1:2] ##Half-Sibling pairs - verified from fishSim vignette
non_POPs <- pairs[pairs$OneTwo == 0,1:2] ##pairs that are not POPs
non_HSPs <- pairs[pairs$TwoTwo == 0,1:2] ##pairs that are not half-sibs

# think of this dataframe as storing only information about the younger fish. Stores ID and birth year only (at present) for every individual in indiv, but renames the ID column to younger so it can be joined with HSPs_tbl below)
youngerbirthyears <- indiv %>%
  dplyr::select(Me, BirthY, Mum, Dad) %>% 
  dplyr::rename("younger" = Me, "Young_sib_birth" = BirthY, "Young_sib_mom" = Mum, "Young_sib_dad" = Dad)

#Create dataframe with IDs of sampled half-sibs and columns named appropriately for joins
HSPs_tbl <- HSPs %>% 
  dplyr::rename(
    Me = Var1,
    younger = Var2) %>%
  mutate_all(as.character)

#Extract the values of indiv for the older individual in each comparison.
#Join all the information from indiv with the IDs of the sampled individuals in HSPs_2_tbl. 
HSPs_2_tbl <- inner_join(indiv, HSPs_tbl, by = "Me") %>%
  left_join(youngerbirthyears, by = "younger")  %>%
  dplyr::rename("Old_sib_birth" = BirthY, "Old_sib_mom" = Mum, "Old_sib_dad" = Dad) %>%   select(c(Old_sib_birth, Young_sib_birth, Old_sib_mom, Young_sib_mom, Old_sib_dad, Young_sib_dad))

#Split HSPs dataframe into MHS and PHS pairs
#Filter out intra-cohort comparisons
mom_positives <- HSPs_2_tbl[HSPs_2_tbl$Old_sib_mom == HSPs_2_tbl$Young_sib_mom,] %>% 
    dplyr::select(Old_sib_birth, Young_sib_birth, Old_sib_mom) %>% 
    dplyr::rename(Mother = Old_sib_mom) %>% #Rename to match column name in model function
    dplyr::filter(Young_sib_birth != Old_sib_birth) %>% 
    plyr::count() %>% 
    dplyr::select(Old_sib_birth, Young_sib_birth, freq, Mother)


dad_positives <- HSPs_2_tbl[HSPs_2_tbl$Old_sib_dad == HSPs_2_tbl$Young_sib_dad,] %>% 
    dplyr::select(Old_sib_birth, Young_sib_birth, Old_sib_dad) %>% 
    dplyr::rename(Father = Old_sib_dad) %>% #Rename to match column name in model function
    dplyr::filter(Young_sib_birth != Old_sib_birth) %>% 
    plyr::count() %>% 
    dplyr::select(Old_sib_birth, Young_sib_birth, freq, Father)

parent_positives <- HSPs_2_tbl[HSPs_2_tbl$Old_sib_dad == HSPs_2_tbl$Young_sib_dad | HSPs_2_tbl$Old_sib_mom == HSPs_2_tbl$Young_sib_mom,] %>% 
  select(Old_sib_birth, Young_sib_birth) %>% 
  filter(Young_sib_birth != Old_sib_birth) %>% 
  plyr::count() %>% 
  select(Old_sib_birth, Young_sib_birth, freq)


nrow(HSPs_2_tbl[HSPs_2_tbl$Old_sib_dad == HSPs_2_tbl$Young_sib_dad & HSPs_2_tbl$Old_sib_mom == HSPs_2_tbl$Young_sib_mom,]) #Should be 0 or there are full sibs

# Rename the columns so that the join functions know what to work with. The column `Me` is still the actual id column. However, we also need the column `younger` so that we can look up the 
#birth years stored in `youngerbirthyears`
non_HSPs_tbl <- non_HSPs %>% 
  dplyr::rename(
    Me = Var1,
    younger = Var2) %>%
  mutate_all(as.character) # convert columns to characters since not all levels of indiv$Me are present in non_POPs_tbl
#head(non_HSPs_tbl)

##Extract the values of indiv for the older individual in each comparison.
#Join all the information from indiv with the IDs of the sampled individuals in non_HSPs_2_tbl. 
non_HSPs_2_tbl <- inner_join(indiv, non_HSPs_tbl, by = "Me") %>%
  left_join(youngerbirthyears, by = "younger")  %>%
  dplyr::rename(Old_sib_birth = BirthY) %>% 
  dplyr::select(c(Old_sib_birth, Young_sib_birth))

#Create table of all negative comparisons, count occurences, and filter intracohort comparisons
mom_negatives = dad_negatives = parent_negatives <- non_HSPs_2_tbl %>% 
    plyr::count() %>% 
    dplyr::filter(Young_sib_birth != Old_sib_birth)

##Define variables and functions
#source("~/R/R_working_dir/CKMR/LemonSharkCKMR/constant_abundance/models/get_P_lemon_HS_fishSim.R")
#source("~/R/R_working_dir/CKMR/LemonSharkCKMR/likelihood_functions/lemon_neg_log_like_HS_Sex_specific.R")

#P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)

#P$P_Mother

#Check P function
#P$P_Mother[19,40,45] #Dimensions are parent birth year, parent capture year, offspring birth year (all of which are specified by n_yrs)

#Fit model
CK_fit1 <- optimx(par=Pars1,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Mother=mom_negatives, Negatives_Father=dad_negatives, Pairs_Mother=mom_positives, Pairs_Father=dad_positives, P_Mother=P_Mother, P_Father=P_Father, n_yrs=n_yrs, t_start=t_start, t_end=t_end)

CK_fit2 <- optimx(par=Pars2,fn=lemon_neg_log_lik_TotalA,hessian=TRUE, method="BFGS", Negatives_Parent=parent_negatives, Pairs_Parent=parent_positives, P_Parent=P_Parent, t_start=t_start, t_end=t_end)

summary(CK_fit1)
exp(CK_fit1[1:2])
exp(CK_fit2[1])

#compute variance covariance matrix - sex-specific model
D1=diag(length(Pars1))*c(exp(CK_fit1$p1[1]),exp(CK_fit1$p2[1])) #derivatives of transformations
VC_trans1 = solve(attr(CK_fit1, "details")["BFGS" ,"nhatend"][[1]])
VC1 = (t(D1)%*%VC_trans1%*%D1) #delta method
SE1=round(sqrt(diag(VC1)),0)


#compute variance covariance matrix - sex-aggregated model
D2=diag(length(Pars2))*exp(CK_fit2$p1[1]) #derivatives of transformations
VC_trans2 = solve(attr(CK_fit2, "details")["BFGS" ,"nhatend"][[1]])
VC2 = (t(D2)%*%VC_trans2%*%D2) #delta method
SE2 = round(sqrt(diag(VC2)),0)

#Rename columns
#colnames(CK_fit)[1:2] <- c("N_est_F", "N_est_M")

#Combine above to make dataframe with truth and estimates side-by-side
#store years from youngest sibling in comparisons to end of study
yrs <- c(min(mom_positives$Young_sib_birth, dad_positives$Young_sib_birth):t_end)

sampled_indvs <- indiv %>% filter(is.na(SampY) == FALSE)

#Set minimum cohort
min_cohort <- n_yrs - maxAge

Mom_truth_mean <- mean(Mom_truth[min_cohort:n_yrs])
Dad_truth_mean <- mean(Dad_truth[min_cohort:n_yrs])
All_truth_mean <- mean(All_truth[min_cohort:n_yrs])

estimates <- data.frame(cbind(round(exp(c(CK_fit1$p1[1], CK_fit1$p2[1], CK_fit2$p1[1])),0)), c(SE1, SE2), c("F", "M", "All"))
estimates <- cbind(estimates, c(Mom_truth_mean, Dad_truth_mean, All_truth_mean))
colnames(estimates) <- c("N_est", "SE", "Sex", "Mean_truth")
#estimates



#Store number of parents detected, observed lambda, and number of samples
metrics <- cbind(c(sum(mom_positives[,3]), sum(dad_positives[,3]), sum(parent_positives[,3])), c(rep(n_samples_total, times=3)))
colnames(metrics) <- c("Parents_detected", "Total_samples")

#-----------------Loop end-----------------------------    
results <- rbind(results, cbind(estimates, metrics))

  print(paste0("finished iteration", iter, " at: ", Sys.time()))
}
 
results <- results %>% 
   mutate(N_est = as.numeric(N_est),
          SE = as.numeric(SE),
          Mean_truth = as.numeric(Mean_truth),
          Parents_detected = as.numeric(Parents_detected),
          Samples = as.numeric(Total_samples)) %>%
    mutate(Relative_bias = round(((N_est - Mean_truth)/Mean_truth)*100,1))
  
  write.table(results, file = paste0("fishSim_AvgN_", n_samples, "_samples_05.10.2021.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
}

#Quick viz of results
library(ggpubr)

#null model results
IBS.null_CNR <- Sys.glob("~/R/R_working_dir/LemonSharkCKMR_GitHub/*.10.2021.csv")
IBS.null_results_CNR <- sapply(IBS.null_CNR, read_csv, simplify=FALSE) %>% 
  bind_rows()

ggplot(data=IBS.null_results_CNR, aes(x = factor(Samples))) + 
  geom_boxplot(aes(y = Relative_bias, fill=Sex)) +
  ylim(-100, 100) +
  geom_hline(yintercept=0, col="black", size=1.25) + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = -20, ymax = 20, alpha = 0.5, col = "red") + 
  labs(x = "Sample size", y = "Relative bias", title = "Relative bias by sample size") +
  scale_fill_brewer(palette="Set2") +
  font("title", size = 10, face = "bold")
