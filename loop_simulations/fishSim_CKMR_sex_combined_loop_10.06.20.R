library(fishSim)
library(foreach)
library(parallel)
library(doParallel)
library(optimx)
library(tidyverse)
library(popbio)

#---------------Set up experiment parameters--------------------
#Constant Abundance - one estimate
t_start = 40 #First year we take samples
t_end = 45 #Last year we take samples
#min_est_cohort <- 42 #First year we are estimating abundance
#max_est_cohort <- t_end #Last year we are estimating abundance
#est_yrs <- min_est_cohort:max_est_cohort #All years for which we are estimating abundance
n_yrs <- t_end #Number of years of study

#Set up the simulation parameters
sim_yrs = c(1:39)
n_samp_yr = t_end-t_start #Number of years being sampled
n_samples_per_yr = 40
n_samples_total <- n_samples_per_yr * n_samp_yr + n_samples_per_yr #because sampling starts in year t_start, we need to add one more year

#makeFounder parameters
osr <- c(0.5, 0.5) #Operating sex ratio
stocks <- c() #stocks (empty for now)
maxAge <- 30 #max age of lemon sharks
Surv_age <- rep(0.85,maxAge) #Set survival for each age to 85%
Surv_age[1] <- 0.6 #Set survival for first year to 60%
Surv_age[2] <- 0.7 #Set survival for second year to 70%
Prop_age <- rep(1:maxAge)
pop <- 5000

for(iage in 2:maxAge)Prop_age[iage]=Prop_age[iage-1]*Surv_age[iage-1]
Prop_age <- Prop_age/sum(Prop_age)  #stable age distribution
survCurv <- Prop_age #Sets probability of founder cohort belonging to each age class

#Set parameters for altMate
batchSize = 3 #average size of brood -- assumes Poisson distribution for fecundity (which is default for altMate). Average fecundity (6) divided by 2 for skipped-breeding
firstBreed <- 12 #If using same knife-edge maturity for males and females
mat_A <- rep(0,maxAge) #creates an empty vector that has the length of n_ages
mat_A[12:maxAge]=1  #Knife-edge maturity for males, beginning at age 12

matCurve <- c(mat_A, rep(max(mat_A), 1000)) #Set female age at maturity
maxClutch <- 18 #Maximum number of offspring one individual can possibly produce
mat_type <- "age"

#Set parameters for mort()
ageMort <- (1 - Surv_age)
ageMort[maxAge+1] <- ageMort[maxAge]
death_type <- "age"

total_abundance <- pop
N_A <- total_abundance
Pars=log(N_A)

#-------------Kinship probabilities - Half-sib-------------------
#Not assessing survival - keeping constant at 0.85
Surv <- 0.85

max_age = 30 #max age of lemon sharks

adult_age <- c(12:30) #Adult age of females -males assumed to be the same
adult_mat <- c(rep(0,11), rep(1,19)) #knife-edge maturity

P_Parent = array(NA,dim=c(n_yrs,n_yrs)) #creates two empty arrays, one for mother and one for father.  Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)

get_P_lemon <- function(Pars,P_Parent,n_yrs,t_start,t_end){
  N_A=exp(Pars[1]) #number of mature females
  #N_F=exp(Pars[1:n_yrs]) #number of mature females
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in max(os_birth+1):n_yrs){
      P_Parent[os_birth, ys_birth] <- (4/N_A)*(Surv^(ys_birth - os_birth)) #N_F is the number of females alive when the younger sibling was born
    }
  }
  return(list(P_Parent=P_Parent)) #return makes sure this is moved out of the loop into the environment
}
#P=get_P_lemon(Pars=Pars,P_Parent,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
#P$P_Parent

#------------------Likelihood function--------------------------
lemon_neg_log_lik <- function(Pars, Negatives_Parent, Pairs_Parent, P_Parent, n_yrs, t_start, t_end){
  
  P=get_P_lemon(Pars=Pars,P_Parent,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
  
  loglik=0
  
  #likelihood contributions for all negative comparisons
  for(irow in 1:nrow(Negatives_Parent)){
    loglik=loglik+Negatives_Parent[irow,3]*log(1-P$P_Parent[Negatives_Parent[irow,1],Negatives_Parent[irow,2]])
  } 
  #likelihood contributions for positive comparisons
  for(irow in 1:nrow(Pairs_Parent)){
    loglik=loglik+Pairs_Parent[irow,3]*log(P$P_Parent[Pairs_Parent[irow,1],Pairs_Parent[irow,2]])
  }
  -loglik
}


iterations <- 100

for(samps in 1:4){
  results <- data.frame(matrix(0, nrow = iterations*2, ncol = 5))
  n_samples <- c(30, 60, 90, 120)[samps]
  
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

truth <- c()

for (y in 1:length(sim_yrs)) {
  indiv <- altMate(indiv, batchSize = batchSize, fecundityDist = "poisson", year = y, type = mat_type, maxClutch = maxClutch, singlePaternity = FALSE, maturityCurve = matCurve, firstBreed = firstBreed)
  indiv <- mort(indiv, type = death_type, year=y, maxAge = maxAge, ageMort = ageMort) 
  indiv <- birthdays(indiv)
  
  #Store true values for each year
  truth[y] <- indiv %>% 
    filter(AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
    summarize(n())
  }

nrow(indiv[is.na(indiv[,6]),]) ## the currently-alive population size.

#Repeats the above loop, assigning offspring, deaths, and birthdays, but this time also draws samples each year and records which animals were sampled in indiv
for (y in c(t_start:t_end)) {
  indiv <- altMate(indiv, batchSize = batchSize, fecundityDist = "poisson", year = y, type = mat_type, maxClutch = maxClutch, singlePaternity = FALSE, maturityCurve = matCurve, firstBreed = firstBreed)
  indiv <- mort(indiv, type = death_type, year=y, maxAge = maxAge, ageMort = ageMort)
  
  indiv <- birthdays(indiv) #Age each individual by one year
  indiv <- capture(indiv, n = n_samples, year = y, fatal = FALSE) #Capture individuals

  #Store true values for each year
  truth[y] <- indiv %>% 
    filter(AgeLast >= firstBreed & is.na(DeathY)==TRUE) %>% 
    summarize(n())
  
}

#Store truth values as vectors (dplyr makes them a list)
truth <- unlist(truth) 

pairs <- findRelativesPar(indiv = indiv, sampled = TRUE)
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
  rename("Old_sib_birth" = BirthY, "Old_sib_mom" = Mum, "Old_sib_dad" = Dad) %>%   select(c(Old_sib_birth, Young_sib_birth)) %>% 
  plyr::count() %>% 
  select(Old_sib_birth, Young_sib_birth, freq) %>% 
  filter(Young_sib_birth != Old_sib_birth)

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
  select(c(Old_sib_birth, Young_sib_birth)) %>% 
  plyr::count() %>% 
  filter(Young_sib_birth != Old_sib_birth)

P=get_P_lemon(Pars=Pars,P_Parent=P_Parent,n_yrs=n_yrs,t_start=t_start,t_end=t_end)


#Fit model
CK_fit <- optimx(par=Pars,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Parent=non_HSPs_2_tbl, Pairs_Parent=HSPs_2_tbl, P_Parent=P_Parent, n_yrs=n_yrs, t_start=t_start, t_end=t_end)

#summary(CK_fit)
#exp(CK_fit[1])

#compute variance covariance matrix
D=diag(length(Pars))*c(exp(CK_fit$p1[1])) #derivatives of transformations
VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
VC = (t(D)%*%VC_trans%*%D) #delta method
SE=round(sqrt(diag(VC)),0)

#Rename columns
#colnames(CK_fit)[1:2] <- c("N_est_F", "N_est_M")

#Combine above to make dataframe with truth and estimates side-by-side
#store years from youngest sibling in comparisons to end of study
yrs <- c(min(HSPs_2_tbl$Young_sib_birth):t_end)

estimates <- data.frame(cbind(t(round(exp(CK_fit[1]),0)), SE))
colnames(estimates) <- c("CKMR_estimate", "SE")
estimates <- cbind(estimates, truth = c(round(mean(truth[yrs]),0)))

#-----------------Loop end-----------------------------
results[iter,] <- c(round(exp(CK_fit$p1[1]),0), round(SE[1],0), round(mean(truth[yrs]),0), sum(HSPs_2_tbl[,3]), n_samples)
colnames(results) <- c("N_est", "SE", "Mean_truth", "Parents_detected", "Samples")

  print(paste0("finished iteration", iter, " at: ", Sys.time()))
}
 results <- results %>% 
    mutate(Relative_bias = round(((N_est - Mean_truth)/Mean_truth)*100,1))
  
  write.table(results, file = paste0("NEW_fishSim_sexes-combined_single_est_", n_samples, "_samples_10.06.2020.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
}
