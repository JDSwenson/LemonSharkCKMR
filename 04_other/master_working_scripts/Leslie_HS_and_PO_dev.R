#To do 01/26/2021:
#Start at line 158
#How can I create a pairwise comparison matrix while maintaining the other columns?

#Send this to Rusers group?
##Post small dataset with like 1:4

#Try to replicate Paul Conn's for loop with an abbreviated dataset using combn.
##If I run combn on the column of birth dates, I may be able to run combn in such a way that one of the two resulting columns from sex and capture year gives the appropriate value.
##Separate so that the appropriate columns of the two dataframes are together
##Rejoin the dataframes (or do this and the last step at once)
###Sort and only maintain the sex column for the older individual
####May need to make key and separate out again ...

#-------------Troubleshooting begin--------------------
ab <- data.frame(c(1, 3, 4, 2, 5), c("a", "b", "c", "d", "e"))
colnames(ab) <- c("col_a", "col_b")
ab
data.frame(t(combn(ab$col_a, m=2)))
data.frame(t(combn(ab$col_b, m=2)))



test <- data.frame(matrix(0,nrow=nrow(ab)*(nrow(ab)-1)/2,ncol=5)) # massive array for pairwise comparisons

Prob_vec = rep(0,nrow(test))
counter=1
#brute force - could probably be improved w/ expand.grid or something. Fills first four columns of Data
for(iind1 in 1:(nrow(ab)-1)){ #iind1 starts at 1 and counts to n_samples-1
  for(iind2 in (iind1+1):nrow(ab)){#iind2 starts at 2 and counts to n_samples
    if(ab[iind1,2]>ab[iind2,2]){ #compare birth years of each combination of samples. If birth year is greater for iind1 i.e. if iind2 was born before iind1.
      #Fill Data array:
      #Data[,1] = adult birth,
      #Data[,2] = adult capture year
      #Data[,3] = young birth
      #Data[,4] = adult sex
      test[counter,1]=ab[iind2,2] #put birth year of iind2 in adult birth year
      test[counter,2]=ab[iind2,1] #put year of capture of iind2 in adult capture year
      test[counter,3]=ab[iind1,2] #put birth year of iind1 in offspring birth year
      test[counter,4]=ab[iind2,3] #put sex of iind2 in adult sex
    }
    else{ #birth year greater for second value - reverse testove
      test[counter,1]=ab[iind1,2]
      test[counter,2]=ab[iind1,1]
      test[counter,3]=ab[iind2,2]
      test[counter,4]=ab[iind1,3]
    }
  }
}
        if(test[counter,4]==0)Prob_vec[counter]=P$P_Mother[test[counter,1],test[counter,2],test[counter,3]] #If the older sample is a female, then add to Prob_vec the probtestility of kinship calculated for the same values for adult birth, adult death, and offspring birth in the prior probtestility array for the mother.
    if(test[counter,4]==1)Prob_vec[counter]=P$P_Father[test[counter,1],test[counter,2],test[counter,3]] #Same with the father
    counter=counter+1
  }


#Troubleshooting end



#----------------------------Script begin---------------------
#rm(list=ls())
#Sex-specific estimates of N;
#Population growth fixed to observed population growth of the whole population.
#set so R doesn't use scientific notation
options("scipen"=100, "digits"=4)

library(foreach)
library(parallel)
library(doParallel)
library(optimx)
#Load individual packages because tidyverse won't load on cluster
library(plyr)
library(dplyr)
library(tidyr)
library(popbio)

set.seed(47)
#----------Simulation parameters-----------
t_start = 120 #First year we take samples
t_end = 125 #Last year we take samples
#min_est_cohort <- 42 #First year we are estimating abundance
#max_est_cohort <- t_end #Last year we are estimating abundance
#est_yrs <- min_est_cohort:max_est_cohort #All years for which we are estimating abundance
n_yrs <- t_end #Number of years of study

samp_yrs = t_start:t_end #Number of years being sampled
#n_samples_per_yr = 75
#n_samples_total <- n_samples_per_yr * n_samp_yr + n_samples_per_yr #because sampling starts in year t_start, we need to add one more year

#----------------Leslie Matrix and makeFounder parameters---------------
s_YOY = 0.6    # YOY survival
s_adult = 0.85  # Adult (mature) survival 

#Change length.out to ramp up from s_YOY to s_adult, and then edit times so the total number of survival rates is 11
s_juv = c(round(seq(s_YOY, s_adult, length.out = 3), 2), rep(s_adult, times = 8))    # Juvenile survival - ramp up from YOY survival to adult survival; with 13 columns corresponding to years 0 (YOY) - 12 (mature adults), we set 11 different survival rates, starting with the YOY survival and moving to adult. The for loop in the function below starts at index 2, so doesn't repeat YOY survival
a_m = 12         # Age at sexual maturity
b_max = 1.5     #Females recruited into the population. Assume equal sex ratio. Calculated from: average fecundity (6) divided by two (for skipped-breeding) divided by two (for females vs males)
max_age = maxAge  <- 30
age_x <- max_age


init_matrix = function(s_j, s_a, a_mat, age_xx, b_maxx){
  A_matrix = matrix(data = 0, nrow = age_xx, ncol = age_xx) # dimension matrix and fill with zeros
  A_matrix[1, (a_mat:age_xx)] = b_maxx # assign birth rates to mature ages in first row of matrix
  A_matrix[2,1] <- s_YOY
  for(ii in 2:(a_mat-1)) A_matrix[ii+1, ii] = s_j[ii] # assign juvenile survival rates
  #A_matrix[age_xx,a_mat] = s_a # adult survival assumed for maturing animals transitioning into plus-group
  for(jj in a_mat:(age_xx-1)) A_matrix[jj+1, jj] = s_adult # assign juvenile
  A_matrix[age_xx, age_xx] = s_a # adult survival assumed for plus-group
#  A_matrix[age_xx, age_xx] = 0 # No plus group
  return(A_matrix)
}

A = init_matrix(s_j = s_juv, s_a = s_adult, a_mat = a_m, 
                age_xx = age_x, b_maxx = b_max) # Call function to initialize matrix A
#View(A)


#Calculate dominant eigenvalue (i.e. population growth rate) from transition matrix
lam <- lambda(A)


#Calculate stable age structure of the population
stable_age = stable.stage(A)
stable_age


#Project population forward in time
Abundance_year0 <- c(20000, rep(0, times = (max_age-1)))

Year1 <- A %*% Abundance_year0

nYears <- n_yrs        # set the number of years to project
TMat <- A     # define the projection matrix
InitAbund <- Abundance_year0    # define the initial abundance

## Perpetuate year 0 individuals through Leslie Matrix
allYears <- matrix(0,nrow=nrow(TMat),ncol=nYears+1)     # build a storage array for all abundances!
allYears[,1] <- InitAbund  # set the year 0 abundance                                    
for(t in 2:(nYears+1)){   # loop through all years
  allYears[,t] <-  TMat %*% allYears[,t-1]
}

#in allYears, columns are years, rows are ages
allYears <- data.frame(allYears)
colnames(allYears) <- paste("Yr_", c(1:(nYears+1)))
allYears <- allYears %>% mutate_all(round, digits=0)

#Sum mature ages for each year
adult_truth <- allYears %>% slice(12:30) %>% 
  summarise_all(sum)


########### Old way of calculating "stable age distribution"#############
#Surv_age <- rep(0.85,maxAge) #Set survival for each age to 85%
#Surv_age[1] <- 0.6 #Set survival for first year to 60%
#Surv_age[2] <- 0.72 #Set survival for second year to 70%
#Prop_age <- rep(1:maxAge)
#for(iage in 2:maxAge)Prop_age[iage]=Prop_age[iage-1]*Surv_age[iage-1]
#Prop_age <- Prop_age/sum(Prop_age)  #stable age distribution ... supposedly
#survCurv <- Prop_age #Sets probability of founder cohort belonging to each age class


mat_m = mat_f = mat_a <- rep(0,maxAge) #creates an empty vector that has the length of n_ages
mat_m[12:maxAge] <- 1  #Knife-edge maturity for males, beginning at age 12
mat_f[12:maxAge] <- 1  #Knife-edge maturity for females, beginning at age 12
mat_a[12:maxAge] <- 1 #Knife-edge maturity for all

iterations <- 10

for(samps in 1:4){
  #results <- data.frame(matrix(0, nrow = iterations*3, ncol = 7))
  results <- NULL
  n_samples <- c(300, 500, 750, 1000)[samps]
  
  set.seed(47)
  for(iter in 1:iterations) {
    #set.seed(iter) #set the seed to the same ten numbers so the results are repeatable.
    index <- (iter*2) + 1 #Index for filling in the matrix at the end of each loop over n_samples
    

Samples <- c() 

prob_m <- N_m/(N_f+N_m)
n_samples_per_yr <- n_samples/(t_end - t_start +1) #add one to account for first year

  
Samples <- data.frame(sample(c(1:max_age), n_samples, replace=TRUE, prob=stable_age)) #birth year - probability of sampling based on Prop_age
colnames(Samples)[1] <- "age"

Samples2 <- Samples %>% mutate(capture_yr = rep(c(t_start:t_end), each=n_samples_per_yr)) %>% 
  mutate(birth_yr = capture_yr - age,
         sex = rbinom(n=n_samples, size=1, prob=prob_m)) %>% 
  select(capture_yr, birth_yr, sex)
  

Samples_HS <- Samples2$birth_yr
min_est_cohort <- min(Samples_HS)

  #Starting parameters
  N_f = N_m <- adult_truth[,min_est_cohort]/2
  N_a <- adult_truth[,min_est_cohort]
  Pars1 <- c(log(N_f),log(N_m)) #For sex-specific model
  Pars2 <- c(log(N_a)) #For sex-aggregated model
      
#Half-sib pairwise comparisons
Data_HS <- t(combn(Samples_HS, m=2))
Data_HS <- plyr::count(t(apply(Data_HS, 1, sort)))
colnames(Data_HS)[1:2] <- c("Old_sib_birth", "Young_sib_birth")


#Remove same-cohort comparisons
Data_HS2 <- Data_HS[which(Data_HS$Old_sib_birth != Data_HS$Young_sib_birth),]

#Parent-offspring pairwise comparisons
head(Samples2)

############# Pick up here 01/26/2021 ############
####NOT SURE THIS IS VALID####
#I think it may not capture all pairwise comparisons
##In an effort to maintain the sex and capture year, I split the dataset but then not all of the pairwise comparisons were completed.

#Split into list of dataframes by sex and capture year so we can create pairwise comparison matrix of birth years
Samples3 <- split(Samples2, Samples2$capture_yr)
Samples4 <- lapply(Samples3, function(x) split(x, x$sex))

pwise <- NULL
pairwise_fun <- function(x){
  for(i in 1:length(x)){
    for(j in 1:length(x[[i]])){
      p <- data.frame(t(combn(x[[i]][[j]]$birth_yr, m=2)))
      p$capture_yr <- x[[i]][[j]]$capture_yr[1]
      p$sex <- x[[i]][[j]]$sex[1]
      pwise <- rbind(pwise, p)
    }
  }
  return(pwise = pwise)
}

pairwise_PO <- pairwise_fun(Samples4)
pairwise_list <- plyr::count(t(apply(pairwise_PO, 1, sort))) %>% 
  rename(sex = x.1, older_birth = x.2, younger_birth = x.3) %>% 
  select(older_birth, younger_birth, sex, freq) %>% 
  group_by(sex) %>% 
  group_split()

#Show that first column of birth years is always greater
#pairwise_PO %>% filter(older_birth > younger_birth)


Samples2 %>% group_by(sex) %>% 
  group_by(capture_yr) %>% 
  expand(.[["birth_yr"]], .[["birth_yr2"]])
  
  expand_grid(age1 = .[["birth_yr"]], age2 = .[["birth_yr"]]) %>% 
  ungroup() #%>% 
  nrow()
  
nrow(expand.grid(Samples2$birth_yr, Samples2$birth_yr))



#For PC
setwd(".")
source("./functions/get_P_lemon_HS.R")
source("./functions/lemon_neg_log_lik_HS_Sex_specific.R")
source("./functions/get_P_lemon_HS_TotalA.R")
source("./functions/lemon_neg_log_lik_HS_TotalA.R")

# #For cluster
# source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/get_P_lemon_HS.R")
# source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/lemon_neg_log_lik_HS_Sex_specific.R")
# source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/get_P_lemon_HS_TotalA.R")
# source("/home/js16a/R/working_directory/CKMR_simulations/scripts/functions/lemon_neg_log_lik_HS_TotalA.R")


#Add probabilities of parentage to columns 4 & 5
for(i in 1:nrow(Data2)){
  Data2$Mom_prob[i] <- P$P_Mother[Data2[i,1],Data2[i,2]]
  Data2$Dad_prob[i] <- P$P_Father[Data2[i,1],Data2[i,2]]
  Data2$Parent_prob[i] <- P_TotalA$P_Parent[Data2[i,1],Data2[i,2]]
}
#Data2

#Randomly assign matches
for(i in 1:nrow(Data2)){
  Data2$Mom_Matches[i] <- sum(rbinom(n=Data2$freq[i], size=1, p=Data2$Mom_prob[i]))
  Data2$Dad_Matches[i] <- sum(rbinom(n=Data2$freq[i], size=1, p=Data2$Dad_prob[i]))
  Data2$Total_Matches[i] <- sum(rbinom(n=Data2$freq[i], size=1, p=Data2$Parent_prob[i]))
}


######################### Format pairwise comparison for CKMR ###########################
##Ultimately, want four dataframes:
#1) Rows contain positive comparisons for mothers
#2) Rows contain positive comparisons for fathers
#3) Rows contain negative comparisons for mothers
#4) Rows contain negtive comparisons for fathers

###################Create dataframes for CKMR model####################
mom_positives <- Data2[which(Data2$Mom_Matches >0), c(1:2,7)]
dad_positives <- Data2[which(Data2$Dad_Matches >0), c(1:2,8)]
dad_negatives = mom_negatives = parent_negatives <- Data2
dad_negatives[,3] <- dad_negatives[,3] - dad_negatives[,8]
dad_negatives <- dad_negatives[,c(1:3)]
mom_negatives[,3] <- mom_negatives[,3] - mom_negatives[,7]
mom_negatives <- mom_negatives[,c(1:3)]
parent_positives <- Data2[which(Data2$Total_Matches >0), c(1:2,9)]
parent_negatives[,3] <- parent_negatives[,3] - parent_negatives[,9]
parent_negatives <- parent_negatives[,c(1:3)]

    
    #Fit model = nlminb version
    # CK_fit <- nlminb(start = Pars, objective = lemon_neg_log_lik, hessian = TRUE, Negatives_Mother=mom_negatives, Negatives_Father=dad_negatives, Pairs_Mother=mom_positives, Pairs_Father=dad_positives, P_Mother=P_Mother, P_Father=P_Father, t_start=t_start, t_end=t_end)
    # 
    # #exp(CK_fit$par)
    # #compute variance covariance matrix - nlminb
    # D=diag(length(Pars))*c(exp(CK_fit$p[1]),exp(CK_fit$p[2])) #derivatives of transformations
    # VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
    # VC = (t(D)%*%VC_trans%*%D) #delta method
    # SE=round(sqrt(diag(VC)),0)    
    
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
    estimates <- cbind(estimates, c(rep(as.numeric(adult_truth[min_est_cohort]/2), times=2), as.numeric(adult_truth[min_est_cohort])))
    colnames(estimates) <- c("CKMR_estimate", "SE", "Sex", "Truth")
    estimates
    
    metrics <- cbind(c(sum(mom_positives[,3]), sum(dad_positives[,3]), sum(parent_positives[,3])), c(rep(lam, times = 3)), c(rep(n_samples, times=3)))
    colnames(metrics) <- c("Parents_detected", "Pop_growth", "Samples")
    
    
       #-----------------Loop end-----------------------------
    results <- rbind(results, cbind(estimates, metrics))
    

    print(paste0("finished iteration", iter, " at: ", Sys.time()))
  }
    
  results <- results %>% 
    mutate(Relative_bias = round(((CKMR_estimate - Truth)/Truth)*100,1))
  
#  results$Relative_bias <- round(((results$N_est - results$Truth)/results$Truth)*100,1)
  
   #results %>% group_by(Sex) %>% 
   #  summarize(median = median(Relative_bias), n = n())
  
  write.table(results, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/Leslie.null_", n_samples, ".samples_01.19.2021_neutral.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
}
