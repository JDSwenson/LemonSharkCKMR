#THIS IS THE MOST RECENT SCRIPT 3/26/2020

######## This code gives a time-series of CKMR abundance estimates over four years of sampling ########

# Formats the data for a half-sibling model
# The code simulates kinship assignment with the same probabilities as the model i.e. it is code to test the null model

#NEXT TO DO
# 1. Standard error estimates aren't working (line 139). Fix this and then run on a loop.

rm(list=ls())
library(tidyverse)
library(optimx)

#################Set up adult life history info ################
# Simulate and assign ages to adults based on life history
# Max age: two options: 
# 37 from https://www.saveourseasmagazine.com/lemon-sharks-old/
# 25 from White et. al. 2014
# Juvenile survival values from Dibattista 2007: 0.55 (mort = 0.45)
# Adult survival value from White et. al. 2014: 0.85 (mort=0.15)
# Litter size from Feldheim et. al. 2002 (avg=6, max=18)
# Age-at-maturity from Brown and Gruber 1988 (11.6 for M, 12.7 for F)
# About 77 juvenile lemon sharks inhabit the nursery at a given time (White et al 2014). So, 10sqrt(N) = 88 samples total

#Set population parameters for sampling
sample_ages <- 4 #specify how many cohorts are represented in the samples, since abundance estimates are given for the cohorts
Surv_age <- rep(0.85, sample_ages) #Set survival for adults
Surv_age[1:2] <- c(rep(0.55,2))
Prop_age <- rep(1, sample_ages) #Initialize variable Prop_age
for(iage in 2:length(Prop_age)) Prop_age[iage]=Prop_age[iage-1]*Surv_age[iage-1] #Set proportion of animals from each cohort surviving to the next age class based on Surv_age
Prop_age=Prop_age/sum(Prop_age)  #stable age distribution given assumed survival i.e. gives proportion of ages in a population based on assumed survival

n_yrs=sample_ages+(sample_ages-1) #Number of years from earliest individual to end of study

#This will change actual years into "study years", so we set up variables for the conversation later
t_start=4 #First year of sampling
t_end=7 #Last year of sampling
study_yrs = c(t_start:t_end)
#length(samp_study_yrs)

####Set up simulation parameters####

#Set Estimated_truth and number of males and females
#n_yrs
total_abundance <- c(3000, 5000, 7000, 9000, 2000) #Can't estimate for yr=1 because abundance is calculated at time of younger sibling's birth, and there is no scenario where the younger sibling was born in yr1. Also unlikely to get estimate for yr=2 because there are relatively few pairwise comparisons where the younger sibling was born in yr2. So length of vector should be n_yrs-2
N_f <- c(total_abundance/2.5)
N_m <- c(total_abundance-N_f)
Pars=c(log(N_f),log(N_m))
est_yrs <- n_yrs-2 #Number of years we're estimating

#Set number of samples according to Estimated_truth
estimated_truth <- max(total_abundance)
n_samples <- round(10*sqrt(estimated_truth), 0)
n_samples_per_yr <- round(n_samples/length(study_yrs),0) #Total number of samples taken over length of study

P_Mother = P_Father = array(NA,dim=c(n_yrs,n_yrs)) #creates two empty arrays, one for mother and one for father.  Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)
source("./time_series/models/get_P_lemon_HS_time_series_liz.R")
source("./likelihood_functions/lemon_neg_log_like_HS.R")

P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
#P$P_Mother
#P$P_Father

#############################Sample population##############################
Samples = matrix(0,n_samples,2) #matrix filled with 0s; each row is a sample(piece of tissue)
prob_m <- N_m/(N_f+N_m)

Samples[,1]=rep(c(t_start:t_end),times=c(rep(n_samples_per_yr, 3),n_samples_per_yr+1)) #what year did I capture you? repeat each year for the number of samples that were taken.
Samples[,2]=Samples[,1]-sample(c(0:(sample_ages-1)),n_samples,replace=TRUE,prob=Prop_age) #birth year -- simulate age based on Prop_age, and subtract from capture year to get birth year. We sample from 0:sample_ages-1 so we can include individuals born in the last year (age 0), but if we enter 0:sample_ages, then we will have one more cohort than we want.
#Samples[,3]=rbinom(n=n_samples, size=1, prob=prob_m) #sex; 0=female -- moved this to below

colnames(Samples) <- c("Capt", "Birth")
#head(Samples)
#tail(Samples)

Data <- t(combn(Samples[,2], m=2))
Data <- plyr::count(t(apply(Data, 1, sort)))
colnames(Data)[1:2] <- c("Old_sib_birth", "Young_sib_birth")

##Add sex to third column, if helpful
#np <- nrow(samp2)
#samp2 <- data.frame(samp2)
#samp2[,3] <- rbinom(n=np, size=1, prob=prob_m) #sex; 0=female
#colnames(samp2)[1:3] <- c("Old_sib", "Young_sib", "Sex")

#Remove same-cohort comparisons
Data2 <- Data[which(Data$Old_sib_birth != Data$Young_sib_birth),]

#Remove unused rows (where young sib birth is yr 2)
Data2 <- Data2[which(Data2$Young_sib_birth>=3),]

#Add probabilities of parentage to columns 4 & 5
for(i in 1:nrow(Data2)){
  Data2$Mom_prob[i] <- P$P_Mother[Data2[i,1],Data2[i,2]]
  Data2$Dad_prob[i] <- P$P_Father[Data2[i,1],Data2[i,2]]
}
Data2

#Randomly assign matches
for(i in 1:nrow(Data2)){
  Data2$Mom_Matches[i] <- sum(runif(Data2[i,3])<Data2[i,4])
  Data2$Dad_Matches[i] <- sum(runif(Data2[i,3])<Data2[i,5])
  }


######################### Format pairwise comparison for CKMR ###########################
##Ultimately, want four dataframes:
#1) Rows contain positive comparisons for mothers
#2) Rows contain positive comparisons for fathers
#3) Rows contain negative comparisons for mothers
#4) Rows contain negtive comparisons for fathers

###################Create dataframes for CKMR model####################
Data_mom_yes <- Data2[which(Data2$Mom_Matches >0), c(1:2,6)]
Data_dad_yes <- Data2[which(Data2$Dad_Matches >0), c(1:2,7)]
Data_dad_no = Data_mom_no <- Data2
Data_dad_no[,3] <- Data_dad_no[,3] - Data_dad_no[,7]
Data_dad_no <- Data_dad_no[,c(1:3)]
Data_mom_no[,3] <- Data_mom_no[,3] - Data_mom_no[,6]
Data_mom_no <- Data_mom_no[,c(1:3)]

  #########################Fit model!##############################
    CK_fit <- optimx(par=Pars,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Mother=Data_mom_no,Negatives_Father=Data_dad_no,Pairs_Mother=Data_mom_yes,Pairs_Father=Data_dad_yes,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)

    #summary(CK_fit)
  
    #Without SE
      #exp(CK_fit[1,1:5])
      #N_f
      #exp(CK_fit[1,6:10])
      #N_m
 
    #Store estimates in vector -- need to change based on number of parameters
      estimates <- c(exp(CK_fit$p1[1]),exp(CK_fit$p2[1]), exp(CK_fit$p3[1]), exp(CK_fit$p4[1]), exp(CK_fit$p5[1]), exp(CK_fit$p6[1]), exp(CK_fit$p7[1]), exp(CK_fit$p8[1]), exp(CK_fit$p9[1]), exp(CK_fit$p10[1]))
      
      #compute variance covariance matrix
    D=diag(length(Pars))*estimates #derivatives of transformations
    
      VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
    VC = (t(D)%*%VC_trans%*%D) #delta method
    SE=sqrt(diag(VC))
N_f
cat(paste0("Estimated/true female abundance: ",c(round(exp(CK_fit[1,1:5]),0)),"/", N_f[1:5], " SE = ",c(round(SE[1:5],0)), "\n"))
cat(paste0("Estimated/true male abundance: ",c(round(exp(CK_fit[1,6:10]),0)),"/", N_m[1:5], " SE = ",c(round(SE[1:5],0)), "\n"))