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
sample_ages <- 7 #specify how many cohorts are represented in the samples, since abundance estimates are given for the cohorts
Surv_age <- rep(0.85, sample_ages) #Set survival for adults
Surv_age[1] <- 0.55 #First year survival at 0.55; rest at 0.85
#Surv_age[1:2] <- c(rep(0.55,2)) #First two years, survival at 0.55
Prop_age <- rep(1, sample_ages) #Initialize variable Prop_age
for(iage in 2:length(Prop_age)) Prop_age[iage]=Prop_age[iage-1]*Surv_age[iage-1] #Set proportion of animals from each cohort surviving to the next age class based on Surv_age
Prop_age=Prop_age/sum(Prop_age)  #stable age distribution given assumed survival i.e. gives proportion of ages in a population based on assumed survival

est_ages <- 4 #Number of cohorts for which we'll estimate abundance
min_est_cohort <- 4 #Earliest year of estimation
n_yrs=est_ages+(est_ages-1) #Number of years from earliest individual to end of study

#This will change actual years into "study years", so we set up variables for the conversation later
t_start=4 #First year of sampling
t_end=7 #Last year of sampling
study_yrs = c(t_start:t_end)
#length(samp_study_yrs)

####Set up simulation parameters####

#Set Estimated_truth and number of males and females
#n_yrs
total_abundance <- c(5000, 5000, 5000, 5000) #Can't estimate for yr=1 because abundance is calculated at time of younger sibling's birth, and there is no scenario where the younger sibling was born in yr1. Also unlikely to get estimate for yr=2 because there are relatively few pairwise comparisons where the younger sibling was born in yr2. So length of vector should be n_yrs-2
N_f <- c(total_abundance/2.5)
N_m <- c(total_abundance-N_f)
Pars=c(log(N_f),log(N_m))
est_yrs <- n_yrs-2 #Number of years we're estimating

#Set number of samples according to Estimated_truth
estimated_truth <- max(total_abundance)
n_samples <- round(10*sqrt(estimated_truth)*2, 0)
#n_samples_per_yr <- round(n_samples/length(study_yrs),0) #Total number of samples taken over length of study

source("./time_series/models/get_P_lemon_HS_time_series_4yrs.R")
source("./likelihood_functions/lemon_neg_log_like_HS.R")

P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
#P$P_Mother
#P$P_Father

##############################Begin simulation##############################
sim_start_time <- Sys.time()
print(paste0("Simulation started at ", Sys.time()))

iterations <- 500 #Set number of iterations to run in the loop
sim_results <- matrix(0, nrow = iterations, ncol = 39)

for(iter in 1:iterations) {

#############################Sample population##############################
Samples <- c()
  #Samples = matrix(0,n_samples,2) #matrix filled with 0s; each row is a sample(piece of tissue)
prob_m <- N_m/(N_f+N_m)

#Samples=sample(c(7:1),n_samples,replace=TRUE) #birth year - probability of sampling is random
Samples=sample(c(7:1),n_samples,replace=TRUE,prob=Prop_age) #birth year - probability of sampling based on Prop_age

#Samples[,1]=rep(c(t_start:t_end),times=c(rep(n_samples_per_yr, 3),n_samples_per_yr+1)) #what year did I capture you? repeat each year for the number of samples that were taken.
#Samples[,2]=sample(c(7:1),n_samples,replace=TRUE,prob=Prop_age) #birth year -- simulate age based on Prop_age, and subtract from capture year to get birth year. We sample from 0:sample_ages-1 so we can include individuals born in the last year (age 0), but if we enter 0:sample_ages, then we will have one more cohort than we want.
#Samples[,3]=rbinom(n=n_samples, size=1, prob=prob_m) #sex; 0=female -- moved this to below

#colnames(Samples) <- c("Capt", "Birth")
#head(Samples)
#tail(Samples)

#Data <- t(combn(Samples[,2], m=2))
Data <- t(combn(Samples, m=2))
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
Data2 <- Data2[which(Data2$Young_sib_birth>=4),]

#Add probabilities of parentage to columns 4 & 5
for(i in 1:nrow(Data2)){
  Data2$Mom_prob[i] <- P$P_Mother[Data2[i,1],Data2[i,2]]
  Data2$Dad_prob[i] <- P$P_Father[Data2[i,1],Data2[i,2]]
}
Data2

#Randomly assign matches
for(i in 1:nrow(Data2)){
  Data2$Mom_Matches_unif[i] <- sum(runif(Data2[i,3])<Data2[i,4])
  Data2$Dad_Matches_unif[i] <- sum(runif(Data2[i,3])<Data2[i,5])
  Data2$Mom_Matches_pois[i] <- sum(rpois(n=Data2$freq[i], lambda=Data2$Mom_prob[i]))
  Data2$Dad_Matches_pois[i] <- sum(rpois(n=Data2$freq[i], lambda=Data2$Dad_prob[i]))
  }
Data2

######################### Format pairwise comparison for CKMR ###########################
##Ultimately, want four dataframes:
#1) Rows contain positive comparisons for mothers
#2) Rows contain positive comparisons for fathers
#3) Rows contain negative comparisons for mothers
#4) Rows contain negtive comparisons for fathers

###################Create dataframes for CKMR model####################
Data_mom_yes <- Data2[which(Data2$Mom_Matches_pois >0), c(1:2,6)]
Data_dad_yes <- Data2[which(Data2$Dad_Matches_pois >0), c(1:2,7)]
Data_dad_no = Data_mom_no <- Data2
Data_dad_no[,3] <- Data_dad_no$freq - Data_dad_no$Dad_Matches_pois
Data_dad_no <- Data_dad_no[,c(1:3)]
Data_mom_no[,3] <- Data_mom_no$freq - Data_mom_no$Mom_Matches_pois
Data_mom_no <- Data_mom_no[,c(1:3)]

####LOOP####
#if(nrow(Data_dad_yes)>=1 & nrow(Data_mom_yes)>=1){
  
  #########################Fit model!##############################
    CK_fit <- optimx(par=Pars,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Mother=Data_mom_no,Negatives_Father=Data_dad_no,Pairs_Mother=Data_mom_yes,Pairs_Father=Data_dad_yes,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)

    #summary(CK_fit)
  
    #Without SE
      #exp(CK_fit[1,1:5])
      #N_f
      #exp(CK_fit[1,6:10])
      #N_m
 
    #Store estimates in vector -- need to change based on number of parameters
      estimates <- c(exp(CK_fit$p1[1]),exp(CK_fit$p2[1]), exp(CK_fit$p3[1]), exp(CK_fit$p4[1]), exp(CK_fit$p5[1]), exp(CK_fit$p6[1]), exp(CK_fit$p7[1]), exp(CK_fit$p8[1]))
      
      #compute variance covariance matrix
    D=diag(length(Pars))*estimates #derivatives of transformations
    
      VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
    VC = (t(D)%*%VC_trans%*%D) #delta method
    SE=sqrt(diag(VC))

  cat(paste("Estimated mature female abundance: ",c(exp(CK_fit[1,1:4])),", SE = ",c(SE[1:4]),"\n"))
  cat(paste("Estimated mature male abundance: ",c(exp(CK_fit[1,5:8])),", SE = ",c(SE[5:8]),"\n"))
    
    sim_results[iter, 1] <- round(exp(CK_fit$p1[1]),0) #Nf
    sim_results[iter, 2] <- round(SE[1],0) #NfSE
    sim_results[iter, 3] <- exp(Pars[1]) #Truth
    sim_results[iter, 4] <- sum(Data_mom_yes[which(Data_mom_yes$Young_sib_birth==4),3]) #Number of HS mothers
    sim_results[iter, 5] <- round(exp(CK_fit$p2[1]),0) #Nf
    sim_results[iter, 6] <- round(SE[2],0) #NfSE
    sim_results[iter, 7] <- exp(Pars[2]) #Truth
    sim_results[iter, 8] <- sum(Data_mom_yes[which(Data_mom_yes$Young_sib_birth==5),3]) #Number of HS mothers
    sim_results[iter, 9] <- round(exp(CK_fit$p3[1]),0) #Nf
    sim_results[iter, 10] <- round(SE[3],0) #NfSE
    sim_results[iter, 11] <- exp(Pars[3]) #Truth
    sim_results[iter, 12] <- sum(Data_mom_yes[which(Data_mom_yes$Young_sib_birth==6),3]) #Number of HS mothers
    sim_results[iter, 13] <- round(exp(CK_fit$p4[1]),0) #Nf
    sim_results[iter, 14] <- round(SE[4],0) #NfSE
    sim_results[iter, 15] <- exp(Pars[4]) #Truth
    sim_results[iter, 16] <- sum(Data_mom_yes[which(Data_mom_yes$Young_sib_birth==7),3]) #Number of HS mothers
    sim_results[iter, 17] <- round(exp(CK_fit$p5[1]),0) #Nm
    sim_results[iter, 18] <- round(SE[5],0) #NmSE
    sim_results[iter, 19] <- exp(Pars[5]) #Truth
    sim_results[iter, 20] <- sum(Data_dad_yes[which(Data_dad_yes$Young_sib_birth==4),3]) #Number of HS fathers
    sim_results[iter, 21] <- round(exp(CK_fit$p6[1]),0) #Nm
    sim_results[iter, 22] <- round(SE[6],0) #NmSE
    sim_results[iter, 23] <- exp(Pars[6]) #Truth
    sim_results[iter, 24] <- sum(Data_dad_yes[which(Data_dad_yes$Young_sib_birth==5),3]) #Number of HS fathers
    sim_results[iter, 25] <- round(exp(CK_fit$p7[1]),0) #Nm
    sim_results[iter, 26] <- round(SE[7],0) #NmSE
    sim_results[iter, 27] <- exp(Pars[7]) #Truth
    sim_results[iter, 28] <- sum(Data_dad_yes[which(Data_dad_yes$Young_sib_birth==6),3]) #Number of HS fathers
    sim_results[iter, 29] <- round(exp(CK_fit$p8[1]),0) #Nm
    sim_results[iter, 30] <- round(SE[8],0) #NmSE
    sim_results[iter, 31] <- exp(Pars[8]) #Truth
    sim_results[iter, 32] <- sum(Data_dad_yes[which(Data_dad_yes$Young_sib_birth==7),3]) #Number of HS fathers
    sim_results[iter, 33] <- sum(Samples==1)
    sim_results[iter, 34] <- sum(Samples==2)
    sim_results[iter, 35] <- sum(Samples==3)
    sim_results[iter, 36] <- sum(Samples==4)
    sim_results[iter, 37] <- sum(Samples==5)
    sim_results[iter, 38] <- sum(Samples==6)
    sim_results[iter, 39] <- sum(Samples==7)
    
    save(CK_fit, file=paste0("../Results/model_objects/Lemon_CKModel_HS_time_series_iteration_", iter))
    #  } else {
    #sim_results[iter, 1]=NA
    #sim_results[iter, 2]=NA
    #sim_results[iter, 3]=NA
    #sim_results[iter, 4]=NA
    #sim_results[iter, 5]=NA
    #sim_results[iter, 6]=NA
    #sim_results[iter, 7]=NA
    #}
    print(paste0("finished iteration", iter, " at: ", Sys.time()))
}

sim_end_time <- Sys.time()
sim_end_time-sim_start_time
#Need to figure out which units for below

#print(paste0("Run time of simulation: ", round(sim_end_time - sim_start_time, 2), ""))
print(paste0("Simulation finished at ", Sys.time()))
colnames(sim_results) <- c("Nf_yr4", "NfSE_yr4", "NfTruth_yr4", "Moms_detected_yr4",
                           "Nf_yr5", "NfSE_yr5", "NfTruth_yr5", "Moms_detected_yr5",
                           "Nf_yr6", "NfSE_yr6", "NfTruth_yr6", "Moms_detected_yr6",
                           "Nf_yr7", "NfSE_yr7", "NfTruth_yr7", "Moms_detected_yr7",
                           "Nm_yr4", "NmSE_yr4", "NmTruth_yr4", "Dads_detected_yr4",
                           "Nm_yr5", "NmSE_yr5", "NmTruth_yr5", "Dads_detected_yr5",
                           "Nm_yr6", "NmSE_yr6", "NmTruth_yr6", "Dads_detected_yr6",
                           "Nm_yr7", "NmSE_yr7", "NmTruth_yr7", "Dads_detected_yr7", "Yr1_samples", "Yr2_samples", "Yr3_samples", "Yr4_samples", "Yr5_samples", "Yr6_samples", "Yr7_samples")
#head(sim_results, 30)

write.table(sim_results, file = paste0("HS_time_series_null_4.30.2020.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
