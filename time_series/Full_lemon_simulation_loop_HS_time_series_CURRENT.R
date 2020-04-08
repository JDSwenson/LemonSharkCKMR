#THIS IS THE MOST RECENT SCRIPT 3/26/2020

#Where I left off
#Just changed the abundance/starting values for the parameters, but haven't tested

#NEXT TO DO
# 1. Need to speed up constructing the pairwise comparison matrix. 
# 2. Is there still an error with the calculation of Standard Error, or was that fixed by changing the starting values?


#To quickly navigate:
#Search for DELETE for sections that can likely be deleted
#Search for VARIABLES for variables that can probably be moved to top of script
#Search for TROUBLESHOOT for areas that may have issues
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
sample_ages <- 4 #specify how many age classes are represented in the samples, since abundance estimates are given for the cohorts
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
total_abundance <- c(3000, 4000, 5000, 6000, 7000, 9000, 2000)
N_f <- c(total_abundance/2.5)
N_m <- c(total_abundance-N_f)
Pars=c(log(N_f),log(N_m))

#Set number of samples according to Estimated_truth
estimated_truth <- max(total_abundance)
n_samples <- round(10*sqrt(estimated_truth), 0)
n_samples_per_yr <- round(n_samples/length(study_yrs),0) #Total number of samples taken over length of study

source("functions/get_P_lemon_HS_time_series.R")
source("functions/lemon_neg_log_like_HS.R")

P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
##############################Begin simulation##############################
sim_start_time <- Sys.time()
print(paste0("Simulation started at ", Sys.time()))

iterations <- 100 #Set number of iterations to run in the loop
sim_results <- matrix(0, nrow = iterations, ncol = 7)
for(iter in 1:iterations) {

#############################Sample population##############################
Samples = matrix(0,n_samples,3) #matrix filled with 0s, 500 rows (n_samples) and 3 columns; each row is a sample(piece of tissue)
prob_m <- N_m/(N_f+N_m)

Samples[,1]=rep(c(t_start:t_end),times=c(rep(n_samples_per_yr, 3),n_samples_per_yr+1)) #what year did I capture you? repeat each year for the number of samples that were taken.
Samples[,2]=Samples[,1]-sample(c(0:(sample_ages-1)),n_samples,replace=TRUE,prob=Prop_age) #birth year -- simulate age based on Prop_age, and subtract from capture year to get birth year. We sample from 0:sample_ages-1 so we can include individuals born in year 7 (age 0), but if we enter 0:sample_ages, then we will have one more age class than we want.
Samples[,3]=rbinom(n=n_samples, size=1, prob=prob_m) #sex; 0=female

colnames(Samples) <- c("Capt", "Birth", "Sex")
head(Samples)
tail(Samples)

######################Set up pairwise comparison matrix####################
#We don't want to compare individuals with themselves (hence -1), and we don't want to do each comparison twice (divide by 2)
Data <- data.frame(matrix(0,nrow=n_samples*(n_samples-1)/2,ncol=4)) # massive array for pairwise comparisons

counter=1
#brute force - could probably be improved w/ expand.grid or something. Fills first four columns of Data

#head(Data)
#View(Samples)

#TROUBLESHOOT - need to make this faster
for(iind1 in 1:(n_samples-1)){ #iind1 starts at 1 and counts to n_samples-1
  for(iind2 in (iind1+1):n_samples){ #iind2 starts at 2 and counts to n_samples
    if(Samples[iind1,2] > Samples[iind2,2]){ #compare birth years of each combination of samples. If birth year is greater for iind1 i.e. if iind2 was born before iind1. Lethal sampling, so iind2 starts at the number after iind1, and there is no need to go back before.
      #Fill Data array:
      #Data[,1] = sib1 (older) birth,
      #Data[,2] = sib2 (younger) birth
      Data[counter,1]=Samples[iind2,2] #put birth year of iind2 in sib1 (older) birth
      Data[counter,2]=Samples[iind1,2] #put birth year of iind2 in sib2 (younger) birth
    }
    else{ #birth year greater for second value (or equal) - reverse above
      Data[counter,1]=Samples[iind1,2]
      Data[counter,2]=Samples[iind2,2]
    }
    Data[counter,3]=P$P_Mother[Data[counter,1],Data[counter,2]] #Add to column 3 the probability that half-sibs of those ages are related via MHS
    Data[counter,4]=P$P_Father[Data[counter,1],Data[counter,2]] #Add to column 3 the probability that half-sibs of those ages are related via PHS    
    counter=counter+1
  }
}

colnames(Data) <- c("Old_sib_birth","Young_sib_birth", "Mother_prob", "Father_prob")

#length(Data[,1])
#head(Data2)

#Remove same-cohort comparisons
Data2 <- Data[which(Data$Old_sib_birth != Data$Young_sib_birth),]

Data2$Mom_Matches <- (runif(nrow(Data2))<Data2[,3]) #Randomly assign HS mothers based on probability
Data2$Dad_Matches <- (runif(nrow(Data2))<Data2[,4]) #Randomly assign HS fathers based on probability
#cat(paste("Number of HS:",sum(Data2[,5], Data2[,6]))) 

######################### Format pairwise comparison for CKMR ###########################
##Ultimately, want four dataframes:
#1) Rows contain positive comparisons for mothers
#2) Rows contain positive comparisons for fathers
#3) Rows contain negative comparisons for mothers
#4) Rows contain negtive comparisons for fathers

###################Create (and band-aid) dataframes for CKMR model####################
Data_mom_yes <- Data2[Data2$Mom_Matches == TRUE, c(1:2)] 
Data_dad_yes <- Data2[Data2$Dad_Matches == TRUE, c(1:2)]
Data_dad_no <- Data2[Data2$Mom_Matches == FALSE, c(1:2)]
Data_mom_no <- Data2[Data2$Dad_Matches == FALSE, c(1:2)]

library(plyr) #for frequencies
Data_dad_no <- plyr::count(Data_dad_no[,1:2]) #count (from package plyr) gives frequencies of each combination of values
Data_mom_no <- plyr::count(Data_mom_no[,1:2])
Data_dad_yes <- plyr::count(Data_dad_yes[,1:2])
Data_mom_yes <- plyr::count(Data_mom_yes[,1:2])
##END of data generation -- we've set up a specific population and age/sex structure, with no growth over time.

####LOOP####
#if(nrow(Data_dad_yes)>=1 & nrow(Data_mom_yes)>=1){
  
  #########################Fit model!##############################
    CK_fit <- optimx(par=Pars,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Mother=Data_mom_no,Negatives_Father=Data_dad_no,Pairs_Mother=Data_mom_yes,Pairs_Father=Data_dad_yes,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)

    summary(CK_fit)
    #compute variance covariance matrix
    D=diag(length(Pars))*c(exp(CK_fit$p1[1]),exp(CK_fit$p2[1]), exp(CK_fit$p3[1]), exp(CK_fit$p4[1]), exp(CK_fit$p5[1]), exp(CK_fit$p6[1]), exp(CK_fit$p7[1]), exp(CK_fit$p8[1]), exp(CK_fit$p9[1]), exp(CK_fit$p10[1]), exp(CK_fit$p11[1]), exp(CK_fit$p12[1]), exp(CK_fit$p13[1]), exp(CK_fit$p14[1])) #derivatives of transformations
    VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
    VC = (t(D)%*%VC_trans%*%D) #delta method
    SE=sqrt(diag(VC))

    cat(paste("Estimated mature female abundance: ",c(exp(CK_fit[1,1:7])),", SE = ",c(SE[1:7]),"\n"))
    cat(paste("Estimated mature male abundance: ",c(exp(CK_fit[1,8:14])),", SE = ",c(SE[8:14]),"\n"))
    
    #cat(paste("Estimated mature female abundance: ",exp(CK_fit$p1),", SE = ",SE[1],"\n"))
    #cat(paste("Estimated mature male abundance: ",exp(CK_fit$p2),", SE = ",SE[2],"\n"))
    #cat(paste("Estimated survival: ", CK_fit$p3, ", SE = ", SE[3], "\n"))

    sim_results[iter, 1] <- round(exp(CK_fit$p1[1]),0) #Nf
    sim_results[iter, 2] <- round(SE[1],0) #NfSE
    sim_results[iter, 3] <- round(exp(CK_fit$p2[1]),0) #Nm
    sim_results[iter, 4] <- round(SE[2],0) #NmSE
    sim_results[iter, 5] <- sum(Data_dad_yes[,3]) #Number of HS fathers
    sim_results[iter, 6] <- sum(Data_mom_yes[,3]) #Number of HS mothers
    sim_results[iter, 7] <- n_samples
    save(CK_fit, file=paste0("models/Lemon_CKModel_HS_sim_", n_samples,"_samples_", iter))
  } else {
    sim_results[iter, 1]=NA
    sim_results[iter, 2]=NA
    sim_results[iter, 3]=NA
    sim_results[iter, 4]=NA
    sim_results[iter, 5]=NA
    sim_results[iter, 6]=NA
    sim_results[iter, 7]=NA
  }
  print(paste0("finished iteration", iter, " at: ", Sys.time()))
#}

sim_end_time <- Sys.time()
sim_end_time-sim_start_time
#Need to figure out which units for below
#print(paste0("Run time of simulation: ", round(sim_end_time - sim_start_time, 2), ""))
  print(paste0("Simulation finished at ", Sys.time()))
colnames(sim_results) <- c("Nf", "NfSE", "Nm", "NmSE", "Fathers_detected", "Mothers_detected", "Total_samples")
#head(sim_results, 30)

write.table(sim_results, file = paste0("Halfsib_sim_", n_samples, "samps_11.30.2019.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
