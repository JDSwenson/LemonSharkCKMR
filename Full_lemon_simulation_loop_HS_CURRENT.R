#THIS IS THE MOST RECENT SCRIPT 9/20/2019
#To quickly navigate:
#Search for DELETE for sections that can likely be deleted
#Search for VARIABLES for variables that can probably be moved to top of script
#Search for Troubleshoot for areas that may have issues
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


##############Set up experiment parameters##############
n_ages=25 #Max age of lemon sharks
n_yrs=29 #Number of years from earliest individual to end of study
Surv_age = rep(0.85,n_ages) #Set survival for each age to 95%
Surv_age[1]=0.85 #Set survival for first year to 85%
Prop_age = rep(1,n_ages) #Set proportion of eage age class to 1

for(iage in 2:n_ages)Prop_age[iage]=Prop_age[iage-1]*Surv_age[iage-1] #Sets proportion of animals from each cohort surviving to iage based on Surv_age
Prop_age=Prop_age/sum(Prop_age)  #stable age distribution given assumed survival i.e. gives proportion of ages in a population based on assumed survival

t_start=26 #First year of sampling
t_end=29 #Last year of sampling
samp_study_yrs = c(26:29)
n_study_yrs <- length(samp_study_yrs)

####Set up simulation parameters####
#Set Estimated_truth to set the sample size
Estimated_truth <- 500
n_samples <- round(10*sqrt(Estimated_truth), 0) #Total number of samples taken over length of study
a_priori_abund <- 500 #For prior estimation - sets initial parameter values, so should match Estimated_truth value
a_priori_surv <- 0.85
iterations <- 100 #Set number of iterations to run in the loop

source("functions/get_P_lemon_HS.R")
source("functions/lemon_neg_log_like_HS.R")

##############################Begin simulation##############################
sim_start_time <- Sys.time()
print(paste0("Simulation N", Estimated_truth, " started at ", Sys.time()))

sim_results <- matrix(0, nrow = iterations, ncol = 6)
#for(iter in 1:iterations) {

#############################Sample population##############################
Samples = matrix(0,n_samples,3)#Make vector called samps that randomly samples the population
Samples[,1]=rep(c(t_start:t_end),each=n_samples/n_study_yrs) #Capture year
Samples[,2]=Samples[,1]-sample(c(1:n_ages),n_samples,replace=TRUE,prob=Prop_age) #birth year
Samples[,3]=round(runif(n_samples)) #sex; 0=female
colnames(Samples) <- c("Indiv_ID", "Birth", "Capt", "Father", "Mother")
#head(Samples)

#Change columns to numeric; changing earlier doesn't work
Samples[,2] <- as.numeric(Samples[,2])
Samples[,3] <- as.numeric(Samples[,3])
#Samples
#head(Juv_ref)
#tail(Samples)


######################Set up pairwise comparison matrix####################
#We don't want to compare individuals with themselves (hence -1), and we don't want to do each comparison twice (divide by 2)
Data <- data.frame(matrix(0,nrow=n_samples*(n_samples-1)/2,ncol=3)) # massive array for pairwise comparisons

counter=1
#brute force - could probably be improved w/ expand.grid or something. Fills first four columns of Data

#head(Data)
#View(Samples)

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
    counter=counter+1
  }
}

colnames(Data) <- c("Old_sib_birth","Young_sib_birth", "Matches")
#length(Data[,1])
#head(Data)

#Remove same-cohort comparisons
Data2 <- Data[which(Data$Old_sib_birth != Data$Young_sib_birth),]
Data_Mother = Data_Father = plyr::count(Data2[,])
colnames(Data_Mother)[4] = colnames(Data_Father)[4] = "No_Matches"
head(Samples)

#Pairwise comparison among sample birth years; if they have the same mother, add 1 to the Matches column and subtract one from the No_Matches column for the appropriate comparison
for(i in 1:(n_samples-1)){ #iind1 starts at 1 and counts to n_samples-1
  for(j in (i+1):n_samples){ #iind2 starts at 2 and counts to n_samples
    if(Samples$Birth[i] != Samples$Birth[j] & Samples$Mother[i] == Samples$Mother[j] & Samples$Mother[i] != "Unknown") {
      Data_Mother[which(Data_Mother$Young_sib_birth == Samples$Birth[j] & Data_Mother$Old_sib_birth == Samples$Birth[i]),3] <- Data_Mother[which(Data_Mother$Young_sib_birth == Samples$Birth[j] & Data_Mother$Old_sib_birth == Samples$Birth[i]),3] +1
      Data_Mother[which(Data_Mother$Young_sib_birth == Samples$Birth[j] & Data_Mother$Old_sib_birth == Samples$Birth[i]),4] <- Data_Mother[which(Data_Mother$Young_sib_birth == Samples$Birth[j] & Data_Mother$Old_sib_birth == Samples$Birth[i]),4] -1
      }
  }
}

#Data_Mother
#Data_Father

#Same for Fathers
for(i in 1:(n_samples-1)){ #iind1 starts at 1 and counts to n_samples-1
  for(j in (i+1):n_samples){ #iind2 starts at 2 and counts to n_samples
    if(Samples$Birth[i] != Samples$Birth[j] & Samples$Father[i] == Samples$Father[j] & Samples$Father[i] != "Unknown") {
      Data_Father[which(Data_Father$Young_sib_birth == Samples$Birth[j] & Data_Father$Old_sib_birth == Samples$Birth[i]),3] <- Data_Father[which(Data_Father$Young_sib_birth == Samples$Birth[j] & Data_Father$Old_sib_birth == Samples$Birth[i]),3] +1
      Data_Father[which(Data_Father$Young_sib_birth == Samples$Birth[j] & Data_Father$Old_sib_birth == Samples$Birth[i]),4] <- Data_Father[which(Data_Father$Young_sib_birth == Samples$Birth[j] & Data_Father$Old_sib_birth == Samples$Birth[i]),4] -1
    }
  }
}
#View(Samples)
#Data_Father

######################### Format pairwise comparison for CKMR ###########################
##Ultimately, want four dataframes:
#1) Rows contain positive comparisons for mothers
#2) Rows contain positive comparisons for fathers
#3) Rows contain negative comparisons for mothers
#4) Rows contain negtive comparisons for fathers

###################Create (and band-aid) dataframes for CKMR model####################
Data_mom_yes <- Data_Mother[Data_Mother$Matches > 0, c(1:3)] 
Data_dad_yes <- Data_Father[Data_Father$Matches > 0, c(1:3)]
Data_dad_no <- Data_Father[Data_Father$No_Matches > 0, c(1:2,4)]
Data_mom_no <- Data_Mother[Data_Mother$No_Matches > 0, c(1:2,4)]

####LOOP####
#if(nrow(Data_dad_yes)>=1 & nrow(Data_mom_yes)>=1){
  
  #########################Fit model!##############################
    CK_fit <- optimx(par=Pars,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Mother=Data_mom_no,Negatives_Father=Data_dad_no,Pairs_Mother=Data_mom_yes,Pairs_Father=Data_dad_yes,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)

    summary(CK_fit)
    #compute variance covariance matrix
    D=diag(length(Pars))*c(exp(CK_fit$p1[1]),exp(CK_fit$p2[1])) #derivatives of transformations
    VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
    VC = (t(D)%*%VC_trans%*%D) #delta method
    SE=sqrt(diag(VC))
    
    cat(paste("Estimated mature female abundance: ",exp(CK_fit$p1),", SE = ",SE[1],"\n"))
    cat(paste("Estimated mature male abundance: ",exp(CK_fit$p2),", SE = ",SE[2],"\n"))
    #cat(paste("Estimated survival: ", CK_fit$p3, ", SE = ", SE[3], "\n"))

    
    
    
    sim_results[iter, 1] <- round(exp(CK_fit$p1[1]),0) #Nf
    sim_results[iter, 2] <- round(SE[1],0) #NfSE
    sim_results[iter, 3] <- round(exp(CK_fit$p2[1]),0) #Nm
    sim_results[iter, 4] <- round(SE[2],0) #NmSE
    sim_results[iter, 5] <- Estimated_truth
    sim_results[iter, 6] <- n_samples
    save(CK_fit, file=paste0("Lemon_CKModel_HS_aprioriN_",Estimated_truth,"_", iter))
  } else {
    sim_results[iter, 1]=NA
    sim_results[iter, 2]=NA
    sim_results[iter, 3]=NA
    sim_results[iter, 4]=NA
    sim_results[iter, 5]=NA
    sim_results[iter, 6]=NA
  }
  print(paste0("finished iteration", iter, " at: ", Sys.time()))
}

sim_end_time <- Sys.time()
sim_end_time-sim_start_time
#Need to figure out which units for below
#print(paste0("Run time of simulation: ", round(sim_end_time - sim_start_time, 2), ""))
  print(paste0("Simulation N", Estimated_truth, " finished at ", Sys.time()))
colnames(sim_results) <- c("Nf", "NfSE", "Nm", "NmSE", "Estimated_truth", "Total_samples")
#head(sim_results, 30)

write.table(sim_results, file = paste0("Lemon_CKMR_loop_nonlethal_aprioriN_",Estimated_truth,".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)