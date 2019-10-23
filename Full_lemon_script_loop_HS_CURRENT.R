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
y_start <- 2004 #Start year of study
y_end <- 2007 #End year of study
yr_1=y_start-n_ages 
yrs <- c(seq(yr_1,y_end))
samp_yrs <- c(y_start:y_end)

#This will change actual years into "study years", so we set up variables for the conversation later
t_start=26 #First year of sampling
t_end=29 #Last year of sampling
study_yrs <- c(1:((y_end+1)-yr_1)) #add the +1 so the first and last years of sampling are counted
samp_study_yrs = c(26:29)
#study_yrs[which(yrs == lems$DOB[1])]

####Set up simulation parameters####
#Set Estimated_truth to set the sample size
Estimated_truth <- 500
n_samples <- round(10*sqrt(Estimated_truth), 0) #Total number of samples taken over length of study
a_priori_abund <- 500 #For prior estimation - sets initial parameter values, so should match Estimated_truth value
a_priori_surv <- 0.85
iterations <- 100 #Set number of iterations to run in the loop

source("functions/get_P_lemon_HS.R")
source("functions/lemon_neg_log_like_HS.R")

####################Read in data####################
lemon_data <- read.csv("Main_lemon_shark.csv", header=TRUE, stringsAsFactors = FALSE)
lemon_data <- lemon_data[which(lemon_data$Island=="BIM"),] #Subset for bimini
lemon_data %>% separate(Tube.Label, sep=";", into = c("Juvenile_ID", "Recapture_1", "Recapture_2", "Recapture_3", "Recapture_4")) -> lemon_data2 #Separate Tube labels into capture history (there is one tube label per sample instance)
colnames(lemon_data2) <- c("PIT_tag","Capture_ID", "Recapture_1", "Recapture_2", "Recapture_3", "Recapture_4", "Capture_Year", "Capture_Date", "Island", "Site", "Sex", "PCL_cm", "TL_cm", "DOB", "Father", "Mother") #relabel columns
#head(lemon_data2)

lemon_ref <- lemon_data2[which(lemon_data2$DOB!=""),] #Subset for juveniles with known birth dates

#head(lemon_ref)
#length(lemon_data2[,1]) #How many total bimini lemon sharks in the dataset
#length(lemon_ref[,1]) #How many bimini lemon sharks with known birth dates

Juv_ref <- lemon_ref[,c(2,7,11,14:16)] #Add columns 15 & 16 to add parent IDs
Juv_ref %>% separate(DOB, sep="/", into=c("Yr1", "Yr2")) -> Juv_ref #Separate uncertain birth years (e.g. 1993/1994) into separate columns
colnames(Juv_ref) <- c("Indiv_ID", "Capture_Year","Sex","Yr1","Yr2","Father","Mother")
#head(Juv_ref)

#Randomly assign year if more than one specified
assign_yr <- function(yr1,yr2) {
  if(is.na(yr2)==TRUE) c <- yr1
  else sample(c(yr1, yr2), 1) -> c
  return(c)
}

DOB = c()
for(i in 1:length(Juv_ref[,1])) {
  assign_yr(Juv_ref$Yr1[i], Juv_ref$Yr2[i]) -> DOB[i]
}

Juv_ref$DOB = DOB #Add assigned birth year to Juv_ref
Juv_ref$Father <- sub("^$", "Unknown", Juv_ref$Father) #Replace empty Father cells with Unknown
Juv_ref$Mother <- sub("^$", "Unknown", Juv_ref$Mother) #Replace empty Mother cells with Unknown
#Remove duplicate entries from Juv_ref
Juv_ref <- Juv_ref[!duplicated(Juv_ref$Indiv_ID),]

#Set up dataframe for population we will sample
# Want to base estimate on animals caught between 2004-2007 (the most heavily sampled years) so subset for animals captured in 2004-2007
lems <- subset(Juv_ref, subset=Capture_Year %in% c(2004:2007), select=c(Indiv_ID, Capture_Year, Father, Mother, DOB, Sex))
#head(Juv_ref)
#head(lems)
#Create vector of birth year for juveniles that reflects birth year relative to study
Juv_study_DOB <- c()
for(i in 1:length(lems[,1])){
  Juv_study_DOB[i] <- study_yrs[which(yrs == lems$DOB[i])]
}

#Create reference dataframe for juveniles from which we will sample
Juv_lems <- as.data.frame(cbind(lems$Indiv_ID, Juv_study_DOB, lems$Capture_Year, lems$Father, lems$Mother), stringsAsFactors = FALSE) 
colnames(Juv_lems) <- c("Indiv", "Juv_study_DOB", "Capture_Year", "Father", "Mother")
#head(Juv_lems)

#Add column for capture year that corresponds to study year (e.g. 26 instead of 2004)
for(i in 1:length(Juv_lems[,1])){ 
  Juv_lems$Capt_study_yr[i] <- samp_study_yrs[which(samp_yrs==Juv_lems$Capture_Year[i])]
}
Juv_lems$Capture_Year <- NULL
##############################Begin simulation##############################
sim_start_time <- Sys.time()
print(paste0("Simulation N", Estimated_truth, " started at ", Sys.time()))

sim_results <- matrix(0, nrow = iterations, ncol = 6)
#for(iter in 1:iterations) {

#############################Sample population##############################
head(Juv_lems)
samps <- c()
#Make vector called samps that randomly samples the population
for(i in 1:4){
  cp_yr = c(2004:2007)[i]
  samps <- c(samps, sample(Juv_lems[,1], size=ceiling(n_samples/length(samp_study_yrs)), replace=FALSE))
}

#Create sample dataframe that has year of capture and year of birth
Samp_birth = c() #Initialize vectors
Samp_cap_yr = c()
Samp_mom = c()
Samp_dad = c()

#Simulate birth year based on study year (i.e. 26 instead of 2004)
for(i in 1:length(samps)){
    Samp_birth[i] <- as.numeric(Juv_lems$Juv_study_DOB[which(Juv_lems==samps[i])])
}

#search for the Indiv ID and extract the associated capture year
for(i in 1:length(samps))Samp_cap_yr[i] = as.numeric(Juv_lems[which(Juv_lems$Indiv==samps[i]),5])

#search for the Indiv ID and extract the associated mother
for(i in 1:length(samps))Samp_mom[i] = Juv_lems[which(Juv_lems$Indiv==samps[i]),4]

#search for the Indiv ID and extract the associated father
for(i in 1:length(samps))Samp_dad[i] = Juv_lems[which(Juv_lems$Indiv==samps[i]),3]

#Combine vectors into dataframe and rename columns
Samples <- as.data.frame(cbind(samps,Samp_birth,Samp_cap_yr, Samp_dad, Samp_mom), stringsAsFactors = FALSE)

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