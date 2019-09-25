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

#Stable age distribution for adult males (We have ages for juveniles)
max_age = 25
m_adult_age <- c(12:25) #Set ages at which males are mature -- used to simulate ages
m_mat <- c(rep(0,11), rep(1,14)) #Set proportion of mature males at each age -- assumes knife-edge maturity
Surv_age_m <- rep(0.85,length(m_adult_age)) #Set survival for adults
Prop_age_m <- rep(1, length(m_adult_age)) #Initialize variable Prop_age_m
for(iage in 2:length(Prop_age_m)) Prop_age_m[iage]=Prop_age_m[iage-1]*Surv_age_m[iage-1] #Set proportion of animals from each cohort surviving to the next age class based on Surv_age
Prop_age_m=Prop_age_m/sum(Prop_age_m)  #stable age distribution given assumed survival i.e. gives proportion of ages in a population based on assumed survival

#Stable age distribution for mature females
f_adult_age <- c(13:25)
f_mat <- c(rep(0,12), rep(1,13))
Surv_age_f <- rep(0.85,length(f_adult_age))
Prop_age_f <- rep(1, length(f_adult_age))
for(iage in 2:length(Prop_age_f)) Prop_age_f[iage]=Prop_age_f[iage-1]*Surv_age_f[iage-1] 
Prop_age_f=Prop_age_f/sum(Prop_age_f)

##############Set up experiment parameters##############
n_ages=25 #Max age of lemon sharks
n_yrs=29 #Number of years from earliest individual to end of study
y_start <- 2004 #Start year of study
y_end <- 2007 #End year of study
yr_1=y_start-max_age 
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
Estimated_truth <- 1000
n_samples <- round(10*sqrt(Estimated_truth), 0) #Total number of samples taken over length of study
a_priori_abund <- 1000 #For prior estimation - sets initial parameter values, so should match Estimated_truth value
a_priori_surv <- 0.85
iterations <- 100 #Set number of iterations to run in the loop

source("functions/get_P_lemon_nonlethal.R")
source("functions/lemon_neg_log_like_nonlethal.R")

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
#Mot_lems <- Mot_lems[!duplicated(Mot_lems$Indiv),]

#Make reference table for adults and frequency of offspring by year
Father_ref <- plyr::count(Juv_ref[,c(6,8)])
Mother_ref <- plyr::count(Juv_ref[,c(7:8)])

#What is the max number of fathers detected in one of the study years?
#Set up dataframe of unique fathers detected during study years
#Father_refmax <- Father_ref[which(Father_ref$Father != "Unknown"),] #Remove unknowns
#Father_refmax <- Father_refmax[which(Father_refmax$DOB %in% samp_yrs),] #subset for study years
#table(Father_refmax[,2])

#What is the max number of mothers detected in one of the study years?
#Set up dataframe of unique mothers detected during study years
Mother_refmax <- Mother_ref[which(Mother_ref$Mother != "Unknown"),] #Remove unknowns
Mother_refmax <- Mother_refmax[which(Mother_refmax$DOB %in% samp_yrs),] #Subset for sample years
#table(Mother_refmax[,2])

#length(Father_ref[,1])
#head(Juv_ref)
#head(Father_ref)
#tail(Father_ref)
#tail(Mother_ref)

####Simulate ages for parents####
Uniq_Father <- unique(Father_ref$Father) #Save unique fathers
Uniq_Mother <- unique(Mother_ref$Mother)

#Simulate ages of mothers

#Set up dataframe for population we will sample
# Want to base estimate on animals caught between 2004-2007 (the most heavily sampled years) so subset for animals captured in 2004-2007
lems <- subset(Juv_ref, subset=Capture_Year %in% c(2004:2007), select=c(Indiv_ID, Capture_Year, Father, Mother, DOB, Sex))
#head(Juv_ref)
#View(lems)

##############################Begin simulation##############################
sim_start_time <- Sys.time()
print(paste0("Simulation N", Estimated_truth, " started at ", Sys.time()))

sim_results <- matrix(0, nrow = iterations, ncol = 6)
for(iter in 1:iterations) {
  
#Create vector of birth year for juveniles that reflects birth year relative to study
Juv_study_DOB <- c()
for(i in 1:length(lems[,1])){
  Juv_study_DOB[i] <- study_yrs[which(yrs == lems$DOB[i])]
}

#Juv_lems is referenced later when we add values from Juv_study_DOB to the sampled individuals
#Juv_lems2 will be combined with father and mother dataframes for sampling
Juv_lems <- as.data.frame(cbind(lems$Indiv_ID, Juv_study_DOB, lems$Capture_Year, lems$Sex), stringsAsFactors = FALSE) 
colnames(Juv_lems) <- c("Indiv", "Juv_study_DOB", "Capture_Year", "Sex")
Juv_lems2 <- Juv_lems[,c(1,3:4)] #Subset for columns also present in Father and Mother dataframes so they can be combined

#Create dataframe for fathers
Fat_lems <- as.data.frame(cbind(lems$Father, lems$Capture_Year, rep("M", length(lems[,1]))), stringsAsFactors = FALSE)

#Create dataframe for mothers
Mot_lems <- as.data.frame(cbind(lems$Mother, lems$Capture_Year, rep("F", length(lems[,1]))), stringsAsFactors = FALSE)

#Name all columns appropriately
colnames(Fat_lems) = colnames(Mot_lems) = colnames(Juv_lems2) = c("Indiv", "Capture_Year", "Sex")

#Add column for capture year that corresponds to study year (e.g. 26 instead of 2004)
for(i in 1:length(Juv_lems2[,1])){ 
  Juv_lems2$Capt_study_yr[i] <- samp_study_yrs[which(samp_yrs==Juv_lems2$Capture_Year[i])]
}
for(i in 1:length(Fat_lems[,1])){
  Fat_lems$Capt_study_yr[i] <- samp_study_yrs[which(samp_yrs==Fat_lems$Capture_Year[i])]
}
for(i in 1:length(Mot_lems[,1])){
  Mot_lems$Capt_study_yr[i] <- samp_study_yrs[which(samp_yrs==Mot_lems$Capture_Year[i])]
}

#head(Fat_lems)
#head(Mot_lems)
#length(Mot_lems[,1])

#Remove duplicate entries (There are no duplicate juveniles - double-checked 9/13/19)
#Troubleshoot - Below will only keep the first instance of each, so there end up being more captures in first year
Mot_lems <- Mot_lems[!duplicated(Mot_lems$Indiv),]
Mot_lems <- Mot_lems[-c(which(Mot_lems[,1]=="Unknown")),] #Remove Unknown individuals
#length(Mot_lems[,1]) #Check number of unique mothers

#Fat_lems$uniqID <- paste0(Fat_lems$Indiv,Fat_lems$DOB)
#length(Fat_lems[,1])
Fat_lems <- Fat_lems[!duplicated(Fat_lems$Indiv),]
Fat_lems <- Fat_lems[-c(which(Fat_lems[,1]=="Unknown")),] #Remove Unknown individuals
#length(Fat_lems[,1])

#Combine data frames
pop_df <- c()
pop_df <- rbind(Juv_lems2, Fat_lems, Mot_lems)
#head(pop_df)
#tail(pop_df)

#length(pop_df[,1]) #Total number of animals in sampled population
#length(Juv_lems[,1]) #How many juveniles included in sample pool
#length(Fat_lems[,1]) #How many fathers in sample pool
#length(Mot_lems[,1]) #How many mothers in sample pool

#############################Sample population##############################
samps=c()

#Make vector called samps that randomly samples the population
for(i in 1:4){
  cp_yr = c(2004:2007)[i]
  samps <- c(samps, sample(pop_df[which(pop_df$Capture_Year==cp_yr),1], size=n_samples/length(samp_study_yrs), replace=FALSE))
}

#length(samps)
#head(samps)
#head(pop_df)

#Create sample dataframe that has year of capture(death) year of birth, and sex
Samp_birth = c() #Initialize vectors
Samp_cap_yr = c()
Samp_sex = c()

#Simulate birth year based on study year (i.e. 26 instead of 2004)
for(i in 1:length(samps)){
  if(samps[i] %in% Juv_lems$Indiv){
    Samp_birth[i] <- as.numeric(Juv_lems$Juv_study_DOB[which(Juv_lems==samps[i])])
  } else if(samps[i] %in% Fat_lems$Indiv){
    Samp_birth[i] <- Fat_lems$Capt_study_yr[which(Fat_lems$Indiv==samps[i])] - sample(m_adult_age, replace=TRUE, prob=Prop_age_m, size=1)
  } else if(samps[i] %in% Mot_lems$Indiv){
    Samp_birth[i] <- Mot_lems$Capt_study_yr[which(Mot_lems$Indiv==samps[i])] - sample(f_adult_age, replace=TRUE, prob=Prop_age_f, size=1)
  }
}

#search pop_df for the Indiv ID and extract the associated capture year
for(i in 1:length(samps))Samp_cap_yr[i] = as.numeric(pop_df[which(pop_df$Indiv==samps[i]),4])
#search pop_df for the Indiv ID and extract the associated sex
for(i in 1:length(samps))Samp_sex[i] = pop_df[which(pop_df$Indiv==samps[i]),3] 

#Combine vectors into dataframe and rename columns
Samples <- as.data.frame(cbind(samps,Samp_birth,Samp_cap_yr, Samp_sex), stringsAsFactors = FALSE)

colnames(Samples) <- c("Indiv_ID", "Birth", "Capt", "Sex")
#head(Samples)

#Change columns to numeric; changing earlier doesn't work
Samples[,2] <- as.numeric(Samples[,2])
Samples[,3] <- as.numeric(Samples[,3])
#Samples
#head(Juv_ref)
#tail(Samples)

# If known, assign father and mother for each sampled individual
for(i in 1:length(Samples[,1])){
  if(Samples[i,1]%in%Juv_ref[,1]) {
    Samples$Father[i] <- Juv_ref[which(Juv_ref[,1]==Samples[i,1]),6]
    Samples$Mother[i] <- Juv_ref[which(Juv_ref[,1]==Samples[i,1]),7]
  } else {
    Samples$Father[i]=NA
    Samples$Mother[i]=NA
  }
}

Samples$Dad_birth=Samples$Dad_capt=Samples$Mom_birth=Samples$Mom_capt=NA #Set birth and death year of parents to NA (so missing data are easier to work with downstream)

#If the sampled animal also is a parent to another sampled animal, assign it a birth and capture date in the columns Dad_birth and Mom_birth
for(i in 1:length(Samples[,1])){
  if(Samples$Indiv_ID[i] %in% Samples$Father){
    Samples$Dad_birth[which(Samples$Father==Samples$Indiv_ID[i])] <- Samples$Birth[i]
    Samples$Dad_capt[which(Samples$Father==Samples$Indiv_ID[i])] <- Samples$Capt[i]
  }
  if(Samples$Indiv_ID[i] %in% Samples$Mother){
    Samples$Mom_birth[which(Samples$Mother==Samples$Indiv_ID[i])] <- Samples$Birth[i]
    Samples$Mom_capt[which(Samples$Mother==Samples$Indiv_ID[i])] <- Samples$Capt[i]
  }
}

#Create binary Mom and Data column that says whether the sampled individual in the row has a parent that was also sampled.
Samples$Mom=Samples$Dad=0 
for(i in 1:length(Samples[,1])){
  if(is.na(Samples$Dad_birth[i])==FALSE) Samples$Dad[i]=1
  if(is.na(Samples$Mom_birth[i])==FALSE) Samples$Mom[i]=1
}

######################Set up pairwise comparison matrix####################
#We don't want to compare individuals with themselves (hence -1), and we don't want to do each comparison twice (divide by 2)
Data <- data.frame(matrix(0,nrow=n_samples*(n_samples-1)/2,ncol=5)) # massive array for pairwise comparisons

counter=1
#brute force - could probably be improved w/ expand.grid or something. Fills first four columns of Data
for(iind1 in 1:(n_samples-1)){ #iind1 starts at 1 and counts to n_samples-1
  for(iind2 in (iind1+1):n_samples){ #iind2 starts at 2 and counts to n_samples
    if(Samples[iind1,2] > Samples[iind2,2]){ #compare birth years of each combination of samples. If birth year is greater for iind1 i.e. if iind2 was born before iind1. Lethal sampling, so iind2 starts at the number after iind1, and there is no need to go back before.
      #Fill Data array:
      #Data[,1] = adult birth,
      #Data[,2] = adult death year
      #Data[,3] = young birth
      #Data[,4] = adult sex
      Data[counter,1]=Samples[iind2,2] #put birth year of iind2 in adult birth year
      Data[counter,2]=Samples[iind2,3] #put year of capture of iind2 in adult death year (remember, lethal sampling)
      Data[counter,3]=Samples[iind1,2] #put birth year of iind1 in offspring birth year
      Data[counter,4]=Samples[iind2,4] #put sex of iind2 in adult sex
    }
    else{ #birth year greater for second value - reverse above
      Data[counter,1]=Samples[iind1,2]
      Data[counter,2]=Samples[iind1,3]
      Data[counter,3]=Samples[iind2,2]
      Data[counter,4]=Samples[iind1,4]
    }
    counter=counter+1
  }
}
colnames(Data) <- c("Adult_birth","Adult_capt","Offspring_birth","Adult_sex", "Matches")
#length(Data[,1])
#head(Data)

######################### Format pairwise comparison for CKMR ###########################
##Ultimately, want four dataframes:
#1) Rows contain positive comparisons for mothers
#2) Rows contain positive comparisons for fathers
#3) Rows contain negative comparisons for mothers
#4) Rows contain negtive comparisons for fathers

dm <- Data[Data$Adult_sex == "F",] #Subset pairwise matrix for mother comparisons
dd <- Data[Data$Adult_sex == "M",] #Subset pairwise matrix for father comparisons

#Create column for unique combinations of adult birth, adult death, and offspring birth year
dm$MasterID <- paste(dm$Adult_birth, dm$Adult_capt, dm$Offspring_birth, sep="-")
dd$MasterID <- paste(dd$Adult_birth, dd$Adult_capt, dd$Offspring_birth, sep="-")
#head(dd)

#Count unique combos of adult birth, adult death, and offspring birth AND Master ID (should be the same)
dm2 <- plyr::count(dm[,c(1:3,6)])
dd2 <- plyr::count(dd[,c(1:3,6)])
colnames(dm2)[5] = colnames(dd2)[5] <- "No_matches" #Change frequency column to No_matches (assume for now all comparisons are no matches - will subtract matches later)
#head(dd2)

#Subset sampled moms and sampled dads from dataframe of all samples
Sampled_Moms <- Samples[Samples$Mom == 1,]
Sampled_Dads <- Samples[Samples$Dad == 1,]
#View(Samples)

#Create dataframe of sampled moms (sm) that includes a column that counts the occurrances of specific combinations of offspring birth, adult birth, and adult death
sm <- Sampled_Moms %>% 
  group_by(Birth, Mom_birth, Mom_capt) %>%  #ID which columns to group by
  summarise(Matches = sum(!is.na(Mom)))

#Create dataframe of sampled dads (sd) that includes a column that counts the occurrances of specific combinations of offspring birth, adult birth, and adult death
sd <- Sampled_Dads %>% 
  group_by(Birth, Dad_birth, Dad_capt) %>%  #ID which columns to group by
  summarise(Matches = sum(!is.na(Dad)))

#Create MasterID column to merge dataframes with matches (sm & sd) with dataframes with all comparisons (dm2, dd2)
sd$MasterID <- paste(sd$Dad_birth, sd$Dad_capt, sd$Birth, sep="-")
sm$MasterID <- paste(sm$Mom_birth, sm$Mom_capt, sm$Birth, sep="-")

all_dads <- merge(dd2, sd[,4:5], by="MasterID", all=TRUE) #Merge father dataframes
all_dads$Matches[is.na(all_dads$Matches)==TRUE] <- 0 #Put 0 in all Matches cells with NA
all_dads$No_matches <- all_dads$No_matches-all_dads$Matches #Subtract matches from no matches

#sum(all_dads$Matches)

all_moms <- merge(dm2, sm[,4:5], by="MasterID", all=TRUE) #Merge mother dataframes
all_moms$Matches[is.na(all_moms$Matches)==TRUE] <- 0 #Put 0 in all Matches cells with NA
all_moms$No_matches <- all_moms$No_matches-all_moms$Matches #Calculate number of no matches

#sum(all_moms$Matches)

#which(is.na(all_moms[,2])==TRUE) #Shouldn't be any
#which(is.na(all_dads[,2])==TRUE) #Shouldn't be any
#which(is.na(sm[,2])==TRUE) #Shouldn't be any
#which(is.na(dm[,2])==TRUE) #Shouldn't be any

###################Create (and band-aid) dataframes for CKMR model####################
Data_mom_yes <- all_moms[all_moms$Matches > 0, c(2:4,6)] 
Data_dad_yes <- all_dads[all_dads$Matches > 0, c(2:4,6)]
Data_dad_no <- all_dads[all_dads$Matches == 0, c(2:4,5)]
Data_mom_no <- all_moms[all_moms$Matches == 0, c(2:4,5)]

if(nrow(Data_dad_yes)>1 & nrow(Data_mom_yes)>1){

  #If adult was born too early to be the parent, then change the birth date so the adult is mature (and alive) when offspring is born.
  for(i in 1:length(Data_dad_yes[,1])){
    if(Data_dad_yes$Offspring_birth[i] - Data_dad_yes$Adult_birth[i] < min(m_adult_age)){
      Data_dad_yes$Adult_birth[i] <- Data_dad_yes$Offspring_birth[i] - min(m_adult_age) #Make sure adult is mature
    } else if (Data_dad_yes$Offspring_birth[i] - Data_dad_yes$Adult_birth[i] > max_age){
      Data_dad_yes$Adult_birth[i] <- Data_dad_yes$Offspring_birth[i] - max_age #Make sure adult is alive
    }
  }

  for(i in 1:length(Data_mom_yes[,1])){
    if(Data_mom_yes$Offspring_birth[i] - Data_mom_yes$Adult_birth[i] < min(f_adult_age)){
      Data_mom_yes$Adult_birth[i] <- Data_mom_yes$Offspring_birth[i] - min(f_adult_age) #Make sure adult is mature
    } else if (Data_mom_yes$Offspring_birth[i] - Data_mom_yes$Adult_birth[i] > max_age){
      Data_mom_yes$Adult_birth[i] <- Data_mom_yes$Offspring_birth[i] - max_age #Make sure adult is alive
    }
  }
  
  #########################Fit model!##############################
    CK_fit <- optimx(par=Pars,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Mother=Data_mom_no,Negatives_Father=Data_dad_no,Pairs_Mother=Data_mom_yes,Pairs_Father=Data_dad_yes,P_Mother=P_Mother,P_Father=P_Father,n_yrs=n_yrs,t_start=t_start,t_end=t_end)
    
    summary(CK_fit)
    #compute variance covariance matrix
    D=diag(length(Pars))*c(exp(CK_fit$p1[1]),exp(CK_fit$p2[1]), CK_fit$p3[1]) #derivatives of transformations
    VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
    VC = (t(D)%*%VC_trans%*%D) #delta method
    SE=sqrt(diag(VC))
    
    #cat(paste("Estimated mature female abundance: ",exp(CK_fit$p1),", SE = ",SE[1],"\n"))
    #cat(paste("Estimated mature male abundance: ",exp(CK_fit$p2),", SE = ",SE[2],"\n"))
    #cat(paste("Estimated survival: ", CK_fit$p3, ", SE = ", SE[3], "\n"))
    
    sim_results[iter, 1] <- round(exp(CK_fit$p1[1]),0) #Nf
    sim_results[iter, 2] <- round(SE[1],0) #NfSE
    sim_results[iter, 3] <- round(exp(CK_fit$p2[1]),0) #Nm
    sim_results[iter, 4] <- round(SE[2],0) #NmSE
    sim_results[iter, 5] <- Estimated_truth
    sim_results[iter, 6] <- n_samples
    save(CK_fit, file=paste0("Lemon_CKModel_nonlethal_aprioriN_",Estimated_truth,"_", iter))
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