#THIS IS THE MOST RECENT SCRIPT 05/19/2020

##NEXT: 
#1) Adapt model for the input dataframes
#2) Fit model


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
# Systematically and exhaustively sampled juvenile lemon sharks from 1995-2000

##############Set up experiment parameters##############
n_ages=25 #Max age of lemon sharks
y_start <- 1995 #Start year of data
y_end <- 2000 #End year of data
total_yrs <- c(seq(y_start,y_end)) #all years represented in dataset
study_yrs <- c(1:length(total_yrs)) #above translated to study years ie beginning at 1
min_est_cohort <- min(study_yrs)+2 #First year we are estimating abundance
max_est_cohort <- max(study_yrs) #Last year we are estimating abundance
est_yrs <- min_est_cohort:max(study_yrs) #All years for which we are estimating abundance
n_yrs <- length(est_yrs) #Number of years for which we are estimating abundance

####Set up simulation parameters####
#Set Estimated_truth to set the sample size
total_abundance <- c(rep(1000, times=n_yrs)) #Can't estimate abundance for every year because HS models only estimate for younger sibling in each pairwise comparison.

N_f <- c(total_abundance/2)
N_m <- c(total_abundance-N_f)
Pars=c(log(N_f),log(N_m))

#iterations <- 100 #Set number of iterations to run in the loop

source("./time_series/models/get_P_lemon_HS_time_series_ind_het.R")
source("./likelihood_functions/lemon_neg_log_like_HS.R")

#################### LOAD AND CLEAN UP DATA ####################
lemon_data <- read.csv("./data/Main_lemon_shark.csv", header=TRUE, stringsAsFactors = FALSE)
lemon_data <- lemon_data[which(lemon_data$Island=="BIM"),] #Subset for bimini
lemon_data %>% separate(Tube.Label, sep=";", into = c("Juvenile_ID", "Recapture_1", "Recapture_2", "Recapture_3", "Recapture_4")) -> lemon_data2 #Separate Tube labels into capture history (there is one tube label per sample instance)
colnames(lemon_data2) <- c("PIT_tag","Capture_ID", "Recapture_1", "Recapture_2", "Recapture_3", "Recapture_4", "Capture_Year", "Capture_Date", "Island", "Site", "Sex", "PCL_cm", "TL_cm", "DOB", "Father", "Mother") #relabel columns
#head(lemon_data2)

lemon_ref <- lemon_data2[which(lemon_data2$DOB!=""),] #Subset for juveniles with known birth dates

#head(lemon_ref)
#nrow(lemon_data2) #How many total bimini lemon sharks in the dataset
#nrow(lemon_ref) #How many bimini lemon sharks with known birth dates
#colnames(lemon_ref)

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

#Now Juv_ref is cleaned and ready to be subsetted and sampled

#############3########### SUBSET AND SAMPLE ###########################
#head(Juv_ref)
table(factor(Juv_ref$DOB))

#Set up dataframe for population we will sample
# Want to base estimate on animals caught between 2004-2007 (the most heavily sampled years) so subset for animals captured in 2004-2007
lems <- subset(Juv_ref, subset=DOB %in% c(y_start:y_end), select=c(Indiv_ID, Capture_Year, Father, Mother, DOB, Sex))

#head(Juv_ref)
#head(lems)

#Create vector of birth year for juveniles that reflects birth year relative to study
Juv_study_DOB <- c()
for(i in 1:length(lems[,1])){
  Juv_study_DOB[i] <- study_yrs[which(total_yrs == lems$DOB[i])]
}

#Create reference dataframe for juveniles from which we will sample
Juv_lems <- as.data.frame(cbind(lems$Indiv_ID, Juv_study_DOB, lems$Capture_Year, lems$Father, lems$Mother), stringsAsFactors = FALSE) 
colnames(Juv_lems) <- c("Indiv", "Juv_study_DOB", "Capture_Year", "Father", "Mother")
head(Juv_lems)
tail(Juv_lems)

##############################Begin simulation##############################
sim_start_time <- Sys.time()
print(paste0("Simulation N", Estimated_truth, " started at ", Sys.time()))

sim_results <- matrix(0, nrow = iterations, ncol = 6)
#for(iter in 1:iterations) {

#############################Sample population##############################
head(Juv_lems)

#samps <- c()
#Make vector called samps that randomly samples the population
#for(i in 1:4){
 # cp_yr = c(2004:2007)[i]
  #samps <- c(samps, sample(Juv_lems[,1], size=ceiling(n_samples/length(samp_study_yrs)), replace=FALSE))
#}

#Create sample dataframe that has year of capture and year of birth
#Samp_birth = c() #Initialize vectors
#Samp_cap_yr = c()
#Samp_mom = c()
#Samp_dad = c()

#Simulate birth year based on study year (i.e. 26 instead of 2004)
#for(i in 1:length(samps)){
#    Samp_birth[i] <- as.numeric(Juv_lems$Juv_study_DOB[which(Juv_lems==samps[i])])
#}

#search for the Indiv ID and extract the associated capture year
#for(i in 1:length(samps))Samp_cap_yr[i] = as.numeric(Juv_lems[which(Juv_lems$Indiv==samps[i]),5])

#search for the Indiv ID and extract the associated mother
#for(i in 1:length(samps))Samp_mom[i] = Juv_lems[which(Juv_lems$Indiv==samps[i]),4]

#search for the Indiv ID and extract the associated father
#for(i in 1:length(samps))Samp_dad[i] = Juv_lems[which(Juv_lems$Indiv==samps[i]),3]

#Combine vectors into dataframe and rename columns
#Samples <- as.data.frame(cbind(samps,Samp_birth,Samp_cap_yr, Samp_dad, Samp_mom), stringsAsFactors = FALSE)

#colnames(Samples) <- c("Indiv_ID", "Birth", "Capt", "Father", "Mother")
#head(Samples)

#Change columns to numeric; changing earlier doesn't work
#Samples[,2] <- as.numeric(Samples[,2])
#Samples[,3] <- as.numeric(Samples[,3])
#Samples
#head(Juv_ref)
#tail(Samples)





#head(Samples)

##Split Juv_lems into dataframes of positive and negative comparisons
Juv_lems$Juv_study_DOB <- as.numeric(Juv_lems$Juv_study_DOB) #change DOB to type numeric
head(Juv_lems, 10)

#################### Create reference tables for individual/mean reproductive output ######################

#Summarize offspring by parent by year
grouped_by_DOB_mom <- Juv_lems %>%
  group_by(Juv_study_DOB, Mother) %>% #Group by DOB and Mother
  filter(Mother!="Unknown") %>% #Filter out unknown mothers
  summarize(pups = n()) #count offspring per mother per year

#Summarize total offspring by year
sum_yr_mom <- grouped_by_DOB_mom %>% 
  summarize(sumyr = sum(pups))

#Calculate mean offspring per year and add to table
sum_table_mom <- grouped_by_DOB_mom %>% 
  summarize(mean_pups = round(mean(pups), 0)) %>%
  left_join(grouped_by_DOB_mom) %>% 
  select(Juv_study_DOB, Mother, pups, mean_pups)

#Now we have a dataframe with the mother, her reproductive output each year, and the mean reproductive output each year
#Let's do the same for the fathers
#Summarize offspring by parent by year
grouped_by_DOB_dad <- Juv_lems %>%
  group_by(Juv_study_DOB, Father) %>% #Group by DOB and Father
  filter(Father!="Unknown") %>% #Filter out unknown Fathers
  summarize(pups = n()) #count offspring per Father per year

#Summarize total offspring by year
sum_yr <- grouped_by_DOB_dad %>% 
  summarize(sumyr = sum(pups))

#Calculate mean offspring per year and add to table
sum_table_dad <- grouped_by_DOB_dad %>% 
  summarize(mean_pups = round(mean(pups), 0)) %>%
  left_join(grouped_by_DOB_dad) %>% 
  select(Juv_study_DOB, Father, pups, mean_pups)

head(sum_table_dad)
tail(sum_table_dad)
nrow(sum_table_dad)

##################Set up pairwise comparison matrices########################
##Create dataframe of positive comparisons for MOTHERS ##
#Group by mother and split into a list of dataframes, each corresponding to the offspring of a mother
grouped_by_mom <- Juv_lems %>% 
  group_by(Mother) %>% 
  #filter(Mother!="Unknown") %>% #Filter out unknown mothers
  group_split() 

#Filter out the dataframes with no matches
group_mom2 <- Filter(function(x) nrow(x) > 1, grouped_by_mom)
View(group_mom2[[1]])

#Function for creating all possible combinations of values and preserving mother
comb_df <- function(x) {
  as.data.frame(cbind(t(combn(x$Juv_study_DOB, m=2)), x$Mother[1]))
    }

#Apply the combn function to every dataframe in the list
group_mom3 <- lapply(group_mom2, comb_df)
head(group_mom3[[1]],)

#combine list of matrices into one matrix
group_mom_df <- do.call(rbind, group_mom3)
head(group_mom_df)
tail(group_mom_df)
str(group_mom_df)

#Count number of times a comparison occurs and convert columns to numeric
group_mom_df <- plyr::count(group_mom_df)
group_mom_df[,1] <- as.numeric(group_mom_df[,1])
group_mom_df[,2] <- as.numeric(group_mom_df[,2])
group_mom_df$freq <- as.numeric(group_mom_df$freq)

#Sort so older/younger sibling in each comparison are in correct column
for(i in 1:nrow(group_mom_df)){
  if(group_mom_df[i,1] > group_mom_df[i,2]){
    group_mom_df$Old_sib_birth[i] <- group_mom_df[i,2]
    group_mom_df$Young_sib_birth[i] <- group_mom_df[i,1]
  } else {
    group_mom_df$Old_sib_birth[i] <- group_mom_df[i,1]
    group_mom_df$Young_sib_birth[i] <- group_mom_df[i,2]
  }
}
group_mom_df[,1:2] <- NULL #Remove unsorted columns

#Check that sorting worked
#nrow(group_mom_df[group_mom_df[,2] > group_mom_df[,1],])
#nrow(group_mom_df[group_mom_df[,1] > group_mom_df[,2],])
#nrow(group_mom_df[group_mom_df$Young_sib_birth > group_mom_df$Old_sib_birth,]) #This total should equal the two above

#
colnames(group_mom_df) <- c("Mother", "Mom_matches", "Old_sib_birth", "Young_sib_birth") #Name columns
group_mom_df <- group_mom_df %>% select(c(Old_sib_birth, Young_sib_birth, Mom_matches, Mother)) #Order columns
head(group_mom_df)

#Remove same cohort comparisons and those we're not interested in
mom_positives <- group_mom_df[which(group_mom_df$Old_sib_birth != group_mom_df$Young_sib_birth),]
mom_positives <-  mom_positives[which(mom_positives$Young_sib_birth >= min_est_cohort),]
head(mom_positives)
tail(mom_positives)

#####REPEAT ABOVE FOR DADS
#Create dataframe of positive comparisons for fatherS - #Group by father and split into a list of dataframes, each corresponding to the offspring of a father
grouped_by_dad <- Juv_lems %>% 
  group_by(Father) %>% 
  #filter(Father!="Unknown") %>% #Filter out unknown Fathers
  group_split() 

#Filter out the dataframes that only have one observed sibling i.e. those with no matches
group_dad2 <- Filter(function(x) nrow(x) > 1, grouped_by_dad)
#View(group_dad2[[1]])

#Function for creating all possible combinations of values and preserving Father
comb_df <- function(x) {
  as.data.frame(cbind(t(combn(x$Juv_study_DOB, m=2)), x$Father[1]))
}

#Apply the combn function to every dataframe in the list
group_dad3 <- lapply(group_dad2, comb_df)
group_dad3[[1]]

#combine list of matrices into one matrix
group_dad_df <- do.call(rbind, group_dad3)
head(group_dad_df)
tail(group_dad_df)
str(group_dad_df)

#Count number of times a comparison occurs and convert columns to numeric
group_dad_df <- plyr::count(group_dad_df)
group_dad_df[,1] <- as.numeric(group_dad_df[,1])
group_dad_df[,2] <- as.numeric(group_dad_df[,2])
group_dad_df$freq <- as.numeric(group_dad_df$freq)

#Sort so older/younger sibling in each comparison are in correct column
for(i in 1:nrow(group_dad_df)){
  if(group_dad_df[i,1] > group_dad_df[i,2]){
    group_dad_df$Old_sib_birth[i] <- group_dad_df[i,2]
    group_dad_df$Young_sib_birth[i] <- group_dad_df[i,1]
  } else {
    group_dad_df$Old_sib_birth[i] <- group_dad_df[i,1]
    group_dad_df$Young_sib_birth[i] <- group_dad_df[i,2]
  }
}
group_dad_df[,1:2] <- NULL #Remove unsorted columns

colnames(group_dad_df) <- c("Father", "dad_matches", "Old_sib_birth", "Young_sib_birth") #Name columns
group_dad_df <- group_dad_df %>% select(c(Old_sib_birth, Young_sib_birth, dad_matches, Father)) #Order columns
head(group_dad_df)

#Remove same cohort comparisons and those we're not interested in
dad_positives <- group_dad_df[which(group_dad_df$Old_sib_birth != group_dad_df$Young_sib_birth),]
dad_positives <-  dad_positives[which(dad_positives$Young_sib_birth >= min_est_cohort),]
head(dad_positives)
tail(dad_positives)

#Create matrix of all possible pairwise comparisons
all_comparisons <- t(combn(as.numeric(Juv_lems$Juv_study_DOB), m=2))
all_comparisons <- plyr::count(t(apply(all_comparisons, 1, sort)))
colnames(all_comparisons)[1:2] <- c("Old_sib_birth", "Young_sib_birth")

#Remove same cohort comparisons
all_comparisons = all_comparisons2 <- all_comparisons[which(all_comparisons$Old_sib_birth != all_comparisons$Young_sib_birth),]

#Before subtracting positive comparisons from total, need to create dataframe of total matches for moms and dads based on years of pairwise comparisons 
(mom_total_matches <- mom_positives %>% 
  group_by(Old_sib_birth, Young_sib_birth) %>% 
  summarize(matches=sum(Mom_matches)))

(dad_total_matches <- dad_positives %>% 
  group_by(Old_sib_birth, Young_sib_birth) %>% 
  summarize(matches=sum(dad_matches)))


#Subtract positive comparisons from all comparisons to make negatives dataframes
all_comparisons <- merge(mom_total_matches, all_comparisons2, by=c("Old_sib_birth", "Young_sib_birth")) %>%
  merge(dad_total_matches, by=c("Old_sib_birth", "Young_sib_birth")) %>% 
  rename(Mom_matches=matches.x, Dad_matches=matches.y)
#head(all_comparisons)

mom_negatives <- all_comparisons[,c(1:2)]
mom_negatives$Mom_negs <- all_comparisons$freq - all_comparisons$Mom_matches
  
dad_negatives <- all_comparisons[,c(1:2)]
dad_negatives$dad_negs <- all_comparisons$freq - all_comparisons$Dad_matches

#Combine positive and negative dataframes for mother and father
mom_all <- mom_positives %>% 
  left_join(mom_negatives, by = c("Old_sib_birth", "Young_sib_birth"))

dad_all <- dad_positives %>% 
  left_join(dad_negatives, by = c("Old_sib_birth", "Young_sib_birth"))

####LOOP####
#if(nrow(Data_dad_yes)>=1 & nrow(Data_mom_yes)>=1){
  
  #########################Fit model!##############################
    CK_fit <- optimx(par=Pars,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Mother=mom_negatives,Negatives_Father=dad_negatives,Pairs_Mother=mom_positives,Pairs_Father=dad_positives, P_Mother=P_Mother, P_Father=P_Father, n_yrs=n_yrs, t_start=t_start, t_end=t_end)

    summary(CK_fit)
    #compute variance covariance matrix
    D=diag(length(Pars))*c(exp(CK_fit$p1[1]),exp(CK_fit$p2[1])) #derivatives of transformations
    VC_trans = solve(attr(CK_fit, "details")["BFGS" ,"nhatend"][[1]])
    VC = (t(D)%*%VC_trans%*%D) #delta method
    SE=sqrt(diag(VC))
    exp(CK_fit[,1:16])
    
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
#}

sim_end_time <- Sys.time()
sim_end_time-sim_start_time
#Need to figure out which units for below
#print(paste0("Run time of simulation: ", round(sim_end_time - sim_start_time, 2), ""))
  print(paste0("Simulation N", Estimated_truth, " finished at ", Sys.time()))
colnames(sim_results) <- c("Nf", "NfSE", "Nm", "NmSE", "Estimated_truth", "Total_samples")
#head(sim_results, 30)

write.table(sim_results, file = paste0("Lemon_CKMR_loop_nonlethal_aprioriN_",Estimated_truth,".csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)


########################Troubleshooting##################################
#Check for immortal lemon sharks (ie years where there shouldn't be matches)
test95 <- subset(lems, subset=DOB == 1995, select=c(Indiv_ID, Capture_Year, Father, Mother, DOB, Sex)) #Subset for sharks borm in 1995
test2011 <- subset(lems, subset=DOB == 2011, select=c(Indiv_ID, Capture_Year, Father, Mother, DOB, Sex)) #subset for sharks born in 2011
test95$Mother[which(test95$Mother %in% test2011$Mother)] #see names of mothers that show up in both
length(unique(test2011$Mother)) #What is the minimum number of mothers we should find for this year
