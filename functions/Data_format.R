rm(list=ls())
library(tidyverse)

lemon_data <- read.csv("Main_lemon_shark.csv", header=TRUE, stringsAsFactors = FALSE)
lemon_data <- lemon_data[which(lemon_data$Island=="BIM"),] #Subset for bimini
lemon_data %>% separate(Tube.Label, sep=";", into = c("Juvenile_ID", "Recapture_1", "Recapture_2", "Recapture_3", "Recapture_4")) -> lemon_data2 #Separate Tube labels into capture history
colnames(lemon_data2) <- c("PIT_tag","Capture_ID", "Recapture_1", "Recapture_2", "Recapture_3", "Recapture_4", "Capture_Year", "Capture_Date", "Island", "Site", "Sex", "PCL_cm", "TL_cm", "DOB", "Father", "Mother") #relabel columns
head(lemon_data2)

lemon_ref <- lemon_data2[which(lemon_data2$DOB!=""),] #Subset for juveniles with known birth dates
lemon_ref2 <- lemon_ref[,c(2,15:16)] #Subset just juvenile and parent IDs
colnames(lemon_ref2) <- c("Juv_ID", "Father", "Mother")

head(lemon_ref)
length(lemon_data2[,1]) #How many total bimini lemon sharks in the dataset
length(lemon_ref[,1]) #How many bimini lemon sharks with known birth dates

Juv_ref <- lemon_ref[,c(2,7,11,14:16)] #Add columns 15 & 16 to add parent IDs
#Separate uncertain birth years (e.g. 1993/1994) into columns
Juv_ref %>% separate(DOB, sep="/", into=c("Yr1", "Yr2")) -> Juv_ref
colnames(Juv_ref) <- c("Indiv_ID", "Capture_Year","Sex","Yr1","Yr2","Father","Mother")
head(Juv_ref)

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
Juv_ref$DOB = DOB #Add stochastic birth date to Juv_ref

# #Save Juvenile birth year in separate vector -- why?
# Juv_birth <- sort(as.numeric(Juv_ref$DOB))
# head(Juv_birth)
# length(Juv_birth)

#Make reference table for adults
#library(plyr) -- tidyverse likely loads this already
Adult_ref <- Juv_ref[,6:8] #Subset for Father ID, Mother ID, and Juv DOB
colnames(Adult_ref) <- c("Father", "Mother", "Juv_DOB")
Adult_ref$Father = sub("^$", "Unknown", Adult_ref$Father) #replace empty parent ID cells with "Unknown"
Adult_ref$Mother = sub("^$", "Unknown", Adult_ref$Mother)
Father_ref <- plyr::count(Adult_ref[,c(1,3)]) #Count number of juveniles parents birthed each year
Mother_ref <- plyr::count(Adult_ref[,c(2,3)])
head(Father_ref)
tail(Mother_ref)

################# Simulate adult info ################
# Simulate and assign ages to adults based on life history
# Max age: two options: 
# 37 from https://www.saveourseasmagazine.com/lemon-sharks-old/
# 25 from White et. al. 2014
# Juvenile survival values from Dibattista 2007: 0.55 (mort = 0.45)
# Adult survival value from White et. al. 2014: 0.85 (mort=0.15)
# Litter size from Feldheim et. al. 2002 (avg=6, max=18)
# Age-at-maturity from Brown and Gruber 1988 (11.6 for M, 12.7 for F)
# About 77 juvenile lemon sharks inhabit the nursery at a given time (White et al 2014). So, 10sqrt(N) = 88 samples total

#Assign age to fathers
max_age = 25
m_adult_age <- c(12:25)
m_mat <- c(rep(0,11), rep(1,14))
Surv_age_m <- rep(0.85,length(m_adult_age))
Prop_age_m <- rep(1, length(m_adult_age))
for(iage in 2:length(Prop_age_m)) Prop_age_m[iage]=Prop_age_m[iage-1]*Surv_age_m[iage-1] #Sets proportion of animals from each cohort surviving to iage based on Surv_age
Prop_age_m=Prop_age_m/sum(Prop_age_m)  #stable age distribution given assumed survival i.e. gives proportion of ages in a population based on assumed survival
Prop_age_m
Uniq_Father <- unique(Father_ref$Father) #Save unique fathers
Father_age <- sample(m_adult_age, replace=TRUE, prob=Prop_age_m, size=length(Uniq_Father))
Father_age <- as.data.frame(cbind(Uniq_Father,Father_age), stringsAsFactors=FALSE)
colnames(Father_age) <- c("Father", "Age")

#Assign age to mothers
f_adult_age <- c(13:25)
f_mat <- c(rep(0,12), rep(1,13))
Surv_age_f <- rep(0.85,length(f_adult_age))
Prop_age_f <- rep(1, length(f_adult_age))
for(iage in 2:length(Prop_age_f)) Prop_age_f[iage]=Prop_age_f[iage-1]*Surv_age_f[iage-1] #Sets proportion of animals from each cohort surviving to iage based on Surv_age
Prop_age_f=Prop_age_f/sum(Prop_age_f)  #stable age distribution given assumed survival i.e. gives proportion of ages in a population based on assumed survival

Uniq_Mother <- unique(Mother_ref$Mother)
Mother_age <- sample(f_adult_age, replace=TRUE, prob=Prop_age_f, size=length(Uniq_Mother))
Mother_age <- as.data.frame(cbind(Uniq_Mother,Mother_age), stringsAsFactors=FALSE)
colnames(Mother_age) <- c("Mother", "Age")
head(Mother_age)

# Want to base estimate on animals caught between 2004-2007 (the most heavily sampled years) 
#Subset for animals captured in 2004-2007
t_start = 2004
t_end=2007
lems <- subset(Juv_ref, subset=Capture_Year %in% c(2004:2007), select=c(Indiv_ID, Capture_Year, Father, Mother, DOB, Sex))
head(lems)

# Assign ages to parents in subset of data
#Mothers
for(i in 1:length(lems[,1])) {
  if(lems$Mother[i]!=""){
    lems$Mother_age[i] <- Mother_age[which(Mother_age$Mother==lems$Mother[i]),2]
  }
}

#Change age into DOB (will be used in CKMR model)
lems$Mother_DOB <- t_start-as.numeric(lems$Mother_age)
min(lems$Mother_DOB)

#Fathers
for(i in 1:length(lems[,1])) {
  if(lems$Father[i]!=""){
    lems$Father_age[i] <- Father_age[which(Father_age$Father==lems$Father[i]),2]
  }
}
lems$Father_DOB <- t_start-as.numeric(lems$Father_age)
head(lems)

#Make separate dfs for Juveniles, Mothers, and Fathers to combine
Juv_lems <- as.data.frame(cbind(lems$Indiv_ID, lems$DOB, lems$Capture_Year, lems$Sex), stringsAsFactors = FALSE)
Fat_lems <- as.data.frame(cbind(lems$Father, lems$Father_DOB, lems$Capture_Year, rep("M", length(lems[,1]))), stringsAsFactors = FALSE)
Mot_lems <- as.data.frame(cbind(lems$Mother, lems$Mother_DOB, lems$Capture_Year, rep("F", length(lems[,1]))), stringsAsFactors = FALSE)
colnames(Fat_lems) = colnames(Juv_lems) = colnames(Mot_lems) = c("Indiv", "DOB", "Capture_Year", "Sex")
head(Fat_lems)

# Remove duplicate parents
Mot_lems$uniqID <- paste0(Mot_lems$Indiv,Mot_lems$DOB)
length(Mot_lems[,1])
Mot_lems <- Mot_lems[!duplicated(Mot_lems$uniqID),]
Mot_lems$uniqID <- NULL
length(Mot_lems[,1])

Fat_lems$uniqID <- paste0(Fat_lems$Indiv,Fat_lems$DOB)
length(Fat_lems[,1])
Fat_lems <- Fat_lems[!duplicated(Fat_lems$uniqID),]
Fat_lems$uniqID <- NULL
length(Fat_lems[,1])

#Combine data frames
pop_df <- rbind(Juv_lems, Fat_lems, Mot_lems)
head(pop_df)
tail(pop_df)

#Want 88 samples total, spread over 4 years, so 22 samples per year
#White estimates 77 juveniles in the population
Estimated_truth <- 77
n_samples <- round(10*sqrt(Estimated_truth), 0)
samp_yrs = 4
samps=c()

for(i in 1:4){
  cp_yr = c(2004:2007)[i]
  samps <- c(samps, sample(pop_df[which(pop_df$Capture_Year==cp_yr),1], size=n_samples/samp_yrs, replace=FALSE))
}
length(samps)
head(samps)
#Want sample dataframe to have year of capture(death) year of birth, and sex
Samp_birth = c()
Samp_cap_yr = c()
Samp_sex = c()

for(i in 1:length(samps))Samp_birth[i] = pop_df[which(pop_df$Indiv==samps[i]),2]
for(i in 1:length(samps))Samp_cap_yr[i] = pop_df[which(pop_df$Indiv==samps[i]),3]
for(i in 1:length(samps))Samp_sex[i] = pop_df[which(pop_df$Indiv==samps[i]),4]
Samples <- as.data.frame(cbind(samps,Samp_birth,Samp_cap_yr, Samp_sex), stringsAsFactors = FALSE)
colnames(Samples) <- c("Indiv_ID", "Birth", "Death", "Sex")
head(Samples)

Samples[,2] <- as.numeric(Samples[,2])
Samples[,3] <- as.numeric(Samples[,3])
head(Samples)
head(Juv_ref)

#Assign father and mother for each sampled individual
for(i in 1:length(Samples[,1])){
  for(j in 1:length(Juv_ref[,1])){
    if(Samples[i,1]==Juv_ref[j,1]){
      Samples$Father[i]=Juv_ref[j,6]
      Samples$Mother[i]=Juv_ref[j,7]
    } #else if(Samples[i,1]==Juv_ref[j,6] | Samples[i,1]==Juv_ref[j,7]){
     #Samples$Offspring[i]=Juv_ref[j,1] 
    #}
  }
}

Samples$Dad_birth=Samples$Dad_death=Samples$Mom_birth=Samples$Mom_death=NA
Samples$Mom=Samples$Dad=0 #To identify which sampels have IDed parents

head(Samples)
#If the sampled animal also is a parent to another sampled animal, assign it a birth and death date
for(i in 1:length(Samples[,1])){
  if(Samples$Indiv_ID[i] %in% Samples$Father){
    Samples$Dad_birth[which(Samples$Father==Samples$Indiv_ID[i])]=Samples$Birth[i]
    Samples$Dad_death[which(Samples$Father==Samples$Indiv_ID[i])]=Samples$Death[i]
  }
  if(Samples$Indiv_ID[i] %in% Samples$Mother){
    Samples$Mom_birth[which(Samples$Mother==Samples$Indiv_ID[i])]=Samples$Birth[i]
    Samples$Mom_death[which(Samples$Mother==Samples$Indiv_ID[i])]=Samples$Death[i]
  }
}
head(Samples)

for(i in 1:length(Samples[,1])){
  if(is.na(Samples$Dad_birth[i])==FALSE) Samples$Dad[i]=1
  if(is.na(Samples$Mom_birth[i])==FALSE) Samples$Mom[i]=1
}

yr_1=1979 #Set first year of study (i.e. earliest year from which we might have samples)

Samples$Birth=Samples$Birth-yr_1
Samples$Death=Samples$Death-yr_1
Samples$Mom_birth=Samples$Mom_birth-yr_1
Samples$Mom_death=Samples$Mom_death-yr_1
Samples$Dad_birth=Samples$Dad_birth-yr_1
Samples$Dad_death=Samples$Dad_death-yr_1
head(Samples)

##Set up pairwise comparison matrix
#We don't want to compare individuals with themselves (hence -1), and we don't want to do each comparison twice (divide by 2)
Data <- data.frame(matrix(0,nrow=n_samples*(n_samples-1)/2,ncol=5)) # massive array for pairwise comparisons

#columns are: adult birth; adult death year; young birth; adult sex; probability of POP
counter=1

#brute force - could probably be improved w/ expand.grid or something. Fills first four columns of Data
for(iind1 in 1:(n_samples-1)){ #iind1 starts at 1 and counts to n_samples-1
  for(iind2 in (iind1+1):n_samples){#iind2 starts at 2 and counts to n_samples
    if(Samples[iind1,2]>Samples[iind2,2]){ #compare birth years of each combination of samples. If birth year is greater for iind1 i.e. if iind2 was born before iind1. Lethal sampling, so iind2 starts at the number after iind1, and there is no need to go back before.
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

colnames(Data) <- c("Adult_birth","Adult_death","Offspring_birth","Adult_sex", "Matches")
length(Data[,1])
head(Data)
head(Samples)

############### Assign parentage ################
head(Samples)
tail(Samples)
head(Data)

Sampled_Moms <- Samples[Samples$Mom == 1,]
Sampled_Dads <- Samples[Samples$Dad == 1,]
#SM2 <- aggregate(Sampled_Moms, by=list("Birth", "Mom_death", "Mom_birth"), FUN = "sum")

Data_mom <- Data[Data$Adult_sex == "F",]
Data_dad <- Data[Data$Adult_sex == "M",]

Data_mom$MasterID <- paste(Data_mom$Adult_birth, Data_mom$Adult_death, Data_mom$Offspring_birth, sep="-")
Data_dad$MasterID <- paste(Data_dad$Adult_birth, Data_dad$Adult_death, Data_dad$Offspring_birth, sep="-")

SM2 <- Sampled_Moms %>% 
  #complete(f_year, month, age) %>%  #Expands each years data to include all months and ages
  group_by(Birth, Mom_birth, Mom_death) %>%  #ID which columns to group by
  summarise(Matches = sum(!is.na(Mom)))
SM2

SD2 <- Sampled_Dads %>% 
  #complete(f_year, month, age) %>%  #Expands each years data to include all months and ages
  group_by(Birth, Dad_birth, Dad_death) %>%  #ID which columns to group by
  summarise(Matches = sum(!is.na(Dad)))
SD2$MasterID <- paste(SD2$Dad_birth, SD2$Dad_death, SD2$Birth, sep="-")
SM2$MasterID <- paste(SM2$Mom_birth, SM2$Mom_death, SM2$Birth, sep="-")

dm3 <- plyr::count(Data_mom[,c(1:3,6)]) #Note to self: change data_mom dataframes to dm
dd3 <- plyr::count(Data_dad[,c(1:3,6)])
colnames(dm3)[5] = colnames(dd3)[5] <- "No_matches"

all_dads <- merge(dd3, SD2[,4:5], by="MasterID", all=TRUE)
all_dads$Matches[is.na(all_dads$Matches)==TRUE] <- 0
all_dads$No_matches <- all_dads$No_matches-all_dads$Matches

sum(all_dads$Matches, na.rm=TRUE)
head(all_dads)

all_moms <- merge(dm3, SM2[,4:5], by="MasterID", all=TRUE)
all_moms$Matches[is.na(all_moms$Matches)==TRUE] <- 0
all_moms$No_matches <- all_moms$No_matches-all_moms$Matches

sum(all_moms$Matches, na.rm=TRUE)
head(all_moms)

Data_mom_yes <- all_moms[all_moms$Matches > 0, c(2:4,6)]
Data_dad_yes <- all_dads[all_dads$Matches > 0, c(2:4,6)]

Data_dad_no <- all_dads[all_dads$Matches == 0, c(2:4,5)]
Data_mom_no <- all_moms[all_moms$Matches == 0, c(2:4,5)]


#head(Data_mom)
#Data_mom2 <- Data_mom[!duplicated(Data_mom),]

#Data_dad2 <- Data_dad[!duplicated(Data_dad),]
# nrow(Data_dad2)
# 
# Data_mom3 <- merge(Data_mom2[,-5],SM2[,4:5], by = "MasterID", all = TRUE)
# Data_mom3$Matches[is.na(Data_mom3$Matches)==TRUE] <- 0
# table(Data_mom3$Matches)
# 
# Data_dad3 <- merge(Data_dad2[,-5],SD2[,4:5], by = "MasterID", all = TRUE)
# Data_dad3$Matches[is.na(Data_dad3$Matches)==TRUE] <- 0
# table(Data_dad3$Matches)



####Old code
# for(i in 1:length(Samples[,1])){
#   if(Samples$Mom[i]==1){
#     Data[which(Data[,1]==Samples$Mom_birth[i] & Data[,2]==Samples$Mom_death[i] & Data[,3]==Samples$Birth[i] & Data[,4]=="F"),5] = Data[which(Data[,1]==Samples$Mom_birth[i] & Data[,2]==Samples$Mom_death[i] & Data[,3]==Samples$Birth[i] & Data[,4]=="F"),5] + 1
#   } 
#   if(Samples$Dad[i]==1){
#     Data[which(Data[,1]==Samples$Dad_birth[i] & Data[,2]==Samples$Dad_death[i] & Data[,3]==Samples$Birth[i] & Data[,4]=="M"),5]=Data[which(Data[,1]==Samples$Dad_birth[i] & Data[,2]==Samples$Dad_death[i] & Data[,3]==Samples$Birth[i] & Data[,4]=="M"),5]+1
#   }
# }
# head(Data)

# Data_dad_no = Data[which(Data[,4]=="M" & Data[,5]==0),1:4]
# Data_dad_yes = Data[which(Data[,4]=="M" & Data[,5]>0),1:4]
# Data_mom_no = Data[which(Data[,4]=="F" & Data[,5]==0),1:4]
# Data_mom_yes = Data[which(Data[,4]=="F" & Data[,5]>0),1:4]
# Data_dad_no=count(Data_dad_no[,1:3])
# Data_mom_no=count(Data_mom_no[,1:3])
# Data_dad_yes=count(Data_dad_yes[,1:3])
# Data_mom_yes=count(Data_mom_yes[,1:3])
# sum(Data_mom_yes$freq)

#In case offspring are born after parents die, band-aid fix to allow the model to run
Data_mom_yes$Adult_death <- ifelse(Data_mom_yes$Offspring_birth > Data_mom_yes$Adult_death, Data_mom_yes$Offspring_birth, Data_mom_yes$Adult_death)

Data_dad_yes$Adult_death <- ifelse(Data_dad_yes$Offspring_birth > Data_dad_yes$Adult_death, Data_dad_yes$Offspring_birth, Data_dad_yes$Adult_death)