library(tidyr)
lemon_data <- read.csv("Main_lemon_shark.csv", header=TRUE, stringsAsFactors = FALSE)
lemon_data <- lemon_data[which(lemon_data$Island=="BIM"),] #Subset for bimini
lemon_data %>% separate(Tube.Label, sep=";", into = c("Juvenile_ID", "Recapture_1", "Recapture_2", "Recapture_3", "Recapture_4")) -> lemon_data2 #Separate Tube labels into capture history
head(lemon_data2)
colnames(lemon_data2) <- c("PIT_tag","Capture_ID", "Recapture_1", "Recapture_2", "Recapture_3", "Recapture_4", "Capture_Year", "Capture_Date", "Island", "Site", "Sex", "PCL_cm", "TL_cm", "DOB", "Father", "Mother")

lemon_ref <- lemon_data2[which(lemon_data2$DOB!=""),]
lemon_ref2 <- lemon_ref[,c(2,15:16)]
colnames(lemon_ref2) <- c("Juv_ID", "Father", "Mother")

head(lemon_ref)
length(lemon_data2[,1])
length(lemon_ref[,1])

Juv_ref <- lemon_ref[,c(2,7,11,14:16)] #Add columns 15 & 16 to add parent IDs
#Separate uncertain birth years (e.g. 1993/1994) into columns
Juv_ref %>% separate(DOB, sep="/", into=c("Yr1", "Yr2")) -> Juv_ref
colnames(Juv_ref) <- c("Indiv_ID", "Capture_Year","Sex","Yr1","Yr2","Father","Mother","DOB")

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
Juv_ref$DOB = DOB

#Save Juvenile birth year in separate vector
Juv_birth <- sort(as.numeric(Juv_ref$DOB))
head(Juv_birth)
length(Juv_birth)

#Make reference table for adults
library(plyr)
Adult_ref <- Juv_ref[,6:8]
colnames(Adult_ref) <- c("Father", "Mother", "Juv_DOB")
Adult_ref$Father = sub("^$", "Unknown", Adult_ref$Father)
Adult_ref$Mother = sub("^$", "Unknown", Adult_ref$Mother)
Father_ref <- count(Adult_ref[,c(1,3)])
Mother_ref <- count(Adult_ref[,c(2,3)])
tail(Father_ref)
tail(Mother_ref)

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
Father_age <- unique(Father_ref$Father)
Fat_age <- sample(m_adult_age, replace=TRUE, prob=Prop_age_m, size=length(Father_age))
Father_age <- as.data.frame(cbind(Father_age,Fat_age), stringsAsFactors=FALSE)
colnames(Father_age) <- c("Father", "Age")

#Assign age to mothers
f_adult_age <- c(13:25)
f_mat <- c(rep(0,12), rep(1,13))
Surv_age_f <- rep(0.85,length(f_adult_age))
Prop_age_f <- rep(1, length(f_adult_age))
for(iage in 2:length(Prop_age_f)) Prop_age_f[iage]=Prop_age_f[iage-1]*Surv_age_f[iage-1] #Sets proportion of animals from each cohort surviving to iage based on Surv_age
Prop_age_f=Prop_age_f/sum(Prop_age_f)  #stable age distribution given assumed survival i.e. gives proportion of ages in a population based on assumed survival

Mother_age <- unique(Mother_ref$Mother)
Mot_age <- sample(f_adult_age, replace=TRUE, prob=Prop_age_f, size=length(Mother_age))
Mother_age <- as.data.frame(cbind(Mother_age,Mot_age), stringsAsFactors=FALSE)
colnames(Mother_age) <- c("Mother", "Age")

# Subset for animals captured in 2004-2007
t_start = 2004
t_end=2007
lems <- subset(Juv_ref, subset=Capture_Year %in% c(2004:2007), select=c(Capture_ID, Capture_Year, Father, Mother, DOB, Sex))
head(lems)
head(Mother_age)
# Assign ages to parents
#Mothers
for(i in 1:length(lems[,1])) {
  if(lems$Mother[i]!=""){
    lems$Mother_age[i] <- Mother_age[which(Mother_age$Mother==lems$Mother[i]),2]
  }
}
lems$Mother_DOB <- t_start-as.numeric(lems$Mother_age)

#Fathers
for(i in 1:length(lems[,1])) {
  if(lems$Father[i]!=""){
    lems$Father_age[i] <- Father_age[which(Father_age$Father==lems$Father[i]),2]
  }
}
lems$Father_DOB <- t_start-as.numeric(lems$Father_age)
head(lems)

#Make separate dfs for Juveniles, Mothers, and Fathers to combine
Juv_lems <- as.data.frame(cbind(lems$Capture_ID, lems$DOB, lems$Capture_Year, lems$Sex), stringsAsFactors = FALSE)
Fat_lems <- as.data.frame(cbind(lems$Father, lems$Father_DOB, lems$Capture_Year, rep("M", length(lems[,1]))), stringsAsFactors = FALSE)
Mot_lems <- as.data.frame(cbind(lems$Mother, lems$Mother_DOB, lems$Capture_Year, rep("F", length(lems[,1]))), stringsAsFactors = FALSE)
colnames(Fat_lems) = colnames(Juv_lems) = colnames(Mot_lems) = c("Indiv", "DOB", "Capture_Year", "Sex")

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
n_samples=88
n_yrs = 4
samps=c()
for(i in 1:4){
  cp_yr = c(2004:2007)[i]
  samps <- c(samps, sample(pop_df[which(pop_df$Capture_Year==cp_yr),1], size=n_samples/n_yrs, replace=FALSE))
}
length(samps)
#Want sample dataframe to have year of capture(death) year of birth, and sex
Samp_birth = c()
Samp_cap_yr = c()
Samp_sex = c()
for(i in 1:length(samps))Samp_birth[i] = pop_df[which(pop_df$Indiv==samps[i]),2]
for(i in 1:length(samps))Samp_cap_yr[i] = pop_df[which(pop_df$Indiv==samps[i]),3]
for(i in 1:length(samps))Samp_sex[i] = pop_df[which(pop_df$Indiv==samps[i]),4]
Samples <- as.data.frame(cbind(samps,Samp_birth,Samp_cap_yr, Samp_sex), stringsAsFactors = FALSE)
colnames(Samples) <- c("Indiv_ID", "Birth", "Death", "Sex")
head(Samples, 15)
Samples[,2] <- as.numeric(Samples[,2])
Samples[,3] <- as.numeric(Samples[,3])
head(Samples)
head(Juv_ref)
Samples[,6:8]=NULL

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
head(Samples)
Samples$Dad_birth=Samples$Dad_death=Samples$Mom_birth=Samples$Mom_death=NA
Samples$Mom=Samples$Dad=0

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

##Set up pairwise comparison matrix
Data = data.frame(matrix(0,nrow=n_samples*(n_samples-1)/2,ncol=5)) # massive array for pairwise comparisons
length(Data[,1])
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
head(Data)

#Assign parentage
head(Samples)
head(Data)

for(i in 1:length(Samples[,1])){
  if(Samples$Mom[i]==1){
    Data[which(Data[,1]==Samples$Mom_birth[i] & Data[,2]==Samples$Mom_death[i] & Data[,3]==Samples$Birth[i] & Data[,4]=="F"),5] = Data[which(Data[,1]==Samples$Mom_birth[i] & Data[,2]==Samples$Mom_death[i] & Data[,3]==Samples$Birth[i] & Data[,4]=="F"),5] + 1
  } 
  if(Samples$Dad[i]==1){
    Data[which(Data[,1]==Samples$Dad_birth[i] & Data[,2]==Samples$Dad_death[i] & Data[,3]==Samples$Birth[i] & Data[,4]=="M"),5]=Data[which(Data[,1]==Samples$Dad_birth[i] & Data[,2]==Samples$Dad_death[i] & Data[,3]==Samples$Birth[i] & Data[,4]=="M"),5]+1
  }
}
head(Data)

Data_dad_no = Data[which(Data[,4]=="M" & Data[,5]==0),1:4]
Data_dad_yes = Data[which(Data[,4]=="M" & Data[,5]>0),1:4]
Data_mom_no = Data[which(Data[,4]=="F" & Data[,5]==0),1:4]
Data_mom_yes = Data[which(Data[,4]=="F" & Data[,5]>0),1:4]

library(plyr) #for frequencies
Data_dad_no=count(Data_dad_no[,1:3])
Data_mom_no=count(Data_mom_no[,1:3])
Data_dad_yes=count(Data_dad_yes[,1:3])
Data_mom_yes=count(Data_mom_yes[,1:3])
