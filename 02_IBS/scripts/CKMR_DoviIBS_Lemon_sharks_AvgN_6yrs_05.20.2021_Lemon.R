rm(list=ls())

#Load packages
library(optimx)
library(tidyverse)
library(popbio)
library(mpmtools)

set.seed(444777)
####---------Simulation parameters---------####
init.pop.size <- 3000 # Initial population size
init.prop.female <- .5 # proportion of the initial population size that is female
repro.age <- 12 #set age of reproductive maturity
max.age = maxAge <- 30 #set the maximum age allowed in the simulation
mating.periodicity <- 2 #number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
num.mates <- c(1:3) # vector of potential number of mates per mating
avg.num.offspring <- 3 # set the average number of offspring per mating (from a poisson distribution)
birth.sex.ratio <- c(.5,.5) # The probability that each baby is F:M - has to add up to 1
Adult.survival <- .9 #Adult survival
juvenile.survival <- .9 #juvenile survival
YOY.survival <- .8 # young of year survival
burn.in <- 40 # number of years to use as simulation burn in period
Num.years <- 50 # The number of years to run in the simulation beyond the burn in
n_yrs = t_end <- burn.in + Num.years

#Decide which years to subtract from n_yrs for sampling and include in the vector below (currently 3:0, which means sample the last four years of the simulation)
sample.years <- c(n_yrs - c(5:0))
sample.size <- 30 #sample size per year


iterations <- 100 #Number of iterations to loop over

#--------------Start simulation loop--------------
for(samps in 1:3){
  #set.seed(47)
  results <- NULL #initialize results array
  sample.size <- c(40, 50, 60)[samps] #To loop over different sample sizes, draw a different number of samples each time
  
  for(iter in 1:iterations) {
  
####---------Set up initial population-----####
init.pop <- data.frame() # create a blank data frame which will become the initial popualation

for(i in 1:init.pop.size){ # Loop that creates the below data for each individual in a population the size of "init.pop.size"
  indv.name <- paste(sample(letters, size = 20, replace = T), collapse="") # generate a random name that is 20 letters long
  age.x <- repro.age # every individual in the initial population is the age at first maturity
  mother.x <- "xxxxx" # The inividuals in the initial population do not have known mothers
  father.x <- "xxxxx" # The inividuals in the initial population do not have known fathers
  birth.year <- -1 # this is a place holder for individuas born within the simulation
  sex <- sample(c('F','M'), size = 1, prob = c(init.prop.female, 1-init.prop.female)) #andomly draw sex based on the prportions set in the parameter section
  repro.cycle <- sample(1:mating.periodicity, size = 1) # randomly assign whenther this individual mother would breed in the even or odd years
  init.vec <- cbind.data.frame(indv.name, birth.year, age.x, mother.x, father.x, sex, repro.cycle) # create a row with the attributes for each individual - using cbind.data.frame allows them to keep some numeric and others characters
  init.pop <- rbind(init.pop, init.vec) # place the individual into the init.pop data frame
}


####----------Breeding----------####

repro.cycle.vec <- rep(1:mating.periodicity, times = 100) # Generate a vector which will be used to determine if it is an even or odd breeding year (or a 1/3 breeding year)



####---------for year 0 breeding-------####
mothers <- which(init.pop$sex=='F' & init.pop$age.x>=repro.age & init.pop$repro.cycle == repro.cycle.vec[1]) #determine which females are available to breed in this year
fathers <- which(init.pop$sex=='M' & init.pop$age.x>=repro.age) # determine which fathers are available to breed in this year

YOY.df <- data.frame() # create an empty data frame to populate with the YOY born in this year
for(j in 1:length(mothers)){ # Loop through all of the available mothers
  num.mates.x <- sample(num.mates, size=1) # determine the number of mates for a given mother this year
  inter.df <- data.frame() # Create a data frame for the offspring born from each mother
  for(h in 1:num.mates.x){ # Loop through each sire that this mother mates with
  num.offspring <- rbinom(n = 1, size = 1, prob = .9) # generate the number of offspring born from this mother/sire pairing
  indv.name <- NA #create a place holder for the random name given to each offspring
  age.x <- 0 # assign a 0 age to each offspring born this year
  mother.x <- init.pop[mothers[j],1] # record the mothers name
  father.x <- init.pop[sample(fathers, size = 1),1]# record the sires name
  birth.year <- 0 # note that these pups were born in the 0 year of the simulation
  sex <- NA # create a place holder for the sexes
  repro.cycle <- sample(1:mating.periodicity, size = 1) # assing the newborn to an eventual breeding cycle group
  inter.df2 <- cbind.data.frame(indv.name, birth.year, age.x, mother.x, father.x, sex, repro.cycle) # Create a row with the attributes of each YOY (at this point only one from each mating pair)
  inter.df2 <- inter.df2[rep(seq_len(nrow(inter.df2)), num.offspring), ] # Create the random number of offspring from each pairing that we set above
  inter.df2$sex <- sample(c('F','M'), size = nrow(inter.df2), prob = birth.sex.ratio, replace = T) #Assign biological sex to each new born based on the sex ration set in the parameter section
  baby.names <- vector() # Create a blank vector for random baby names
  for(w in 1:nrow(inter.df2)){ # create a random name for each newborn from this mating pair
  name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
  baby.names <- c(baby.names,name1)  
  }
  if(nrow(inter.df2)==0){next} #if there were no offspring from this mating pair, skips to the next mating pair (if you dont include this you will get an error in the loop)
  inter.df2$indv.name <- baby.names #add the baby names to the data frame
  inter.df <- rbind(inter.df, inter.df2) #add the new borns from this mating pair to the other from the same mother
  }
  YOY.df <- rbind(YOY.df,inter.df) #add all of the offspring from each mother to a data frame of all the YOY for the year
}

year.end.pop.0 <- NULL

year.end.pop.0 <- rbind(init.pop, YOY.df) #Combine the YOY data with the other individuals present this year

#Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
year.end.pop.0$Survival <- ifelse(year.end.pop.0$age.x==0, sample(c("S","M"), size=length(which(year.end.pop.0$age.x==0)), prob=c(YOY.survival, 1-YOY.survival), replace=T),
                                  ifelse(year.end.pop.0$age.x<repro.age, sample(c("S","M"), size=length(which(year.end.pop.0$age.x<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                        ifelse(year.end.pop.0$age.x<max.age, sample(c("S","M"), size=length(which(year.end.pop.0$age.x<max.age)), prob=c(Adult.survival, 1-Adult.survival),replace=T), "M")))


loopy.pop <- year.end.pop.0

####-------Loop for all other breeding--------####
pop.size <- data.frame()

loopy.list <- list() #Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year

for(v in 1:(burn.in + Num.years)){ #loop through all of the years in the simulation - the burn in and the years that matter

data1 <- loopy.pop[loopy.pop$Survival =="S", -8] #Bring in the data from the previous iteration, but only include those that survive
data1$age.x <- data1$age.x+1 # increase each individuals age by one for the new year - happy birthday survivors!

mothers <- which(data1$sex=='F' & data1$age.x>=repro.age & data1$repro.cycle == repro.cycle.vec[v+1])  #determine which females are available to breed in this year; this is an index
fathers <- which(data1$sex=='M' & data1$age.x>=repro.age)  #determine which males are available to breed in this year

YOY.df <- data.frame() # create an empty data frame to populate with the YOY born in this year

for(j in 1:length(mothers)){ # Loop through all of the available mothers
  num.mates.x <- sample(num.mates, size=1) # determine the number of mates for a given mother this year
  inter.df <- data.frame() # Create a data frame for the offspring born from each mother
  
  for(h in 1:num.mates.x){ # Loop through each sire with whom this mother mates 
    num.offspring <- rbinom(n = 1, size = 1, prob = .9) # generate the number of offspring born from this mother/sire pairing
    indv.name <- NA #create a place holder for the random name given to each offspring
    age.x <- 0 # assign a 0 age to each offspring born this year
    mother.x <- data1[mothers[j],1] # record the mothers name
    father.x <- data1[sample(fathers, size = 1),1] # record the fathers name
    birth.year <- v # assign the birth year in the simulation
    sex <- NA # create a place holder for the biological sex of each offspring
    repro.cycle <- sample(1:mating.periodicity, size = 1) # assign the newborn to an eventual breeding cycle group
    inter.df2 <- cbind.data.frame(indv.name, birth.year, age.x, mother.x, father.x, sex, repro.cycle) # Create a row with the attributes of each YOY (at this point only one from each mating pair)
    inter.df2 <- inter.df2[rep(seq_len(nrow(inter.df2)), num.offspring), ] # Create the random number of offspring from each pairing that we set above
    inter.df2$sex <- sample(c('F','M'), size = nrow(inter.df2), prob = birth.sex.ratio, replace = T) #Assign biological sex to each new born based on the sex ratio set in the parameter section
    baby.names <- vector() # Create a blank vector for random baby names
    
    for(w in 1:nrow(inter.df2)){ # create a random name for each newborn from this mating pair
      name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
      baby.names <- c(baby.names,name1)  
    }
    if(nrow(inter.df2)==0){next} #if there were no offspring from this mating pair, skips to the next mating pair (if you dont include this you will get an error in the loop)
    inter.df2$indv.name <- baby.names # add the baby names to the data frame
    inter.df <- rbind(inter.df, inter.df2) #add the new borns from this mating pair to the other from the same mother
  }
  YOY.df <- rbind(YOY.df,inter.df) #add all of the offspring from each mother to a data frame of all the YOY for the year
}


loopy.pop <- rbind(data1, YOY.df) #Combine the YOY data with the other individuals present this year. YOY go at the bottom.

#Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
loopy.pop$Survival <- ifelse(loopy.pop$age.x==0, sample(c("S","M"), size=length(which(loopy.pop$age.x==0)), prob=c(YOY.survival, 1-YOY.survival), replace=T),
                                  ifelse(loopy.pop$age.x<repro.age, sample(c("S","M"), size=length(which(loopy.pop$age.x<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                         ifelse(loopy.pop$age.x<max.age, sample(c("S","M"), size=length(which(loopy.pop$age.x<max.age)), prob=c(Adult.survival, 1-Adult.survival),replace=T), "M")))

#assign(paste("year.end.pop.", v, sep=""),loopy.pop) # save the current year's population data as an object

loopy.list[[v]] <- loopy.pop #Save the current year's population data as a list element, where the index corresponds to the year

print(paste("year", v, "N= ", nrow(loopy.pop) , sep=" ")) # print the simulation year and the population size in the R console so they can be oberved
pop.size.vec <- cbind.data.frame(year=v, population_size=nrow(loopy.pop), Male.adult.pop = nrow(loopy.pop[loopy.pop$sex == "M" & loopy.pop$age.x >= repro.age,]), Female.adult.pop = nrow(loopy.pop[loopy.pop$sex == "F" & loopy.pop$age.x >= repro.age,]))
pop.size <- rbind(pop.size, pop.size.vec)
}

#Label the list elements with the year
names(loopy.list) <- paste0("year.end.pop.", seq(1:(burn.in + Num.years)))


####---------quick analysis of population growth (lambda)-----####
lambda <- vector()
for(l in 2:nrow(pop.size)){
  lambda.1 <- pop.size$population_size[l]/pop.size$population_size[l-1]
  lambda <- c(lambda,lambda.1)
}
lambda <- c(NA, lambda)
pop.size$Lambda <- lambda
mean(pop.size$Lambda, na.rm=T) #Mean Lambda

hist(lambda[30:90])

####---------Checking population parameters-------####
nrow(YOY.df)/length(mothers) #Average fecundity for last year; remeber they're skipped breeding

nrow(loopy.pop[loopy.pop$Survival=='S' & loopy.pop$age.x>=repro.age, ])/nrow(loopy.pop[loopy.pop$age.x>=repro.age,]) #survival of adults
nrow(loopy.pop[loopy.pop$Survival=='S' & loopy.pop$age.x>0 & loopy.pop$age.x<repro.age, ])/nrow(loopy.pop[loopy.pop$age.x>0 & loopy.pop$age.x<repro.age,]) #survival of juveniles
nrow(loopy.pop[loopy.pop$Survival=='S' & loopy.pop$age.x==0, ])/nrow(loopy.pop[loopy.pop$age.x==0,]) #survival of YOY

survival_vec <- NULL
for(i in 1:length(loopy.list)){
  survival_vec[i] <- nrow(loopy.list[[i]][loopy.list[[i]]$Survival=='S' & loopy.list[[i]]$age.x>=repro.age, ])/nrow(loopy.list[[i]][loopy.list[[i]]$age.x>=repro.age, ])
}

nrow(loopy.pop[loopy.pop$Survival=='M' & loopy.pop$age.x>=12, ])+
  nrow(loopy.pop[loopy.pop$Survival=='M' & loopy.pop$age.x>0 & loopy.pop$age.x<12, ])+
  nrow(loopy.pop[loopy.pop$Survival=='M' & loopy.pop$age.x==0, ])/nrow(loopy.pop[loopy.pop$age.x==0,])

nrow(YOY.df)
#The number of rows in a dataframe are the number of individuals in the population

#The initial adult population size
init.adult_pop.size <- loopy.list[[1]] %>% filter(age.x >= repro.age) %>% 
  nrow()

#Set initial parameters for model based on initial abundance of males and females
f.init = m.init <- init.adult_pop.size/2
Pars <- c(log(f.init), log(m.init)) #Pars1 is for the sex-specific model
Pars2 <- log(init.adult_pop.size) #Pars2 is for the sex-aggregated model


####------------------------Collect samples---------------------####

#Initialize dataframes
sample.df_all.info <- NULL
sample.df_temp <- NULL

#Sample population during last six years and make dataframe of samples with all metadata
for(i in sample.years){
  sample.df_temp <- loopy.list[[i]] %>% dplyr::slice_sample(n = sample.size)
  sample.df_all.info <- rbind(sample.df_all.info, sample.df_temp)
}


#Keep just one instance of each individual (to avoid self-recapture) and sort by birth year so when we make the pairwise comparison matrix, Ind_1 is always older than Ind_2 (bc it comes first)
sample.df_all.info <- sample.df_all.info %>% distinct(indv.name, .keep_all = TRUE) %>% 
  dplyr::arrange(birth.year, desc()) 


####-------------Construct pairwise comparison matrix--------------####

#Create dataframe of pairwise comparisons with just individual IDs
pairwise.df <- data.frame(t(combn(sample.df_all.info$indv.name, m=2)))
colnames(pairwise.df) <- c("Ind_1", "indv.name") #Rename columns so they can easily be joined

head(pairwise.df)

#Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
Ind1_birthyears <- sample.df_all.info %>%
  select(indv.name, birth.year, age.x, mother.x, father.x) %>% #select relevant columns only
  dplyr::rename("Ind_1" = indv.name, "Ind_1_birth" = birth.year, "Ind_1_age" = age.x, "Ind_1_mom" = mother.x, "Ind_1_dad" = father.x) #Rename columns for join and also so younger sib birth year and parents are distinguishable from older sib data when joined below.

#Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
#This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.
pairwise.df_all.info <- pairwise.df %>% left_join(Ind1_birthyears, by = "Ind_1") %>% 
  left_join(sample.df_all.info, by = "indv.name") %>% 
  dplyr::rename("Ind_2" = indv.name, "Ind_2_birth" = birth.year, "Ind_2_age" = age.x, "Ind_2_mom" = mother.x, "Ind_2_dad" = father.x) %>% 
  select(Ind_1, Ind_1_birth, Ind_1_age, Ind_2, Ind_2_birth, Ind_2_age, Ind_1_mom, Ind_2_mom, Ind_1_dad, Ind_2_dad) %>% 
  filter(Ind_1_birth != Ind_2_birth)

#Check that there are no repeat comparisons -- compare number of distinct comparisons to total number of comparisons.
#Should return TRUE
pairwise.df_all.info %>% distinct(Ind_1, Ind_2) %>% nrow() == nrow(pairwise.df_all.info)


#Extract positive HS comparisons and exclude make sure there are no full sibs
positives <- pairwise.df_all.info %>% filter(Ind_1_mom == Ind_2_mom | Ind_1_dad == Ind_2_dad)

#Remove full sibs -- adjust JDS
positives %>% filter(Ind_1_mom == Ind_2_mom & Ind_1_dad == Ind_2_dad)

#Second filter to check for self-recaptures
self <- positives %>% filter(Ind_1 == Ind_2)

nrow(self) #They should have been filtered earlier so should be 0


####----------------Split dataframes into final form for model----------####
#Sex-specific
mom_positives <- positives %>% filter(Ind_1_mom == Ind_2_mom) %>% 
select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

dad_positives <- positives %>% filter(Ind_1_dad == Ind_2_dad) %>% 
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

#Sex-aggregated
parent_positives <- positives %>%
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()


#Make dataframes for negative comparisons
#Sex-specific
mom_negatives <- pairwise.df_all.info %>% filter(Ind_1_mom != Ind_2_mom & Ind_1_birth != Ind_2_birth) %>% #filter for same cohort is repetitive
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

dad_negatives <- pairwise.df_all.info %>% filter(Ind_1_dad != Ind_2_dad & Ind_1_birth != Ind_2_birth) %>% 
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

#Sex-aggregated
parent_negatives <- pairwise.df_all.info %>% filter(Ind_1_dad != Ind_2_dad & Ind_1_mom != Ind_2_mom & Ind_1_birth != Ind_2_birth) %>% 
  select(Ind_1_birth, Ind_2_birth) %>% 
  plyr::count()

#-------------Kinship probabilities - Half-sib-------------------
min_cohort <- n_yrs-max.age #set first year for calculating mean (arbitrary)

m_adult_age = f_adult_age <- c(repro.age:max.age) #Set ages at which males and females are mature. Called by kinship probability function.

pop_growth_all_mean <- mean(pop.size$Lambda[1:nrow(pop.size)], na.rm=T)

###Used in model###
surv <- Adult.survival #Set value
#surv <- mean_adult_surv #Observed value

lam <- mean(pop.size$Lambda[min_cohort:n_yrs], na.rm=T) #Mean Lambda over years of estimation


####Sex-specific####

P_Mother = P_Father = array(NA,dim=c(n_yrs,n_yrs)) #creates two empty arrays, one for mother and one for father.  Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)

get_P_lemon <- function(Pars,P_Mother,P_Father,t_start,t_end){
  N_F=exp(Pars[1]) #number of mature females (assume time constant)
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      if((ys_birth - os_birth) <= ((maxAge) - repro.age)){
        P_Mother[os_birth, ys_birth] <- (surv^(ys_birth - os_birth))/N_F
      } else P_Mother[os_birth, ys_birth] <- 0
    }
  }
  N_M=exp(Pars[2]) #number of mature males (time constant) - (total reproductive output from males)
  for(os_birth in 1:(n_yrs-1)){  #> = after, < = before
    for(ys_birth in (os_birth+1):n_yrs){
      if((ys_birth - os_birth) <= ((maxAge) - repro.age)){
        P_Father[os_birth,ys_birth] <- (surv^(ys_birth - os_birth))/N_M
      } else P_Father[os_birth,ys_birth] <- 0
    }
  }
  return(list(P_Mother=P_Mother, P_Father=P_Father)) #return makes sure this is moved out of the loop into the environment
}


####Sex-aggregated####
#Set up empty array that will be filled with function below
P_Parent = array(NA, dim=c(n_yrs,n_yrs)) #Dimensions are older sib birth year and younger sib birth year (all of which are specified by n_yrs)

#CKMR model: populate array with kinship probabilities based on birth years
get_P_lemon_TotalA <- function(Pars2,P_Parent,t_start,t_end){
  N_A=exp(Pars2[1]) #number of mature adults
  
  for(os_birth in 1:(n_yrs-1)){  #Loop over possible ages of older sibs
    for(ys_birth in max(os_birth+1):n_yrs){ #Loop over possible ages of younger sibs
      if((ys_birth - os_birth) <= ((maxAge) - repro.age)){ #If the adult could have been mature in the birth year of both individual, then fill the array with the appropriate probability of kinship
        
        #Probability of kinship based on birth year
        #See Hillary et al (2018) equation (3)
        P_Parent[os_birth, ys_birth] <- (4/N_A)*(surv^(ys_birth - os_birth))
      } else P_Parent[os_birth, ys_birth] <- 0 #If it's not possible, set kinship probability to 0
    }
  }
  return(list(P_Parent=P_Parent)) #return makes sure this is moved out of the loop into the environment
}

#------------------Likelihood functions--------------------------
#####Sex-specific#####
lemon_neg_log_lik <- function(Pars,Negatives_Mother,Negatives_Father,Pairs_Mother,Pairs_Father,P_Mother,P_Father,t_start,t_end){
  
  P=get_P_lemon(Pars=Pars,P_Mother=P_Mother,P_Father=P_Father,t_start=t_start,t_end=t_end)
  
  loglik=0
  
  #likelihood contributions for all negative comparisons
  for(irow in 1:nrow(Negatives_Mother)){
    loglik=loglik+Negatives_Mother[irow,3]*log(1-P$P_Mother[Negatives_Mother[irow,1],Negatives_Mother[irow,2]])
  } 
  for(irow in 1:nrow(Negatives_Father)){
    loglik=loglik+Negatives_Father[irow,3]*log(1-P$P_Father[Negatives_Father[irow,1],Negatives_Father[irow,2]])
  }  
  #likelihood contributions for positive comparisons
  for(irow in 1:nrow(Pairs_Mother)){
    loglik=loglik+Pairs_Mother[irow,3]*log(P$P_Mother[Pairs_Mother[irow,1],Pairs_Mother[irow,2]])
  }
  for(irow in 1:nrow(Pairs_Father)){
    loglik=loglik+Pairs_Father[irow,3]*log(P$P_Father[Pairs_Father[irow,1],Pairs_Father[irow,2]])
  }  
  -loglik
}

####Sex-aggregated####
lemon_neg_log_lik_TotalA <- function(Pars2, Negatives_Parent, Pairs_Parent, P_Parent, t_start, t_end){
  
  P_TotalA=get_P_lemon_TotalA(Pars2=Pars2,P_Parent,t_start=t_start,t_end=t_end)
  
  loglik=0
  
  #likelihood contributions for all negative comparisons
  for(irow in 1:nrow(Negatives_Parent)){
    loglik=loglik+Negatives_Parent[irow,3]*log(1-P_TotalA$P_Parent[Negatives_Parent[irow,1],Negatives_Parent[irow,2]])
  } 
  #likelihood contributions for positive comparisons
  for(irow in 1:nrow(Pairs_Parent)){
    loglik=loglik+Pairs_Parent[irow,3]*log(P_TotalA$P_Parent[Pairs_Parent[irow,1],Pairs_Parent[irow,2]])
  }
  -loglik
}

#Fit model - optimx version
CK_fit1 <- optimx(par=Pars,fn=lemon_neg_log_lik,hessian=TRUE, method="BFGS", Negatives_Mother=mom_negatives, Negatives_Father=dad_negatives, Pairs_Mother=mom_positives, Pairs_Father=dad_positives, P_Mother=P_Mother, P_Father=P_Father, t_start=t_start, t_end=t_end)

CK_fit2 <- optimx(par=Pars2,fn=lemon_neg_log_lik_TotalA,hessian=TRUE, method="BFGS", Negatives_Parent=parent_negatives, Pairs_Parent=parent_positives, P_Parent=P_Parent, t_start=t_start, t_end=t_end)

#summary(CK_fit1)
#exp(CK_fit1[1:2])

####----#Need to check and work on script for calculating SE below----####

#compute variance covariance matrix - optimx
D1=diag(length(Pars))*c(exp(CK_fit1$p1[1]),exp(CK_fit1$p2[1])) #derivatives of transformations
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
yrs <- c(min_cohort:t_end)

#Extract true values from year of estimation (ie min_cohort)
Mom_truth <- round(mean(pop.size$Female.adult.pop[min_cohort:n_yrs]),0)
Dad_truth <- round(mean(pop.size$Male.adult.pop[min_cohort:n_yrs]), 0)
Adult_truth <- round(Mom_truth + Dad_truth, 0)

#Create dataframe of estimates and truth
estimates <- data.frame(cbind(round(exp(c(CK_fit1$p1[1], CK_fit1$p2[1], CK_fit2$p1[1])),0)), c(SE1, SE2), c("F", "M", "All"))
estimates <- cbind(estimates, c(Mom_truth, Dad_truth, Adult_truth))
colnames(estimates) <- c("CKMR_estimate", "SE", "Sex", "Mean_truth")
estimates

#Extract more metrics that can help with troubleshooting
total_samples <- sample.size * length(sample.years)
pop_size_mean <- round(mean(pop.size$population_size[min_cohort:n_yrs]),0)

metrics <- cbind(c(sum(mom_positives[,3]), sum(dad_positives[,3]), sum(parent_positives[,3])), 
                 c(rep(lam, times = 3)), 
                 c(rep(total_samples, times=3)),
                 c(rep(pop_size_mean, times=3)))
colnames(metrics) <- c("Parents_detected", "Pop_growth", "Total_samples", "Pop_size_mean")

#-----------------Loop end-----------------------------    
#Bind results from previous iterations with current iteration
results <- rbind(results, cbind(estimates, metrics))


print(paste0("finished iteration", iter, " at: ", Sys.time()))
}

#Calculate relative bias for all estimates
results <- results %>% 
  mutate(Relative_bias = round(((CKMR_estimate - Mean_truth)/Mean_truth)*100,1))

 results %>% group_by(Sex) %>% 
   dplyr::summarize(median = median(Relative_bias), n = n())

write.table(results, file = paste0("~/R/R_working_dir/LemonSharkCKMR_GitHub/02_IBS/Dovi_IBS_model_validation/Lemon_sharks/results/length_of_sampling/six_year_sampling/Dovi_AvgN_", total_samples, "_samples_05.17.2021_",length(sample.years), "yrs.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

#write.table(age_dist, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/fishSim_age.distributions_", total_samples, ".samples_02.10.2021_ages.correct_age.dist.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

#write.table(survival_at_age, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/fishSim_survival.at.age_", total_samples, ".samples_02.10.2021_ages.correct_surv.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
}

####-------Quick viz of results--------####
library(ggpubr)

#Box plot of relative bias
ggplot(data=results, aes(x=factor(Total_samples))) +
  geom_boxplot(aes(y=Relative_bias, fill=Sex)) +
  ylim(-100, 100) +
  geom_hline(yintercept=0, col="black", size=1.25) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  labs(x="Sample size", y="Relative bias", title="Relative Bias by sample size") +
  scale_fill_brewer(palette="Set2") +
  font("title", size = 10, face = "bold")
