#### LEMON SHARKS: DOVI'S BS MODEL

#Load packages
#devtools::install_github("BruceKendall/mpmtools")

library(optimx)
library(tidyverse) # safe to ignore conflicts with filter() and lag()
#library(popbio) #not needed
#library(mpmtools) #not needed
library(ggpubr)
library(rjags)
library(R2jags)
library(jagsUI)
library(Rlab)
library(postpack)

#######################################################################
### DATA-GENERATING MODEL

rm(list=ls())

# Note on sequencing: Births happen at beginning of each year, followed by deaths 
# (i.e. a female who dies in year 10 can still give birth and have surviving pups in year 10)

####---------Simulation parameters---------####

init.adult.pop.size <- 1000 # CHANGED FROM 3000; Initial adult population size
init.prop.female <- 0.5 # proportion of the initial population size that is female
birth.sex.ratio <- c(0.5,0.5) # The probability that each baby is F:M - has to add up to 1

YOY.survival <- 0.7 # CHANGED FROM 0.8; young of year survival
juvenile.survival <- 0.8 # CHANGED FROM 0.9; juvenile survival
Adult.survival <- 0.825 # CHANGED FROM 0.9; Adult survival

repro.age <- 12 # set age of reproductive maturity
max.age <- maxAge <- 50 # CHANGED FROM 30; set the maximum age allowed in the simulation

mating.periodicity <- 1 # CHANGED FROM 2; number of years between mating; assigned to an individual 
                        # and sticks with them through their life. So they're either a one or two year breeder.
num.mates <- c(1:3) # CHANGED FROM c(1:3); vector of potential number of mates per mating
#avg.num.offspring <- 3 # NOT USED? CHANGED FROM 3; set the average number of offspring per mating (from a poisson distribution) #JDS Q

f <- (1-Adult.survival)/(YOY.survival * juvenile.survival^11) # adult fecundity at equilibrium if no age truncation #JDS Q
ff <- f/init.prop.female * mating.periodicity/mean(num.mates) # female fecundity per breeding cycle
ff

# stable age distribution - JDS Q
props <- rep(NA, max.age+1)
props[1] <- f
props[2] <- f * YOY.survival
for (y in 3:(repro.age+1)) props[y] <- props[y-1] * juvenile.survival
#props[repro.age+1] <- props[repro.age] * juvenile.survival + Adult.survival
for (y in (repro.age+2):(max.age+1)) props[y] <- props[y-1] * Adult.survival

prop.Adult <- sum(props[(repro.age+1):(max.age+1)])/sum(props)

Nages <- round(props[-1] * init.adult.pop.size) 
init.pop.size <- sum(Nages) # all ages except YOYs

burn.in <- 40 # number of years to use as simulation burn in period
Num.years <- 50 # The number of years to run in the simulation beyond the burn in
n_yrs = t_end <- burn.in + Num.years

iterations <- 100 # CHANGED FROM 100; Number of iterations to loop over
#rseeds <- sample(1:1000000,iterations)
load("rseeds.rda")

#set.seed(885767)

####---------Sampling parameters---------####

#Decide which years to subtract from n_yrs for sampling and include in the vector below 
# (currently 5:0, which means sample the last 6 years of the simulation)
sample.years <- c(n_yrs - c(5:0))
sample.size <- 50 #sample size per year

#--------------Start simulation loop--------------

### MOVED SAMPLING BELOW, SO EXTRACT VARIOUS SAMPLE SIZES FROM THE SAME POPULATION

results <- NULL #initialize results array

for(iter in 1:iterations) {
  set.seed(rseeds[iter])
  
  ####---------Set up initial population-----####
  
  init.pop <- data.frame() # create a blank data frame which will become the initial population
  
  age.x <- NULL
  for(y in 1:length(Nages)) age.x <- c(age.x, rep(y, Nages[y]))
  
  for(i in 1:init.pop.size) { # Loop that creates the below data for each individual in a population the size of "init.pop.size"
    
    indv.name <- paste(sample(letters, size = 20, replace = T), collapse="") # generate a random name that is 20 letters long
    
    #age.x <- repro.age # every individual in the initial population is the age at first maturity ### WHY? WOULDN'T IT BE BETTER TO START WITH ~ STABLE AGE DISTRIBUTION?
    
    mother.x <- "xxxxx" # The individuals in the initial population do not have known mothers
    father.x <- "xxxxx" # The individuals in the initial population do not have known fathers
    birth.year <- -1 # this is a place holder for individuals born within the simulation
    sex <- sample(c('F','M'), size = 1, prob = c(init.prop.female, 1-init.prop.female)) # randomly draw sex based on the prportions set in the parameter section
    repro.cycle <- sample(1:mating.periodicity, size = 1) # randomly assign whether this individual mother would breed in the even or odd years
    init.vec <- cbind.data.frame(indv.name, birth.year, age.x[i], mother.x, father.x, sex, repro.cycle) # CHANGED THIS TO AGE.X[I]; create a row with the attributes for each individual - using cbind.data.frame allows them to keep some numeric and others characters
    init.pop <- rbind(init.pop, init.vec) # place the individual into the init.pop data frame
    
  }
  
  head(init.pop)
  names(init.pop)[3] <- "age.x"
  summary(as.factor(init.pop$age.x))
  summary(as.factor(init.pop$sex))
  
  ####----------Breeding----------####
  
  repro.cycle.vec <- rep(1:mating.periodicity, times = 100) # Generate a vector which will be used to determine if it is an even or odd breeding year (or a 1/3 breeding year)
  
  ####---------for year 0 breeding-------####
  
  mothers <- which(init.pop$sex=='F' & init.pop$age.x>=repro.age & init.pop$repro.cycle == repro.cycle.vec[1]) # determine which females are available to breed in this year
  fathers <- which(init.pop$sex=='M' & init.pop$age.x>=repro.age) # determine which fathers are available to breed in this year
  
  YOY.df <- data.frame() # create an empty data frame to populate with the YOY born in this year
  
  for(j in 1:length(mothers)) { # Loop through all of the available mothers
    
    num.mates.x <- sample(num.mates, size=1) # determine the number of mates for a given mother this year
    inter.df <- data.frame() # Create a data frame for the offspring born from each mother
    
    for(h in 1:num.mates.x) { # Loop through each sire that this mother mates with
      
      num.offspring <- rpois(1, ff) # CHANGED FROM rbinom(n = 1, size = 1, prob = ff); generate the number of offspring born from this mother/sire pairing #### WHY IS THIS BINOMIAL NOT, SAY, POISSON
      indv.name <- NA # create a place holder for the random name given to each offspring
      age.x <- 0 # assign a 0 age to each offspring born this year
      mother.x <- init.pop[mothers[j],1] # record the mothers name
      father.x <- init.pop[sample(fathers, size = 1),1] # record the sires name
      birth.year <- 0 # note that these pups were born in the 0 year of the simulation
      sex <- NA # create a place holder for the sexes
      repro.cycle <- sample(1:mating.periodicity, size = 1) # assign the newborn to an eventual breeding cycle group
      inter.df2 <- cbind.data.frame(indv.name, birth.year, age.x, mother.x, father.x, sex, repro.cycle) # Create a row with the attributes of each YOY (at this point only one from each mating pair)
      inter.df2 <- inter.df2[rep(seq_len(nrow(inter.df2)), num.offspring), ] # Create the random number of offspring from each pairing that we set above
      inter.df2$sex <- sample(c('F','M'), size = nrow(inter.df2), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
      baby.names <- vector() # Create a blank vector for random baby names
      
      for(w in 1:nrow(inter.df2)) { # create a random name for each newborn from this mating pair
        name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
        baby.names <- c(baby.names,name1)  
      } # end loop over name
      
      if(nrow(inter.df2)==0){next} # if there were no offspring from this mating pair, skips to the next mating pair (if you don't include this you will get an error in the loop)
      inter.df2$indv.name <- baby.names # add the baby names to the data frame
      inter.df <- rbind(inter.df, inter.df2) # add the newborns from this mating pair to the other from the same mother
    } # end loop over mates
    
    YOY.df <- rbind(YOY.df,inter.df) #add all of the offspring from each mother to a data frame of all the YOY for the year
  } # end loop over mothers
  
  year.end.pop.0 <- NULL
  
  year.end.pop.0 <- rbind(init.pop, YOY.df) #Combine the YOY data with the other individuals present this year
  
  #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
  year.end.pop.0$Survival <- ifelse(year.end.pop.0$age.x==0, sample(c("S","M"), size=length(which(year.end.pop.0$age.x==0)), prob=c(YOY.survival, 1-YOY.survival), replace=T),
                                    ifelse(year.end.pop.0$age.x<repro.age, sample(c("S","M"), size=length(which(year.end.pop.0$age.x<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                           ifelse(year.end.pop.0$age.x<max.age, sample(c("S","M"), size=length(which(year.end.pop.0$age.x<max.age)), prob=c(Adult.survival, 1-Adult.survival),replace=T), "M")))
  
  
  loopy.pop <- year.end.pop.0
  
  # At end of year 0 ...
  print(paste("year 0", "N_mothers=", length(mothers), "N_pups=", nrow(YOY.df), "N_deaths=", sum(loopy.pop$Survival=="M"), "N= ", nrow(loopy.pop[loopy.pop$Survival=="S",]) , sep=" ")) # print the simulation year and the population size in the R console so they can be oberved
  
  ####-------Loop for all other breeding years--------####
  
  pop.size <- data.frame()
  
  loopy.list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year
  
  for(v in 1:(burn.in + Num.years)){ #loop through all of the years in the simulation - the burn in and the years that matter
    
    data1 <- loopy.pop[loopy.pop$Survival =="S", -8] #Bring in the data from the previous iteration, but only include those that survive
    data1$age.x <- data1$age.x+1 # increase each individuals age by one for the new year - happy birthday survivors!
    
    mothers <- which(data1$sex=='F' & data1$age.x>=repro.age & data1$repro.cycle == repro.cycle.vec[v+1])  #determine which females are available to breed in this year; this is an index
    fathers <- which(data1$sex=='M' & data1$age.x>=repro.age)  #determine which males are available to breed in this year
    
    YOY.df <- data.frame() # create an empty data frame to populate with the YOY born in this year
    
    for(j in 1:length(mothers)) { # Loop through all of the available mothers
      
      num.mates.x <- sample(num.mates, size=1) # determine the number of mates for a given mother this year
      inter.df <- data.frame() # Create a data frame for the offspring born from each mother
      
      for(h in 1:num.mates.x) { # Loop through each sire with whom this mother mates 
        
        num.offspring <- rpois(1, ff) # CHANGED FROM rbinom(n = 1, size = 1, prob = ff); generate the number of offspring born from this mother/sire pairing ### AGAIN, WHY BINOMIAL???
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
      } # end loop over mates
      
      YOY.df <- rbind(YOY.df,inter.df) #add all of the offspring from each mother to a data frame of all the YOY for the year
    } # end loop over mothers
    
    loopy.pop <- rbind(data1, YOY.df) #Combine the YOY data with the other individuals present this year. YOY go at the bottom.
    
    #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
    loopy.pop$Survival <- ifelse(loopy.pop$age.x==0, sample(c("S","M"), size=length(which(loopy.pop$age.x==0)), prob=c(YOY.survival, 1-YOY.survival), replace=T),
                                 ifelse(loopy.pop$age.x<repro.age, sample(c("S","M"), size=length(which(loopy.pop$age.x<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                        ifelse(loopy.pop$age.x<max.age, sample(c("S","M"), size=length(which(loopy.pop$age.x<max.age)), prob=c(Adult.survival, 1-Adult.survival),replace=T), "M")))
    
    #assign(paste("year.end.pop.", v, sep=""),loopy.pop) # save the current year's population data as an object
    
    loopy.list[[v]] <- loopy.pop # Save the current year's population data as a list element, where the index corresponds to the year
    
    #  print(paste("year", v, "N= ", nrow(loopy.pop) , sep=" ")) # print the simulation year and the population size in the R console so they can be observed
    print(paste("year", v, "N_mothers=", length(mothers), "N_pups=", nrow(YOY.df), "N_deaths=", sum(loopy.pop$Survival=="M"), "N= ", nrow(loopy.pop[loopy.pop$Survival=="S",]) , sep=" ")) # CHANGED THIS
    
    ### CHANGED THIS TO ONLY COUNT SURVIVORS - JDS Q
    pop.size.vec <- cbind.data.frame(year=v, population_size=nrow(loopy.pop[loopy.pop$Survival=="S",]), # 
                                     Male.adult.pop = nrow(loopy.pop[loopy.pop$sex == "M" & loopy.pop$age.x >10 & loopy.pop$Survival=="S",]), # 
                                     Female.adult.pop = nrow(loopy.pop[loopy.pop$sex == "F" & loopy.pop$age.x >10 & loopy.pop$Survival=="S",])) # 
    pop.size <- rbind(pop.size, pop.size.vec)
    
  } # end loop over sim years
  
  #Label the list elements with the year
  names(loopy.list) <- paste0("year.end.pop.", seq(1:(burn.in + Num.years)))
  
  ##===============================================================================================
  ####---------quick analysis of population growth (lambda)-----####
  
  plot(pop.size$population_size, pch=19, ylim=c(0.9*min(pop.size$population_size), 1.1*max(pop.size$population_size)))
  
  #} # end loop over iters
  
  lambda <- vector()
  for(l in 2:nrow(pop.size)){ 
    lambda.1 <- pop.size$population_size[l]/pop.size$population_size[l-1]
    lambda <- c(lambda,lambda.1)
  }
  
  lambda <- c(NA, lambda)
  pop.size$Lambda <- lambda
  
  # plot(lambda[(burn.in+1):n_yrs], pch=19)
  # abline(h=1, lty=3)
  
  mean_lam <- mean(pop.size$Lambda[(burn.in+1):n_yrs], na.rm=T) # mean Lambda
  sd(pop.size$Lambda[(burn.in+1):n_yrs], na.rm=T) # sd Lambda
  
  ####---------Checking population parameters-------####
  
  nrow(YOY.df)/length(mothers) # Average fecundity for last year; remember whether they've skipped breeding
  
  # Adult.survival <- 0.9 # Adult survival
  # juvenile.survival <- 0.9 # juvenile survival
  # YOY.survival <- 0.8 # young of year survival
  
  nrow(loopy.pop[loopy.pop$Survival=='S' & loopy.pop$age.x>=repro.age, ])/nrow(loopy.pop[loopy.pop$age.x>=repro.age,]) #survival of adults
  nrow(loopy.pop[loopy.pop$Survival=='S' & loopy.pop$age.x>0 & loopy.pop$age.x<repro.age, ])/nrow(loopy.pop[loopy.pop$age.x>0 & loopy.pop$age.x<repro.age,]) #survival of juveniles
  nrow(loopy.pop[loopy.pop$Survival=='S' & loopy.pop$age.x==0, ])/nrow(loopy.pop[loopy.pop$age.x==0,]) #survival of YOY
  
  # survival_vec <- NULL
  # for(i in 1:length(loopy.list)){
  #   survival_vec[i] <- nrow(loopy.list[[i]][loopy.list[[i]]$Survival=='S' & loopy.list[[i]]$age.x>=repro.age, ])/nrow(loopy.list[[i]][loopy.list[[i]]$age.x>=repro.age, ])
  # }
  
  # nrow(loopy.pop[loopy.pop$Survival=='M' & loopy.pop$age.x>=12, ])+
  #   nrow(loopy.pop[loopy.pop$Survival=='M' & loopy.pop$age.x>0 & loopy.pop$age.x<12, ])+
  #   nrow(loopy.pop[loopy.pop$Survival=='M' & loopy.pop$age.x==0, ])/nrow(loopy.pop[loopy.pop$age.x==0,]) ### WHAT IS THIS? ratio of non-pup deaths to all pups
  # 
  # nrow(YOY.df) # the number of rows in a dataframe are the number of individuals in the population
  
  ####---------Starting parameters for estimation model------####
  
  #The initial adult population size ### HOW DO WE KNOW THIS???
  init.adult_pop.size <- loopy.list[[1]] %>% filter(age.x >= repro.age) %>%
    nrow()
  init.adult_pop.size
  
  # Set initial parameters for model based on initial abundance of males and females
  f.init <- m.init <- init.adult_pop.size/2
  Pars <- c(log(f.init), log(m.init)) #Pars1 is for the sex-specific model
  Pars2 <- log(init.adult_pop.size) #Pars2 is for the sex-aggregated model
  
  #####################################################################################
  ### DATA SAMPLING
  
  for(samps in 1:3){
  #try({  
     
    sample.size <- c(40, 50, 60)[samps] #To loop over different sample sizes, draw a different number of samples each time
    
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
    pairwise.df <- data.frame(t(combn(sample.df_all.info$indv.name, m=2))) # generate all combinations of the elements of x taken m at a time.
    colnames(pairwise.df) <- c("Ind_1", "indv.name") #Rename columns so they can easily be joined
    head(pairwise.df)
    
    #Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
    Ind1_birthyears <- sample.df_all.info %>%
      select(indv.name, birth.year, age.x, mother.x, father.x) %>% #select relevant columns only
      dplyr::rename("Ind_1" = indv.name, "Ind_1_birth" = birth.year, "Ind_1_age" = age.x, "Ind_1_mom" = mother.x, "Ind_1_dad" = father.x) 
    #Rename columns for join and also so younger sib birth year and parents are distinguishable from older sib data when joined below.
    
    #Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
    #This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
    pairwise.df_all.info <- pairwise.df %>% left_join(Ind1_birthyears, by = "Ind_1") %>% 
      left_join(sample.df_all.info, by = "indv.name") %>% 
      dplyr::rename("Ind_2" = indv.name, "Ind_2_birth" = birth.year, "Ind_2_age" = age.x, "Ind_2_mom" = mother.x, "Ind_2_dad" = father.x) %>% 
      select(Ind_1, Ind_1_birth, Ind_1_age, Ind_2, Ind_2_birth, Ind_2_age, Ind_1_mom, Ind_2_mom, Ind_1_dad, Ind_2_dad) %>% 
      filter(Ind_1_birth != Ind_2_birth)
    
    head(pairwise.df_all.info )
    
    #Check that there are no repeat comparisons -- compare number of distinct comparisons to total number of comparisons.
    #Should return TRUE
    pairwise.df_all.info %>% distinct(Ind_1, Ind_2) %>% nrow() == nrow(pairwise.df_all.info)
    nrow(pairwise.df_all.info)
    
    #Extract positive half-sib comparisons
    positives <- pairwise.df_all.info %>% filter(Ind_1_mom == Ind_2_mom | Ind_1_dad == Ind_2_dad) 
    nrow(positives)
    
    #Remove full sibs -- adjust JDS
    positives <- positives %>% filter(Ind_1_mom != Ind_2_mom | Ind_1_dad != Ind_2_dad) # EDITED CODE TO REMOVE FULL SIBS 
    nrow(positives)
    
    #Second filter to check for self-recaptures
    self <- positives %>% filter(Ind_1 == Ind_2)
    nrow(self) #They should have been filtered earlier so should be 0
    
    ####----------------Split dataframes into final form for model----------####
    
    #Sex-specific
    mom_positives <- positives %>% filter(Ind_1_mom == Ind_2_mom) %>% 
      select(Ind_1_birth, Ind_2_birth) %>% 
      plyr::count() 
    
    dad_positives <- positives %>% filter(Ind_1_dad == Ind_2_dad)  %>%
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
    
    #-------------Kinship probabilities - Half-sib-------------------
    min_cohort <- n_yrs-40 # CHANGED THIS FROM MAX.AGE; set first year for calculating mean (arbitrary)
    
    m_adult_age <- f_adult_age <- c(repro.age:max.age) # Set ages at which males and females are mature. Called by kinship probability function.
    
    pop_growth_all_mean <- mean(pop.size$Lambda[1:nrow(pop.size)], na.rm=T)
    
    ###Used in model###
    surv <- Adult.survival #Set value ### HOW DO WE KNOW THIS?
    #surv <- mean_adult_surv #Observed value
    
    lam <- mean(pop.size$Lambda[min_cohort:n_yrs], na.rm=T) # mean Lambda over years of estimation ### HOW DO WE KNOW THIS?
    
    #### Fit Bayesian model ####
    
    ################ STEP 1: PREPARE DATA #################
    mom_comps <- mom_positives %>% 
      rename(yes = freq) %>% 
      full_join(mom_negatives, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
      rename(no = freq) %>% 
      mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>% 
      mutate(all = yes + no)
    
    dad_comps <- dad_positives %>% 
      rename(yes = freq) %>% 
      full_join(dad_negatives, by = c("Ind_1_birth", "Ind_2_birth")) %>% 
      rename(no = freq) %>% 
      mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>% 
      mutate(all = yes + no)
    
    #head(mom_comps)
    #head(dad_comps)
    
    #Need priors for:
    #Number of adults (uninformative)
    #maybe eventually survival (beta -- conjugate prior for binomial)
    
    #Define data
    jags_data = list(
      #Moms
      MHSP = mom_comps$yes,
      mom_n_comps = mom_comps$all,
      mom_ys_birth = mom_comps$Ind_2_birth,
      mom_os_birth = mom_comps$Ind_1_birth,
      mom_yrs = nrow(mom_comps),
      
      #Dads
      FHSP = dad_comps$yes,
      dad_n_comps = dad_comps$all,
      dad_ys_birth = dad_comps$Ind_2_birth,
      dad_os_birth = dad_comps$Ind_1_birth,
      dad_yrs = nrow(dad_comps),
      
      #Fix other potential parameters
      surv = surv,
      lam = lam,
      min_cohort = min_cohort
    )
    
    #Define initial values for priors
    
    ################### STEP 2: SPECIFY JAGS MODEL CODE ##################
    ##CKMR code
    HS_model = function(){
      #PRIORS
      Nf ~ dnorm(0, 1.0E-6) #Uninformative prior for female abundance
      Nm ~ dnorm(0, 1.0E-6) ##Uninformative prior for male abundance
      #Likelihood
      for(i in 1:mom_yrs){
        MHSP[i] ~ dbin((surv^(mom_ys_birth[i] - mom_os_birth[i]))/(Nf*lam^(mom_ys_birth[i]-min_cohort)), mom_n_comps[i])
      }
      for(j in 1:dad_yrs){
        FHSP[j] ~ dbin((surv^(dad_ys_birth[j] - dad_os_birth[j]))/(Nm*lam^(dad_ys_birth[j]-min_cohort)), dad_n_comps[j])
      }
    }
    
    jags_file = paste0("~/R/working_directory/LemonSharkCKMR/02_IBS/Dovi_IBS_model_validation/Lemon_sharks/results/basic_bayesian/models/HS_model_iteration_", iter, ".txt")
    write_model(HS_model, jags_file)
    
    
    ########### STEP 3: SPECIFY INITIAL VALUES ##################
    jags_inits = function(nc) {
      inits = list()
      for(c in 1:nc){
        inits[[c]] = list(
          Nf = rnorm(1, mean = 500, sd = 100),
          Nm = rnorm(1, mean = 500, sd = 100)
        )
      }
      return(inits)
    }
    
    ########## STEP 4: SET NODES TO MONITOR ################
    jags_params = c("Nf", "Nm")
    
    
    ########### STEP 5: SET MCMC DIMENSIONS ################
    jags_dims = c(
      ni = 5000,  # number of post-burn-in samples per chain
      nb = 5000,  # number of burn-in samples
      nt = 1,     # thinning rate
      nc = 2      # number of chains
    )
    
    
    ##### STEP 6: RUN THE MODEL WITH JAGS #####
    
    post = jagsUI::jags.basic(data = jags_data, #If using postpack from AFS workshop
                              
                              #post = rjags::jags(data = jags_data, #If wanting to use other diagnostics
                              model.file = jags_file,
                              inits = jags_inits(jags_dims["nc"]),
                              parameters.to.save = jags_params,
                              n.adapt = 1000,
                              n.iter = sum(jags_dims[c("ni", "nb")]),
                              n.thin = jags_dims["nt"],
                              n.burnin = jags_dims["nb"],
                              n.chains = jags_dims["nc"],
                              parallel = F
    )
    
    ##### STEP 7: CONVERGENCE DIAGNOSTICS #####
    # view convergence diagnostic summaries for all monitored nodes
    model_summary <- data.frame(t(post_summ(post, jags_params, Rhat = T, neff = T)))

##### Compile and report results ####
    #Combine above to make dataframe with truth and estimates side-by-side
    #store years from youngest sibling in comparisons to end of study
    yrs <- c(min_cohort:t_end)
    
    #Extract true values from year of estimation (ie min_cohort)
    pop.size <- pop.size %>% mutate(Total.adult.pop = Male.adult.pop + Female.adult.pop)
    Mom_truth <- round(pop.size$Female.adult.pop[min_cohort],0)
    Dad_truth <- round(pop.size$Male.adult.pop[min_cohort], 0)
    Adult_truth <- round(pop.size$Total.adult.pop[min_cohort], 0)
    Mom_min <- min(pop.size$Female.adult.pop[min_cohort:n_yrs])
    Mom_max <- max(pop.size$Female.adult.pop[min_cohort:n_yrs])
    Dad_min <- min(pop.size$Male.adult.pop[min_cohort:n_yrs])
    Dad_max <- max(pop.size$Male.adult.pop[min_cohort:n_yrs])
    Adult_min <- min(pop.size$Total.adult.pop[min_cohort:n_yrs])
    Adult_max <- max(pop.size$Total.adult.pop[min_cohort:n_yrs])
    
    
    #Create dataframe of estimates and truth
    estimates <- model_summary %>% rownames_to_column(var = "parameter") %>% 
      mutate(min = c(Mom_min, Dad_min), max = c(Mom_max, Dad_max), truth = c(Mom_truth, Dad_truth)) %>%
      as_tibble() %>% 
      rename(`50` = X50., `2.5` = X2.5., `97.5` = X97.5.)
    
    #Extract more metrics that can help with troubleshooting
    total_samples <- sample.size * length(sample.years)
    pop_size_mean <- round(mean(pop.size$population_size[min_cohort:n_yrs]),0)
    
    metrics <- cbind(c(sum(mom_positives[,3]), sum(dad_positives[,3])), 
                     c(rep(lam, times = 2)), 
                     c(rep(pop_growth_all_mean, times = 2)),
                     c(rep(total_samples, times=2)),
                     c(rep(pop_size_mean, times=2)),
                     c(rep(iter, times = 2)))
    colnames(metrics) <- c("parents_detected", "pop_growth_est_yrs", "pop_growth_all_yrs", "total_samples", "pop_size_mean", "iteration")
    
    #-----------------Loop end-----------------------------    
    #Bind results from previous iterations with current iteration
    results <- rbind(results, cbind(estimates, metrics))
  
#  }) # end try clause    
  } # end loop over sample sizes
  
  print(paste0("finished iteration", iter, " at: ", Sys.time()))
} # end loop over iterations

##################################################################################
### SAVE AND CHECK RESULTS

#Calculate relative bias for all estimates

results2 <- results %>% 
  mutate(relative_bias = round(((mean - truth)/truth)*100,1)) %>% # CHANGED TABLE NAME SO CAN BUILD & CHECK RESULTS ITERATIVELY
  mutate(in_interval = ifelse(`2.5` < mean & mean < `97.5`, "Y", "N")) 

 results2 %>% group_by(total_samples, parameter) %>% 
   dplyr::summarize(median = median(relative_bias), n = n())

#Home computer: Dell Precision
write.table(results2, file = paste0("~/R/working_directory/LemonSharkCKMR/02_IBS/Dovi_IBS_model_validation/Lemon_sharks/results/basic_bayesian/simple_model_validation.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

#MGHPCC
#write.table(results2, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/Dovi_lambdaModel_06_22.2021_neutralPopGrowth.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

#write.table(age_dist, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/fishSim_age.distributions_", total_samples, ".samples_02.10.2021_ages.correct_age.dist.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)

#write.table(survival_at_age, file = paste0("/home/js16a/R/working_directory/CKMR_simulations/results/fishSim_survival.at.age_", total_samples, ".samples_02.10.2021_ages.correct_surv.csv"), sep=",", dec=".", qmethod="double", row.names=FALSE)
#}

####-------Quick viz of results--------####

#Box plot of relative bias
ggplot(data=results2, aes(x=factor(Total_samples))) +
  geom_boxplot(aes(y=Relative_bias, fill=Sex)) +
  ylim(-100, 100) +
  geom_hline(yintercept=0, col="black", size=1.25) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  labs(x="Sample size", y="Relative bias", title="Relative Bias by sample size") +
  scale_fill_brewer(palette="Set2") +
  font("title", size = 10, face = "bold")

##################################################################################