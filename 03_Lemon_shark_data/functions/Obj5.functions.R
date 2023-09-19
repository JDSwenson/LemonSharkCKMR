#Note that mating periodicity only applies to females, even though males are assigned it as well. And the repro.cycle column determines which year the animal breeds in, not the frequency. So females can have a 1 or 2 for biennial breeding. And males can be anything, but it doesn't actually come into play.

simulate.pop <- function(init.pop.size, init.prop.female, Nages, mating.periodicity, repro.age, YOY.survival, juvenile.survival, Adult.survival, max.age, num.mates, ff, burn.in, Num.years){
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
    repro.strategy <- sample(c("conformist", "non-conformist"), size = 1, prob = c(1 - non.conformists, non.conformists))
    init.vec <- cbind.data.frame(indv.name, birth.year, age.x[i], mother.x, father.x, sex, repro.cycle, repro.strategy) # CHANGED THIS TO AGE.X[I]; create a row with the attributes for each individual - using cbind.data.frame allows them to keep some numeric and others characters
    init.pop <- rbind(init.pop, init.vec) # place the individual into the init.pop data frame
    
  }
  
  head(init.pop)
  names(init.pop)[3] <- "age.x"
  summary(as.factor(init.pop$age.x))
  summary(as.factor(init.pop$sex))
  
  ####----------Breeding----------####
  
  repro.cycle.vec <- rep(1:mating.periodicity, times = 100) # Generate a vector which will be used to determine if it is an even or odd breeding year (or a 1/3 breeding year)
  
  ####---------for year 0 breeding-------####
  
  #------------------------------------Conformist mothers------------------------------------#
  mothers <- which(init.pop$sex=='F' & init.pop$age.x>=repro.age & init.pop$repro.cycle == repro.cycle.vec[1] & init.pop$repro.strategy == "conformist") # determine which females are available to breed in this year
  
  #-----------------------------------Non-conformist mothers------------------------------------#
  #--------------If allowing all annual breeders to breed at age 12
  mothers.2 <- which(init.pop$sex=='F' & init.pop$age.x>=repro.age & init.pop$repro.strategy == "non-conformist")
  mothers <- c(mothers, mothers.2)
  
  #-------------If allowing half of annual breeders to breed at age 12
  #For annual breeders that are the age-at-maturity (12), only allow 50% to breed. Annual breeders older than this age always breed. This means that half of biennial breeders and half of annual breeders give birth for the first time at repro.age, and the other half give birth for the first time the following year. This effectively makes the age at reproduction 12.5 for all animals.
  # mothers.2 <- which(init.pop$sex=='F' & init.pop$age.x>repro.age & init.pop$repro.strategy == "non-conformist")
  # mothers.3 <- which(init.pop$sex=='F' & init.pop$age.x==repro.age & init.pop$repro.strategy == "non-conformist" & init.pop$repro.cycle == repro.cycle.vec[1]) #Only time non-conformists rely on the repro.cycle.
  # #Combine conformist and non-conformist mothers
  # mothers <- c(mothers, mothers.2, mothers.3)
  
  
  #-----------------------------------Fathers--------------------------------------------#
  fathers <- which(init.pop$sex=='M' & init.pop$age.x>=repro.age) # determine which fathers are available to breed in this year
  
  YOY.df <- data.frame() # create an empty data frame to populate with the YOY born in this year
  
  for(j in 1:length(mothers)) { # Loop through all of the available mothers
    
    num.mates.x <- sample(num.mates, size=1) # determine the number of mates for a given mother this year
    inter.df <- data.frame() # Create a data frame for the offspring born from each mother
    
    for(h in 1:num.mates.x) { # Loop through each sire that this mother mates with
      
      #num.offspring <- rpois(1, ff) # CHANGED FROM rbinom(n = 1, size = 1, prob = ff); generate the number of offspring born from this mother/sire pairing
      
      #(liz insert) ====
      xpup <- ceiling(ff)- ff
      num.offspring <- ifelse(runif(1)<xpup, trunc(ff), ceiling(ff))
      #(end liz insert) ====
      
      indv.name <- NA # create a place holder for the random name given to each offspring
      age.x <- 0 # assign a 0 age to each offspring born this year
      mother.x <- init.pop[mothers[j],1] # record the mothers name
      father.x <- init.pop[sample(fathers, size = 1),1] # record the sires name
      birth.year <- 0 # note that these pups were born in the 0 year of the simulation
      sex <- NA # create a place holder for the sexes
      
      repro.cycle <- sample(1:mating.periodicity, size = 1) # assign the newborn to an eventual breeding cycle group
      repro.strategy <- sample(c("conformist", "non-conformist"), size = 1, prob = c(1 - non.conformists, non.conformists)) #assign the newborn to be a conformist or non-conformist if a female
      
      inter.df2 <- cbind.data.frame(indv.name, birth.year, age.x, mother.x, father.x, sex, repro.cycle, repro.strategy) # Create a row with the attributes of each YOY (at this point only one from each mating pair)
      inter.df2 <- inter.df2[rep(seq_len(nrow(inter.df2)), num.offspring), ] # Create the random number of offspring from each pairing that we set above
      inter.df2$sex <- sample(c('F','M'), size = nrow(inter.df2), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ratio set in the parameter section
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
  year.end.pop.0$Survival <- ifelse(year.end.pop.0$age.x==0, 
                                    sample(c("S","M"), 
                                           size=length(which(year.end.pop.0$age.x==0)), 
                                           prob=c(YOY.survival, 1-YOY.survival), 
                                           replace=T),
                                    ifelse(year.end.pop.0$age.x<repro.age, 
                                           sample(c("S","M"), 
                                                  size=length(which(year.end.pop.0$age.x<repro.age)), 
                                                  prob=c(juvenile.survival, 1-juvenile.survival),
                                                  replace=T),
                                           ifelse(year.end.pop.0$age.x<max.age, 
                                                  sample(c("S","M"), 
                                                         size=length(which(year.end.pop.0$age.x<max.age)), 
                                                         prob=c(Adult.survival, 1-Adult.survival),
                                                         replace=T), "M")))
  
  
  loopy.pop <- year.end.pop.0
  
  # At end of year 0 ...
  print(paste("year 0", "N_mothers=", length(mothers), "N_pups=", nrow(YOY.df), "N_deaths=", sum(loopy.pop$Survival=="M"), "N= ", nrow(loopy.pop[loopy.pop$Survival=="S",]) , sep=" ")) # print the simulation year and the population size in the R console so they can be oberved
  
  ####-------Loop for all other breeding years--------####
  
  pop.size <- data.frame()
  
  loopy.list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year
  parents.tibble <- tibble()
  moms.temp = dads.temp <- NULL
  
  
  
  
  
  for(v in 1:(burn.in + Num.years)){ #loop through all of the years in the simulation - the burn in and the years that matter
    
    #Population will increase starting in year 70
    if(v == n_yrs - 20){
      
      (Ma <- -log(Adult.survival)) #Mortality of adults
      (Mj <- -log(juvenile.survival)) #Mortality of juveniles
      (Sa <- exp(-Ma)) #Survival of adults
      (Sj <- exp(-Mj)) #Survival of juveniles
      
      Mx <- -0.05 #Extra mortality (make negative so it reduces mortality)
      (Sa <- exp(-Ma - Mx)) #Survival of adults
      (Sj <- exp(-Mj - Mx)) #Survival of juveniles
      
      Adult.survival <- Sa
      juvenile.survival <- Sj
      
      cat(paste0("Adult survival is ", round(Adult.survival, 3), ";\nJuvenile survival is ", round(juvenile.survival, 3), ";\nPopulation will increase in size"))
      
    }
    
    #Population will stabilize for five years
    if(v == n_yrs - 10){
      
      Adult.survival <- 0.825
      juvenile.survival <- 0.8
      
      cat(paste0("Adult survival is ", round(Adult.survival, 3), ";\nJuvenile survival is ", round(juvenile.survival, 3), ";\nPopulation will stabilize"))
      
    }
    
    if(v == n_yrs - 5){
      
      (Ma <- -log(Adult.survival)) #Mortality of adults
      (Mj <- -log(juvenile.survival)) #Mortality of juveniles
      (Sa <- exp(-Ma)) #Survival of adults
      (Sj <- exp(-Mj)) #Survival of juveniles
      
      Mx <- 0.05 #Extra mortality
      (Sa <- exp(-Ma - Mx)) #Survival of adults
      (Sj <- exp(-Mj - Mx)) #Survival of juveniles
      
      Adult.survival <- Sa
      juvenile.survival <- Sj
      
      cat(paste0("Adult survival is ", round(Adult.survival, 3), ";\nJuvenile survival is ", round(juvenile.survival, 3), ";\nPopulation will now decrease in size"))
      
    }
    
    data1 <- loopy.pop[loopy.pop$Survival =="S", -9] #Bring in the data from the previous iteration, but only include those that survive (and leave the column of survival out)
    data1$age.x <- data1$age.x+1 # increase each individuals age by one for the new year - happy birthday survivors!
    
    
    
    #------------------------------------Conformist mothers------------------------------------#
    
    
    #Add all mothers that are meant to breed in the year, plus a small proportion that are not
    mothers <- which(data1$sex=='F' & data1$age.x>=repro.age & data1$repro.cycle == repro.cycle.vec[v+1] & data1$repro.strategy == "conformist") #determine which females are available to breed in this year; this is an index
    
    #-----------------------------------Switch cycles------------------------#
    
    #Remove a percent of on-cycle breeders from the breeding pool this year
    if(exists("percent.skip_on.cycle") == TRUE){
      if(percent.skip_on.cycle > 0){
        #Randomly select from the vector of on-cycle mothers some that will not breed this year
        mothers.skip_on.cycle <- sample(mothers, size = length(mothers)*percent.skip_on.cycle)
        
        #Remove the non-breeding on-cycle mothers from the vector of breeding mothers for this year
        mothers <- c(mothers[!mothers %in% mothers.skip_on.cycle])
      }
    }
    
    
    #Include a percentage of off-cycle breeders in the breeding pool this year
    if(exists("percent.breed_off.cycle") == TRUE){ 
      if(percent.breed_off.cycle > 0){
        #Which potential mothers are off cycle?
        mothers.all_off.cycle <- which(data1$sex=='F' & data1$age.x>=repro.age & data1$repro.cycle == repro.cycle.vec[v] & data1$repro.strategy == "conformist")
        
        #Randomly select from the vector of off-cycle mothers the percent specified by percent.breed_off.cycle
        mothers.breed_off.cycle <- sample(mothers.all_off.cycle, size = length(mothers.all_off.cycle)*percent.breed_off.cycle)
        
        #Add the off cycle mothers to the vector of breeding mothers for this year
        mothers <- sort(c(mothers, mothers.breed_off.cycle))
      }
    }
    
    #-----------------------------------Non-conformist mothers------------------------------------#
    #------------------If allowing all annual breeders to breed at age 12
    mothers.2 <- which(data1$sex=='F' & data1$age.x>=repro.age & data1$repro.strategy == "non-conformist") #All non-conformists greater than age 12 will breed; assume annual breeders
    mothers <- c(mothers, mothers.2)
    
    #-----------------If allowing half of annual breeders to breed at age 12
    #For annual breeders that are the age-at-maturity (12), only allow 50% to breed. Annual breeders older than this age always breed. This means that half of biennial breeders and half of annual breeders give birth for the first time at repro.age, and the other half give birth for the first time the following year. This effectively makes the age at reproduction 12.5 for all animals.
    #  mothers.2 <- which(data1$sex=='F' & data1$age.x>repro.age & data1$repro.strategy == "non-conformist") #All non-conformists greater than age 12 will breed; assume annual breeders
    #  mothers.3 <- which(data1$sex=='F' & data1$age.x==repro.age & data1$repro.strategy == "non-conformist" & data1$repro.cycle == repro.cycle.vec[v+1]) #Only time non-conformists rely on the repro.cycle. 
    # #Combine conformist and non-conformist mothers
    # mothers <- c(mothers, mothers.2, mothers.3)
    
    
    
    
    
    #     mothers.off.yr <- which(data1$sex=='F' & data1$age.x>=repro.age & data1$repro.cycle != repro.cycle.vec[v+1])
    #   mom.off.add <- sample(mothers.off.yr, 
    #                         size = length(mothers.off.yr)*non.conformists)
    # 
    # mothers <- c(mothers, mom.off.add)
    
    #-----------------------------------Fathers--------------------------------------------#
    fathers <- which(data1$sex=='M' & data1$age.x>=repro.age)  #determine which males are available to breed in this year
    
    YOY.df <- data.frame() # create an empty data frame to populate with the YOY born in this year
    
    for(j in 1:length(mothers)) { # Loop through all of the available mothers
      
      num.mates.x <- sample(num.mates, size=1) # determine the number of mates for a given mother this year
      inter.df <- data.frame() # Create a data frame for the offspring born from each mother
      
      for(h in 1:num.mates.x) { # Loop through each sire with whom this mother mates 
        
        #      num.offspring <- rpois(1, ff) #generate the number of offspring born from this mother/sire pairing
        #Liz insert =====
        xpup <- ceiling(ff)- ff
        num.offspring <- ifelse(runif(1)<xpup, trunc(ff), ceiling(ff))
        #End Liz insert =====
        
        indv.name <- NA #create a place holder for the random name given to each offspring
        age.x <- 0 # assign a 0 age to each offspring born this year
        mother.x <- data1[mothers[j],1] # record the mothers name
        father.x <- data1[sample(fathers, size = 1),1] # record the fathers name
        birth.year <- v # assign the birth year in the simulation
        sex <- NA # create a place holder for the biological sex of each offspring
        repro.cycle <- sample(1:mating.periodicity, size = 1) # assign the newborn to an eventual breeding cycle group
        
        inter.df2 <- cbind.data.frame(indv.name, birth.year, age.x, mother.x, father.x, sex, repro.cycle, repro.strategy) # Create a row with the attributes of each YOY (at this point only one from each mating pair)
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
      
      YOY.df <- rbind(YOY.df,inter.df) %>%  #add all of the offspring from each mother to a data frame of all the YOY for the year
        mutate(repro.strategy = ifelse(sex == "F", 
                                       sample(c("conformist", "non-conformist"), 
                                              n(),
                                              prob = c(1 - non.conformists, non.conformists),
                                              replace = TRUE), 
                                       "conformist")) #If they're a female, assign conformity/non-conformity at the specified probability
      
    } # end loop over mothers
    
    loopy.pop <- rbind(data1, YOY.df) #Combine the YOY data with the other individuals present this year. YOY go at the bottom.
    
    #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
    loopy.pop$Survival <- ifelse(loopy.pop$age.x==0, 
                                 sample(c("S","M"), 
                                        size=length(which(loopy.pop$age.x==0)), 
                                        prob=c(YOY.survival, 1-YOY.survival), 
                                        replace=T),
                                 ifelse(loopy.pop$age.x<repro.age, 
                                        sample(c("S","M"), 
                                               size=length(which(loopy.pop$age.x<repro.age)), 
                                               prob=c(juvenile.survival, 1-juvenile.survival),
                                               replace=T),
                                        ifelse(loopy.pop$age.x<max.age, 
                                               sample(c("S","M"), 
                                                      size=length(which(loopy.pop$age.x<max.age)), 
                                                      prob=c(Adult.survival, 1-Adult.survival),
                                                      replace=T), "M")))
    
    #assign(paste("year.end.pop.", v, sep=""),loopy.pop) # save the current year's population data as an object
    
    loopy.list[[v]] <- loopy.pop # Save the current year's population data as a list element, where the index corresponds to the year
    
    moms.temp <- YOY.df %>% group_by(mother.x) %>% 
      summarize(num.off = n()) %>% 
      rename(parent = mother.x) %>% 
      mutate(year = v, parent.sex = "mother")
    
    dads.temp <- YOY.df %>% group_by(father.x) %>% 
      summarize(num.off = n()) %>% 
      rename(parent = father.x) %>% 
      mutate(year = v, parent.sex = "father")
    
    parents.tibble <- rbind(parents.tibble, moms.temp, dads.temp)
    
    #  print(paste("year", v, "N= ", nrow(loopy.pop) , sep=" ")) # print the simulation year and the population size in the R console so they can be observed
    cat(paste("\nyear", v, 
              "\nN_potential_mothers=", length(mothers), 
              "N_actual_mothers=", nrow(moms.temp), 
              "Percent female breeders=", round(nrow(moms.temp)/length(mothers)*100,0),
              "\nN_potential_fathers=", length(fathers), 
              "N_actual_fathers=", nrow(dads.temp), 
              "Percent male breeders=", round(nrow(dads.temp)/length(fathers)*100, 0),
              "\nN_pups=", nrow(YOY.df), 
              "\nN_deaths=", sum(loopy.pop$Survival=="M"), 
              "\nTotal N= ", nrow(loopy.pop[loopy.pop$Survival=="S",]) , sep=" ")) # CHANGED THIS
    
    ### CHANGED THIS TO ONLY COUNT SURVIVORS - JDS Q
    pop.size.vec <- cbind.data.frame(year=v, 
                                     population_size=nrow(data1), # 
                                     Male.adult.pop = nrow(data1[data1$sex == "M" & data1$age.x >= repro.age,]), # 
                                     Female.adult.pop = nrow(data1[data1$sex == "F" & data1$age.x >= repro.age,]), #difference between this and number of mothers only if reproductive cycle != 1
                                     Num.mothers = nrow(moms.temp), #The actual number that reproduced this year
                                     Num.fathers = nrow(dads.temp)) #The actual number that reproduced this year
    pop.size <- rbind(pop.size, pop.size.vec)
    
  } # end loop over sim years
  
  #Label the list elements with the year
  names(loopy.list) <- paste0("year.end.pop.", seq(1:(burn.in + Num.years)), "_iteration_", iter)
  
  return(invisible(list(loopy.list, pop.size, parents.tibble)))
}





#Find aunt|uncle/niece|nephew pairs
find_aunt.uncle_nephew.niece_pairs <- function(input.sample.info = sample.df_all.info_temp1, loopy.list = loopy.list){
  
  #Make dataframe of sampled individuals parents and the sampled indv's birth year (so we can subset from loopy.list). Preserve column names for join
  #Would have to insert this into the sampling scheme loop ... 
  sampled.indvs.rents.df <- input.sample.info %>% 
    dplyr::select(mother.x, father.x, birth.year) %>% 
    arrange(birth.year)
  
  #Make a vector with all unique parents
  rents.vec <- unique(c(sampled.indvs.rents.df$mother.x, sampled.indvs.rents.df$father.x))
  length(rents.vec) #How many unique parents for sampled individuals?
  
  #Make a dataframe of all loopy.list dataframes from the oldest sampled indv's birth year until present
  b.y <- min(sampled.indvs.rents.df$birth.year)
  
  #Bind rows of loopy.list together so we can grab the info from all parents in the dataset
  #Confirmed it's doing what we want it to
  loopy.list.merged.df <- bind_rows(loopy.list[b.y:n_yrs]) %>%
    distinct(indv.name, .keep_all = TRUE) %>% #Only keep one copy of each individual
    dplyr::filter(indv.name %in% rents.vec) %>% #Only keep individuals that are parents of sampled individuals
    as_tibble() %>% 
    rename(grandmother = mother.x, grandfather = father.x) #The parents of the sampled indvidual's parents are the grandparents of the sampled individuals
  
  #Split into separate dfs for maternal and paternal lineages. Keep the sampled indv as the reference point.
  mothers.grandparents.df <- loopy.list.merged.df %>% dplyr::filter(sex == "F") %>%
    dplyr::select(mother = indv.name, 
                  mom.grandmother = grandmother, 
                  mom.grandfather = grandfather)
  
  fathers.grandparents.df <- loopy.list.merged.df %>% dplyr::filter(sex == "M") %>%
    dplyr::select(father = indv.name,
                  dad.grandmother = grandmother,
                  dad.grandfather = grandfather)
  
  #Add grandparents to each sample dataframe
  sample.info_w_grandparents <- input.sample.info %>% 
    as_tibble() %>% 
    rename(mother = mother.x, father = father.x) %>% 
    distinct(indv.name, .keep_all = TRUE) %>% #Don't want individuals that were captured more than once, or they will be duplicated below
    inner_join(mothers.grandparents.df, by = "mother") %>% #add maternal grandmother info to the df, joining by mom
    inner_join(fathers.grandparents.df, by = "father") %>% #add paternal grandmother info to the df, joining by dad
    dplyr::select(indv.name, 
                  mother, 
                  father, 
                  mom.grandmother, 
                  mom.grandfather, 
                  dad.grandmother, 
                  dad.grandfather, 
                  birth.year, 
                  sex, 
                  capture.year) %>% #Keep most columns but remove the x from mother and father
    mutate(both.rents = paste0(mother, father), #Reduce parents to one variable
           both.mom.grandparents = paste0(mom.grandmother, mom.grandfather),
           both.dad.grandparents = paste0(dad.grandmother, dad.grandfather)) 
  
  
  ### 09/05/2023: REMOVED the below filter, which keeps just one indv from each litter.
  #%>% group_by(mother, father, birth.year) %>% #Filter for full sibs 
  #    slice_sample(n = 1) %>% #Keep one sibling. The pairwise comparison script just keeps one indv from each litter, so we can do that here too. #QUESTION for JDS: is this appropriate? I think it is ... 
  #    ungroup()
  
  #Q for Ben or Anthony: is it good that we're only keeping one indv from each litter? Or can we easily integrate full sibs into the likelihood?
  #Check how many full sibs there are after filtering (should be 0)
  # sample.info_w_grandparents %>% 
  #   group_by(mother, father, birth.year) %>%
  #   summarize(n=n()) %>%
  #   dplyr::filter(n>1)
  
  #Save each combos of parents as vector objects to ease code later
  both.parents.vec <- sample.info_w_grandparents$both.rents %>% unique()
  both.mom.grandparents.vec <- sample.info_w_grandparents$both.mom.grandparents %>% unique()
  both.dad.grandparents.vec<- sample.info_w_grandparents$both.dad.grandparents %>% unique()
  
  #Save dataframe of sampled aunt/uncles
  sampled_aunts_or_uncs.df <- sample.info_w_grandparents %>% dplyr::filter(both.rents %in%  both.mom.grandparents.vec | both.rents %in% both.dad.grandparents.vec) %>% #If this individual's parents are the grandparents of another sampled individual, then this indv is an aunt or uncle to another sampled indv.
    mutate(shared.relation = both.rents) %>% #Save the shared relation for joins
    dplyr::select(aunt.unc = indv.name, #Rename columns for joins
                  aunt.unc_mother = mother, 
                  aunt.unc_father = father, 
                  shared.relation,
                  aunt.unc_birth.year = birth.year,
                  aunt.unc_sex = sex)
  
  #Save dataframe of sampled nieces/nephews. Easier to split maternal and paternal and then join. First focus on maternal line.
  sampled_nieces_or_nephews_maternal.df <- sample.info_w_grandparents %>% dplyr::filter(both.mom.grandparents %in% both.parents.vec) %>% #If the indvidual's maternal grandparents are the parents to another sampled individual, then this is a niece/nephew
    mutate(shared.relation = both.mom.grandparents,
           lineage = "maternal") %>% #Save the shared relation for joins and for double-checking code
    dplyr::select(niece.nephew = indv.name, #Rename columns for joins
                  lineage,
                  shared.relation,
                  niece.nephew_mother = mother,
                  niece.nephew_father = father,
                  niece.nephew_maternal.grandmother= mom.grandmother, 
                  niece.nephew_maternal.grandfather = mom.grandfather,
                  niece.nephew_birth.year = birth.year,
                  niece.nephew_sex = sex)
  
  #Isolate paternal nieces and nephews, then join with maternal niece/nephews
  sampled_nieces_or_nephews_paternal.df <- sample.info_w_grandparents %>% 
    dplyr::filter(both.dad.grandparents %in% both.parents.vec) %>%  #If the individual's paternal grandparents are the parents to another sampled individual, then this is a niece/nephew
    mutate(shared.relation = both.dad.grandparents,
           lineage = "paternal") %>% #Save shared relation for join
    dplyr::select(niece.nephew = indv.name, #Rename for join
                  lineage,
                  shared.relation,
                  niece.nephew_mother = mother,
                  niece.nephew_father = father,
                  niece.nephew_paternal.grandmother = dad.grandmother, 
                  niece.nephew_paternal.grandfather = dad.grandfather,
                  niece.nephew_birth.year = birth.year,
                  niece.nephew_sex = sex) 
  
  all_sampled_nieces_or_nephews.df <- sampled_nieces_or_nephews_maternal.df %>% 
    full_join(sampled_nieces_or_nephews_paternal.df, by = c("niece.nephew", "shared.relation", "niece.nephew_birth.year", "niece.nephew_sex", "lineage", "niece.nephew_mother", "niece.nephew_father")) #Join with maternal nieces and nephews. If related through grandparents in the maternal line, then the paternal line will be NA and vice versa.
  
  #Combine niece/nephew dataframe with aunt/uncle dataframe
  aunt.unc_niece.nephew_pw.comps.all <- sampled_aunts_or_uncs.df %>% 
    full_join(all_sampled_nieces_or_nephews.df, by = c("shared.relation")) %>% #This is the key join. It makes a dataframe where we join by the shared relation, and then we should have the aunt/uncle and niece/nephew name and information in one dataframe, with the positives linked up.
    dplyr::filter(aunt.unc != niece.nephew_mother & aunt.unc != niece.nephew_father) %>% #Remove POPs
    mutate(aunt.unc_relation = as.character(ifelse(aunt.unc_sex == "F", "aunt", 
                                                   ifelse(aunt.unc_sex == "M", "uncle", NA))),
           niece.nephew_relation = as.character(ifelse(niece.nephew_sex == "F", "niece",
                                                       ifelse(niece.nephew_sex == "M", "nephew", NA))))
}



draw.samples <- function(target.samples = target.samples){
  
  for(i in sample.years){
    
    sample.df_temp <- NULL #Df with sample info for EACH sampling scheme for EACH year. So, reset at the beginning of each sampling year.
    
    num.YOY <- loopy.list[[i]] %>%
      dplyr::filter(age.x == 0) %>% 
      nrow()
    
    sample.size <- num.YOY * (sample.prop/100)
    
    #Sample YOY only for half-sib analysis
    sample.df_temp <- loopy.list[[i]] %>% mutate(capture.year = i) %>%
      dplyr::filter(age.x == 0) %>%
      dplyr::slice_sample(n = sample.size) %>%  #Sample each year WITHOUT replacement
      as_tibble()
    
    #Combine samples from all years
    sample.df_all.info_temp1 <- rbind(sample.df_all.info_temp1, sample.df_temp) #Df has info from EACH sampling scheme for ALL years
    
  } #End loop over sample years
  
  
  
  #Identify aunt|uncle / niece|nephew pairs. Will save to be imported later.
  aunt.unc_niece.nephew_pw.comps.all_temp1 <- find_aunt.uncle_nephew.niece_pairs(input.sample.info = sample.df_all.info_temp1, loopy.list = loopy.list) #Looks for aunt|uncle / nephew|niece pairs from EACH sampling scheme for EACH sample proportion over ALL sampling years
  
  return(invisible(list(sample.df_all.info_temp1, 
                        aunt.unc_niece.nephew_pw.comps.all_temp1,
                        sample.size)))
}



calculate.ref.year <- function(sample.info = sample.info){
  
  #Calculate the reference year for calculations of survival (which will vary depending on the samples)
  ref.year <- sample.info %>% dplyr::filter(age.x < repro.age) %>%
    arrange(birth.year) %>%
    slice_min(birth.year) %>%
    distinct(birth.year) %>%
    pull(birth.year)
  
}



####--------------------Split individual recaptures into first and last occasion------------------####
split.dups <- function(samples){
  
  #Double-checked code on 08/06/2023
  
  #Identify duplicates
  dups.first.capture <- samples %>% group_by(indv.name) %>% 
    filter(n() > 1) %>% 
    distinct(indv.name, .keep_all = TRUE) %>% 
    ungroup() 
  
  dups.later.capture <- samples %>% group_by(indv.name) %>% 
    filter(n() > 1) %>% 
    anti_join(dups.first.capture, by = c("indv.name", "capture.year")) %>% 
    ungroup()
  
  #Save the first capture instance of recaptured individuals and sort so older individual comes first (for pairwise comparison matrix)
  samples.save.first <- samples %>% anti_join(dups.later.capture, by = c("indv.name", "capture.year")) %>% 
    dplyr::arrange(birth.year) %>% 
    as_tibble()
  
  #Save the last capture instance of recaptured individuals and sort so older individual comes first (for pairwise comparison matrix)
  samples.save.second <- samples %>% anti_join(dups.first.capture, by = c("indv.name", "capture.year")) %>%
    dplyr::arrange(birth.year) %>% 
    as_tibble()
  
  return(list(samples.save.first, samples.save.second))  
}

####---------------Filter samples -------------####
#Filters for full sibs, but also identifies appropriate sample sets for each offspring birth year (for PO) and each younger sib birth year (for half-sib)
filter.samples <- function(samples){

  #Double-checked code on 08/06/2023
  
  #------------------Filter 1: full siblings----------------------
    NoFullSibs.df <- samples %>% distinct(mother.x, father.x, .keep_all = TRUE) %>% #CHANGE THIS if wanting to use just one indv per litter
    as_tibble() #If there is more than one individual with the same mother AND father, then only keep one. 
    total.sibs <- nrow(samples) - nrow(NoFullSibs.df) #Calculate number of full sibs
      
    dif.cohort.sibs <- samples %>% group_by(mother.x, father.x) %>%
      filter(n() > 1) %>% 
      distinct(birth.year) %>%
      ungroup() %>% 
      group_by(mother.x, father.x) %>% 
      summarize(distinct.birth.years = n()) %>% 
      filter(distinct.birth.years > 1) %>% 
      nrow()

    #The vector of years for which we need to split samples into potential parents and offspring i.e. offspring birth years.
  OffBirth.years <- NoFullSibs.df %>% 
    dplyr::filter(age.x <= repro.age) %>% #CHANGE FROM age.x == 0 (not sure why it was that?)
    distinct(birth.year) %>% 
    arrange(birth.year) %>%
    pull() 
  
  #Initialize lists to save output
  PO.samps.list <- list()
  HS.samps.list <- list()
  parents <- list()
  offspring <- list()
  

  #---------Filter 2: potential parents and offspring for each year------
  #For each offspring birth year, create a dataframe of potential parents i.e. sampled individuals that are old enough to have been a parent in that year.
  for(y in 1:length(OffBirth.years)){
    OffBirth.year <- OffBirth.years[y] #year of focus for this iteration
    
    #Dataframe of reproductively mature individuals during offspring birth year
    parents[[y]] <- NoFullSibs.df %>%
      mutate(age.in.OffBirth.year = age.x + (OffBirth.year - capture.year)) %>% 
      filter(age.in.OffBirth.year >= repro.age) %>% 
      mutate(relation = "parent")
    
    #Dataframe of offspring born in this year
    offspring[[y]] <- NoFullSibs.df %>% 
      mutate(age.in.OffBirth.year = age.x + (OffBirth.year - capture.year)) %>% #unnecessary line, I think
      filter(birth.year == OffBirth.year) %>%
      mutate(relation = "offspring")
      
    #Add dataframe of all potential offspring and parents to a list where each dataframe corresponds to an offspring birth year
    PO.samps.list[[y]] <- rbind(parents[[y]], offspring[[y]]) %>% 
      dplyr::select(indv.name, 
                    birth.year, 
                    capture.year,
                    age.in.OffBirth.year,
                    mother.x, 
                    father.x,
                    sex,
                    relation)
    
    #Not sure if this is needed, but add dataframe of potential half-sib samples for the year to a list where each dataframe corresponds to a specific birth year
    # HS.samps.list[[y]] <- NoFullSibs.df %>%
    #   filter(birth.year <= OffBirth.year) #Filters for offspring that were born before or in Offbirth.year
  }
  
  #Name list elements to correspond to sample year
  names(PO.samps.list) <- paste0("PO.samples.year_", c(OffBirth.years))
  #names(HS.samps.list) <- paste0("HS.samples.year_", c(OffBirth.years))
  
  #Not sure if this is needed, but create dataframes of samples for PO and HS comparisons. Think only need this for half-sibs
#  PO.samps.df <- bind_rows(PO.samps.list) #Confirmed that it's the correct number of rows
  HS.samps.df <- NoFullSibs.df %>% 
    dplyr::filter(age.x < repro.age) #Only use juveniles for half-sib model
    
  print(paste0("There are ", total.sibs, " pairs of full siblings in the dataset, ", dif.cohort.sibs, " of these pairs were born in different years."))
  return(list(PO.samps.list, HS.samps.df, NoFullSibs.df, total.sibs, dif.cohort.sibs)) #Return list of possible parents and offspring for each year, and dataframe of potential half-sibs
  
}



####--------------Build pairwise comparison matrices------------####
#Input is the filtered list from above for POs, and the filtered dataframe for HS

#Double-checked code 08/06/2023

build.pairwise <- function(filtered.samples.HS.df){

  HS.mom.pairwise.list <- list()
  HS.dad.pairwise.list <- list()
  
#--------------Half sibling--------------------------
filtered.samples.HS.df <- filtered.samples.HS.df %>% dplyr::arrange(birth.year) 
  #Arrange so older sib always comes first in pairwise comparison matrix
pairwise.df.HS <- data.frame(t(combn(filtered.samples.HS.df$indv.name, m=2))) %>% 
  as_tibble() # generate all combinations of the elements of x, taken m at a time.
colnames(pairwise.df.HS) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined

#Create dataframe of pairwise comparisons with just individual IDs
head(pairwise.df.HS)
#head(pairwise.df)

#Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
OlderSib_birthyears.HS <- filtered.samples.HS.df %>%
  dplyr::select(older.sib = indv.name, 
         older.sib.birth = birth.year, 
         older.sib.age = age.x, 
         older.sib.mom = mother.x,
         older.sib.dad = father.x)

#Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
#This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
pairwise.df_all.info.HS <- pairwise.df.HS %>% left_join(OlderSib_birthyears.HS, by = "older.sib") %>% 
  as_tibble() %>%  
  left_join(filtered.samples.HS.df, by = "indv.name") %>% 
  dplyr::rename("younger.sib" = indv.name, 
                "younger.sib.birth" = birth.year, 
                "younger.sib.age" = age.x, 
                "younger.sib.mom" = mother.x, 
                "younger.sib.dad" = father.x) %>%
  dplyr::select(older.sib, 
         older.sib.birth, 
         older.sib.age, 
         younger.sib, 
         younger.sib.birth, 
         younger.sib.age, 
         older.sib.mom, 
         younger.sib.mom, 
         older.sib.dad, 
         younger.sib.dad)

pairwise.df_HS.filt <- pairwise.df_all.info.HS %>% 
  dplyr::filter(older.sib.birth != younger.sib.birth) %>% #Filter intra-cohort comparisons (for now)
  dplyr::filter((younger.sib.birth - older.sib.birth) <= (max.age - repro.age))
  

#Extract positive half-sib comparisons
positives.HS <- pairwise.df_HS.filt %>% filter(older.sib.mom == younger.sib.mom | older.sib.dad == younger.sib.dad) %>% 
  mutate(shared.parent = ifelse(older.sib.mom == younger.sib.mom, "mother", "father"))
#nrow(positives)



####----------------Split dataframes into final form for model----------####
#Sex-specific half-sib
mom_positives.HS <- positives.HS %>% filter(shared.parent == "mother") %>% 
  dplyr::select(older.sib.birth, younger.sib.birth) %>% 
  plyr::count() %>% 
  rename(yes = freq)

dad_positives.HS <- positives.HS %>% filter(shared.parent == "father")  %>%
  dplyr::select(older.sib.birth, younger.sib.birth) %>%
  plyr::count() %>% 
  rename(yes = freq)


#-------- Make dataframes for negative comparisons -- Discount possibility of finding half-sibs when there are full sibs ------
 hsp.negs <- pairwise.df_HS.filt %>%
   dplyr::filter(younger.sib.mom != older.sib.mom & younger.sib.dad != older.sib.dad) %>%
   dplyr::select(older.sib.birth, younger.sib.birth) %>%
   plyr::count() %>%
   as_tibble()

 mom_comps.HS <- hsp.negs %>%
   dplyr::rename(no = freq) %>%
   left_join(mom_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>%
   mutate(yes = replace_na(yes, 0)) %>%
   mutate(year_gap = younger.sib.birth - older.sib.birth) %>%
   mutate(all = yes + no)

 dad_comps.HS <- hsp.negs %>%
   dplyr::rename(no = freq) %>%
   left_join(dad_positives.HS, by = c("older.sib.birth", "younger.sib.birth")) %>%
   mutate(yes = replace_na(yes, 0)) %>%
   mutate(year_gap = younger.sib.birth - older.sib.birth) %>%
   mutate(all = yes + no)


HS.mom.pairwise.df <- mom_comps.HS %>%
  dplyr::rename(ref.year = younger.sib.birth,
                mort.yrs = year_gap) %>% 
  dplyr::select(ref.year, all, yes, mort.yrs) %>% 
  #dplyr::filter(mort.yrs < (2*repro.age)) %>% #Added 04/01/2022 to remove comparisons that could be confused for grand-parents
  mutate(type = "HS",
         parent = "mother")

HS.dad.pairwise.df <- dad_comps.HS %>%
  dplyr::rename(ref.year = younger.sib.birth,
                mort.yrs = year_gap) %>% 
  dplyr::select(ref.year, all, yes, mort.yrs) %>% 
  #dplyr::filter(mort.yrs < (2*repro.age)) %>% #Added 04/01/2022 to remove comparisons that could be confused for grand-parents
  mutate(type = "HS",
         parent = "father")

#SWITCHED ref.year and estimation.year in pop growth calculation
#Combine PO and HS dataframes for each sex
mom_comps.all <- HS.mom.pairwise.df %>% 
  mutate(pop.growth.yrs = ref.year - estimation.year) %>% 
  arrange(desc(ref.year), mort.yrs) %>% 
  mutate(BI = ifelse(mort.yrs %% mating.periodicity == 0, "on", "off"))

dad_comps.all <- HS.dad.pairwise.df %>% 
  mutate(pop.growth.yrs = ref.year - estimation.year) %>% 
  arrange(desc(ref.year), mort.yrs)

return(list(mom_comps.all, dad_comps.all, positives.HS))
}


#----------------Calculate number of HSPs for downsampling---------------
downsample <- function(mom_comps.all, dad_comps.all, HS.samps.df, PO.samps.list){

  #Downsample for HSPs
#Calculate proportion of samples born in each birth year
  HS.samps.props <- HS.samps.df %>% group_by(birth.year) %>%
    summarize(num = n()) %>%
  ungroup() %>%
  mutate(prop = num/sum(num))

  #Save vector of proportions corresponding to the correct index in the HS sample dataframe
HS.props.vec <- HS.samps.df %>% left_join(HS.samps.props, by = "birth.year") %>%
  pull(prop)

  #Calculate proportion of positive HS comparisons for mom
HS_comps.mom.yes <- mom_comps.all %>% dplyr::filter(type == "HS") %>%
  summarize(yes = sum(yes)) %>% 
  pull(yes)

HS_comps.mom.all <- mom_comps.all %>% dplyr::filter(type == "HS") %>%
  summarize(all = sum(all)) %>% 
  pull(all)

HS_prop.mom.yes <- HS_comps.mom.yes/HS_comps.mom.all #Calculate proportion of positive comparisons

#Calculate proportion of positive HS comparisons for dad
HS_comps.dad.yes <- dad_comps.all %>% dplyr::filter(type == "HS") %>%
  summarize(yes = sum(yes)) %>% 
  pull(yes)

HS_comps.dad.all <- dad_comps.all %>% dplyr::filter(type == "HS") %>%
  summarize(all = sum(all)) %>% 
  pull(all)

HS_prop.dad.yes <- HS_comps.dad.yes/HS_comps.dad.all

HS_prop.mean.yes <- mean(c(HS_prop.mom.yes, HS_prop.dad.yes)) 
HS_target.comps = round(max.HSPs/HS_prop.mean.yes, 0)
HS_target.samples <- round(sqrt(HS_target.comps), 0)

HS.samps.df.down <- NULL #Initialize so I don't get an error
if(HS_comps.mom.yes + HS_comps.dad.yes > max.HSPs){
  HS.samps.df.down <- HS.samps.df %>% slice_sample(n = HS_target.samples, weight_by = HS.props.vec)
} else{
  HS.samps.df.down <- HS.samps.df
}


#-----------------------Downsample for POPs-------------------------------------#
PO_comps.mom.yes <- mom_comps.all %>% dplyr::filter(type == "PO") %>%
  summarize(yes = sum(yes)) %>% 
  pull(yes)

PO_comps.mom.all <- mom_comps.all %>% dplyr::filter(type == "PO") %>%
  summarize(all = sum(all)) %>% 
  pull(all)

PO_prop.mom.yes <- PO_comps.mom.yes/PO_comps.mom.all #Calculate proportion of positive comparisons

#Calculate proportion of positive PO comparisons for dad
PO_comps.dad.yes <- dad_comps.all %>% dplyr::filter(type == "PO") %>%
  summarize(yes = sum(yes)) %>% 
  pull(yes)

PO_comps.dad.all <- dad_comps.all %>% dplyr::filter(type == "PO") %>%
  summarize(all = sum(all)) %>% 
  pull(all)

PO_prop.dad.yes <- PO_comps.dad.yes/PO_comps.dad.all

PO.all.yes <- PO_comps.dad.yes + PO_comps.mom.yes

#This doesn't work for POPs the same way it does for half sibs. So, just streamline it and cut number of samples down proportionally
# PO_prop.mean.yes <- mean(c(PO_prop.mom.yes, PO_prop.dad.yes)) 
# PO_target.comps = round(max.POPs/PO_prop.mean.yes, 0)
# PO_target.samples <- round(sqrt(PO_target.comps), 0)
# PO_target.samples.per.yr <- round(PO_target.samples/length(sample.years), 0)

#sample.size.rents.reduced <- (max.POPs/PO.all.yes)*sample.size.rents
#sample.size.juvs.reduced <- (max.POPs/PO.all.yes)*sample.size.juvs
#PO.down.prop <- max.POPs/PO.all.yes

PO.samps.list.down <- list()
if(PO_comps.mom.yes + PO_comps.dad.yes > max.POPs){

  #---------Filter 2: potential parents and offspring for each year------
#For each offspring birth year, create a dataframe of potential parents i.e. sampled individuals that are old enough to have been a parent in that year.
for(y in 1:length(PO.samps.list)){
  
  off.down <- PO.samps.list[[y]] %>% dplyr::filter(relation == "offspring") %>% 
    slice_sample(prop = PO.down.prop)
  
  rents.down <- PO.samps.list[[y]] %>% dplyr::filter(relation == "parent") %>% 
    slice_sample(prop = PO.down.prop)
  
  PO.samps.list.down[[y]] <- rbind(off.down, rents.down)
  
  } #End creation of PO sample list

names(PO.samps.list.down) <- paste0("PO.samples.year_", c(sample.years))

} else {PO.samps.list.down <- PO.samps.list} #End PO downsample

return(list(HS.samps.df.down,
            PO.samps.list.down))
}



calc.psi <- function(loopy.list){
  #Calculate breeding interval
  #Initialize sample dataframes
  BI.df <- NULL
  BI.df_temp <- NULL
  
  ref.years <- c(50:n_yrs)
  
  for(i in ref.years){ #Extract all YOY for each year sampled
    BI.df_temp <- loopy.list[[i]] %>% mutate(capture.year = i) %>% 
      dplyr::filter(age.x == 0)
    
    #Combine all
    BI.df <- bind_rows(BI.df, BI.df_temp)
  }
  
  #Extract unique mothers for each birth year
  BI.df <- BI.df %>% as_tibble() %>% 
    arrange(birth.year) %>% 
    distinct(mother.x, birth.year, .keep_all = TRUE)
  
  #-----------------Identify all relatives in population since estimation year and find breeding interval--------------------#
  pairwise.BI.df <- data.frame(t(combn(BI.df$indv.name, m=2)))  # generate all combinations of the elements of x, taken m at a time.
  colnames(pairwise.BI.df) <- c("older.sib", "indv.name") #Rename columns so they can easily be joined
  
  #Create dataframe of pairwise comparisons with just individual IDs
  head(BI.df)
  #head(pairwise.df)
  
  #Create dataframe that will be used to extract the birth years for the younger fish from each pairwise comparison using joins.
  OlderSib_birthyears.BI <- BI.df %>%
    dplyr::select(older.sib = indv.name, 
           older.sib.birth = birth.year, 
           older.sib.age = age.x, 
           older.sib.mom = mother.x,
           older.sib.dad = father.x)
  
  #Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
  #This is the main pairwise comparison matrix with all (relevant) comparisons and individual data.###
  pairwise.BI.all.info <- pairwise.BI.df %>% left_join(OlderSib_birthyears.BI, by = "older.sib") %>% 
    as_tibble() %>%  
    left_join(BI.df, by = "indv.name") %>% 
    dplyr::rename("younger.sib" = indv.name, 
                  "younger.sib.birth" = birth.year, 
                  "younger.sib.age" = age.x, 
                  "younger.sib.mom" = mother.x, 
                  "younger.sib.dad" = father.x) %>%
    dplyr::select(older.sib, 
                  older.sib.birth, 
                  older.sib.age, 
                  younger.sib, 
                  younger.sib.birth, 
                  younger.sib.age, 
                  older.sib.mom, 
                  younger.sib.mom, 
                  older.sib.dad, 
                  younger.sib.dad) %>% 
    dplyr::filter(older.sib.birth != younger.sib.birth)
  
  
  #Extract all positive maternal half-sib comparisons from sampled years
  positives.mom.BI <- pairwise.BI.all.info %>% filter(older.sib.mom == younger.sib.mom) %>% 
    mutate(yr.gap = younger.sib.birth - older.sib.birth) %>% 
    mutate(BI = ifelse(yr.gap %% mating.periodicity == 0, "on", "off")) %>% #If the year gap between ha;f-sibs is divisible by 2, then call the breeding interval (BI) "even"
    dplyr::distinct(older.sib.birth, younger.sib.birth, younger.sib.mom, .keep_all = TRUE) #Duplicative of earlier; making sure cohort size isn't biasing things here
  
  #Everything matches up re: positives and 
  positives.mom.BI %>% group_by(BI) %>% summarize(n()) #How many positive comparisons for 
  
  skipped.moms <- positives.mom.BI %>% dplyr::filter(BI == "on") %>% 
    dplyr::distinct(younger.sib.mom) %>% 
    pull(younger.sib.mom) #Should be the same as the older sib mom
  
  annual.moms <- positives.mom.BI %>% dplyr::filter(BI == "odd") %>% 
    dplyr::distinct(younger.sib.mom) %>% 
    pull(younger.sib.mom) #Should be the same as the older sib mom
  
  #Which moms ONLY reproduced in even years?
  skipped.only.moms <- skipped.moms[which(!skipped.moms %in% annual.moms)]
  all.moms <- unique(positives.mom.BI$younger.sib.mom)
  
  #Percent of individuals that only breed when years are evenly spaced
  (psi.truth <- length(skipped.only.moms)/length(all.moms))
  
}







#----------------Calculate expected matches-------------------
calc.Exp <- function(mom_comps.all, dad_comps.all){
#-----------Calculate C and expected HSPs-------------------
#Calculate Exp(R) for MHSPs
mom_C.HS <- mom_comps.all %>% dplyr::filter(type == "HS") #Filter for HS relationships

#Summarize number of distinct mortality years
mom_age.gap.denom <- mom_C.HS %>% distinct(mort.yrs) %>% 
  summarize(s = sum(mort.yrs)) %>% 
  pull(s)

#What's the true abundance?
# mom_N.truth <- pop.size.tibble %>% filter(year == estimation.year) %>% 
#   distinct(Num.mothers) %>% 
#   pull(Num.mothers)

#Number of comparisons for each year gap (mort.yrs)
mom_comps.yr.df <- mom_C.HS %>% group_by(mort.yrs) %>% 
  summarize(comps.yr = sum(all))

#Calculate expected MHSPs 
mom_C.HS.Exp <- mom_C.HS %>% 
  group_by(mort.yrs) %>% 
  summarize(cross.cohort.instances = n()) %>% 
  mutate(cum.surv = adult.survival ^ mort.yrs) %>% 
  mutate(Ct_y = cross.cohort.instances * cum.surv) %>% 
  mutate(C_yr = Ct_y/mom_age.gap.denom) %>% 
  left_join(mom_comps.yr.df, by = "mort.yrs") %>% 
  mutate(Exp.R = (cum.surv*comps.yr)/Nf.truth)

mom_weighted.mean.C <- sum(mom_C.HS.Exp$C_yr)
mom_total.comps <- sum(mom_comps.yr.df$comps.yr)

#(mom.Exp.HS <- round((mom_weighted.mean.C * mom_total.comps)/mom_N.truth, 0))
mom.Exp.HS <- round(sum(mom_C.HS.Exp$Exp.R), 0)


#Calculate Exp(R) for MPOPs
mom_C.PO <- mom_comps.all %>% dplyr::filter(type == "PO")

#Need to use Na to get the correct number of expected POs 
# Nf.truth <- pop.size.tibble %>% filter(year == 85) %>% 
#   distinct(Female.adult.pop) %>% 
#   pull(Female.adult.pop)

mom.Exp.PO <- round(sum(mom_C.PO$all)/Nf.truth, 0)

#---------------------Dad--------------------------------
#Calculate Exp(R) for MHSPs
dad_C.HS <- dad_comps.all %>% dplyr::filter(type == "HS")

#Summarize umber of distinct mortality years
dad_age.gap.denom <- dad_C.HS %>% distinct(mort.yrs) %>% 
  summarize(s = sum(mort.yrs)) %>% 
  pull(s)

#What's the true abundance?
# dad_N.truth <- pop.size.tibble %>% filter(year == 85) %>% 
#   distinct(Male.adult.pop) %>% 
#   pull(Male.adult.pop)

#Number of comparisons for each year gap (mort.yrs)
dad_comps.yr.df <- dad_C.HS %>% group_by(mort.yrs) %>% 
  summarize(comps.yr = sum(all))

#Calculate expected FHSPs 
dad_C.HS.Exp <- dad_C.HS %>% 
  group_by(mort.yrs) %>% 
  summarize(cross.cohort.instances = n()) %>% 
  mutate(cum.surv = adult.survival ^ mort.yrs) %>% 
  mutate(Ct_y = cross.cohort.instances * cum.surv) %>% 
  mutate(C_yr = Ct_y/dad_age.gap.denom) %>% 
  left_join(dad_comps.yr.df, by = "mort.yrs") %>% 
  mutate(Exp.R = (cum.surv*comps.yr)/Nm.truth)

dad_weighted.mean.C <- sum(dad_C.HS.Exp$C_yr)
dad_total.comps <- sum(dad_comps.yr.df$comps.yr)

#(dad.Exp.HS <- round((dad_weighted.mean.C * dad_total.comps)/dad_N.truth, 0))
dad.Exp.HS <- round(sum(dad_C.HS.Exp$Exp.R), 0) #Sum expected HAPs per year


#Calculate Exp(R) for MPOPs
dad_C.PO <- dad_comps.all %>% dplyr::filter(type == "PO")

dad.Exp.PO <- round(sum(dad_C.PO$all)/Nm.truth, 0)

return(list(mom.Exp.HS, mom.Exp.PO, dad.Exp.HS, dad.Exp.PO))
}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}


#Need to save aunt.unc_niece.nephew_pw.comps.all as RDS object in Population simulation script.Then import into model estimation script, and use function below to reduce the dataframe to charlatan HSPs.

#GOT EM
identify_imposters <- function(imposters.df = imposters.df){

   charlatan_HSPs <- imposters.df %>%
     dplyr::filter(iteration == iter,
                   sampling.scheme == s.scheme,
                   sample.prop == sample.proportion) %>%
     mutate(shared.parent = as.character(ifelse(lineage == "maternal", "mother", 
                                   ifelse(lineage == "paternal", "father", NA))),
            older.sib.age = as.numeric(NA), #Need this column for row bind, but it's not used for any calculation
            younger.sib.age = as.numeric(NA) #Need this column for row bind, but it's not used for any calculation
     ) %>% 
     dplyr::select(older.sib = aunt.unc,
                   older.sib.birth = aunt.unc_birth.year,
                   older.sib.age,
                   younger.sib = niece.nephew,
                   younger.sib.birth.year = niece.nephew_birth.year,
                   younger.sib.age,
                   older.sib.mom = aunt.unc_mother,
                   younger.sib.mom = aunt.unc_mother, #Not true, but need to make sure we get the right number of positives and negatives, so need this to be the same for the aunt|uncle and niece|nephew.
                   older.sib.dad = aunt.unc_father,
                   younger.sib.dad = aunt.unc_father, #Not true, but need to make sure we get the right number of positives and negatives, so need this to be the same for the aunt|uncle and niece|nephew.
                   older.sib.birth = aunt.unc_birth.year,
                   younger.sib.birth = niece.nephew_birth.year,
                   shared.parent)
   
   if(filter.decoys == "yes"){
     
     #See what would happen if we filtered them by age difference
     charlatan_HSPs <- charlatan_HSPs %>% dplyr::filter(younger.sib.birth - older.sib.birth < year.gap.threshold)

   }
   
  return(charlatan_HSPs)

}






calc.truth <- function(truth.df = truth.df){
  #If there's no lambda in the model, then the true value of abundance is the mean over the estimated years. Here, we make this correction for the scenarios that do not include lambda in the model; for those that do, the values in truth.df are correct.
  if(scenario == "scenario_1_model.validation" | 
     scenario == "scenario_1.2.1" | 
     scenario == "scenario_1.2.2" | 
     scenario == "scenario_1.2.3"){
    
    #Make truth equal to mean over years
    #est.year.calibrate is the first year of data aka T0. This is different than the object estimation.year, which is what we loop over and change.
    Nf.truth <- pop.size.tibble %>% dplyr::filter(year >= est.year.calibrate) %>% 
      summarize(mean.female.pop = mean(Female.adult.pop)) %>% 
      pull(mean.female.pop)
    
    Nfb.truth <- pop.size.tibble %>% dplyr::filter(year >= est.year.calibrate) %>% 
      summarize(mean.mothers = mean(Num.mothers)) %>% 
      pull(mean.mothers)
    
    Nm.truth <- pop.size.tibble %>% dplyr::filter(year >= est.year.calibrate) %>% 
      summarize(mean.male.pop = mean(Male.adult.pop)) %>% 
      pull(mean.male.pop)
    
    Nmb.truth <- pop.size.tibble %>% dplyr::filter(year >= est.year.calibrate) %>% 
      summarize(mean.fathers = mean(Num.fathers)) %>% 
      pull(mean.fathers)
    
    truth.iter <- truth.iter %>% mutate(all.truth = ifelse(parameter == "Nf", Nf.truth, 
                                                           ifelse(parameter == "Nfb1" | parameter == "Nfb2", Nfb.truth, 
                                                                  ifelse(parameter == "Nm", Nm.truth, 
                                                                         ifelse(parameter == "Nmb1" | parameter == "Nmb2", Nmb.truth, all.truth)))))
    cat("\nUsing average abundance for truth\n")
  } else cat("\nUsing year-specific truth, rather than average.\n")
  
  
  
  (truth.iter_all.params <- truth.df %>% dplyr::filter(sampling.scheme == s.scheme,
                                                       iteration == iter,
                                                       sample.prop == sample.proportion) %>% 
      mutate(estimation.year = estimation.year) %>% 
      dplyr::mutate(T0 = ref.yr + 1) %>% 
      dplyr::select(-c(surv_min, surv_max, population.growth)) %>% 
      bind_rows(truth.iter) %>% 
      dplyr::select(parameter, all.truth, estimation.year, T0, iteration, sampling.scheme, sample.prop, seed) %>% #Change order of columns
      dplyr::arrange(parameter))
  
  
  results.temp <- model.summary2 %>% left_join(truth.iter_all.params, by = c("parameter", "iteration", "seed")) %>% 
    dplyr::select(parameter,
                  Q50,
                  all.truth,
                  HPD2.5,
                  HPD97.5,
                  mean,
                  sd,
                  Q2.5,
                  Q97.5,
                  Rhat,
                  estimation.year,
                  T0,
                  sampling.scheme,
                  iteration,
                  neff,
                  sample.prop,
                  seed) %>% 
    mutate(estimation.sim = ifelse(est == 1, "T0", 
                                   ifelse(est == 2, "T0-10",
                                          ifelse(est == 3, "present-5",
                                                 ifelse(est == 4, "present", NA)))))

  
  if(HS.only == "yes"){
    
    results.temp <- results.temp %>% 
      mutate(HSPs_detected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.HSPs, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.HSPs, mom.HSPs + dad.HSPs)),
             HSPs_expected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.Exp.HS, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.Exp.HS, mom.Exp.HS + dad.Exp.HS)),
             POPs_detected = NA,
             POPs_expected = NA,
             scenario = scenario)
    
  } else if(HS.only != "yes"){
    
    results.temp <- results.temp %>% 
      mutate(HSPs_detected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.HSPs, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.HSPs, mom.HSPs + dad.HSPs)),
             HSPs_expected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.Exp.HS, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.Exp.HS, mom.Exp.HS + dad.Exp.HS)),
             POPs_detected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.POPs, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.POPs, mom.POPs + dad.POPs)),
             POPs_expected = ifelse(parameter %in% c("Nf", "Nfb1", "Nfb2", "psi"), mom.Exp.PO, 
                                    ifelse(parameter %in% c("Nm", "Nmb1", "Nmb2"), dad.Exp.PO, mom.Exp.PO + dad.Exp.PO)),
             scenario = scenario)
    
  }
  
  return(results.temp)
}
