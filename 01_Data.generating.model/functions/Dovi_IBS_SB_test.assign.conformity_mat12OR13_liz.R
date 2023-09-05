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
  
  if(population.growth == "lambda.extreme" & v >= n_yrs-10){
    Adult.survival <- Sa
    juvenile.survival <- Sj
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





###------------------ Create function to identify aunt|uncle/niece|nephew pairs----------------
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
                  capture.year, 
                  sampling.scheme) %>% #Keep most columns but remove the x from mother and father
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
                  aunt.unc_sex = sex,
                  sampling.scheme)
  
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
                  niece.nephew_sex = sex,
                  sampling.scheme)
  
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
                  niece.nephew_sex = sex,
                  sampling.scheme) 
  
  all_sampled_nieces_or_nephews.df <- sampled_nieces_or_nephews_maternal.df %>% 
    full_join(sampled_nieces_or_nephews_paternal.df, by = c("niece.nephew", "shared.relation", "niece.nephew_birth.year", "niece.nephew_sex", "sampling.scheme", "lineage", "niece.nephew_mother", "niece.nephew_father")) #Join with maternal nieces and nephews. If related through grandparents in the maternal line, then the paternal line will be NA and vice versa.
  
  #Combine niece/nephew dataframe with aunt/uncle dataframe
  aunt.unc_niece.nephew_pw.comps.all <- sampled_aunts_or_uncs.df %>% 
    full_join(all_sampled_nieces_or_nephews.df, by = c("shared.relation", "sampling.scheme")) %>% #This is the key join. It makes a dataframe where we join by the shared relation, and then we should have the aunt/uncle and niece/nephew name and information in one dataframe, with the positives linked up.
    dplyr::filter(aunt.unc != niece.nephew_mother & aunt.unc != niece.nephew_father) %>% #Remove POPs
    mutate(aunt.unc_relation = as.character(ifelse(aunt.unc_sex == "F", "aunt", 
                                      ifelse(aunt.unc_sex == "M", "uncle", NA))),
           niece.nephew_relation = as.character(ifelse(niece.nephew_sex == "F", "niece",
                                          ifelse(niece.nephew_sex == "M", "nephew", NA))))
}



draw.samples <- function(target.samples = target.samples){

  for(i in sample.years){
  
    sample.df_temp <- NULL #Df with sample info for EACH sampling scheme for EACH year. So, reset at the beginning of each sampling year.
    
  #Extract the relevant row from the pop size dataframe
  pop.size.yr <- pop.size.tibble %>% dplyr::filter(year == i)
  
  sample.size <- pop.size.yr %>%
    mutate(sample.size = round(population_size*(sample.prop/100), 0)) %>%
    pull(sample.size)
  
  if(target.samples == "target.YOY"){ #If targeting YOY for juvenile samples
    
    #Sample YOY only for half-sib analysis
    sample.df_temp <- loopy.list[[i]] %>% mutate(capture.year = i,
                                                 sampling.scheme = target.samples) %>%
      dplyr::filter(age.x == 0) %>%
      dplyr::slice_sample(n = sample.size) %>%  #Sample each year WITHOUT replacement
      as_tibble()
    
  } else if(target.samples == "sample.all.juvenile.ages"){ #If sampling juveniles
    sample.df_temp <- loopy.list[[i]] %>% dplyr::filter(age.x < repro.age & age.x > 0) %>% 
      mutate(capture.year = i,
             sampling.scheme = target.samples) %>%
      dplyr::slice_sample(n = sample.size) %>% #Sample each year WITHOUT replacement (doesn't affect cross-year sampling since it's in a loop)
      as_tibble()
    
  } else if(target.samples == "sample.ALL.ages"){ #If sampling all ages
    
    sample.df_temp <- loopy.list[[i]] %>% mutate(capture.year = i,
                                                 sampling.scheme = target.samples) %>%
      dplyr::slice_sample(n = sample.size) %>% #Sample each year WITHOUT replacement (doesn't affect cross-year sampling since it's in a loop)
      as_tibble()
    
    
  } #End conditional statement and sampling for this year/sampling scheme
  
  #Combine samples from all years
  sample.df_all.info_temp1 <- rbind(sample.df_all.info_temp1, sample.df_temp) #Df has info from EACH sampling scheme for ALL years
  
  } #End loop over sample years
 
  
  
  #Identify aunt|uncle / niece|nephew pairs. Will save to be imported later.
  aunt.unc_niece.nephew_pw.comps.all_temp1 <- find_aunt.uncle_nephew.niece_pairs(input.sample.info = sample.df_all.info_temp1, loopy.list = loopy.list) #Looks for aunt|uncle / nephew|niece pairs from EACH sampling scheme for EACH sample proportion over ALL sampling years
   
  return(invisible(list(sample.df_all.info_temp1, 
                        aunt.unc_niece.nephew_pw.comps.all_temp1,
                        sample.size)))
}

calculate.ref.year <- function(){
  
  #Calculate the reference year for calculations of survival (which will vary depending on the samples)
  ref.year <- sample.df_all.info_temp1 %>% dplyr::filter(age.x < repro.age) %>%
    arrange(birth.year) %>%
    slice_min(birth.year) %>%
    distinct(birth.year) %>%
    pull(birth.year)
  
  ref.temp <- tibble(sampling.scheme = target.samples, 
                     sample.prop = sample.prop,
                     ref.yr = ref.year,
                     iteration = iter)
  
  ref.tibble <- bind_rows(ref.tibble, ref.temp)
  
}