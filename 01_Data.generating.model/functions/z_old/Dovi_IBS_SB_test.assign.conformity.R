simulate.pop <- function(init.pop.size, init.prop.female, Nages, mating.periodicity, repro.age, YOY.survival, juvenile.survival, Adult.survival, max.age, num.mates, ff, burn.in, Num.years){

  ####---------Set up initial population-----####
init.pop <- data.frame() # create a blank data frame which will become the initial population

age.x <- NULL
for(y in 1:length(Nages)) age.x <- c(age.x, rep(y, Nages[y]))

repro.cycle.vec_init <- c(rep(NA, length.out = length(2:repro.age)), #Use 2 bc we don't want to include age 0 and we don't want to double count the age at reproduction. This will give us the right number for our code below to correctly draw from the vector.
                          rep(c(1:mating.periodicity), length.out = length(repro.age:max.age)))

for(i in 1:init.pop.size) { # Loop that creates the below data for each individual in a population the size of "init.pop.size"
  
  indv.name <- paste(sample(letters, size = 20, replace = T), collapse="") # generate a random name that is 20 letters long
  
  #age.x <- repro.age # every individual in the initial population is the age at first maturity ### WHY? WOULDN'T IT BE BETTER TO START WITH ~ STABLE AGE DISTRIBUTION?
  
  mother.x <- "xxxxx" # The individuals in the initial population do not have known mothers
  father.x <- "xxxxx" # The individuals in the initial population do not have known fathers
  birth.year <- -1 # this is a place holder for individuals born within the simulation
  sex <- sample(c('F','M'), size = 1, prob = c(init.prop.female, 1-init.prop.female)) # randomly draw sex based on the prportions set in the parameter section
  
  init.vec <- cbind.data.frame(indv.name, birth.year, age.x[i], mother.x, father.x, sex) # CHANGED THIS TO AGE.X[I]; create a row with the attributes for each individual - using cbind.data.frame allows them to keep some numeric and others characters
  init.pop <- rbind(init.pop, init.vec) # place the individual into the init.pop data frame
  
#  print(paste0("Finished iteration ", i, " of ", init.pop.size))
}

head(init.pop)
names(init.pop)[3] <- "age.x"
summary(as.factor(init.pop$age.x))
summary(as.factor(init.pop$sex))
n.mothers_init <- init.pop %>% dplyr::filter(age.x >= repro.age & sex == "F") %>% 
  nrow()

#Confirmed that this assigns conformists and non-conformists at the appropriate ratio; focus on age.x == 0 bc it accounts for having one less element in the vector.
init.pop <- init.pop %>% mutate(repro.cycle = ifelse(age.x == 0, 
                                                     NA,
                                                     repro.cycle.vec_init[age.x]),
                                repro.strategy = ifelse(sex == "F" & age.x >= repro.age,
                                                        sample(c("conformist", "non-conformist"),
                                                               size = n.mothers_init,
                                                               prob = c(1 - non.conformists, non.conformists),
                                                               replace = TRUE),
                                                        NA))


####----------Breeding----------####
repro.cycle.vec <- rep(1:mating.periodicity, length.out = (n_yrs+1)) #+1 for year 0 i.e. initial population
#repro.cycle.vec <- rep(1:mating.periodicity, times = 100) # Generate a vector which will be used to determine if it is an even or odd breeding year (or a 1/3 breeding year)

####---------for year 0 breeding-------####
#------------------------------------Conformist mothers------------------------------------#
mothers <- which(init.pop$sex=='F' & init.pop$age.x >= repro.age & init.pop$repro.cycle == repro.cycle.vec[1] & init.pop$repro.strategy == "conformist") # determine which females are available to breed in this year

#-----------------------------------Non-conformist mothers------------------------------------#
mothers.2 <- which(init.pop$sex == 'F' & init.pop$age.x >= repro.age & init.pop$repro.strategy == "non-conformist") #Assume for now that non-conformist means annual breeding

#Combine conformist and non-conformist mothers
mothers <- sort(c(mothers, mothers.2))

#-----------------------------------Fathers--------------------------------------------#
fathers <- which(init.pop$sex=='M' & init.pop$age.x>=repro.age) # determine which fathers are available to breed in this year



YOY.df <- data.frame() # create an empty data frame to populate with the YOY born in this year

for(j in 1:length(mothers)) { # Loop through all of the available mothers
  num.mates.x <- sample(num.mates, size=1) # determine the number of mates for a given mother this year
  inter.df <- data.frame() # Create a data frame for the offspring born from each mother
  
  for(h in 1:num.mates.x) { # Loop through each sire that this mother mates with
    
    num.offspring <- rpois(1, ff) #generate the number of offspring born from this mother/sire pairing
    indv.name <- NA # create a place holder for the random name given to each offspring
    age.x <- 0 # assign a 0 age to each offspring born this year
    mother.x <- init.pop[mothers[j],1] # record the mothers name
    father.x <- init.pop[sample(fathers, size = 1),1] # record the sires name
    birth.year <- 0 # note that these pups were born in the 0 year of the simulation
    sex <- NA # create a place holder for the sexes
    
    repro.cycle <- NA # assign the newborn to an eventual breeding cycle group
    repro.strategy <- NA #assign the newborn to be a conformist or non-conformist if a female
      
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
  
  #Re-assign fecundity to produce a population with slight positive or negative growth.
  if(v > burn.in & exists("ff.shift") == TRUE){
    ff <- ff.shift
  }
  
  #Re-assign survival to produce substantial population decline.
  if(population.growth == "lambda.extreme" & v >= n_yrs-10){
    Adult.survival <- Sa
    juvenile.survival <- Sj
  }
  
  
  data1 <- loopy.pop[loopy.pop$Survival =="S", -9] #Bring in the data from the previous iteration, but only include those that survive (and leave the column of survival out)
  data1$age.x <- data1$age.x+1 # increase each individuals age by one for the new year - happy birthday survivors!
  
data1 <- data1 %>% mutate(repro.cycle = ifelse(sex == "F" & 
                                          age.x == repro.age &
                                          is.na(repro.cycle) == TRUE, 
                                        repro.cycle.vec[v+1],
                                        repro.cycle),
                               repro.strategy = ifelse(sex == "F" & 
                                                         age.x == repro.age &
                                                         is.na(repro.strategy) == TRUE,
                                                       sample(c("conformist", "non-conformist"),
                                                              size = n.mothers_init,
                                                              prob = c(1 - non.conformists, non.conformists),
                                                              replace = TRUE),
                                                       repro.strategy))

  
  #------------------------------------Conformist mothers------------------------------------#
  #Add all mothers that are meant to breed in the year, plus a small proportion that are not
  mothers <- which(data1$sex=='F' & data1$age.x>=repro.age & data1$repro.cycle == repro.cycle.vec[v+1] & data1$repro.strategy == "conformist") #determine which females are available to breed in this year; this is an index

  
  #-----------------------------------Non-conformist mothers------------------------------------#
  mothers.2 <- which(data1$sex=='F' & data1$age.x>=repro.age & data1$repro.strategy == "non-conformist") #Assume for now that non-conformist means annual breeding
  
  #Combine conformist and non-conformist mothers
  mothers <- c(mothers, mothers.2)
  
  #-----------------------------------Fathers--------------------------------------------#
  fathers <- which(data1$sex=='M' & data1$age.x>=repro.age)  #determine which males are available to breed in this year
  
  YOY.df <- data.frame() # create an empty data frame to populate with the YOY born in this year
  
  for(j in 1:length(mothers)) { # Loop through all of the available mothers
    
    num.mates.x <- sample(num.mates, size=1) # determine the number of mates for a given mother this year
    inter.df <- data.frame() # Create a data frame for the offspring born from each mother
    
    for(h in 1:num.mates.x) { # Loop through each sire with whom this mother mates 
      
      num.offspring <- rpois(1, ff) #generate the number of offspring born from this mother/sire pairing
      indv.name <- NA #create a place holder for the random name given to each offspring
      age.x <- 0 # assign a 0 age to each offspring born this year
      mother.x <- data1[mothers[j],1] # record the mothers name
      
      father.x <- data1[sample(fathers, size = 1),1] # record the fathers name
      birth.year <- v # assign the birth year in the simulation
      sex <- NA # create a place holder for the biological sex of each offspring
      repro.cycle <- NA
      repro.strategy <- NA

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
      #inter.df2$Survival <- NA
      inter.df <- rbind(inter.df, inter.df2) #add the new borns from this mating pair to the other from the same mother
        } # end loop over mates
    
    YOY.df <- rbind(YOY.df,inter.df) #add all of the offspring from each mother to a data frame of all the YOY for the year
  } # end loop over mothers
  
  loopy.pop <- bind_rows(data1, YOY.df) #Combine the YOY data with the other individuals present this year. YOY go at the bottom. Use "bind_rows" bc there are different numbers of columns.
  
  #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section.
  
  loopy.pop$Survival <- ifelse(loopy.pop$age.x==0, 
                               sample(c("S","M"), 
                                      size=length(which(loopy.pop$age.x==0)), 
                                      prob=c(YOY.survival, 1-YOY.survival), 
                                      replace=T),
                               ifelse(loopy.pop$age.x > 0 &
                                        loopy.pop$age.x < repro.age, 
                                      sample(c("S","M"), 
                                             size=length(which(loopy.pop$age.x > 0 & 
                                                                 loopy.pop$age.x < repro.age)), 
                                             prob=c(juvenile.survival, 1-juvenile.survival),
                                             replace=T),
                                      ifelse(loopy.pop$age.x >= repro.age &
                                               loopy.pop$age.x < max.age, 
                                             sample(c("S","M"), 
                                                    size=length(which(loopy.pop$age.x >= repro.age &
                                                                        loopy.pop$age.x<max.age)),
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
                                   Num.fathers = nrow(dads.temp), #The actual number that reproduced this year
                                   Num.YOY = nrow(YOY.df))
  pop.size <- rbind(pop.size, pop.size.vec)
  
} # end loop over sim years

#Label the list elements with the year
names(loopy.list) <- paste0("year.end.pop.", seq(1:(burn.in + Num.years)), "_iteration_", iter)

return(invisible(list(loopy.list, pop.size, parents.tibble)))
}