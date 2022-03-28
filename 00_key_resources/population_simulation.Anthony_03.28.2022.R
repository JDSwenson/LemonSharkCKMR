####-------------ANTHONY SEVEQUE, MISSISSIPPI STATE UNIVERSITY, 2022--------------####
####-------------FORWARD IN TIME POPULATION SIMULATION AND KINSHIP ANALYSIS FOR CKMR USE IN LARGE TERRESTRIAL MAMMALS--------------####

####----------A1 PROMISCUITY / RANDOM SAMPLING----------####

a1.promiscuity.random <- function(init.pop.size,
                            init.prop.female,
                            repro.age, 
                            num.mates,
                            birth.sex.ratio, 
                            num.years,
                            yoy.survival,
                            juvenile.survival,
                            adult.survival,
                            max.age, 
                            harvest.size){

###---Set up your initial population---###
  
init.pop <- data.frame() # create a blank data frame which will become the initial population

# Loop that creates the below data for each individual in a population the size of "init.pop.size"
for(i in 1:init.pop.size){
  
  indv.name <- paste(sample(letters, size = 20, replace = T), collapse="") # generate a random name that is 20 letters long
  mother <- "xxxxx" # The individuals in the initial population do not have known mothers
  father <- "xxxxx" # The individuals in the initial population do not have known fathers
  birth.year <- -1 # this is a place holder for individuals born within the simulation
  sex <- sample(c('F','M'), size = 1, prob = c(init.prop.female, 1-init.prop.female)) # randomly draw sex based on the proportions set in the parameter section
  sex.bias <- ifelse(sex == "F", 1, 5) #Sex bias variable for sampling probabilities (weights)
  age <- sample(1:25, size = 1) #Randomly assign an age between 1 and 25 years old
  init.vec <- cbind.data.frame(indv.name, birth.year, age, mother, father, sex, sex.bias) # Create a row with the attributes for each individual - using cbind.data.frame allows them to keep some numeric and others characters
  init.pop <- rbind(init.pop, init.vec) # place the individual into the init.pop data frame
}

# Have a look at your population
head(init.pop)

###---Breeding---###

###---Year 0 breeding---###

mothers <- which(init.pop$sex=='F' & init.pop$age>=repro.age) # Determine which females are available to breed in this year
fathers <- which(init.pop$sex=='M' & init.pop$age>=repro.age) # Determine which fathers are available to breed in this year

YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year

# Loop through all the available mothers
for(j in 1:length(mothers)){
  
  inter.df <- data.frame() # Create a data frame for the offspring born from each mother
  
  num.offspring <- round(rnorm(1, 2.61, 0.82)) # Generate the number of offspring born from this mother/male pairing. Data from black bears den checks in UP.
  num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
  indv.name <- NA # Create a place holder for the random name given to each offspring
  age <- 0 # Assign a 0 age to each offspring born this year
  mother <- init.pop[mothers[j],1] # record the mothers name
  father <- init.pop[sample(fathers, size = 1, replace = FALSE),1] # randomly select and record a male. Replace means the male can be selected more than one time.
  birth.year <- 0 # Note that these pups were born in the 0 year of the simulation
  sex <- NA # Create a place holder for the sexes
  sex.bias <- NA
  inter.df <- cbind.data.frame(indv.name, birth.year, age, mother, father, sex, sex.bias) # Create a row with the attributes of each YOY (at this point only one from each mating pair)
  inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring), ] # Create the random number of offspring from each pairing that we set above
  inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
  inter.df$sex.bias <- ifelse(inter.df$sex == "F", 1, 5) #Sex bias variable for sampling probabilities (weights)
  baby.names <- vector() # Create a blank vector for random baby names
    
  for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
    
    name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
    baby.names <- c(baby.names,name1)  
  } # end loop over name
    
  if(nrow(inter.df)==0){next} # if there were no offspring from this mating pair, skips to the next mating pair (if you don't include this you will get an error in the loop)
  
    inter.df$indv.name <- baby.names # add the baby names to the data frame
 

  YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from each mother to a data frame of all the YOY for the year
} # end loop over mothers

year.end.pop.0 <- NULL
year.end.pop.0 <- rbind(init.pop, YOY.df) #Combine the YOY data with the other individuals present this year

#Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
year.end.pop.0$Survival <- ifelse(year.end.pop.0$age==0, sample(c("S","M"), size=length(which(year.end.pop.0$age==0)), prob=c(yoy.survival, 1-yoy.survival), replace=T),
                                  ifelse(year.end.pop.0$age<repro.age, sample(c("S","M"), size=length(which(year.end.pop.0$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                         ifelse(year.end.pop.0$age<max.age, sample(c("S","M"), size=length(which(year.end.pop.0$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M")))


# Lethal sampling of individuals. Sample individuals that are already naturally dead (M), so that overall survival isn't impacted by harvest size (easier to find stable age structure)

# This is based on Michigan's black bear management strategy, where young, immature, individuals can be hunted.
number.adults = nrow(year.end.pop.0[year.end.pop.0$age >= 1,]) # Size of this year's adult & juvenile population (before survival draw).

sample.pop <- year.end.pop.0 %>%
  subset (Survival == "M") %>% # Only harvest dead individuals
  subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
  slice_sample(n = ((harvest.size*number.adults)/100), replace = FALSE) # Harvest 10% of this year's total adult population

loopy.pop <- year.end.pop.0

# At end of year 0
print(paste("year 0", "N_mothers=", length(mothers), "N_fathers=", length(fathers), "N_newborns=", nrow(YOY.df), "N_deaths=", sum(loopy.pop$Survival=="M"), "N_survivors= ", nrow(loopy.pop[loopy.pop$Survival=="S",]) , sep=" ")) # print the simulation year and the population size in the R console so they can be observed

# Note: not all males or females are fathers or mothers
# Mothers and fathers may have contributed to reproduction this year and die afterwards (since survival is estimated at a late stage)


###---Loop for all other breeding years---###

pop.size <- data.frame()

loopy.list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year

sample.list <- list() # Make list to store dataframe of sampled individual for each year

for(v in 1:(num.years)){ #loop through all of the years in the simulation
  
  data1 <- loopy.pop[loopy.pop$Survival =="S", -8] # Bring in the data from the previous iteration, but only include those that survive
  data1$age <- data1$age+1 # Increase each individuals age by one for the new year
  
  mothers <- which(data1$sex=='F' & data1$age>=repro.age)  # Determine which females are available to breed in this year
  fathers <- which(data1$sex=='M' & data1$age>=repro.age)  # Determine which males are available to breed in this year
  
  juveniles <- which (data1$age >= 1 & data1$age < repro.age)
  
  YOY.df <- data.frame() # create an empty data frame to populate with the YOY born in this year
  
  for(j in 1:length(mothers)) { # Loop through all of the available mothers
    
    inter.df <- data.frame() # Create a data frame for the offspring born from each mother
    
    num.offspring <- round(rnorm(1, 2.61, 0.82)) # Generate the number of offspring born from this mother/male pairing. Data from black bears den checks in UP.
    num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
    indv.name <- NA # Create a place holder for the random name given to each offspring
    age <- 0 # Assign a 0 age to each offspring born this year
    mother <- data1[mothers[j],1] # record the mothers name
    father <- data1[sample(fathers, size = 1, replace = FALSE),1] # record the male name
    birth.year <- v # Note that these pups were born in the 0 year of the simulation
    sex <- NA # Create a place holder for the sexes
    sex.bias <- NA
    inter.df <- cbind.data.frame(indv.name, birth.year, age, mother, father, sex, sex.bias) # Create a row with the attributes of each YOY (at this point only one from each mating pair)
    inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring), ] # Create the random number of offspring from each pairing that we set above
    inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
    inter.df$sex.bias <- ifelse(inter.df$sex == "F", 1, 5) #Sex bias variable for sampling probabilities (weights)
    baby.names <- vector() # Create a blank vector for random baby names
      
      for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
        name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
        baby.names <- c(baby.names,name1)  
      }
      
      if(nrow(inter.df)==0){next} #if there were no offspring from this mating pair, skips to the next mating pair (if you dont include this you will get an error in the loop)
      inter.df$indv.name <- baby.names # add the baby names to the data frame

    YOY.df <- rbind(YOY.df,inter.df) #add all of the offspring from each mother to a data frame of all the YOY for the year
  } # end loop over mothers
  
  loopy.pop <- rbind(data1, YOY.df) #Combine the YOY data with the other individuals present this year. YOY go at the bottom.
  
  #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
  loopy.pop$Survival <- ifelse(loopy.pop$age==0, sample(c("S","M"), size=length(which(loopy.pop$age==0)), prob=c(yoy.survival, 1-yoy.survival), replace=T),
                               ifelse(loopy.pop$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                      ifelse(loopy.pop$age<max.age, sample(c("S","M"), size=length(which(loopy.pop$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 

  # Lethal sampling of individuals.
  
  number.adults = nrow(loopy.pop[loopy.pop$age >= 1,]) # Size of this year's juveniles and adult population
  
  sample.pop <- loopy.pop %>%
    subset (Survival == "M") %>% # Only harvest dead individuals
    subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
    slice_sample(n = ((harvest.size*number.adults)/100), replace = FALSE) # Harvest 10% of this year's total adult population
  
  sample.list[[v]] <- sample.pop
  
  #assign(paste("year.end.pop.", v, sep=""),loopy.pop) # save the current year's population data as an object
   
  loopy.list[[v]] <- loopy.pop # Save the current year's population data as a list element, where the index corresponds to the year
   
  # Print the simulation year and the population size in the R console so they can be observed
  print(paste("year", v, "N_mothers=", length(mothers), "N_fathers=", length(fathers), "N_juveniles=", length(juveniles), "N_newborns=", nrow(YOY.df), "N_deaths=", sum(loopy.pop$Survival=="M", loopy.pop$Survival=="H"), "N_survivors= ", nrow(loopy.pop[loopy.pop$Survival=="S",]) , sep=" "))
   
  ### Count survivors for each year v
  pop.size.vec <- cbind.data.frame(year=v, population_size=nrow(loopy.pop[loopy.pop$Survival=="S",]), #
                                   Nm_preharvest = length (fathers), # Includes males that may have reproduced this year but were subsequently killed (harvest or natural)
                                   Nf_preharvest = length(mothers), # Includes females that may have reproduced this year but were subsequently killed (harvest or natural)
                                    Nm_postharvest = nrow(loopy.pop[loopy.pop$sex == "M" & loopy.pop$age >= repro.age & loopy.pop$Survival=="S",]), # Adult males alive at the end of the year
                                    Nf_postharvest = nrow(loopy.pop[loopy.pop$sex == "F" & loopy.pop$age >= repro.age & loopy.pop$Survival=="S",])) # Adult females alive at the end of the year
  pop.size <- rbind(pop.size, pop.size.vec)
   
} # end loop over sim years

# Label the list elements with the year
names(loopy.list) <- paste0("year.end.pop.", seq(1:num.years))

return(invisible(list(loopy.list, sample.list, pop.size)))
}


####----------A2 PROMISCUITY / SEX BIAS SAMPLING----------####

a2.promiscuity.sexbias <- function(init.pop.size,
                                  init.prop.female,
                                  repro.age, 
                                  num.mates,
                                  birth.sex.ratio, 
                                  num.years,
                                  yoy.survival,
                                  juvenile.survival,
                                  adult.survival,
                                  max.age, 
                                  harvest.size){
  
  ###---Set up your initial population---###
  
  init.pop <- data.frame() # create a blank data frame which will become the initial population
  
  # Loop that creates the below data for each individual in a population the size of "init.pop.size"
  for(i in 1:init.pop.size){
    
    indv.name <- paste(sample(letters, size = 20, replace = T), collapse="") # generate a random name that is 20 letters long
    mother <- "xxxxx" # The individuals in the initial population do not have known mothers
    father <- "xxxxx" # The individuals in the initial population do not have known fathers
    birth.year <- -1 # this is a place holder for individuals born within the simulation
    sex <- sample(c('F','M'), size = 1, prob = c(init.prop.female, 1-init.prop.female)) # randomly draw sex based on the proportions set in the parameter section
    sex.bias <- ifelse(sex == "F", 1, 5) #Sex bias variable for sampling probabilities (weights)
    age <- sample(1:25, size = 1) #Randomly assign an age between 1 and 25 years old
    init.vec <- cbind.data.frame(indv.name, birth.year, age, mother, father, sex, sex.bias) # Create a row with the attributes for each individual - using cbind.data.frame allows them to keep some numeric and others characters
    init.pop <- rbind(init.pop, init.vec) # place the individual into the init.pop data frame
  }
  
  # Have a look at your population
  head(init.pop)
  
  ###---Breeding---###
  
  ###---Year 0 breeding---###
  
  mothers <- which(init.pop$sex=='F' & init.pop$age>=repro.age) # Determine which females are available to breed in this year
  fathers <- which(init.pop$sex=='M' & init.pop$age>=repro.age) # Determine which fathers are available to breed in this year
  
  YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
  
  # Loop through all the available mothers
  for(j in 1:length(mothers)){
    
    inter.df <- data.frame() # Create a data frame for the offspring born from each mother
   
    num.offspring <- round(rnorm(1, 2.61, 0.82)) # Generate the number of offspring born from this mother/male pairing.
    num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
    indv.name <- NA # Create a place holder for the random name given to each offspring
    age <- 0 # Assign a 0 age to each offspring born this year
    mother <- init.pop[mothers[j],1] # record the mothers name
    father <- init.pop[sample(fathers, size = 1, replace = TRUE),1] # randomly select and record a male
    birth.year <- 0 # Note that these pups were born in the 0 year of the simulation
    sex <- NA # Create a place holder for the sexes
    sex.bias <- NA
    inter.df <- cbind.data.frame(indv.name, birth.year, age, mother, father, sex, sex.bias) # Create a row with the attributes of each YOY (at this point only one from each mating pair)
    inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring), ] # Create the random number of offspring from each pairing that we set above
    inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
    inter.df$sex.bias <- ifelse(inter.df$sex == "F", 1, 5) #Sex bias variable for sampling probabilities (weights)
    baby.names <- vector() # Create a blank vector for random baby names
    
    for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
      
      name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
      baby.names <- c(baby.names,name1)  
    } # end loop over name
    
    if(nrow(inter.df)==0){next} # if there were no offspring from this mating pair, skips to the next mating pair (if you don't include this you will get an error in the loop)
    
    inter.df$indv.name <- baby.names # add the baby names to the data frame
    
    YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from each mother to a data frame of all the YOY for the year
  } # end loop over mothers
  
  year.end.pop.0 <- NULL
  year.end.pop.0 <- rbind(init.pop, YOY.df) #Combine the YOY data with the other individuals present this year
  
  #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
  year.end.pop.0$Survival <- ifelse(year.end.pop.0$age==0, sample(c("S","M"), size=length(which(year.end.pop.0$age==0)), prob=c(yoy.survival, 1-yoy.survival), replace=T),
                                    ifelse(year.end.pop.0$age<repro.age, sample(c("S","M"), size=length(which(year.end.pop.0$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                           ifelse(year.end.pop.0$age<max.age, sample(c("S","M"), size=length(which(year.end.pop.0$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M")))
  
  
  # Lethal sampling of individuals. Sample individuals that are already naturally dead (M), so that overall survival isn't impacted by harvest size (easier to find stable age structure)
  
  # This is based on Michigan's black bear management strategy, where young, immature, individuals can be hunted.
  number.adults = nrow(year.end.pop.0[year.end.pop.0$age >= 1,]) # Size of this year's adult & juvenile population (before survival draw).
  
  sample.pop <- year.end.pop.0 %>%
    subset (Survival == "M") %>% # Only harvest dead individuals
    subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
    slice_sample(n = ((harvest.size*number.adults)/100), weight_by = sex.bias, replace = FALSE) # Harvest 10% of this year's total adult population
  
  loopy.pop <- year.end.pop.0
  
  # At end of year 0
  print(paste("year 0", "N_mothers=", length(mothers), "N_fathers=", length(fathers), "N_pups=", nrow(YOY.df), "N_deaths=", sum(loopy.pop$Survival=="M", loopy.pop$Survival=="H"), "N_survivors= ", nrow(loopy.pop[loopy.pop$Survival=="S",]) , sep=" ")) # print the simulation year and the population size in the R console so they can be observed

  # Note: not all males or females are fathers or mothers
  # Mothers and fathers may have contributed to reproduction this year and die afterwards (since survival is estimated at a late stage)
  
  
  ###---Loop for all other breeding years---###
  
  pop.size <- data.frame()
  
  loopy.list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year
  
  sample.list <- list() # Make list to store dataframe of sampled individual for each year
  
  for(v in 1:(num.years)){ #loop through all of the years in the simulation
    
    data1 <- loopy.pop[loopy.pop$Survival =="S", -8] # Bring in the data from the previous iteration, but only include those that survive
    data1$age <- data1$age+1 # Increase each individuals age by one for the new year
    
    mothers <- which(data1$sex=='F' & data1$age>=repro.age)  # Determine which females are available to breed in this year
    fathers <- which(data1$sex=='M' & data1$age>=repro.age)  # Determine which males are available to breed in this year
    
    YOY.df <- data.frame() # create an empty data frame to populate with the YOY born in this year
    
    for(j in 1:length(mothers)) { # Loop through all of the available mothers
      
      inter.df <- data.frame() # Create a data frame for the offspring born from each mother
      
      num.offspring <- round(rnorm(1, 2.61, 0.82)) # Generate the number of offspring born from this mother/male pairing.
      num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
      indv.name <- NA # Create a place holder for the random name given to each offspring
      age <- 0 # Assign a 0 age to each offspring born this year
      mother <- data1[mothers[j],1] # record the mothers name
      father <- data1[sample(fathers, size = 1, replace = TRUE),1] # record the male name
      birth.year <- v # Note that these pups were born in the 0 year of the simulation
      sex <- NA # Create a place holder for the sexes
      sex.bias <- NA
      inter.df <- cbind.data.frame(indv.name, birth.year, age, mother, father, sex, sex.bias) # Create a row with the attributes of each YOY (at this point only one from each mating pair)
      inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring), ] # Create the random number of offspring from each pairing that we set above
      inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
      inter.df$sex.bias <- ifelse(inter.df$sex == "F", 1, 5) #Sex bias variable for sampling probabilities (weights)
      baby.names <- vector() # Create a blank vector for random baby names
      
      for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
        name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
        baby.names <- c(baby.names,name1)  
      }
      
      if(nrow(inter.df)==0){next} #if there were no offspring from this mating pair, skips to the next mating pair (if you dont include this you will get an error in the loop)
      inter.df$indv.name <- baby.names # add the baby names to the data frame
      
      YOY.df <- rbind(YOY.df,inter.df) #add all of the offspring from each mother to a data frame of all the YOY for the year
    } # end loop over mothers
    
    loopy.pop <- rbind(data1, YOY.df) #Combine the YOY data with the other individuals present this year. YOY go at the bottom.
    
    #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
    loopy.pop$Survival <- ifelse(loopy.pop$age==0, sample(c("S","M"), size=length(which(loopy.pop$age==0)), prob=c(yoy.survival, 1-yoy.survival), replace=T),
                                 ifelse(loopy.pop$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                        ifelse(loopy.pop$age<max.age, sample(c("S","M"), size=length(which(loopy.pop$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
    
    # Lethal sampling of individuals. Sample individuals that are already naturally dead (M), so that overall survival isn't impacted by harvest size (easier to find stable age structure)
    
    # This is based on Michigan's black bear management strategy, where young, immature, individuals can be hunted.
    number.adults = nrow(year.end.pop.0[year.end.pop.0$age >= 1,]) # Size of this year's adult & juvenile population (before survival draw).
    
    sample.pop <- loopy.pop %>%
      subset (Survival == "M") %>% # Only harvest dead individuals
      subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
      slice_sample(n = ((harvest.size*number.adults)/100), weight_by = sex.bias, replace = FALSE) # Harvest 10% of this year's total adult population
    
    sample.list[[v]] <- sample.pop
    
    #assign(paste("year.end.pop.", v, sep=""),loopy.pop) # save the current year's population data as an object
    
    loopy.list[[v]] <- loopy.pop # Save the current year's population data as a list element, where the index corresponds to the year
    
    # Print the simulation year and the population size in the R console so they can be observed
    print(paste("year", v, "N_mothers=", length(mothers), "N_fathers=", length(fathers), "N_pups=", nrow(YOY.df), "N_deaths=", sum(loopy.pop$Survival=="M", loopy.pop$Survival=="H"), "N_survivors= ", nrow(loopy.pop[loopy.pop$Survival=="S",]) , sep=" "))
    
    ### Count survivors for each year v
    pop.size.vec <- cbind.data.frame(year=v, population_size=nrow(loopy.pop[loopy.pop$Survival=="S",]), #
                                     Nm_preharvest = length (fathers), # Includes males that may have reproduced this year but were subsequently killed (harvest or natural)
                                     Nf_preharvest = length(mothers), # Includes females that may have reproduced this year but were subsequently killed (harvest or natural)
                                     Nm_postharvest = nrow(loopy.pop[loopy.pop$sex == "M" & loopy.pop$age >= repro.age & loopy.pop$Survival=="S",]), # Adult males alive at the end of the year
                                     Nf_postharvest = nrow(loopy.pop[loopy.pop$sex == "F" & loopy.pop$age >= repro.age & loopy.pop$Survival=="S",])) # Adult females alive at the end of the year
    pop.size <- rbind(pop.size, pop.size.vec)
    
  } # end loop over sim years
  
  # Label the list elements with the year
  names(loopy.list) <- paste0("year.end.pop.", seq(1:num.years))
  
  return(invisible(list(loopy.list, sample.list, pop.size)))
}

####----------B1. POLYGYNY / RANDOM SAMPLING----------####

b1.polygyny.random <- function(init.pop.size,
                                  init.prop.female,
                                  repro.age, 
                                  num.mates,
                                  birth.sex.ratio, 
                                  num.years,
                                  yoy.survival,
                                  juvenile.survival,
                                  adult.survival,
                                  max.age, 
                                  harvest.size){
  
  ###---Set up your initial population---###
  
  init.pop <- data.frame() # create a blank data frame which will become the initial population
  
  # Loop that creates the below data for each individual in a population the size of "init.pop.size"
  for(i in 1:init.pop.size){
    
    indv.name <- paste(sample(letters, size = 20, replace = T), collapse="") # generate a random name that is 20 letters long
    mother <- "xxxxx" # The individuals in the initial population do not have known mothers
    father <- "xxxxx" # The individuals in the initial population do not have known fathers
    birth.year <- -1 # this is a place holder for individuals born within the simulation
    sex <- sample(c('F','M'), size = 1, prob = c(init.prop.female, 1-init.prop.female)) # randomly draw sex based on the proportions set in the parameter section
    sex.bias <- ifelse(sex == "F", 1, 5) # Sex bias variable for sampling probabilities (weights)
    male.repro <- rbeta(1,0.5,1) # Create "probability of reproduction" for males. Change values of rbeta to create different distributions.
    age <- sample(1:25, size = 1) #Randomly assign an age between 1 and 25 years old
    init.vec <- cbind.data.frame(indv.name, birth.year, age, mother, father, sex, sex.bias, male.repro) # Create a row with the attributes for each individual - using cbind.data.frame allows them to keep some numeric and others characters
    init.pop <- rbind(init.pop, init.vec) # place the individual into the init.pop data frame
  }
  
  # Have a look at your population
  head(init.pop)
  
  ###---Breeding---###
  
  ###---Year 0 breeding---###
  
  mothers <- which(init.pop$sex=='F' & init.pop$age>=repro.age) # Determine which females are available to breed in this year
  fathers <- which(init.pop$sex=='M' & init.pop$age>=repro.age) # Determine which fathers are available to breed in this year
  
  YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
  
  # Loop through all the available mothers
  for(j in 1:length(mothers)){
    
    inter.df <- data.frame() # Create a data frame for the offspring born from each mother
    num.offspring <- round(rnorm(1, 2.5, 0.5)) # Generate the number of offspring born from this mother/male pairing.
    indv.name <- NA # Create a place holder for the random name given to each offspring
    age <- 0 # Assign a 0 age to each offspring born this year
    mother <- init.pop[mothers[j],1] # record the mothers name
    
    #Loop for mate selection
    for (k in 1:length(fathers)){
      dice.roll <- runif(1,0,1) #Create random value between 0 and 1
      potential.father <- init.pop[sample(fathers, size = 1, replace = TRUE), 1:8] # randomly select a male
      ifelse(dice.roll <= potential.father$male.repro, # Simulate a dice roll that will select, or not, the male
             father <- potential.father[,1],
             father <- ("next try")) # Cannot break a loop within an ifelse, so you have to do that in two separate steps.
      if (father != "next try"){ # If mate has been found, then break the loop
        break()
      }}
    
    birth.year <- 0 # Note that these pups were born in the 0 year of the simulation
    sex <- NA # Create a place holder for the sexes
    sex.bias <- NA
    male.repro <- NA
    inter.df <- cbind.data.frame(indv.name, birth.year, age, mother, father, sex, sex.bias, male.repro) # Create a row with the attributes of each YOY (at this point only one from each mating pair)
    inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring), ] # Create the random number of offspring from each pairing that we set above
    inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
    inter.df$sex.bias <- ifelse(inter.df$sex == "F", 1, 5) #Sex bias variable for sampling probabilities (weights)
    inter.df$male.repro <- rbeta(1,0.5,1) # Create "probability of reproduction" for males. Change values of rbeta to create different distributions.
    baby.names <- vector() # Create a blank vector for random baby names
    
    for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
      
      name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
      baby.names <- c(baby.names,name1)  
    } # end loop over name
    
    if(nrow(inter.df)==0){next} # if there were no offspring from this mating pair, skips to the next mating pair (if you don't include this you will get an error in the loop)
    
    inter.df$indv.name <- baby.names # add the baby names to the data frame
    
    
    YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from each mother to a data frame of all the YOY for the year
  } # end loop over mothers
  
  year.end.pop.0 <- NULL
  year.end.pop.0 <- rbind(init.pop, YOY.df) #Combine the YOY data with the other individuals present this year
  
  #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
  year.end.pop.0$Survival <- ifelse(year.end.pop.0$age==0, sample(c("S","M"), size=length(which(year.end.pop.0$age==0)), prob=c(yoy.survival, 1-yoy.survival), replace=T),
                                    ifelse(year.end.pop.0$age<repro.age, sample(c("S","M"), size=length(which(year.end.pop.0$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                           ifelse(year.end.pop.0$age<max.age, sample(c("S","M"), size=length(which(year.end.pop.0$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M")))
  
  # Lethal sampling of individuals. No need to save them for the first year of the simulation, but it prevents the first year of simulation from "jumping high" if no sampling at all.
  sample.pop <- slice_sample(year.end.pop.0, prop = harvest.size, weight_by = sex.bias, replace = FALSE) 
  
  # Assign all sampled individuals as "harvested", so that they are considered a dead for the next year
  year.end.pop.0 <- within (year.end.pop.0, {
    Survival[year.end.pop.0$indv.name %in% sample.pop$indv.name] = "H"
  })
  
  loopy.pop <- year.end.pop.0
  
  # At end of year 0
  print(paste("year 0", "N_mothers=", length(mothers), "N_fathers=", length(fathers), "N_pups=", nrow(YOY.df), "N_deaths=", sum(loopy.pop$Survival=="M", loopy.pop$Survival=="H"), "N_survivors= ", nrow(loopy.pop[loopy.pop$Survival=="S",]) , sep=" ")) # print the simulation year and the population size in the R console so they can be observed
  
  # Note: not all males or females are fathers or mothers
  # Mothers and fathers may have contributed to reproduction this year and die afterwards (since survival is estimated at a late stage)

    ###---Loop for all other breeding years---###
  
  pop.size <- data.frame()
  
  loopy.list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year
  
  sample.list <- list() # Make list to store dataframe of sampled individual for each year
  
  for(v in 1:(num.years)){ #loop through all of the years in the simulation
    
    data1 <- loopy.pop[loopy.pop$Survival =="S", -9] # Bring in the data from the previous iteration, but only include those that survive
    data1$age <- data1$age+1 # Increase each individuals age by one for the new year
    
    mothers <- which(data1$sex=='F' & data1$age>=repro.age)  # Determine which females are available to breed in this year
    fathers <- which(data1$sex=='M' & data1$age>=repro.age)  # Determine which males are available to breed in this year
    
    YOY.df <- data.frame() # create an empty data frame to populate with the YOY born in this year
    
    for(j in 1:length(mothers)) { # Loop through all of the available mothers
      
      inter.df <- data.frame() # Create a data frame for the offspring born from each mother
      
      num.offspring <- round(rnorm(1, 2.5, 0.5)) # Generate the number of offspring born from this mother/male pairing.
      indv.name <- NA # Create a place holder for the random name given to each offspring
      age <- 0 # Assign a 0 age to each offspring born this year
      mother <- data1[mothers[j],1] # record the mothers name
      
      #Loop for mate selection
      for (k in 1:length(fathers)){
        dice.roll <- runif(1,0,1) #Create random value between 0 and 1
        potential.father <- init.pop[sample(fathers, size = 1, replace = TRUE), 1:8] # randomly select a male
        ifelse(dice.roll <= potential.father$male.repro, # Simulate a dice roll that will select, or not, the male
               father <- potential.father[,1],
               father <- ("next try")) # Cannot break a loop within an ifelse, so you have to do that in two separate steps.
        if (father != "next try"){ # If mate has been found, then break the loop
          break()
        }}
      
      
      # This error may pop up randomly (at any given year? changes constantly) and kill the loop, need to troubleshoot why
      
      # Error in `$<-.data.frame`(`*tmp*`, "male.repro", value = 0.259101062327585) : 
      # replacement has 1 row, data has 0
      
      birth.year <- v # Note that these pups were born in the 0 year of the simulation
      sex <- NA # Create a place holder for the sexes
      sex.bias <- NA
      male.repro <- NA
      inter.df <- cbind.data.frame(indv.name, birth.year, age, mother, father, sex, sex.bias, male.repro) # Create a row with the attributes of each YOY (at this point only one from each mating pair)
      inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring), ] # Create the random number of offspring from each pairing that we set above
      inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
      inter.df$sex.bias <- ifelse(inter.df$sex == "F", 1, 5) #Sex bias variable for sampling probabilities (weights)
      inter.df$male.repro <- rbeta(1,0.5,1) # Create "probability of reproduction" for males. Change values of rbeta to create different distributions.
      baby.names <- vector() # Create a blank vector for random baby names
      
      for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
        name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
        baby.names <- c(baby.names,name1)  
      }
      
      if(nrow(inter.df)==0){next} #if there were no offspring from this mating pair, skips to the next mating pair (if you dont include this you will get an error in the loop)
      inter.df$indv.name <- baby.names # add the baby names to the data frame
      
      YOY.df <- rbind(YOY.df,inter.df) #add all of the offspring from each mother to a data frame of all the YOY for the year
    } # end loop over mothers
    
    loopy.pop <- rbind(data1, YOY.df) #Combine the YOY data with the other individuals present this year. YOY go at the bottom.
    
    # Lethal sampling of individuals.
    
    sample.pop <- slice_sample(loopy.pop, prop = harvest.size, replace = FALSE) 
    
    sample.list[[v]] <- sample.pop
    
    #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
    loopy.pop$Survival <- ifelse(loopy.pop$age==0, sample(c("S","M"), size=length(which(loopy.pop$age==0)), prob=c(yoy.survival, 1-yoy.survival), replace=T),
                                 ifelse(loopy.pop$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                        ifelse(loopy.pop$age<max.age, sample(c("S","M"), size=length(which(loopy.pop$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
    
    # Assign all sampled individuals as "harvested", so that they are considered a dead for the next year
    loopy.pop <- within (loopy.pop, {
      Survival[loopy.pop$indv.name %in% sample.pop$indv.name] = "H"
    })
    
    # To check with Dana: right now, we are also harvesting individuals that would be naturally dead the same year. 
    # So, natural mortality and harvest do not fully "add up", as they borrow from each other. 
    # If we want to separate both processes, above line of code needs changing.
    
    #assign(paste("year.end.pop.", v, sep=""),loopy.pop) # save the current year's population data as an object
    
    loopy.list[[v]] <- loopy.pop # Save the current year's population data as a list element, where the index corresponds to the year
    
    # Print the simulation year and the population size in the R console so they can be observed
    print(paste("year", v, "N_mothers=", length(mothers), "N_fathers=", length(fathers), "N_pups=", nrow(YOY.df), "N_deaths=", sum(loopy.pop$Survival=="M", loopy.pop$Survival=="H"), "N_survivors= ", nrow(loopy.pop[loopy.pop$Survival=="S",]) , sep=" "))
    
    ### Count survivors for each year v
    pop.size.vec <- cbind.data.frame(year=v, population_size=nrow(loopy.pop[loopy.pop$Survival=="S",]), # 
                                     Male.adult.pop = nrow(loopy.pop[loopy.pop$sex == "M" & loopy.pop$age >= repro.age & loopy.pop$Survival=="S",]), # 
                                     Female.adult.pop = nrow(loopy.pop[loopy.pop$sex == "F" & loopy.pop$age >= repro.age & loopy.pop$Survival=="S",])) # 
    pop.size <- rbind(pop.size, pop.size.vec)
    
  } # end loop over sim years
  
  # Label the list elements with the year
  names(loopy.list) <- paste0("year.end.pop.", seq(1:num.years))
  
  return(invisible(list(loopy.list, sample.list, pop.size)))
}

####----------B2. POLYGYNY / SEX BIAS SAMPLING----------####
####----------C1. SERIAL MONOGAMY / RANDOM SAMPLING----------####

# When a male is selected, he can no longer reproduce this year (so remove him from fathers)
# Change how many males a female can mate with?
# Promiscuity = male can reproduce with more than one female - too similar to serial monogamy? (ie whats the proba that a male will be selected more than once in promiscuity)

####----------C2. SERIAL MONOGAMY / SEX BIAS SAMPLING----------####
