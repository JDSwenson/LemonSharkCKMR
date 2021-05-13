####---------Parameters---------####
init.pop.size <- 2000 # Initial population size
init.prop.female <- .5 # proportion of the initial population size that is female
repro.age <- 12 #set age of reproductive maturity
max.age <- 30 #set the maximum age allowed in the simulation
mating.periodicity <- 2 #number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
num.mates <- c(1:3) # vector of potential number of mates per mating
avg.num.offspring <- 3 # set the average number of offspring per mating (from a poisson distribution)
birth.sex.ratio <- c(.5,.5) # The probability that each baby is F:M - has to add up to 1
Adult.survival <- .85 #Adult survival
juvenile.survival <- .85 #juvenile survival
YOY.survival <- .85 # young of year survival
burn.in <- 20 # number of years to use as simulation burn in period
Num.years <- 30 # The number of years to run in the simulation beyond the burn in

#Saves a dataframe for each year of the simulation so I can go back and sample specific years



####---------Set up initial population-----####
init.pop <- data.frame() # create a bank data frame which will become the initial popualation

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

YOY.df <- data.frame() # create an empty data frame to populate with the YOY bor in this year
for(j in 1:length(mothers)){ # Loop through all of the available mothers
  num.mates.x <- sample(num.mates, size=1) # determine the number of mates for a given mother this year
  inter.df <- data.frame() # Create a data frame for the offspring born from each mother
  for(h in 1:num.mates.x){ # Loop through each sire that this mother mates with
  num.offspring <- rpois(1, lambda = avg.num.offspring) # generate the number of offspring born from this mother/sire pairing
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

year.end.pop.0 <- rbind(init.pop, YOY.df) #Combine the YOY data with the other individuals present this year

#Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
year.end.pop.0$Survival <- ifelse(year.end.pop.0$age.x==0, sample(c("S","M"), size=length(which(year.end.pop.0$age.x==0)), prob=c(YOY.survival, 1-YOY.survival), replace=T),
                                  ifelse(year.end.pop.0$age.x<repro.age, sample(c("S","M"), size=length(which(year.end.pop.0$age.x<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                        ifelse(year.end.pop.0$age.x<max.age, sample(c("S","M"), size=length(which(year.end.pop.0$age.x<max.age)), prob=c(Adult.survival, 1-Adult.survival),replace=T), "M")))


loopy.pop <- year.end.pop.0
####-------Loop for all other breeding--------####
pop.size <- data.frame()
for(v in 1:(burn.in + Num.years)){ #loop through all of the years in the simulation - the burn in and the years that matter

data1 <- loopy.pop[loopy.pop$Survival =="S", -8] #Bring in the data from the previous iteration, but only include those that survive
data1$age.x <- data1$age.x+1 # increase each individuals age by one for the new year - happy birthday survivors!

mothers <- which(data1$sex=='F' & data1$age.x>=repro.age & data1$repro.cycle == repro.cycle.vec[v+1])  #determine which females are available to breed in this year
fathers <- which(data1$sex=='M' & data1$age.x>=repro.age)  #determine which males are available to breed in this year

YOY.df <- data.frame() # create an empty data frame to populate with the YOY bor in this year
for(j in 1:length(mothers)){ # Loop through all of the available mothers
  num.mates.x <- sample(num.mates, size=1) # determine the number of mates for a given mother this year
  inter.df <- data.frame() # Create a data frame for the offspring born from each mother
  for(h in 1:num.mates.x){ # Loop through each sire with whom this mother mates 
    num.offspring <- rpois(1, lambda = avg.num.offspring) # generate the number of offspring born from this mother/sire pairing
    indv.name <- NA #create a place holder for the random name given to each offspring
    age.x <- 0 # assign a 0 age to each offspring born this year
    mother.x <- data1[mothers[j],1] # record the mothers name
    father.x <- data1[sample(fathers, size = 1),1] # record the fathers name
    birth.year <- v # assign the birth year in the simulation
    sex <- NA # create a place holder for the biological sex of each offspring
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
    inter.df2$indv.name <- baby.names # add the baby names to the data frame
    inter.df <- rbind(inter.df, inter.df2) #add the new borns from this mating pair to the other from the same mother
  }
  YOY.df <- rbind(YOY.df,inter.df) #add all of the offspring from each mother to a data frame of all the YOY for the year
}


loopy.pop <- rbind(data1, YOY.df) #Combine the YOY data with the other individuals present this year

#Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
loopy.pop$Survival <- ifelse(loopy.pop$age.x==0, sample(c("S","M"), size=length(which(loopy.pop$age.x==0)), prob=c(YOY.survival, 1-YOY.survival), replace=T),
                                  ifelse(loopy.pop$age.x<repro.age, sample(c("S","M"), size=length(which(loopy.pop$age.x<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                         ifelse(loopy.pop$age.x<max.age, sample(c("S","M"), size=length(which(loopy.pop$age.x<max.age)), prob=c(Adult.survival, 1-Adult.survival),replace=T), "M")))
assign(paste("year.end.pop.", v, sep=""),loopy.pop) # save the current year's population data as an object
print(paste("year", v, "N= ", nrow(loopy.pop) , sep=" ")) # print the simulation year and the population size in the R console so they can be oberved
pop.size.vec <- cbind.data.frame(year=v, population_size=nrow(loopy.pop), Male.adult.pop = nrow(loopy.pop[loopy.pop$sex == "M" & loopy.pop$age.x >= repro.age,]), Female.adult.pop = nrow(loopy.pop[loopy.pop$sex == "F" & loopy.pop$age.x >= repro.age,]))
pop.size <- rbind(pop.size, pop.size.vec)
}

####---------quick analysis of population growth (lambda)-----####
lambda <- vector()
for(l in 2:nrow(pop.size)){
  lambda.1 <- pop.size$population_size[l]/pop.size$population_size[l-1]
  lambda <- c(lambda,lambda.1)
}
lambda <- c(NA, lambda)
pop.size$Lambda <- lambda
mean(pop.size$Lambda, na.rm=T)


####---------Checking population parameters-------####
nrow(YOY.df)/length(mothers) #Average fecundity

nrow(loopy.pop[loopy.pop$Survival=='S' & loopy.pop$age.x>=12, ])/nrow(loopy.pop[loopy.pop$age.x>=12,]) 
nrow(loopy.pop[loopy.pop$Survival=='S' & loopy.pop$age.x>0 & loopy.pop$age.x<12, ])/nrow(loopy.pop[loopy.pop$age.x>0 & loopy.pop$age.x<12,])
nrow(loopy.pop[loopy.pop$Survival=='S' & loopy.pop$age.x==0, ])/nrow(loopy.pop[loopy.pop$age.x==0,])

nrow(loopy.pop[loopy.pop$Survival=='M' & loopy.pop$age.x>=12, ])+
  nrow(loopy.pop[loopy.pop$Survival=='M' & loopy.pop$age.x>0 & loopy.pop$age.x<12, ])+
  nrow(loopy.pop[loopy.pop$Survival=='M' & loopy.pop$age.x==0, ])/nrow(loopy.pop[loopy.pop$age.x==0,])

nrow(YOY.df)

#I can go through the output dataframes and say "do they have the same mother AND father? Let Dovi know if I need his help with the pairwise comparison dataframe setup.
#Once I've formatted the data and fit the CKMR model, I can add parameters to the top of the script to make it more customizable.
#The number of rows in a dataframe are the number of individuals in the population