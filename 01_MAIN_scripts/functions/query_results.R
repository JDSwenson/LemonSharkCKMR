####------------- Examine simulation results ---------------####
plot(pop.size$population_size, pch=19, ylim=c(0.9*min(pop.size$population_size), 1.1*max(pop.size$population_size))) #Plot population size through time

#nrow(YOY.df)/length(mothers) # Average fecundity for last year; remember whether they've skipped breeding

#Make dataframe of abundance for each year
pop.size <- pop.size %>% mutate(Total.adult.pop = Male.adult.pop + Female.adult.pop)

#Calculate population growth for whole population
total.lambda <- NULL
for(l in 2:nrow(pop.size)){ 
  total.lambda.1 <- pop.size$population_size[l]/pop.size$population_size[l-1]
  total.lambda <- c(total.lambda, total.lambda.1)
}

#Calculate population growth for adults only
adult.lambda <- NULL
for(l in 2:nrow(pop.size)){ 
  adult.lambda.1 <- pop.size$Total.adult.pop[l]/pop.size$Total.adult.pop[l-1]
  adult.lambda <- c(adult.lambda, adult.lambda.1)
}

#Add NA to first element of population growth vectors
total.lambda <- c(NA, total.lambda)
adult.lambda <- c(NA, adult.lambda)

#Add population growth per year to pop.size dataframe
pop.size$total.lambda <- total.lambda
pop.size$adult.lambda <- adult.lambda

#plot(total.lambda[(burn.in+1):n_yrs], pch=19)
#abline(h=1, lty=3)

(mean.total.lam <- mean(pop.size$total.lambda[(burn.in+1):n_yrs], na.rm=T)) # mean population growth for whole population
sd(pop.size$total.lambda[(burn.in+1):n_yrs], na.rm=T) # sd Lambda


#Calculate SURVIVAL for each year
sVec <- NULL #Make empty vector to save yearly survival rates

#Store annual survival of adults
for(yr in 2:length(loopy.list)){
  sVec[yr] <- length(which(loopy.list[[yr]]$Survival=='S' & loopy.list[[yr]]$age.x>=repro.age))/length(which(loopy.list[[yr]]$age.x>=repro.age))
}
sVec <- c(NA, sVec) #Add NA to survival vector for first year of simulation

length(which(loopy.list[[n_yrs]]$Survival=='S' & loopy.list[[n_yrs]]$age.x>0 & loopy.list[[n_yrs]]$age.x<repro.age))/length(which(loopy.list[[n_yrs]]$age.x>0 & loopy.list[[n_yrs]]$age.x<repro.age)) #survival of juveniles

length(which(loopy.list[[n_yrs]]$Survival=='S' & loopy.list[[n_yrs]]$age.x==0))/length(which(loopy.list[[n_yrs]]$age.x==0)) #survival of YOY

