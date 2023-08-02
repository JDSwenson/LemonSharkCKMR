####------------- Examine simulation results ---------------####
plot(pop.size.tibble$population_size, pch=19, ylim=c(0.9*min(pop.size.tibble$population_size), 1.1*max(pop.size.tibble$population_size))) #Plot population size through time

#nrow(YOY.df)/length(mothers) # Average fecundity for last year; remember whether they've skipped breeding

#Make dataframe of abundance for each year
pop.size.tibble <- pop.size.tibble %>% mutate(Total.adult.pop = Male.adult.pop + Female.adult.pop)

#Calculate population growth for whole population
total.lambda <- NULL
for(l in 2:nrow(pop.size.tibble)){ 
  total.lambda.1 <- pop.size.tibble$population_size[l]/pop.size.tibble$population_size[l-1]
  total.lambda <- c(total.lambda, total.lambda.1)
}

#Calculate population growth for adults only
adult.lambda <- NULL
for(l in 2:nrow(pop.size.tibble)){ 
  adult.lambda.1 <- pop.size.tibble$Total.adult.pop[l]/pop.size.tibble$Total.adult.pop[l-1]
  adult.lambda <- c(adult.lambda, adult.lambda.1)
}

#Add NA to first element of population growth vectors
total.lambda <- c(NA, total.lambda)
adult.lambda <- c(NA, adult.lambda)

#Add population growth per year to pop.size dataframe
pop.size.tibble$total.lambda <- total.lambda
pop.size.tibble$adult.lambda <- adult.lambda

#plot(total.lambda[(burn.in+1):n_yrs], pch=19)
#abline(h=1, lty=3)


#Calculate SURVIVAL for each year
sVec <- NULL #Make empty vector to save yearly survival rates

#Store annual survival of adults
for(yr in 2:length(loopy.list)){
  sVec[yr] <- length(which(loopy.list[[yr]]$Survival=='S' & loopy.list[[yr]]$age.x>=repro.age))/length(which(loopy.list[[yr]]$age.x>=repro.age))
}
# <- c(NA, sVec) #Add NA to survival vector for first year of simulation

length(which(loopy.list[[n_yrs]]$Survival=='S' & loopy.list[[n_yrs]]$age.x>0 & loopy.list[[n_yrs]]$age.x<repro.age))/length(which(loopy.list[[n_yrs]]$age.x>0 & loopy.list[[n_yrs]]$age.x<repro.age)) #survival of juveniles

length(which(loopy.list[[n_yrs]]$Survival=='S' & loopy.list[[n_yrs]]$age.x==0))/length(which(loopy.list[[n_yrs]]$age.x==0)) #survival of YOY


#Calculate ERRO and expected number of HSPs
