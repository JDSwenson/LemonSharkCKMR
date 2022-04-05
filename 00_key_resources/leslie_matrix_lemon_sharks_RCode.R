library(tidyverse)
library(mpmtools)
library(popbio)

batchSize <- 3.5
fb <- batchSize/2 #Batch size for Leslie Matrix - females only
maxAge <- 50
YOY.sx <- 0.6
Juv.sx <- 0.8
Adult.sx <- 0.9
repro.age <- 12

Leslie_input <- data.frame(
  x = c(0:maxAge), #age
  sx = c(YOY.sx, rep(Juv.sx, times = (repro.age-1)), rep(Adult.sx, times = maxAge - repro.age), 0), #survival
  mx = c(rep(0, times = repro.age), rep(fb, times = (maxAge+1) - repro.age)) #age-specific birth rates (female proportion of the population only)
)

A1_pre <- make_Leslie_matrix(Leslie_input)
#View(A1)
A1_post <- pre_to_post(Amat = A1_pre, S0 = .6)
#View(A1_post)

#Calculate dominant eigenvalue (i.e. population growth rate) from transition matrix
lambda1(A1_pre)
lambda1(A1_post)
#stable_age <- mpmtools::stable_stage(A1_post)

stable_age_pre <- mpmtools::stable_stage(A1_pre) #fishSim seems to use a pre-breeding census for the age distribution

stable_age_post <- mpmtools::stable_stage(A1_post)

#Calculate stable age structure of the population
# stable_age = stable.stage(A)
# stable_age
survCurv <- stable_age #Sets probability of founder cohort belonging to each age class

#--------------Population simulation-------------------
#Project population forward in time
Abundance_year0 <- c(20000, rep(0, times = (max_age-1)))

Year1 <- A %*% Abundance_year0

nYears <- n_yrs        # set the number of years to project
TMat <- A     # define the projection matrix
InitAbund <- Abundance_year0    # define the initial abundance

## NOTE: the code below can be re-used without modification:
allYears <- matrix(0,nrow=nrow(TMat),ncol=nYears+1)     # build a storage array for all abundances!
allYears[,1] <- InitAbund  # set the year 0 abundance                                    
for(t in 2:(nYears+1)){   # loop through all years
  allYears[,t] <-  TMat %*% allYears[,t-1]
}

#in allYears, columns are years, rows are ages
allYears <- data.frame(allYears)
colnames(allYears) <- paste("Yr_", c(1:(nYears+1)))
allYears <- allYears %>% mutate_all(round, digits=0)

#Sum mature ages for each year
adult_truth <- allYears %>% slice(12:30) %>% 
  summarise_all(sum)




#------------Table 2 from Waples and Feutry----------------
batchSize <- 1
fb <- batchSize/2 #Batch size for Leslie Matrix - females only
maxAge <- 10
YOY.sx <- 0.7
Juv.sx <- 0.7
Adult.sx <- 0.7
repro.age <- 3

Leslie_input <- data.frame(
  x = c(0:maxAge), #age
  sx = c(YOY.sx, rep(Juv.sx, times = (repro.age-1)), rep(Adult.sx, times = maxAge - repro.age), 0), #survival
  mx = c(rep(0, times = repro.age), rep(fb, times = (maxAge+1) - repro.age)) #age-specific birth rates (female proportion of the population only)
)

A1_pre <- make_Leslie_matrix(Leslie_input)
#View(A1)
A1_post <- pre_to_post(Amat = A1_pre, S0 = .7)
#View(A1_post)

#Calculate dominant eigenvalue (i.e. population growth rate) from transition matrix
popbio::lambda(A1_pre)
popbio::stable.stage(A1_pre)
rv <- popbio::reproductive.value(A1_pre)

rv.stand <- c(rv[repro.age:maxAge]/sum(rv[repro.age:maxAge]))

