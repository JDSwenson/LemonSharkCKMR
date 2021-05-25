#rm(list=ls())
#Sex-specific estimates of N;
#Population growth fixed to observed population growth of the whole population.
#set so R doesn't use scientific notation
options("scipen"=100, "digits"=4)

library(foreach)
library(parallel)
library(doParallel)
library(optimx)
#Load individual packages because tidyverse won't load on cluster
library(plyr)
library(dplyr)
library(tidyr)
library(popbio)
library(mpmtools)

set.seed(47)

#-----------------Leslie Matrix 2--------------------
batchSize <- 1
fb <- batchSize/2 #Batch size for Leslie Matrix - females only
maxAge <- 20
repro.age <- 7
YOY.surv <- .76
adult.surv <- 0.9 #adult survival
juv.surv <- seq(YOY.surv, adult.surv, length.out = repro.age)

#Starts with age 0 individuals
Leslie_input <- data.frame(
  x = c(0:maxAge), #age
  sx = c(juv.surv, rep(adult.surv, times = maxAge-repro.age), 0), #survival
  mx = c(rep(0, times = repro.age), rep(fb, times = maxAge-repro.age+1)) #age-specific birth rates (female proportion of the population only)
)

A1_pre <- make_Leslie_matrix(Leslie_input)
#View(A1)
A1_post <- pre_to_post(Amat = A1_pre, S0 = YOY.surv)
#View(A1_post)

#Calculate dominant eigenvalue (i.e. population growth rate) - should be the same
(lam <- lambda1(A1_pre))
lambda1(A1_post)
#stable_age <- mpmtools::stable_stage(A1_post)

stable_age <- mpmtools::stable_stage(A1_pre) #fishSim seems to use a pre-breeding census for the age distribution

#Calculate stable age structure of the population
# stable_age = stable.stage(A)
# stable_age
survCurv <- stable_age #Sets probability of founder cohort belonging to each age class

#--------------Population simulation-------------------
#Project population forward in time
Abundance_year0 <- c(20000, rep(0, times = (maxAge-1)))

Year1 <- A1_pre %*% Abundance_year0

nYears <- n_yrs        # set the number of years to project
TMat <- A1_pre     # define the projection matrix
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
