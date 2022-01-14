rm(list=ls())
#Sex-specific estimates of N;
#Population growth fixed to observed population growth of the whole population.
#set so R doesn't use scientific notation
options("scipen"=100, "digits"=4)

library(tidyverse)
library(ggpubr)
library(MASS)
library(popbio)
#devtools::install_github("BruceKendall/mpmtools")
library(mpmtools)

#-----------------Leslie Matrix--------------------
cv <- c(0.05, 0.1)

maxAge <- 50
repro.age <- 12
juv.ages <- repro.age - 1

ages <- c(0:maxAge)
adult.ages <- length(ages) - (juv.ages + 1)
YOY.survival <- 0.7 # CHANGED FROM 0.8; young of year survival
juvenile.survival <- 0.8 # CHANGED FROM 0.9; juvenile survival
mean.survival <- 0.825
surv.sd <- cv * mean.survival


mating.periodicity <- 1 # CHANGED FROM 2; number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
num.mates <- c(1:3) # CHANGED FROM c(1:3); vector of potential number of mates per mating
f <- (1-mean.survival)/(YOY.survival * juvenile.survival^11) # adult fecundity at equilibrium if no age truncation
init.prop.female = 0.5
ff <- f/init.prop.female * mating.periodicity/mean(num.mates) # female fecundity per breeding cycle
ff

mean.fecundity <- ff
fec.sd <- cv * mean.fecundity #Fecundity per breeding cycle i.e. ff

corr.vec <- c(0.25, 0, -0.25, -0.5, -0.75)

vec.mean <- c(mean.fecundity, mean.survival)   #vector of the means

n.draws <- 100

lambda.df <- NULL

#Loop through different correlation values
for(c in 1:length(corr.vec)){
  correlation <- corr.vec[c]

#Loop through different cv values for survival
  for(s in 1:length(surv.sd)){
    #mean.survival <- surv.vec[s]
    sd.survival <- surv.sd[s]
    surv.cv <- sd.survival/mean.survival

  #Loop through different cv values for fecundity
    for(fec in 1:length(fec.sd)){
      #mean.fecundity <- fec.vec[fec]
      sd.fecundity <- fec.sd[fec]
      fec.cv <- sd.fecundity/mean.fecundity
      
      #create sd matrix
      vec.sd <- c(sd.fecundity, sd.survival)
      mat.sd <- diag(vec.sd)  
      
      
      #create correlation matrix (assuming correlation = -0.5) ====
      corr.fec.surv <- correlation
      corr.mat <- matrix(c(1, corr.fec.surv, corr.fec.surv, 1), nrow=2, ncol=2, byrow=T)
      # get covariance matrix
      covar.mat <-mat.sd %*% corr.mat %*% mat.sd
      
      fecSurv.draw <- mvrnorm(n=n.draws, mu=vec.mean, Sigma=covar.mat) %>% 
        as_tibble() %>% 
        rename(fec = V1, surv = V2) %>% 
        mutate(case = correlation)
      
      
#End with a dataframe of 100 values for survival and fecundity at each cv for each parameter

#Now, loop over the different values of survival and fecundity from the dataframes above and run a Leslie matrix with these values each time to extract lambda

lambda.temp.df <- NULL

#Loop over all values of survival
for(k in 1:nrow(fecSurv.draw)){

  batchSize <-  fecSurv.draw[1] * mean(num.mates) #Adult fecundity -- average fecundity per breeding (2.91) x average number of mates (2)
  
  #Dataframe of fecundity (female offspring only)
  fb <- batchSize/2 #Adult fecundity - female offspring only (assume equal sexes)
  
  #Dataframe of survival
  as <- fecSurv.draw[2] # CHANGED FROM 0.9; Adult survival
  
    adult.survival <- as[[1]][k]
    fecundity <- fb[[1]][k]
  
    
#Input to Leslie matrix
survival.vec <- c(YOY.survival, rep(juvenile.survival, times = juv.ages), rep(adult.survival, times = adult.ages - 1), 0)
fecund.vec <- c(rep(0, times = repro.age), rep(fecundity, times = maxAge - juv.ages))

Leslie_input <- data.frame(
  x = c(0:maxAge), #age
  sx = survival.vec, #survival
  mx = fecund.vec
)

A1_pre <- make_Leslie_matrix(Leslie_input)
#View(A1)
#A1_post <- pre_to_post(Amat = A1_pre, S0 = YOY.survival)
#View(A1_post)

#Calculate dominant eigenvalue (i.e. lambda) from transition matrix
lambda1 <- as_tibble(lambda1(A1_pre)) %>% 
  rename(lambda = value)

#lambda1(A1_post)
#stable_age <- mpmtools::stable_stage(A1_post)
lambda.temp.df <- lambda1 %>% 
  mutate(mean.survival = adult.survival, mean.fecundity = fecundity, correlation = correlation, surv.cv = surv.cv, fec.cv = fec.cv, iteration = k)

#Store lambda dataframe
lambda.df <- rbind(lambda.df, lambda.temp.df)

print(paste0("Finished iteration ", k, " with fecundity cv ", fec.cv, ", survival sv ", surv.cv, " and correlation ", correlation, "."))
    }
  }
}
  }

lambda.df %>% arrange(lambda) %>% 
  View()

lambda.plot.df <- lambda.df %>% 
  mutate(key = paste0("corr_", correlation, "_surv.cv_", surv.cv, "_fec.cv_", fec.cv))

low.lambda <- min(lambda.df$lambda)
high.lambda <- max(lambda.df$lambda)
mean.lambda <- mean(c(low.lambda, high.lambda))
sd.lambda <- sd(lambda.df$lambda)

#Seems like reasonable values for lambda span 0.95 - 1.05
#Examine different bounds of lambda
# lambda.plot.df <- rnorm(n = 10000, mean = mean.lambda, sd = sd.lambda) %>% 
#   as_tibble() %>% 
#   rename(lambda = value)


lambda.plot.df %>% filter(correlation == -0.75) %>%
  ggplot(aes(x=lambda, col = key)) +
  geom_density(alpha=0.6)


ggplot(fecSurv.draw, aes(x=as.numeric(fec), y=as.numeric(surv), col=case)) +
  geom_point()

ggplot(lambda.plot.df, aes(x=lambda)) +
  geom_density(alpha=0.6)



#Irrelevant
#--------------Population simulation-------------------
#Project population forward in time
# Abundance_year0 <- c(20000, rep(0, times = (max_age-1)))
# 
# Year1 <- A %*% Abundance_year0
# 
# nYears <- n_yrs        # set the number of years to project
# TMat <- A     # define the projection matrix
# InitAbund <- Abundance_year0    # define the initial abundance
# 
# ## NOTE: the code below can be re-used without modification:
# allYears <- matrix(0,nrow=nrow(TMat),ncol=nYears+1)     # build a storage array for all abundances!
# allYears[,1] <- InitAbund  # set the year 0 abundance                                    
# for(t in 2:(nYears+1)){   # loop through all years
#   allYears[,t] <-  TMat %*% allYears[,t-1]
# }
# 
# #in allYears, columns are years, rows are ages
# allYears <- data.frame(allYears)
# colnames(allYears) <- paste("Yr_", c(1:(nYears+1)))
# allYears <- allYears %>% mutate_all(round, digits=0)
# 
# #Sum mature ages for each year
# adult_truth <- allYears %>% slice(12:30) %>% 
#   summarise_all(sum)