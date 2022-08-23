#Examine lambda values with a Leslie matrix

#set so R doesn't use scientific notation
options("scipen"=100, "digits"=4)
library(tidyverse)
library(popbio)
library(ggpubr)
library(MASS)
#devtools::install_github("BruceKendall/mpmtools") #In case you don't have mpmtools
library(mpmtools)

rm(list=ls())

#source("./Objective.3_life_history_sensitivity/functions/Obj3.functions.R") #If there's a way to translate CV to the beta distribution and still take correlated draws, then will need these functions

#-------------Age and survival settings--------------------------
YOY.survival <- 0.7 #young of year survival
juvenile.survival <- 0.8 # juvenile survival
adult.survival <- 0.825 # Adult survival
repro.age <- 12 # set age of reproductive maturity
max.age <- maxAge <- 50 #set the maximum age allowed in the simulation
juv.ages <- repro.age - 1 #Years of being a juvenile
ages <- c(0:maxAge) #Total ages
adult.ages <- length(ages) - (juv.ages + 1)

#-----------------------Fecundity and reproduction settings---------------------------
mating.periodicity <- 2 #number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
non.conformists <- 0.05 #proportion of off-year breeders to randomly include off their breeding cycle - want to change this to non.conformists
num.mates <- c(1:3) #vector of potential number of mates per mating
f <- (1-adult.survival)/(YOY.survival * juvenile.survival^11) # adult fecundity at equilibrium if no age truncation

#Fecundity for Leslie matrix
leslie.fecundity <- f #Divide by 2 for biennial breeding; divide by 2 again for females only



#-----------------CV and correlations------------------------
cv.vec <- c(seq(from = 0.05, to = 0.40, by = 0.01))
corr.vec <- c(0, -0.25, -0.5)
vec.mean <- c(YOY.survival, juvenile.survival, adult.survival, leslie.fecundity)   #vector of the means
n.draws <- 100 #Number of draws from a multivariate normal distribution
#Calculate parameters for beta distribution from mean and variance for survival


lambda.df <- NULL #Initialize dataframe to store lambda values
lambda.temp.df <- NULL

####Start loop over CV and correlations####
for(corr in 1:length(corr.vec)){
  correlation <- corr.vec[corr]
  
  for(js in 1:length(cv.vec)){
    YOY.juv_survival.cv <- cv.vec[js]

    YOY_survival.sd <- YOY.juv_survival.cv * YOY.survival
    juvenile_survival.sd <- YOY.juv_survival.cv * juvenile.survival
   
    for(as in 1:length(cv.vec)){
     adult_survival.cv <- cv.vec[as]
     adult_survival.sd <- adult_survival.cv * adult.survival
     
     for(fs in 1:length(cv.vec)){
       fecundity.cv <- cv.vec[fs]
       fecundity.sd <- fecundity.cv * leslie.fecundity
       
       #------------Convert survival means and cvs to beta distribution----------------
       #YOY
       # YOY.surv.betaParams <- estBetaParams(YOY.survival, YOY.survival.sd^2)
       # YOY.surv.alpha <- YOY.surv.betaParams[[1]]
       # YOY.surv.beta <- YOY.surv.betaParams[[2]]
       # 
       # #juveniles
       # juvenile.surv.betaParams <- estBetaParams(juvenile.survival, juvenile.survival.sd^2)
       # juvenile.surv.alpha <- juvenile.surv.betaParams[[1]]
       # juvenile.surv.beta <- juvenile.surv.betaParams[[2]]
       # 
       # #adults
       # adult.surv.betaParams <- estBetaParams(adult.survival, adult.survival.sd^2)
       # adult.surv.alpha <- adult.surv.betaParams[[1]]
       # adult.surv.beta <- adult.surv.betaParams[[2]]
       
       
       #----------Draw values from multivariate normal distribution------------#
       vec.sd <- c(YOY_survival.sd, juvenile_survival.sd, adult_survival.sd, fecundity.sd)
       matrix.sd <- diag(vec.sd)
       
       #Columns are in the same order as the rows
          corr.matrix <- matrix(c(1, 0, 0, 0, #row is YOY.survival
                               0, 1, 0, 0, #row is juvenile.survival
                               0, 0, 1, correlation, #row is adult.survival
                               0, 0, correlation, 1), #row is fecundity
                             nrow=4, ncol=4, byrow=T)

          # get covariance matrix
          covar.matrix <- matrix.sd %*% corr.matrix %*% matrix.sd
          colnames(covar.matrix) <- rownames(covar.matrix) <- c("YOY_survival", "juvenile_survival", "adult_survival", "fecundity")


                  fecSurv.draw <- mvrnorm(n=n.draws, mu=vec.mean, Sigma=covar.matrix) %>%
                    as_tibble() %>%
                    mutate(survival.fecundity_correlation = factor(correlation),
                           YOY_survival.cv = factor(YOY.juv.survival.cv),
                           juvenile_survival.cv = factor(YOY.juv.survival.cv),
                           adult_survival.cv = factor(adult.survival.cv),
                           fecund.cv = factor(fecundity.cv)) %>%
                    mutate(YOY_survival = ifelse(YOY_survival >=1, 0.99, YOY_survival),
                           juvenile_survival = ifelse(juvenile_survival >=1, 0.99, juvenile_survival), 
                           adult_survival = ifelse(adult_survival >=1, 0.99, adult_survival)) %>%  #Make sure that survival is never greater than 1
                    dplyr::select(YOY_survival, YOY_survival.cv,
                                  juvenile_survival, juvenile_survival.cv,
                                  adult_survival, adult_survival.cv,
                                  fecundity, fecund.cv) %>% 
                    dplyr::filter(YOY_survival < juvenile_survival & juvenile_survival < adult_survival)
                  
                  for(l in 1:nrow(fecSurv.draw)){
                    #-----------------Leslie Matrix parameters--------------------
                    #Prep Leslie matrix input
                    survival.draw.vec <- c(fecSurv.draw$YOY_survival[l], 
                                           rep(fecSurv.draw$juvenile_survival[l], times = juv.ages), 
                                           rep(fecSurv.draw$adult_survival[l], times = adult.ages - 1), 0)
                    
                    fecundity.draw.vec <- c(rep(0, times = repro.age), 
                                         rep(fecSurv.draw$fecundity[l], times = maxAge - juv.ages))
                    
                    #Create dataframe for input to Leslie matrix
                    Leslie_input <- data.frame(
                      x = c(0:maxAge), #age
                      sx = survival.draw.vec, #survival
                      mx = fecundity.draw.vec)
                    
                    #Make Leslie matrix
                    A1_pre <- make_Leslie_matrix(Leslie_input)
                    
                    #Calculate dominant eigenvalue (i.e. lambda) from transition matrix
                    (lambda1 <- as_tibble(lambda1(A1_pre)) %>% 
                        rename(lambda = value))
                    
                    lambda.temp.df <- lambda1 %>%
                      mutate(YOY.cv = factor(YOY.juv_survival.cv),
                             juvenile.cv = factor(YOY.juv_survival.cv),
                             adult.cv = factor(adult_survival.cv),
                             fecund.cv = factor(fecundity.cv)) %>%
                      mutate(YOY_survival = fecSurv.draw$YOY_survival[l],
                             juvenile_survival = fecSurv.draw$juvenile_survival[l], 
                             adult_survival = fecSurv.draw$adult_survival[l],
                             fecundity = fecSurv.draw$fecundity[l],
                             survival.fecundity_correlation = factor(correlation),
                             case = paste0(YOY.cv, "_", juvenile.cv, "_", adult.cv, "_fecCV_", fecund.cv, "_corr_", correlation),
                             iteration = l) %>% 
                      dplyr::select(lambda,
                                    YOY_survival, YOY.cv,
                                    juvenile_survival, juvenile.cv,
                                    adult_survival, adult.cv,
                                    fecundity, fecund.cv,
                                    survival.fecundity_correlation,
                                    case,
                                    iteration)
                    
                    #Store lambda dataframe
                    lambda.df <- rbind(lambda.df, lambda.temp.df)
                  } #Finish loop over fecSurv.draw
                  print(paste0("Working on fecundity CV ",fecundity.cv, "; adult survival ", adult_survival.cv, "; YOY survival ", YOY.juv_survival.cv, "; correlation ", correlation))
     } #Finish loop over fecundity CV
    } #Finish loop over adult survival cv
  } #Finish loop over juvenile/YOY survival CV
  
  print(paste0("Finished with correlation ", corr, " of ", length(corr.vec)))
} #Finish loop over correlations

saveRDS(lambda.df, file = "G://My Drive/Personal_Drive/R/CKMR/Objective.3_life_history_sensitivity/lambda.df")

lambda.df.summ <- lambda.df %>% group_by(YOY.cv, juvenile.cv, adult.cv, fecund.cv, survival.fecundity_correlation) %>%
  dplyr::filter(survival.fecundity_correlation == -0.25) %>% 
  summarize(lambda.sd = sd(lambda, na.rm = TRUE),
            lambda.mean = mean(lambda, na.rm = TRUE))

lambda.df.summ_filt <- lambda.df %>% inner_join(lambda.df.summ, by = c("YOY.cv", "juvenile.cv", "adult.cv", "fecund.cv", "survival.fecundity_correlation")) %>%
  mutate(z_score = (lambda - lambda.mean)/lambda.sd) %>% 
  dplyr::filter(z_score < 3 & z_score > -3) %>% 
  group_by(YOY.cv, juvenile.cv, adult.cv, fecund.cv, survival.fecundity_correlation) %>% 
  summarize(lambda.min = min(lambda),
            lambda.max = max(lambda))

#If wanting to convert to a post-breeding census Leslie matrix
A1_post <- pre_to_post(Amat = A1_pre, S0 = YOY.survival)
#View(A1_post)