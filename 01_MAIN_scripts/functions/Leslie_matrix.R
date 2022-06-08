rm(list=ls())

#set so R doesn't use scientific notation
options("scipen"=100, "digits"=4)

library(tidyverse)
library(ggpubr)
library(MASS)
library(popbio)
#devtools::install_github("BruceKendall/mpmtools")
library(mpmtools)

#The goal here is to find a reasonable prior for lambda by setting a CV on survival and another CV on fecundity and playing around with the correlation. I can either run simulations using several of these values, or I can pick one and justify it in the text.

#-----------------Leslie Matrix parameters--------------------
#First, set the common values that will be used in the Leslie matrix
cv <- c(0.05, 0.1)

#-------------------Survival---------------------#
#Keep these values here
maxAge <- maxAge
repro.age <- repro.age
juv.ages <- repro.age - 1 #Years of being a juvenile
ages <- c(0:maxAge) #Total ages
adult.ages <- length(ages) - (juv.ages + 1)
YOY.survival <- YOY.survival # young of year survival
juvenile.survival <- juvenile.survival # juvenile survival
mean.survival <- Adult.survival #adult survival


#-------------------Fecundity--------------------#
mating.periodicity <- mating.periodicity # number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
num.mates <- num.mates # vector of potential number of mates per mating
f <- f # adult fecundity at equilibrium if no age truncation
init.prop.female = init.prop.female
ff <- ff # female fecundity per breeding cycle
mean.fecundity <- ff


#------------------Correlation-------------------#
corr.vec <- c(0.25, 0, -0.25, -0.5, -0.75)
vec.mean <- c(mean.fecundity, mean.survival)   #vector of the means
n.draws <- 100 #Number of draws from a multivariate normal distribution

#------------------Loop through different correlation values----------------
lambda.df <- NULL #Initialize dataframe to store lambda values

for(c in 1:length(corr.vec)){
  correlation <- corr.vec[c]

#Loop through different cv values for survival
  for(s in 1:length(cv)){
    #mean.survival <- surv.vec[s]
    surv.cv <- cv[s]
    sd.survival <- surv.cv * mean.survival #survival standard deviation
    

  #Loop through different cv values for fecundity
    for(fec in 1:length(cv)){
      #mean.fecundity <- fec.vec[fec]
      fec.cv <- cv[fec]
      sd.fecundity <- fec.cv * mean.fecundity #Fecundity per breeding cycle i.e. ff
      
      
      #create sd matrix
      vec.sd <- c(sd.fecundity, sd.survival)
      mat.sd <- diag(vec.sd)  
      
      
      #create correlation matrix (assuming correlation = -0.5) ====
      corr.fec.surv <- correlation
      corr.mat <- matrix(c(1, corr.fec.surv, corr.fec.surv, 1), nrow=2, ncol=2, byrow=T)
      # get covariance matrix
      covar.mat <- mat.sd %*% corr.mat %*% mat.sd
      
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
  mutate(mean.survival = adult.survival, mean.fecundity = fecundity, correlation = correlation, scv = surv.cv, fcv = fec.cv, iteration = k)

#Store lambda dataframe
lambda.df <- rbind(lambda.df, lambda.temp.df)

print(paste0("Finished iteration ", k, " with fecundity cv ", fec.cv, ", survival sv ", surv.cv, " and correlation ", correlation, "."))
    }
  }
}
  } #End loop

#----------------Analyze and plot output----------------------

#---------Summarize output---------------#
#View lambda values
lambda.df %>% arrange(lambda) %>% 
  View()

#Extract range of lambda values
low.lambda <- min(lambda.df$lambda)
high.lambda <- max(lambda.df$lambda)
mean.lambda <- mean(lambda.df$lambda)
sd.lambda <- sd(lambda.df$lambda)

#Check mean values of lambda when survival and fecundity cv are the same
lambda.df %>% dplyr::filter(scv == fcv) %>% 
  group_by(correlation, scv) %>% 
  summarize(mean(lambda))


#Check mean values of lambda when survival and fecundity cv are not the same
lambda.df %>% dplyr::filter(scv != fcv) %>% 
  group_by(correlation, scv, fcv) %>% 
  summarize(mean(lambda))


#---------Plot output---------------
#Make dataframe for plotting
lambda.plot.df <- lambda.df %>% 
  mutate(key = paste0("Surv.CV: ", scv, "; Fec.CV: ", fcv))

densPlot.list1 <- list()
densPlot.list2 <- list()

for(c in 1:length(corr.vec)){
  corr <- corr.vec[c]
  
  densPlot.list1[[c]] <- lambda.plot.df %>% dplyr::filter(scv == fcv & correlation == corr) %>%
  ggplot(aes(x=lambda, fill = factor(scv))) +
  geom_density(alpha=0.4) + 
  xlim(0.9, 1.1) + 
    labs(fill = "CV") +
  ggtitle(paste0("Correlation = ", corr))

  
  densPlot.list2[[c]] <- lambda.plot.df %>% dplyr::filter(scv != fcv & correlation == corr) %>%
    ggplot(aes(x=lambda, fill = key)) +
    geom_density(alpha=0.4) + 
    xlim(0.9, 1.1) + 
    labs(fill = "CV") +
    ggtitle(paste0("Correlation = ", corr))
}


#Save to pdf
today <- format(Sys.Date(), "%d%b%Y")
prior_plots_location <- "G://My Drive/Personal_Drive/R/CKMR/Sensitivity.tests/Prior_sensitivity/figures/"
prior_purpose <- "test_CV_on_fecundity_and_survival"
priorPlots.file <- paste0(prior_plots_location, prior_purpose, ".pdf")

pdf(file = priorPlots.file) #Open pdf

#Print to pdf
ggarrange(plotlist = densPlot.list1, common.legend = TRUE, legend = "right") %>% 
  annotate_figure(top = text_grob("Fecundity CV = Survival CV", face = "bold", size = 14))
ggarrange(plotlist = densPlot.list2, common.legend = TRUE) %>% 
  annotate_figure(top = text_grob("Fecundity CV != Survival CV", face = "bold", size = 14))

#Close pdf
dev.off()


ggplot(fecSurv.draw, aes(x=as.numeric(fec), y=as.numeric(surv), col=case)) +
  geom_point()




####---------------------Sensitivity and elasticity------------------------------####
#-------------------Survival---------------------#
maxAge <- 50 
repro.age <- 12
juv.ages <- repro.age - 1 #Years of being a juvenile
ages <- c(0:maxAge) #Total ages
adult.ages <- length(ages) - (juv.ages + 1)
YOY.survival <- 0.7 # young of year survival
juvenile.survival <- 0.8 # juvenile survival
adult.survival <- 0.825 #adult survival


#-------------------Fecundity--------------------#
mating.periodicity <- 1 # number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
num.mates <- c(1:3) # vector of potential number of mates per mating
f <- (1-mean.survival)/(YOY.survival * juvenile.survival^11) # adult fecundity at equilibrium if no age truncation
init.prop.female = 0.5
ff <- f/init.prop.female * mating.periodicity/mean(num.mates) # female fecundity per breeding cycle
fecundity <- ff


#Input to Leslie matrix
survival.vec <- c(YOY.survival, rep(juvenile.survival, times = juv.ages), rep(adult.survival, times = adult.ages - 1), 0)
fecund.vec <- c(rep(0, times = repro.age), rep(fecundity, times = maxAge - juv.ages))

Leslie_input <- data.frame(
  x = c(0:maxAge), #age
  sx = survival.vec, #survival
  mx = fecund.vec
)

A1_pre <- make_Leslie_matrix(Leslie_input)

sensitivity(A1_pre)

#----------------------Load life history info----------------------
elh <- read_csv("00_key_resources/Elasmobranch_life_history.csv") %>% 
  mutate(litter.size = as.numeric(litter.size)) %>% #Read in file with elasmobranch life history traits; convert litter size column from character to numeric
  mutate(annual.fecundity = litter.size/periodicity) %>% 
mutate(age.mat.bin = ifelse(age.at.maturity <= 3, "early", 
                                    ifelse(age.at.maturity > 3 & age.at.maturity <= 7, "moderate", 
                                           "late"))) %>% 
  mutate(fecundity.bin = ifelse(annual.fecundity <= 5, "low", 
                                  ifelse(annual.fecundity > 5 & annual.fecundity <= 10, "moderate"
                                         , "high"))) %>%
  as.data.frame()

elh2 <- elh %>% dplyr::select(Species, `Common name`, age.at.maturity, litter.size, periodicity, annual.fecundity, fecundity.bin, age.mat.bin)

elh2 %>% group_by(fecundity.bin) %>% 
  summarize(n())

elh2 %>% group_by(age.mat.bin, fecundity.bin) %>% #View numbers in each category
  summarize(n())

elh2 %>% filter(annual.fecundity < 2)

elh2 %>% filter(`Common name` == "whale shark")

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