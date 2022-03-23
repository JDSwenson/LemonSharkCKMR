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
f <- (1-adult.survival)/(YOY.survival * juvenile.survival^11) # adult fecundity at equilibrium if no age truncation
init.prop.female = 0.5
ff <- f/init.prop.female * mating.periodicity/mean(num.mates) # female fecundity per breeding cycle

batchSize <-  ff * mean(num.mates) #Adult fecundity -- average fecundity per breeding (2.91) x average number of mates (2)
        
#Dataframe of fecundity (female offspring only)
fecundity <- batchSize/2 #Adult fecundity - female offspring only (assume equal sexes)
        
        
#Input to Leslie matrix
survival.vec <- c(YOY.survival, rep(juvenile.survival, times = juv.ages), rep(adult.survival, times = adult.ages - 1), 0)

fecund.vec <- c(rep(0, times = repro.age), rep(fecundity, times = maxAge - juv.ages))
        
Leslie_input <- data.frame(
  x = c(0:maxAge), #age
  sx = survival.vec, #survival
  mx = fecund.vec
  )
        
A1_pre <- make_Leslie_matrix(Leslie_input)
A1_post <- pre_to_post(Amat = A1_pre, S0 = YOY.survival)
#View(A1_post)

View(A1_pre)
        

##Start here on 3/15/2022##
#Sensitivity and elasticity
#Examine which values are the highest in the output matrices; these will be the most sensitive values.
#I can convert to a dataframe/tibble, and convert to long format so that each value is associated with a value and parameter
#Then, I can make a vector corresponding to the different values of fecundity
#Make another vector corresponding to survival



####----------------Elasticity - Pre matrix------------------------####
A1.e <- elasticity(A1_post) %>% 
  as_tibble()

#Add age to column name
colnames(A1.e) <- paste0("age_", seq(from = 1, to = maxAge, by = 1))

#Convert first row (i.e. fecundity) to its own tidy dataframe for joining later
A1.e_fec.e <- A1.e[1,] %>% pivot_longer(cols = starts_with("age"), names_to = "age", values_to = "fecundity")

#Convert the rest of the elasticity matrix (i.e. survival) to its own tidy dataframe and join with the fecundity values
A1.e_tidy <- A1.e[-1,] %>% pivot_longer(cols = starts_with("age"), names_to = "age", values_to = "survival") %>% 
  dplyr::filter(survival > 0) %>% #Loads of survival values have 0s in the matrix
  full_join(A1.e_fec.e, by = "age") %>% 
  replace_na(list(survival = 0, fecundity = 0))

A1.e_tidy$age <- as.character(A1.e_tidy$age)
A1.e_tidy$age <- factor(A1.e_tidy$age, levels = unique(A1.e_tidy$age))

A1.e_tidy2 <- A1.e_tidy %>% pivot_longer(cols = !age,
                           names_to = "parameter",
                           values_to = "elasticity")

#Plot elasticity values
A1.e_tidy2 %>% ggplot(aes(x = age, y = elasticity, group = parameter)) + 
  geom_point(aes(colour = parameter)) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = -.15)) +
  labs(title = "Elasticity")



#----------Elasticity: post matrix
A1.ePost <- elasticity(A1_post) %>% 
  as_tibble()

#Add age to column name
colnames(A1.ePost) <- paste0("age_", seq(from = 1, to = maxAge, by = 1))

#Convert first row (i.e. fecundity) to its own tidy dataframe for joining later
A1.ePost_fec.e <- A1.ePost[1,] %>% pivot_longer(cols = starts_with("age"), names_to = "age", values_to = "fecundity")

#Convert the rest of the elasticity matrix (i.e. survival) to its own tidy dataframe and join with the fecundity values
A1.ePost_tidy <- A1.ePost[-1,] %>% pivot_longer(cols = starts_with("age"), names_to = "age", values_to = "survival") %>% 
  dplyr::filter(survival > 0) %>% #Loads of survival values have 0s in the matrix
  full_join(A1.ePost_fec.e, by = "age") %>% 
  replace_na(list(survival = 0, fecundity = 0))

A1.ePost_tidy$age <- as.character(A1.ePost_tidy$age)
A1.ePost_tidy$age <- factor(A1.ePost_tidy$age, levels = unique(A1.ePost_tidy$age))

A1.ePost_tidy2 <- A1.ePost_tidy %>% pivot_longer(cols = !age,
                                         names_to = "parameter",
                                         values_to = "elasticity")

#Plot elasticity values
A1.ePost_tidy2 %>% ggplot(aes(x = age, y = elasticity, group = parameter)) + 
  geom_point(aes(colour = parameter)) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = -.15)) +
  labs(title = "Elasticity")



####--------------------------------Sensitivity Pre-matrix------------------------####
#Sensitivity
A1.s <- sensitivity(A1_pre, zero = TRUE) %>% 
  as_tibble()

#Add age to column name
colnames(A1.s) <- paste0("age_", seq(from = 1, to = maxAge, by = 1))

#Convert first row (i.e. fecundity) to its own tidy dataframe for joining later
A1.s_fec.s <- A1.s[1,] %>% pivot_longer(cols = starts_with("age"), names_to = "age", values_to = "fecundity")

#Convert the rest of the elasticity matrix (i.e. survival) to its own tidy dataframe and join with the fecundity values
A1.s_tidy <- A1.s[-1,] %>% pivot_longer(cols = starts_with("age"), names_to = "age", values_to = "survival") %>% 
  dplyr::filter(survival > 0) %>% #Loads of survival values have 0s in the matrix
  full_join(A1.s_fec.s, by = "age") %>% 
  replace_na(list(survival = 0, fecundity = 0))

A1.s_tidy$age <- as.character(A1.s_tidy$age)
A1.s_tidy$age <- factor(A1.s_tidy$age, levels = unique(A1.s_tidy$age))

A1.s_tidy2 <- A1.s_tidy %>% pivot_longer(cols = !age,
                                         names_to = "parameter",
                                         values_to = "sensitivity")

#Plot sensitivity values
A1.s_tidy2 %>% ggplot(aes(x = age, y = sensitivity, group = parameter)) + 
  geom_point(aes(colour = parameter)) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = -.15)) +
  labs(title = "Sensitivity")




























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

