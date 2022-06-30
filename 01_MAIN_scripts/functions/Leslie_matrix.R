#Examine lambda values with a Leslie matrix

#set so R doesn't use scientific notation
options("scipen"=100, "digits"=4)
library(tidyverse)
library(popbio)
#devtools::install_github("BruceKendall/mpmtools") #In case you don't have mpmtools
library(mpmtools)

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
f <- (1-Adult.survival)/(YOY.survival * juvenile.survival^11) # adult fecundity at equilibrium if no age truncation
ff <- f/init.prop.female * mating.periodicity/mean(num.mates) # female fecundity per breeding cycle
ff
ff <- ff*(1-non.conformists) #Change female fecundity per breeding cycle to account for non-conformists
ff

#Fecundity for Leslie matrix
leslie.fecundity <- (ff/2)/2 #Divide by 2 for biennial breeding; divide by 2 again for females only


#-----------------Leslie Matrix parameters--------------------
#Prep Leslie matrix input
survival.vec <- c(YOY.survival, rep(juvenile.survival, times = juv.ages), rep(adult.survival, times = adult.ages - 1), 0)
fecund.vec <- c(rep(0, times = repro.age), rep(leslie.fecundity, times = maxAge - juv.ages))

#Create dataframe for input to Leslie matrix
Leslie_input <- data.frame(
  x = c(0:maxAge), #age
  sx = survival.vec, #survival
  mx = fecund.vec
)

#Make Leslie matrix
A1_pre <- make_Leslie_matrix(Leslie_input)

#Calculate dominant eigenvalue (i.e. lambda) from transition matrix
(lambda1 <- as_tibble(lambda1(A1_pre)) %>% 
  rename(lambda = value))


#If wanting to convert to a post-breeding census Leslie matrix
A1_post <- pre_to_post(Amat = A1_pre, S0 = YOY.survival)
#View(A1_post)
