#Running a leslie matrix assuming 

#set so R doesn't use scientific notation
options("scipen"=100, "digits"=4)

#CV is the ratio between the standard deviation and the mean

#The goal here is to find a reasonable prior for lambda by setting a CV on survival and another CV on fecundity and playing around with the correlation. I can either run simulations using several of these values, or I can pick one and justify it in the text.

#-------------------Survival---------------------#
#Keep these values here
# maxAge <- maxAge
# repro.age <- repro.age
juv.ages <- repro.age - 1 #Years of being a juvenile
ages <- c(0:maxAge) #Total ages
adult.ages <- length(ages) - (juv.ages + 1)
YOY.survival <- YOY.survival # young of year survival
juvenile.survival <- juvenile.survival # juvenile survival
Adult.survival <- Adult.survival #adult survival


#-------------------Fecundity--------------------#
# mating.periodicity <- mating.periodicity # number of years between mating; assigned to an individual and sticks with them through their life. So they're either a one or two year breeder.
# num.mates <- num.mates # vector of potential number of mates per mating
# f <- f # adult fecundity at equilibrium if no age truncation
# init.prop.female = init.prop.female
# ff <- ff # female fecundity per breeding cycle
leslie.fecundity <- leslie.fecundity


#------------------Correlation-------------------#
#Established in main script
#corr.vec <- c(0.25, 0, -0.25, -0.5, -0.75)
#n.draws <- 100 #Number of draws from a multivariate normal distribution
#Needs to be run after above
vec.mean <- c(leslie.fecundity, YOY.survival, juvenile.survival, Adult.survival)   #vector of the means

#------------------Loop through different correlation values----------------
#Want to account for survival and fecundity potentially compensating for one another
lambda.df <- NULL #Initialize dataframe to store lambda values

#Set CV and SD for survival
YOY.surv.cv = YOY.surv.cv
YOY.survival.sd <- YOY.surv.cv * YOY.survival

adult.surv.cv = adult.surv.cv
adult.survival.sd <- surv.cv * Adult.survival #survival standard deviation


#Set CV and SD for fecundity
fec.cv <- fec.cv
fecundity.sd <- fec.cv * leslie.fecundity #Fecundity per breeding cycle i.e. ff

#create matrix of SD values
vec.sd <- c(fecundity.sd, YOY.survival.sd, adult.survival.sd, adult.survival.sd) #Assume the same cv and sd for juvenile and adult survival

mat.sd <- diag(vec.sd)  

  corr.mat <- matrix(c(1, fec.surv.corr, fec.surv.corr, fec.surv.corr,
                       fec.surv.corr, 1, survival.corr, survival.corr,
                       fec.surv.corr, survival.corr, 1, survival.corr,
                       fec.surv.corr, survival.corr, survival.corr, 1), nrow=4, ncol=4, byrow=T)
  # get covariance matrix
  covar.mat <- mat.sd %*% corr.mat %*% mat.sd
      
  fecSurv.draw <- mvrnorm(n=n.draws, mu=vec.mean, Sigma=covar.mat) %>% 
    as_tibble() %>% 
    rename(fec = V1, YOY.surv = V2, juv.surv = V3, adult.surv = V4) %>% 
    mutate(case = correlation) %>% 
    mutate(YOY.surv = ifelse(YOY.surv >=1, 0.99, YOY.surv),
           juv.surv = ifelse(juv.surv >=1, 0.99, juv.surv),
           adult.surv = ifelse(adult.surv >= 1, 0.99, adult.surv)) #Make sure that survival is never greater than 1
      
fecSurv.draw <- fecSurv.draw %>% dplyr::filter(juv.surv > YOY.surv & adult.surv > juv.surv)      
#End with a dataframe of 100 values for survival and fecundity at each cv for each parameter

      
#Now, loop over the different values of survival and fecundity from the dataframes above and run a Leslie matrix with these values each time to extract lambda
lambda.temp.df =lambda.df <- NULL

#Loop over all values of survival
for(k in 1:nrow(fecSurv.draw)){

  batchSize <-  fecSurv.draw[1] * mean(num.mates) #Adult fecundity -- average fecundity per breeding (2.91) x average number of mates (2)
  
  #Dataframe of fecundity (female offspring only)
  fb <- batchSize/2 #Adult fecundity - female offspring only (assume equal sexes)
  
  #Dataframe of survival
  ys <- fecSurv.draw[2]
  js <- fecSurv.draw[3]
  as <- fecSurv.draw[4] # CHANGED FROM 0.9; Adult survival
  
  YOY.survival <- ys[[1]][k]
  juvenile.survival <- js[[1]][k]
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
(lambda1 <- as_tibble(lambda1(A1_pre)) %>% 
  rename(lambda = value))

#lambda1(A1_post)
#stable_age <- mpmtools::stable_stage(A1_post)
lambda.temp.df <- lambda1 %>% 
  mutate(Adult.survival = adult.survival, leslie.fecundity = fecundity, correlation = correlation, scv = surv.cv, fcv = fec.cv, iteration = k)

#Store lambda dataframe
lambda.df <- rbind(lambda.df, lambda.temp.df)

}
print(paste0("Finished iteration ", k, " with fecundity cv ", fec.cv, ", survival cv ", surv.cv, " and correlation ", correlation, "."))


#Extract range of lambda values
(low.lambda <- min(lambda.df$lambda))
(high.lambda <- max(lambda.df$lambda))
(mean.lambda <- mean(lambda.df$lambda))
(lambda.sd <- sd(lambda.df$lambda))